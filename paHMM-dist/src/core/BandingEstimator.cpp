//==============================================================================
// Pair-HMM phylogenetic tree estimator
// 
// Copyright (c) 2015-2019 Marcin Bogusz.
//               2020 Mazen Mardini for library wrapper compatibility.
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses>.
//==============================================================================


#include "core/BandingEstimator.hpp"
#include "models/GTRModel.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"
#include "models/NegativeBinomialGapModel.hpp"
#include "hmm/DpMatrixFull.hpp"

namespace EBC
{


BandingEstimator::BandingEstimator(Definitions::AlgorithmType at, Sequences* inputSeqs, Definitions::ModelType model ,std::vector<double> indel_params,
        std::vector<double> subst_params, Definitions::OptimizationType /*ot*/, unsigned int rateCategories, double alpha, GuideTree* g) :
                inputSequences(inputSeqs), gt(g), algorithm(at), gammaRateCategories(rateCategories),
                /*hmms(pairCount), bands(pairCount),*/ pairCount(inputSequences->getPairCount()), divergenceTimes(pairCount, NAN)
{
	//Banding estimator means banding enabled!

	DEBUG("Starting Banding Estimator");
	maths = new Maths();
	dict = inputSequences->getDictionary();

	//Helper models
	if (model == Definitions::ModelType::GTR)
	{
		substModel = new GTRModel(dict, maths,gammaRateCategories);
	}
	else if (model == Definitions::ModelType::HKY85)
	{
		substModel = new HKY85Model(dict, maths,gammaRateCategories);
	}
	else if (model >= Definitions::ModelType::LG)
	{
		switch(model){
			case Definitions::ModelType::LG :
				substModel = new AminoacidSubstitutionModel(dict, maths,gammaRateCategories,Definitions::aaLgModel);
				DEBUG("Using LG model");
			break;
			case Definitions::ModelType::JTT :
				substModel = new AminoacidSubstitutionModel(dict, maths,gammaRateCategories,Definitions::aaJttModel);
				DEBUG("Using JTT model");
			break;
			case Definitions::ModelType::WAG :
				substModel = new AminoacidSubstitutionModel(dict, maths,gammaRateCategories,Definitions::aaWagModel);
				DEBUG("Using WAG model");
			break;
			case Definitions::ModelType::GTR :
			  throw HmmException("Not implemented");
			  
			case Definitions::ModelType::HKY85 :
			  throw HmmException("Not implemented");
					}
	}

	indelModel = new NegativeBinomialGapModel();

	estimateSubstitutionParams = false;
	estimateIndelParams = false;
	estimateAlpha = false;

	//pairwise comparison mode
	modelParams = new OptimizedModelParameters(substModel, indelModel,2, 1, estimateSubstitutionParams,
			estimateIndelParams, estimateAlpha, true, maths);

	modelParams->boundDivergenceBasedOnLambda(indel_params[0]);

	if(!estimateIndelParams)
		modelParams->setUserIndelParams(indel_params);
	if(!estimateSubstitutionParams)
		modelParams->setUserSubstParams(subst_params);
	modelParams->setAlpha(alpha);

	substModel->setObservedFrequencies(inputSequences->getElementFrequencies());
	if (estimateSubstitutionParams == false)
	{
		//set parameters and calculate the model
		substModel->setAlpha(modelParams->getAlpha());
		substModel->setParameters(modelParams->getSubstParameters());
		substModel->calculateModel();
	}

	if (estimateIndelParams == false)
	{
            //set parameters and calculate the model
            indelModel->setParameters(modelParams->getIndelParameters());
	}


    //EvolutionaryPairHMM *hmm;

    numopt = new BrentOptimizer(modelParams, nullptr);

}

BandingEstimator::~BandingEstimator()
{
  delete numopt;
  delete modelParams;
  delete maths;
  delete indelModel;
  delete substModel;
}

void BandingEstimator::optimizePairByPair()
{
	for(unsigned int i =0; i< pairCount; i++)
	{
        optimizePair(i);
	}

	INFO("Optimized divergence times:");
    INFO(this->divergenceTimes);
}

double BandingEstimator::optimizePair(int i)
{
    if (!std::isnan(this->divergenceTimes[i])) {
        return this->divergenceTimes[i];
    }

    EvolutionaryPairHMM* hmm;
    Band* band;
    DistanceMatrix* dm = gt->getDistanceMatrix();
    PairHmmCalculationWrapper* wrapper = new PairHmmCalculationWrapper();
    double result;

    DEBUG("Optimizing distance for pair #" << i);
    std::pair<unsigned int, unsigned int> idxs = inputSequences->getPairOfSequenceIndices(i);
    INFO("Running pairwise calculator for sequence id " << idxs.first << " and " << idxs.second
            << " ,number " << i+1 <<" out of " << pairCount << " pairs" );
    BandCalculator* bc = new BandCalculator(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
            substModel, indelModel, dm->getDistance(idxs.first,idxs.second));
    band = bc->getBand();
    if (algorithm == Definitions::AlgorithmType::Viterbi)
    {
        DEBUG("Creating Viterbi algorithm to optimize the pairwise divergence time...");
        hmm = new ViterbiPairHMM(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
                substModel, indelModel, Definitions::DpMatrixType::Full, band);
    }
    else if (algorithm == Definitions::AlgorithmType::Forward)
    {
        DEBUG("Creating forward algorithm to optimize the pairwise divergence time...");
        hmm = new ForwardPairHMM(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
                substModel, indelModel, Definitions::DpMatrixType::Full, band);
    }

    //hmm->setDivergenceTimeAndCalculateModels(modelParams->getDivergenceTime(0)); //zero as there's only one pair!

    //LikelihoodSurfacePlotter lsp;
    //lsp.setTargetHMM(hmm);
    //lsp.getLikelihoodSurface();

    wrapper->setTargetHMM(hmm);
    DUMP("Set model parameter in the hmm...");
    wrapper->setModelParameters(modelParams);
    modelParams->setUserDivergenceParams({bc->getClosestDistance()});
    numopt->setTarget(wrapper);
    numopt->setAccuracy(bc->getBrentAccuracy());
    numopt->setBounds(bc->getLeftBound(), bc->getRightBound() < 0 ? modelParams->divergenceBound : bc->getRightBound());

    result = numopt->optimize() * -1.0;
    DEBUG("Likelihood after pairwise optimization: " << result);
    if (result <= (Definitions::minMatrixLikelihood /2.0))
    {
        DEBUG("Optimization failed for pair #" << i << " Zero probability FWD");
        band->output();
        dynamic_cast<DpMatrixFull*>(hmm->M->getDpMatrix())->outputValuesWithBands(band->getMatchBand() ,band->getInsertBand(),band->getDeleteBand(),'|', '-');
        dynamic_cast<DpMatrixFull*>(hmm->X->getDpMatrix())->outputValuesWithBands(band->getInsertBand(),band->getMatchBand() ,band->getDeleteBand(),'\\', '-');
        dynamic_cast<DpMatrixFull*>(hmm->Y->getDpMatrix())->outputValuesWithBands(band->getDeleteBand(),band->getMatchBand() ,band->getInsertBand(),'\\', '|');
    }

    delete band;
    delete bc;
    delete hmm;

    return modelParams->getDivergenceTime(0);
}


double BandingEstimator::runIteration()
{
    double result = 0;
    //double tmp;
    return result;
}

void BandingEstimator::outputDistanceMatrix(stringstream& ss)
{
	unsigned int count, pairCount;
	count = this->inputSequences->getSequenceCount();
	pairCount = this->inputSequences->getPairCount();

	ss << "\t" << this->inputSequences->getSequenceCount() << endl;

	for (unsigned int i = 0; i< count; i++)
	{
		ss << "S" << i << " ";
		for(unsigned int j=0; j< count; j++)
		{
			ss << this->modelParams->getDistanceBetween(i,j) << " ";
		}
		ss << endl;
	}
}

} /* namespace EBC */


