//==============================================================================
// Pair-HMM phylogenetic tree estimator
// 
// Copyright (c) 2015 Marcin Bogusz.
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

#include "core/CommandReader.hpp"
#include "core/Sequences.hpp"
#include "core/HmmException.hpp"
#include "core/BandingEstimator.hpp"
#include "core/BioNJ.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "heuristics/ModelEstimator.hpp"
#include <array>
#include <ctime>


#include "core/OptimizedModelParameters.hpp"
#include "hmm/ForwardPairHMM.hpp"
#include "hmm/ViterbiPairHMM.hpp"


using namespace std;
using namespace EBC;

int main(int argc, char ** argv) {


	cout << fixed << setprecision(8);
	cerr << fixed << setprecision(8);

	try
	{
		CommandReader* cmdReader = new CommandReader(argc, argv);
		ofstream distfile;

		INFO("Reading input sequences...");
		IParser* parser = cmdReader->getParser();

		//Remove gaps if the user provides a MSA file
		bool removeGaps = true;

		Sequences* inputSeqs = new Sequences(parser, cmdReader->getSequenceType(),removeGaps);

		INFO("Creating Model Parameters heuristics...");

		cout << "Estimating evolutionary model parameters..." << endl;

		ModelEstimator* tme = new ModelEstimator(inputSeqs, cmdReader->getModelType(),
				cmdReader->getOptimizationType(), cmdReader->getCategories(), cmdReader->getAlpha(),
				cmdReader->estimateAlpha());

		vector<double> indelParams;
		vector<double> substParams;
		double alpha = cmdReader->getAlpha();

		substParams = tme->getSubstitutionParameters();
		indelParams = tme->getIndelParameters();

		if(cmdReader->estimateAlpha()){
			alpha = tme->getAlpha();
		}


		try{
			substParams = cmdReader->getSubstParams();
		}
		//do nothing - if exception, this means no user-specified params
		catch(HmmException& pe){
			substParams = tme->getSubstitutionParameters();
		}

		try{
			indelParams = cmdReader->getIndelParams();
		}
		//do nothing - if exception, this means no user-specified params
		catch(HmmException& pe){
			indelParams = tme->getIndelParameters();
		}

		cout << "Estimating pairwise distances..." << endl;

		BandingEstimator* be = new BandingEstimator(Definitions::AlgorithmType::Forward, inputSeqs, cmdReader->getModelType() ,indelParams,
				substParams, cmdReader->getOptimizationType(), cmdReader->getCategories(),alpha, tme->getGuideTree());
		be->optimizePairByPair();


		auto distances = be->getOptimizedTimes();
		auto seqCount =  inputSeqs->getSequenceCount();



		//output distance matrix
		distfile.open((string(cmdReader->getInputFileName()).append(Definitions::distMatExt)).c_str(),ios::out);
		distfile << inputSeqs->getSequenceCount() << endl;
		for (unsigned int seqId = 0; seqId < seqCount; seqId++){
			distfile << inputSeqs->getSequenceName(seqId) << "        ";
			for(unsigned int j = 0; j<seqId; j++)
			{

				distfile << " " << distances[(seqId - j - 1) + (j*seqCount) - (((1+j)/2.0)*(j*1.0))];
			}
			distfile << endl;
		}
		distfile.close();

	}
	catch(HmmException& pe)
	{
		ERROR(pe.what());
	}
	catch(exception &ex)
	{
		ERROR(ex.what());
	}

	return 0;
}
