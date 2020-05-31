//==============================================================================
// Pair-HMM phylogenetic tree estimator
//
// Copyright (c) 2015-2019 Marcin Bogusz. 2020 Mazen Mardini.
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
#include "core/Sequences.hpp"
#include "core/BrentOptimizer.hpp"
#include "core/Definitions.hpp"
#include "heuristics/ModelEstimator.hpp"
#include "StreamParser.hpp"

#include <sstream>
#include <iomanip>

#include "cpahmm.h"

using namespace std;
using namespace EBC;

int test_cpp_api()
{
    // REQUIRED INPUT:
    Definitions::ModelType model = Definitions::ModelType::JTT;

    double indel_NB_probability = 0.0;
    double indel_rate = 0.0;

    // Let ModelEstimator::getIndelParameters() decide
    bool estimate_indel_params = true;

    double hky85_param = 0.0;

    double gtr_param1 = 0.0;
    double gtr_param2 = 0.0;
    double gtr_param3 = 0.0;
    double gtr_param4 = 0.0;
    double gtr_param5 = 0.0;

    // Let ModelEstimator::getSubstitutionParameters() decide
    bool estimate_model_params = true;

    double alpha = 0.5;

    // Let ModelEstimator::getAlpha() decide
    bool estimate_alpha = true;

    unsigned int gamma_rate_categories = 4;

    stringstream inputStream(R"SEQ(>H0
ENVVDDTSDRPTICQKWNTTSAAISKYDFLSFYPHYRPASVETFLNLLLK
>H4
ENVVDDKSDRPTICQKWNATSAAISKYNFLEFYPHVRTASVEMFLNLLLK
>H21
SPATQSSKDDALLSMAATVGEASLDKRSHIFSFPSMHVRTVTSDLSGLAF
>H26
SSLTQSSKDDEILSMIAIVGDACIDWRSHIVSFSYIHVLTVTSNLSGINF
>H35
SKASQENKTDQLLKRDAIVGEACIDKKKHNFGYKSVRVRSVTTNLAGLAF
)SEQ");

    // AUXILIARY VALUES:
    Definitions::SequenceType sequenceType;
    if (model == Definitions::ModelType::GTR || model == Definitions::ModelType::HKY85) {
        sequenceType = Definitions::SequenceType::Nucleotide;
    } else {
        sequenceType = Definitions::SequenceType::Aminoacid;
    }

    // SURPRESS UNWANTED MESSAGES:
    FileLogger::InfoLogger().active = false;
    FileLogger::WarningLogger().active = false;
    FileLogger::DebugLogger().active = false;

    // CALCULATE:
    vector<double> indelParams;
    vector<double> substParams;

    StreamParser parser(inputStream);
    Sequences inputSeqs(&parser, sequenceType, true);

    // Caluclate (among other things) a rough estimate of O(n*m) distances.
    // Optimizing this is a challange for another day...
    ModelEstimator tme(&inputSeqs, model, Definitions::OptimizationType::BFGS,
                       gamma_rate_categories, alpha, estimate_alpha);

    if(estimate_alpha){
        alpha = tme.getAlpha();
    }

    if (estimate_indel_params) {
        indelParams = tme.getIndelParameters();
    } else {
        indelParams = {indel_NB_probability, indel_rate};
    }

    if (estimate_model_params) {
        substParams = tme.getSubstitutionParameters();
    } else {
        if (model == Definitions::ModelType::HKY85) {
            substParams = {hky85_param};
        } else if (model == Definitions::ModelType::GTR) {
            substParams = {gtr_param1, gtr_param2, gtr_param3, gtr_param4, gtr_param5};
        }
    }

    // Estimating pairwise distances...
    BandingEstimator be(Definitions::AlgorithmType::Forward, &inputSeqs, model ,indelParams,
                        substParams, Definitions::OptimizationType::BFGS, gamma_rate_categories,
                        alpha, tme.getGuideTree());

    be.optimizePairByPair();
    auto distances = be.getOptimizedTimes();
    auto seqCount = inputSeqs.getSequenceCount();

    // OUTPUT DISTANCE MATRIX:
    ostringstream distfile;
    distfile << seqCount << "\n";
    for (unsigned int seqId = 0; seqId < seqCount; seqId++)
    {
        distfile << inputSeqs.getSequenceName(seqId) << "        ";
        for(unsigned int j = 0; j<seqId; j++)
        {
            distfile << " " << distances[(seqId - j - 1) + (j*seqCount) - (((1+j)/2)*j)];
        }
        distfile << "\n";
    }

    cout << distfile.str() << endl;

 /*
    cout << fixed << setprecision(8);
    cerr << fixed << setprecision(8);

    try
    {
    }
    catch(HmmException& pe)
    {
        ERROR(pe.what());
    }
    catch(exception &ex)
    {
        ERROR(ex.what());
    }
*/
    return 0;
}

int test_c_api()
{
    EBCBandingEstimator *be = ebc_be_create();
    /*ebc_be_set_input(be, R"SEQ(>H0
ENVVDDTSDRPTICQKWNTTSAAISKYDFLSFYPHYRPASVETFLNLLLK
>H4
ENVVDDKSDRPTICQKWNATSAAISKYNFLEFYPHVRTASVEMFLNLLLK
>H21
SPATQSSKDDALLSMAATVGEASLDKRSHIFSFPSMHVRTVTSDLSGLAF
>H26
SSLTQSSKDDEILSMIAIVGDACIDWRSHIVSFSYIHVLTVTSNLSGINF
>H35
SKASQENKTDQLLKRDAIVGEACIDKKKHNFGYKSVRVRSVTTNLAGLAF
)SEQ");*/
    ebc_be_set_input_from_file(be, "nonsense50x70.fasta");

    EBCSequences *seq = ebc_be_execute_jtt_model(be);
    cout << ebc_be_last_error_msg(be) << endl;
    unsigned int seqCount = ebc_seq_count(seq);

    cout << seqCount << endl;
    for (unsigned int seqId1 = 0; seqId1 < seqCount; seqId1++)
    {
        cout << ebc_seq_get_name(seq, seqId1) << "        ";
        for (unsigned int seqId2 = 0; seqId2 < seqId1; seqId2++)
        {
            cout << ebc_seq_get_distance(seq, seqId1, seqId2) << " ";
        }
        cout << endl;
    }

    // Clean up
    cout.flush();
    ebc_seq_free(seq);
    ebc_be_free(be);

    // Exit
    return 0;
}


int main(int /*argc*/, char ** /*argv*/)
{
    return test_c_api();
    //return test_cpp_api();
}
