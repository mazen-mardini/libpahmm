//==============================================================================
// libpahmm - library for paHMM-Tree, a phylogenetic tree estimator
//
// Copyright (c) 2020 Mazen Mardini.
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

#include "cpahmm.h"
#include "cpahmm_p.h"
#include <cstring>
#include <sstream>

#include "core/BandingEstimator.hpp"
#include "core/Sequences.hpp"
#include "core/Definitions.hpp"
#include "heuristics/ModelEstimator.hpp"
#include "StreamParser.hpp"

using namespace std;
using namespace EBC;

extern "C"
{

#define EBC_BE_DEFAULTS_INDEL_NB_PROBABILITY 0.0
#define EBC_BE_DEFAULTS_INDEL_RATE 0.0
#define EBC_BE_DEFAULTS_ALPHA 0.5
#define EBC_BE_DEFAULTS_GAMMA_RATE_CATEGORIES 4


EBCBandingEstimator *ebc_be_create()
{
    EBCBandingEstimator *be = new EBCBandingEstimator;
    be->_parser = nullptr;
    be->_error = new HmmException;
    be->indel_NB_probability = EBC_BE_DEFAULTS_INDEL_NB_PROBABILITY;
    be->indel_rate = EBC_BE_DEFAULTS_INDEL_RATE;
    be->alpha = EBC_BE_DEFAULTS_ALPHA;
    be->gamma_rate_categories = EBC_BE_DEFAULTS_GAMMA_RATE_CATEGORIES;
    be->estimate_indel_params = true;
    be->estimate_alpha = true;
    be->estimate_categories = true;

    return be;
}

void ebc_be_free(EBCBandingEstimator *be)
{
    if (!be) {
        return;
    }

    if (be->_parser) {
        delete reinterpret_cast<EBC::StreamParser *>(be->_parser);
    }

    if (be->_error) {
        delete reinterpret_cast<HmmException *>(be->_error);
    }

    delete be;
}

const char *ebc_be_last_error_msg(EBCBandingEstimator *be)
{
    return reinterpret_cast<HmmException *>(be->_error)->what();
}

EBCSequences *ebc_be_execute_gtr_model(EBCBandingEstimator *be, double param1, double param2,
                                       double param3, double param4, double param5)
{
    try {
        return ebc_seq_create(be, Definitions::ModelType::GTR, false, 5,
                              param1, param2, param3, param4, param5);
    } catch (HmmException& e) {
        *reinterpret_cast<HmmException *>(be->_error) = e;
        return nullptr;
    }
}

EBCSequences *ebc_be_execute_gtr_modelv2(EBCBandingEstimator *be)
{
    try {
        return ebc_seq_create(be, Definitions::ModelType::GTR, true, 0);
    } catch (HmmException& e) {
        *reinterpret_cast<HmmException *>(be->_error) = e;
        return nullptr;
    }
}

EBCSequences *ebc_be_execute_hky85_model(EBCBandingEstimator *be, double param)
{
    try {
        return ebc_seq_create(be, Definitions::ModelType::HKY85, false, 1, param);
    } catch (HmmException& e) {
        *reinterpret_cast<HmmException *>(be->_error) = e;
        return nullptr;
    }
}

EBCSequences *ebc_be_execute_hky85_modelv2(EBCBandingEstimator *be)
{
    try {
        return ebc_seq_create(be, Definitions::ModelType::HKY85, true, 0);
    } catch (HmmException& e) {
        *reinterpret_cast<HmmException *>(be->_error) = e;
        return nullptr;
    }
}

EBCSequences *ebc_be_execute_jtt_model(EBCBandingEstimator *be)
{
    try {
        return ebc_seq_create(be, Definitions::ModelType::JTT, true, 0);
    } catch (HmmException& e) {
        *reinterpret_cast<HmmException *>(be->_error) = e;
        return nullptr;
    }
}

EBCSequences *ebc_be_execute_lg_model(EBCBandingEstimator *be)
{
    try {
        return ebc_seq_create(be, Definitions::ModelType::LG, true, 0);
    } catch (HmmException& e) {
        *reinterpret_cast<HmmException *>(be->_error) = e;
        return nullptr;
    }
}

EBCSequences *ebc_be_execute_wag_model(EBCBandingEstimator *be)
{
    try {
        return ebc_seq_create(be, Definitions::ModelType::WAG, true, 0);
    } catch (HmmException& e) {
        *reinterpret_cast<HmmException *>(be->_error) = e;
        return nullptr;
    }
}

void ebc_be_set_indel_parameters(EBCBandingEstimator *be, double NB_probability, double rate)
{
    be->estimate_indel_params = false;
    be->indel_NB_probability = NB_probability;
    be->indel_rate = rate;
}

void ebc_be_unset_indel_parameters(EBCBandingEstimator *be)
{
    be->estimate_indel_params = true;
    be->indel_NB_probability = EBC_BE_DEFAULTS_INDEL_NB_PROBABILITY;
    be->indel_rate = EBC_BE_DEFAULTS_INDEL_RATE;
}

void ebc_be_set_alpha(EBCBandingEstimator *be, double alpha)
{
    be->estimate_alpha = false;
    be->alpha = alpha;
}

void ebc_be_unset_alpha(EBCBandingEstimator *be)
{
    be->estimate_alpha = true;
    be->alpha = EBC_BE_DEFAULTS_ALPHA;
}

void ebc_be_set_categories(EBCBandingEstimator *be, unsigned int categories)
{
    be->estimate_categories = false;
    be->gamma_rate_categories = categories;
}

void ebc_be_unset_categories(EBCBandingEstimator *be)
{
    be->estimate_categories = true;
    be->gamma_rate_categories = EBC_BE_DEFAULTS_GAMMA_RATE_CATEGORIES;
}

bool ebc_be_set_input(EBCBandingEstimator *be, const char *fasta)
{
    if (be->_parser) {
        delete reinterpret_cast<EBC::StreamParser *>(be->_parser);
        be->_parser = nullptr;
    }

    stringstream inputStream(fasta);

    try {
        be->_parser = new StreamParser(inputStream);
    } catch (HmmException &e) {
        *reinterpret_cast<HmmException *>(be->_error) = e;
        return false;
    }

    return true;
}

bool ebc_be_set_input_from_file(EBCBandingEstimator *be, const char *file_name)
{
    if (be->_parser) {
        delete reinterpret_cast<EBC::StreamParser *>(be->_parser);
        be->_parser = nullptr;
    }

    ifstream inputStream(file_name, ios::in);

    try {
        be->_parser = new StreamParser(inputStream);
    } catch (HmmException &e) {
        *reinterpret_cast<HmmException *>(be->_error) = e;
        return false;
    }

    return true;
}

void ebc_seq_free(EBCSequences *seq)
{
    if (!seq) {
        return;
    }

    delete reinterpret_cast<EBC::Sequences *>(seq->_sequences);
    delete reinterpret_cast<EBC::ModelEstimator *>(seq->_modelEstimator);
    delete reinterpret_cast<EBC::BandingEstimator *>(seq->_bandingEstimator);
    delete seq;
}

unsigned int ebc_seq_count(EBCSequences *seq)
{
    if (!seq) {
        return 0;
    }

    return reinterpret_cast<EBC::Sequences *>(seq->_sequences)->getSequenceCount();
}

double ebc_seq_get_distance(EBCSequences *seq, unsigned int seq_id1, unsigned int seq_id2)
{
    EBC::BandingEstimator * be = reinterpret_cast<EBC::BandingEstimator *>(seq->_bandingEstimator);
    EBC::Sequences * sequences = reinterpret_cast<EBC::Sequences *>(seq->_sequences);
    unsigned int size = sequences->getSequenceCount();

    if (seq_id1 > seq_id2) {
        // Swap if they are in an inconvinient order for the calculations below.
        unsigned int tmp = seq_id1;
        seq_id1 = seq_id2;
        seq_id2 = tmp;
    } else if (seq_id1 == seq_id2) {
        return 0.0;
    }

    /*
     * The distances stored inside paHMM are the elements in an upper-triangular
     * distance matrix (excluding the diagonal). They are actually stored in a list,
     * starting from the first row of the distance matrix, from left to right.
     *
     * We therefore have to figure out which element in the list corresponds to the
     * given sequence numbers (IDs).
     *
     * The following formula gives us the right index: ((2s - 3)*i - i^2)/2 + j - 1
     * where s is the total number of sequences and (i,j) is a position on the
     * distance matrix's upper-triangular part.
     */
    return be->optimizePair(((2*size-3) * seq_id1 - seq_id1*seq_id1)/2 + seq_id2 - 1);
}

double ebc_seq_get_distance_from_names(EBCSequences *seq,
                                       const char *seq_name1, const char *seq_name2)
{
    return ebc_seq_get_distance(seq,
                reinterpret_cast<EBC::Sequences *>(seq->_sequences)
                                ->getSequenceIdFromCString(seq_name1),
                reinterpret_cast<EBC::Sequences *>(seq->_sequences)
                                ->getSequenceIdFromCString(seq_name2));
}

const char *ebc_seq_get_name(EBCSequences *seq, unsigned int seq_id)
{
    return reinterpret_cast<EBC::Sequences *>(seq->_sequences)->getSequenceName(seq_id).c_str();
}

const char *ebc_seq_get_sequence(EBCSequences *seq, unsigned int seq_id)
{
    EBC::Sequences * sequences = reinterpret_cast<EBC::Sequences *>(seq->_sequences);

    if (seq_id >= sequences->getSequenceCount()) {
        return nullptr;
    }

    return sequences->getRawSequenceAt(seq_id).c_str();
}

const char *ebc_seq_get_sequence_from_name(EBCSequences *seq, const char *seq_name)
{
    return ebc_seq_get_sequence(seq, reinterpret_cast<EBC::Sequences *>(seq->_sequences)
                                ->getSequenceIdFromCString(seq_name));
}

}

