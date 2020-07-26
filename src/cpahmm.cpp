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

EBCBandingEstimator *ebc_be_create()
{
    EBCBandingEstimator *be = new EBCBandingEstimator;
    be->_parser = nullptr;
    be->_error = nullptr;
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
    if (!be || !be->_error) {
        return nullptr;
    }

    return reinterpret_cast<HmmException *>(be->_error)->what();
}

EBCSequences *ebc_be_execute_gtr_model(EBCBandingEstimator *be, double param1, double param2,
                                       double param3, double param4, double param5)
{
    try {
        if (!be) {
            return nullptr;
        }

        EBCSequences *sequences = ebc_seq_create(be, Definitions::ModelType::GTR, false, 5,
                                                 param1, param2, param3, param4, param5);
        ebc_be_unset_error(be);
        return sequences;
    } catch (HmmException& e) {
        ebc_be_set_error(be, e);
        return nullptr;
    }
}

EBCSequences *ebc_be_execute_gtr_modelv2(EBCBandingEstimator *be)
{
    try {
        if (!be) {
            return nullptr;
        }

        EBCSequences *sequences = ebc_seq_create(be, Definitions::ModelType::GTR, true, 0);
        ebc_be_unset_error(be);
        return sequences;
    } catch (HmmException& e) {
        ebc_be_set_error(be, e);
        return nullptr;
    }
}

EBCSequences *ebc_be_execute_hky85_model(EBCBandingEstimator *be, double param)
{
    try {
        if (!be) {
            return nullptr;
        }

        EBCSequences *sequences = ebc_seq_create(be, Definitions::ModelType::HKY85, false, 1, param);
        ebc_be_unset_error(be);
        return sequences;
    } catch (HmmException& e) {
        ebc_be_set_error(be, e);
        return nullptr;
    }
}

EBCSequences *ebc_be_execute_hky85_modelv2(EBCBandingEstimator *be)
{
    try {
        if (!be) {
            return nullptr;
        }

        EBCSequences *sequences = ebc_seq_create(be, Definitions::ModelType::HKY85, true, 0);
        ebc_be_unset_error(be);
        return sequences;
    } catch (HmmException& e) {
        ebc_be_set_error(be, e);
        return nullptr;
    }
}

EBCSequences *ebc_be_execute_jtt_model(EBCBandingEstimator *be)
{
    try {
        if (!be) {
            return nullptr;
        }

        EBCSequences *sequences = ebc_seq_create(be, Definitions::ModelType::JTT, true, 0);
        ebc_be_unset_error(be);
        return sequences;
    } catch (HmmException& e) {
        ebc_be_set_error(be, e);
        return nullptr;
    }
}

EBCSequences *ebc_be_execute_lg_model(EBCBandingEstimator *be)
{
    try {
        if (!be) {
            return nullptr;
        }

        EBCSequences *sequences = ebc_seq_create(be, Definitions::ModelType::LG, true, 0);
        ebc_be_unset_error(be);
        return sequences;
    } catch (HmmException& e) {
        ebc_be_set_error(be, e);
        return nullptr;
    }
}

EBCSequences *ebc_be_execute_wag_model(EBCBandingEstimator *be)
{
    try {
        if (!be) {
            return nullptr;
        }

        EBCSequences *sequences = ebc_seq_create(be, Definitions::ModelType::WAG, true, 0);
        ebc_be_unset_error(be);
        return sequences;
    } catch (HmmException& e) {
        ebc_be_set_error(be, e);
        return nullptr;
    }
}

void ebc_be_set_indel_parameters(EBCBandingEstimator *be, double NB_probability, double rate)
{
    if (!be) {
        return;
    }

    be->estimate_indel_params = false;
    be->indel_NB_probability = NB_probability;
    be->indel_rate = rate;

    ebc_be_unset_error(be);
}

void ebc_be_unset_indel_parameters(EBCBandingEstimator *be)
{
    if (!be) {
        return;
    }

    be->estimate_indel_params = true;
    be->indel_NB_probability = EBC_BE_DEFAULTS_INDEL_NB_PROBABILITY;
    be->indel_rate = EBC_BE_DEFAULTS_INDEL_RATE;

    ebc_be_unset_error(be);
}

void ebc_be_set_alpha(EBCBandingEstimator *be, double alpha)
{
    if (!be) {
        return;
    }

    be->estimate_alpha = false;
    be->alpha = alpha;

    ebc_be_unset_error(be);
}

void ebc_be_unset_alpha(EBCBandingEstimator *be)
{
    if (!be) {
        return;
    }

    be->estimate_alpha = true;
    be->alpha = EBC_BE_DEFAULTS_ALPHA;

    ebc_be_unset_error(be);
}

void ebc_be_set_categories(EBCBandingEstimator *be, unsigned int categories)
{
    if (!be) {
        return;
    }

    be->estimate_categories = false;
    be->gamma_rate_categories = categories;

    ebc_be_unset_error(be);
}

void ebc_be_unset_categories(EBCBandingEstimator *be)
{
    if (!be) {
        return;
    }

    be->estimate_categories = true;
    be->gamma_rate_categories = EBC_BE_DEFAULTS_GAMMA_RATE_CATEGORIES;

    ebc_be_unset_error(be);
}

bool ebc_be_set_input(EBCBandingEstimator *be, const char *fasta)
{
    if (!be) {
        return false;
    }

    if (be->_parser) {
        delete reinterpret_cast<EBC::StreamParser *>(be->_parser);
        be->_parser = nullptr;
    }

    stringstream inputStream(fasta);

    try {
        be->_parser = new StreamParser(inputStream);
    } catch (HmmException &e) {
        ebc_be_set_error(be, e);
        return false;
    }

    ebc_be_unset_error(be);
    return true;
}

bool ebc_be_set_input_from_file(EBCBandingEstimator *be, const char *file_name)
{
    if (!be) {
        return false;
    }

    if (be->_parser) {
        delete reinterpret_cast<EBC::StreamParser *>(be->_parser);
        be->_parser = nullptr;
    }

    ifstream inputStream(file_name, ios::in);

    try {
        be->_parser = new StreamParser(inputStream);
    } catch (HmmException &e) {
        ebc_be_set_error(be, e);
        return false;
    }

    ebc_be_unset_error(be);
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
    if (!seq) {
        return NAN;
    }

    EBC::BandingEstimator * be = reinterpret_cast<EBC::BandingEstimator *>(seq->_bandingEstimator);
    EBC::Sequences * sequences = reinterpret_cast<EBC::Sequences *>(seq->_sequences);
    unsigned int size = sequences->getSequenceCount();

    if (seq_id1 >= size) {
        ebc_seq_set_error(seq, string("Sequence with ID ") + to_string(seq_id1) + " not found.");
        return NAN;
    }

    if (seq_id2 >= size) {
        ebc_seq_set_error(seq, string("Sequence with ID ") + to_string(seq_id2) + " not found.");
        return NAN;
    }

    if (seq_id1 > seq_id2) {
        // Swap if they are in an inconvenient order for the calculations below.
        unsigned int tmp = seq_id1;
        seq_id1 = seq_id2;
        seq_id2 = tmp;
    } else if (seq_id1 == seq_id2) {
        ebc_seq_unset_error(seq);
        return 0.0;
    }

    double distance;

    try {
        /*
         * The distances stored inside paHMM are the elements in an upper-triangular
         * distance matrix (excluding the diagonal). They are actually stored in a list,
         * starting from the first row of the distance matrix, from left to right.
         *
         * We have to figure out which element in the list correspond to the distance
         * between the given sequences.
         *
         * The following formula gives us the right index: ((2s - 3)*i - i^2)/2 + j - 1
         * where s is the total number of sequences and (i,j) is a position on the
         * distance matrix's upper-triangular part.
         */
        distance = be->optimizePair(((2*size-3) * seq_id1 - seq_id1*seq_id1)/2 + seq_id2 - 1);
    }  catch (HmmException &error) {
        ebc_seq_set_error(seq, error);
        return NAN;
    }

    ebc_seq_unset_error(seq);
    return distance;
}

double ebc_seq_get_distance_from_names(EBCSequences *seq,
                                       const char *seq_name1, const char *seq_name2)
{
    if (!seq) {
        return NAN;
    }

    int seq_id1, seq_id2;

    try {
        seq_id1 = reinterpret_cast<EBC::Sequences *>(seq->_sequences)
                ->getSequenceIdFromCString(seq_name1);
        seq_id2 = reinterpret_cast<EBC::Sequences *>(seq->_sequences)
                ->getSequenceIdFromCString(seq_name2);
    } catch (HmmException &error) {
        ebc_seq_set_error(seq, error);
        return NAN;
    }

    return ebc_seq_get_distance(seq, seq_id1, seq_id2);
}

const char *ebc_seq_get_name(EBCSequences *seq, unsigned int seq_id)
{
    if (!seq) {
        return nullptr;
    }

    EBC::Sequences *sequences = reinterpret_cast<EBC::Sequences *>(seq->_sequences);

    if (seq_id >= sequences->getSequenceCount()) {
        ebc_seq_set_error(seq, string("Sequence with ID ") + to_string(seq_id) + " not found.");
        return nullptr;
    }

    ebc_seq_unset_error(seq);
    return sequences->getSequenceName(seq_id).c_str();
}

const char *ebc_seq_get_sequence(EBCSequences *seq, unsigned int seq_id)
{
    if (!seq) {
        return nullptr;
    }

    EBC::Sequences * sequences = reinterpret_cast<EBC::Sequences *>(seq->_sequences);

    if (seq_id >= sequences->getSequenceCount()) {
        ebc_seq_set_error(seq, string("Sequence with ID ") + to_string(seq_id) + " not found.");
        return nullptr;
    }

    ebc_seq_unset_error(seq);
    return sequences->getRawSequenceAt(seq_id).c_str();
}

const char *ebc_seq_get_sequence_from_name(EBCSequences *seq, const char *seq_name)
{
    if (!seq) {
        return nullptr;
    }

    int seq_id;

    try {
        seq_id = reinterpret_cast<EBC::Sequences *>(seq->_sequences)
                ->getSequenceIdFromCString(seq_name);
    } catch (HmmException &error) {
        ebc_seq_set_error(seq, error);
        return nullptr;
    }

    return ebc_seq_get_sequence(seq, seq_id);
}

}

