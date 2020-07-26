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

#include "cpahmm_p.h"
#include "cpahmm.h"

#include <sstream>
#include <cstdarg>
#include "core/BandingEstimator.hpp"
#include "core/Sequences.hpp"
#include "heuristics/ModelEstimator.hpp"
#include "StreamParser.hpp"

using namespace std;
using namespace EBC;

EBCSequences *ebc_seq_create(EBCBandingEstimator *be, Definitions::ModelType model,
                             bool estimate_model_params, int model_param_count, ...)
{
    if (!be) {
        return nullptr;
    }

    va_list args;
    va_start(args, model_param_count);

    EBCSequences *seq = new EBCSequences;

    StreamParser *parser;
    if (be->_parser) {
        parser = reinterpret_cast<StreamParser *>(be->_parser);
    } else {
        stringstream emptyStream;
        parser = new StreamParser(emptyStream);
    }

    Definitions::SequenceType sequenceType;
    if (model == Definitions::ModelType::GTR || model == Definitions::ModelType::HKY85) {
        sequenceType = Definitions::SequenceType::Nucleotide;
    } else {
        sequenceType = Definitions::SequenceType::Aminoacid;
    }
    seq->sequenceType = sequenceType;

    Sequences* inputSeqs = new Sequences(parser, sequenceType, true);
    seq->_sequences = inputSeqs;

    ModelEstimator *tme = new ModelEstimator(inputSeqs, model, Definitions::OptimizationType::BFGS,
                                             be->gamma_rate_categories, be->alpha, be->estimate_alpha);
    seq->_modelEstimator = tme;

    if(be->estimate_alpha){
        be->alpha = tme->getAlpha();
    }

    vector<double> indelParams;
    vector<double> substParams;

    if (be->estimate_indel_params) {
        indelParams = tme->getIndelParameters();
    } else {
        indelParams = {be->indel_NB_probability, be->indel_rate};
    }

    if (estimate_model_params) {
        substParams = tme->getSubstitutionParameters();
    } else {
        if (model == Definitions::ModelType::HKY85) {
            substParams = {va_arg(args, double)};
        } else if (model == Definitions::ModelType::GTR) {
            substParams = {va_arg(args, double),
                           va_arg(args, double),
                           va_arg(args, double),
                           va_arg(args, double),
                           va_arg(args, double)};
        }
    }

    va_end(args);

    BandingEstimator* bandingEstimator =
            new BandingEstimator(Definitions::AlgorithmType::Forward, inputSeqs, model,
                                 indelParams, substParams, Definitions::OptimizationType::BFGS,
                                 be->gamma_rate_categories, be->alpha, tme->getGuideTree());
    seq->_bandingEstimator = bandingEstimator;
    seq->_ebcBandingEstimator = be;

    return seq;
}

void ebc_be_set_error(EBCBandingEstimator *be, const string &message)
{
    if (!be) {
        return;
    }

    HmmException *& error = reinterpret_cast<HmmException *&>(be->_error);

    if (!error) {
        error = new HmmException(message);
    } else {
        *error = HmmException(message);
    }
}

void ebc_be_set_error(EBCBandingEstimator *be, const HmmException &exception)
{
    if (!be) {
        return;
    }

    HmmException *& error = reinterpret_cast<HmmException *&>(be->_error);

    if (!error) {
        error = new HmmException(exception);
    } else {
        *error = exception;
    }
}

void ebc_be_unset_error(EBCBandingEstimator *be)
{
    if (!be) {
        return;
    }

    HmmException *& error = reinterpret_cast<HmmException *&>(be->_error);

    if (error) {
        delete error;
        error = nullptr;
    }
}

void ebc_seq_set_error(EBCSequences *seq, const string &message)
{
    if (seq) {
        return ebc_be_set_error(seq->_ebcBandingEstimator, message);
    }
}

void ebc_seq_set_error(EBCSequences *seq, const HmmException &exception)
{
    if (seq) {
        return ebc_be_set_error(seq->_ebcBandingEstimator, exception);
    }
}

void ebc_seq_unset_error(EBCSequences *seq)
{
    if (seq) {
        return ebc_be_unset_error(seq->_ebcBandingEstimator);
    }
}
