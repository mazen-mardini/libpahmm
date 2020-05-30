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

#ifndef CPAHMM_P_H
#define CPAHMM_P_H

#include "core/Definitions.hpp"

struct EBCSequences;
struct EBCBandingEstimator;

/*
 * Creates a Sequence, executes a model, and stores the result all at once!
 */
EBCSequences *ebc_seq_create(EBCBandingEstimator *be, EBC::Definitions::ModelType model,
                             bool estimate_model_params, int model_param_count, ...);

#endif // CPAHMM_P_H
