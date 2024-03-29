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


#ifndef FORWARDPAIRHMM_HPP_
#define FORWARDPAIRHMM_HPP_

#include "hmm/EvolutionaryPairHMM.hpp"


namespace EBC
{

class ForwardPairHMM: public EBC::EvolutionaryPairHMM
{
friend class BackwardPairHMM;
protected:

	vector<double> userIndelParameters;
	vector<double> userSubstParameters;

public:
	ForwardPairHMM(vector<SequenceElement*>* s1, vector<SequenceElement*>* s2,
			SubstitutionModelBase* smdl, IndelModel* imdl,
			Definitions::DpMatrixType mt, Band* bandObj = nullptr, bool useEquilibriumProbabilities = true);

	virtual ~ForwardPairHMM();

	double runAlgorithm();

};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
