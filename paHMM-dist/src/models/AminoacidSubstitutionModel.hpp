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


#ifndef AA_MODEL_H_
#define AA_MODEL_H_

#include "core/Dictionary.hpp"
#include "core/Definitions.hpp"
#include "models/SubstitutionModelBase.hpp"
#include "core/Maths.hpp"
#include <cmath>

namespace EBC
{

class AminoacidSubstitutionModel : public EBC::SubstitutionModelBase
{
protected:

	bool eigenDecomposed;

	double maxRate;

public:

	AminoacidSubstitutionModel(Dictionary*, Maths*, unsigned int, Definitions::aaModelDefinition);

	virtual ~AminoacidSubstitutionModel();

	void calculateModel();

	void summarize();

	void setParameters(const vector<double>&){}

	//void setObservedFrequencies(double* observedFrequencies) {};
};

} /* namespace EBC */
#endif /* MODEL_H_ */
