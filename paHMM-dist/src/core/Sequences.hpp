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


#ifndef SEQUENCES_H_
#define SEQUENCES_H_

#include <vector>
#include "core/Definitions.hpp"
#include "core/IParser.hpp"
#include "core/Dictionary.hpp"
#include "core/HmmException.hpp"
#include "core/SequenceElement.hpp"
#include "core/Definitions.hpp"
#include <array>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <string_view>


using namespace std;

namespace EBC
{

class Sequences
{
private:
    IParser *iParser;
    vector<string>* rawSequences;

    vector<string>* sequenceNames;
    unordered_map<string_view, unsigned int> sequenceNamesToIds;

    vector<vector<SequenceElement*> > translatedSequences;
    vector<std::pair<unsigned int, unsigned int> > pairs;
    vector<std::pair<unsigned int, unsigned int> >::iterator pairIterator;

    unsigned int sequenceCount;
    double* observedFrequencies;

    bool removeGaps;

    Dictionary *dict;

public:

    //Input from file or console
    Sequences(IParser*, Definitions::SequenceType, bool fixedAlignment=false) /* throw (HmmException&) */;

    virtual ~Sequences();

    inline Dictionary *getDictionary()
    {
        return dict;
    }

    //Get then in a dictionary order eg. T C A G, check the definitions in constants
    //TODO - change to A T C G
    inline double* getElementFrequencies()
    {
        if(not observedFrequencies)
            calculateObservedFrequencies();
        //DEBUGV(observedFrequencies,4);
        return observedFrequencies;
    }

    inline double* getElementFrequencies(array<unsigned int, 3>& triplet);

    //void getSequencePair(vector<SequenceElement> s1, vector<SequenceElement> s2 );
    inline vector<SequenceElement*>* getSequencesAt(unsigned int pos){
            return &translatedSequences[pos];
    }

    inline unsigned int getPairCount()
    {
        unsigned int ct = static_cast<unsigned int>(translatedSequences.size());
        return (ct*(ct-1))/2;
    }

    inline unsigned int getSequenceCount()
    {
        return static_cast<unsigned int>(translatedSequences.size());
    }

    inline string& getSequenceName(unsigned int pos)
    {
        return (*sequenceNames)[pos];
    }

    inline unsigned int getSequenceId(const string& seqname)
    {
        try {
            return sequenceNamesToIds[seqname];
        } catch (out_of_range &e) {
            throw HmmException(seqname + " not found");
        }
    }

    inline unsigned int getSequenceIdFromCString(const char *seqname)
    {
        try {
            return sequenceNamesToIds[seqname];
        } catch (out_of_range &e) {
            throw HmmException(string(seqname) + " not found");
        }
    }

    inline string& getRawSequenceAt(unsigned int pos)
    {
        return (*rawSequences)[pos];
    }

    std::pair<unsigned int, unsigned int> getPairOfSequenceIndices(unsigned int idx)
    {
        return pairs[idx];
    }
private:

    void calculateObservedFrequencies();
    inline void buildDictionary(Definitions::SequenceType);

};

} /* namespace EBC */



#endif /* SEQUENCES_H_ */
