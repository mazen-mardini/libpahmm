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


#include "core/Dictionary.hpp"
#include "core/Definitions.hpp"
#include "core/HmmException.hpp"
#include <algorithm>

namespace EBC
{

void Dictionary::setAlphabet(char dict[], unsigned short size)
{
	//this->alphabet.reserve(size+1);

	//includes gap
	alphabet.append(dict,size+1);

    for(unsigned char i=0; i<=size; i++)
	{
		//this->alphabet.push_back(string(1,dict[i]));
		unsigned char* idxptr = new unsigned char[1];
		idxptr[0] = i;
		SequenceElement* sel  = new SequenceElement(i==gapId, i, idxptr, alphabet[i]);
        translator[static_cast<size_t>(alphabet[i])] = sel;
        translator[static_cast<size_t>(tolower(alphabet[i]))] = sel;
	}

	//alphabet size does not include gap e.g. size is 4 for nucleotides
	this->alphabetSize = size;
}

void Dictionary::addFastaClasses(const map<char, vector<char>>& classmap)
{
    unsigned char currId = static_cast<unsigned char>(this->alphabetSize + 1);

	for(auto const &fcls : classmap){
		unsigned char* ids = new unsigned char[fcls.second.size()];
		for(unsigned int i = 0; i < fcls.second.size(); i++){
            ids[i] = (translator[static_cast<size_t>(fcls.second[i])])->getMatrixIndex();
		}

        SequenceElement* sel = new SequenceElement(false, currId, ids, fcls.first,
                                                   static_cast<unsigned short>(fcls.second.size()));
        translator[static_cast<size_t>(fcls.first)] = sel;
        translator[static_cast<size_t>(tolower(fcls.first))] = sel;
        alphabet.append(static_cast<size_t>(fcls.first), static_cast<char>(currId));
		currId++;
	}

}

void Dictionary::outputAlphabet()
{
	cout << "Model dictionary: " << endl;

	cout << alphabet << endl;
}

char Dictionary::getSymbolAt(unsigned char index)
{
    return alphabet[index];
}


Dictionary::~Dictionary()
{}

unsigned char Dictionary::getSymbolIndex(char symbol)
{
    // Yes, if symbol doesn't exist everything will crash and burn. ;)
    return (translator[static_cast<size_t>(symbol)])->getMatrixIndex();
}

inline SequenceElement* Dictionary::getSequenceElement(char symbol)
{
    SequenceElement* se = nullptr;

    if (static_cast<size_t>(symbol) < translator.size()) {
        se = translator.at(static_cast<size_t>(symbol));
    }

    if (!se) {
        string message = "Symbol not found in the dictionary: ";
        message += symbol;
        throw HmmException(message);
    }

    return se;
}

vector<SequenceElement*> Dictionary::translate(const string& sequence, bool removeGaps)
{

    vector<SequenceElement*> translatedVector;
    translatedVector.reserve(sequence.size());

    for(auto it = sequence.cbegin(); it < sequence.cend(); it++)
    {
		auto el = getSequenceElement(*it);

        if(el->isIsGap() && removeGaps) {
			continue;
        }

        translatedVector.push_back(el);
	}

	return translatedVector;
}

const char Dictionary::nucleotides[] = {'T', 'C', 'A', 'G', '-'};
//const char Dictionary::nucleotides[] = {'A', 'C', 'G', 'T'};
const char Dictionary::aminoacids[] = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','-'};

const char Dictionary::gapChar = '-';

const map<char,vector<char> > Dictionary::nucFastaClasses = {
		{'R',{'A','G'}},
		{'Y',{'C','T'}},
		{'K',{'G','T'}},
		{'M',{'A','C'}},
		{'S',{'C','G'}},
		{'W',{'A','T'}},
		{'B',{'C','G','T'}},
		{'D',{'A','G','T'}},
		{'H',{'A','C','T'}},
		{'V',{'A','C','G'}},
		{'N',{'A','C','G','T'}}

};

const map<char,vector<char> > Dictionary::aaFastaClasses = {
		{'B',{'D','N'}},
		{'J',{'L','I'}},
		{'Z',{'E','Q'}},
		{'X',{'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'}}
};

unsigned short Dictionary::getAlphabetSize()
{
	return this->alphabetSize;
}

NucleotideDictionary::NucleotideDictionary()
{
	gapId = 4;
    setAlphabet(const_cast<char*>(Dictionary::nucleotides), 4);
    handleTUequivalence();
    addFastaClasses(Dictionary::nucFastaClasses);

}

void NucleotideDictionary::handleTUequivalence()
{
    translator['U'] = getSequenceElement('T');
    translator['u'] = getSequenceElement('t');
}

NucleotideDictionary *NucleotideDictionary::instance = nullptr;

NucleotideDictionary *NucleotideDictionary::getDictionary()
{
    if (not NucleotideDictionary::instance)
        NucleotideDictionary::instance = new NucleotideDictionary;

    return NucleotideDictionary::instance;
}

AminoacidDictionary::AminoacidDictionary()
{
	gapId = 20;
    setAlphabet(const_cast<char*>(Dictionary::aminoacids), 20);
    addFastaClasses(Dictionary::aaFastaClasses);
}

AminoacidDictionary *AminoacidDictionary::instance = nullptr;

AminoacidDictionary *AminoacidDictionary::getDictionary()
{
    if (not AminoacidDictionary::instance)
        AminoacidDictionary::instance = new AminoacidDictionary;

    return AminoacidDictionary::instance;
}



}

