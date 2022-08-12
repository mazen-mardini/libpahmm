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


#ifndef DICTIONARY_H_
#define DICTIONARY_H_

#include <vector>
#include <string>
#include <array>
#include <map>
#include <cmath>
#include "core/SequenceElement.hpp"
#include "core/Definitions.hpp"

using namespace std;

namespace EBC
{
	class Dictionary
	{
	protected:
		unsigned short alphabetSize;
		unsigned char gapId;
		string alphabet;
        array<SequenceElement*, 256> translator = {nullptr};

	public:
        virtual ~Dictionary();

		static const char nucleotides[5];
		static const char aminoacids[21];
		static const char gapChar;
        static const map<char, vector<char>> nucFastaClasses;
        static const map<char, vector<char>> aaFastaClasses;

		//virtual short getSymbolIndex(string &symbol);
		virtual unsigned char getSymbolIndex(char symbol);
        virtual vector<SequenceElement*> translate(const string &sequence, bool removeGaps = false);
		virtual unsigned short getAlphabetSize();

        inline SequenceElement* getSequenceElement(char symbol);

		virtual char getSymbolAt(unsigned char i);

		virtual void outputAlphabet();

		inline unsigned char getGapID()
		{
			return gapId;
		}

        static constexpr size_t maxAlphabetSize = 20;

    protected:
        virtual void setAlphabet(char alphabet[], unsigned short size);
        void addFastaClasses(const map<char, vector<char>>& classmap);

	};

    /*
     * NucleotideDictionary and AminoacidDictionary defined as singletons to
     * prevent the same dictionary from being rebuilt and save a few CPU cycles.
     */

	class NucleotideDictionary : public Dictionary
	{
        NucleotideDictionary();
        void handleTUequivalence();

        static NucleotideDictionary *instance;

	public:
        static NucleotideDictionary *getDictionary();
	};

	class AminoacidDictionary : public Dictionary
	{
        AminoacidDictionary();

        static AminoacidDictionary *instance;

	public:
        static AminoacidDictionary *getDictionary();
	};
}


#endif /* DICTIONARY_H_ */
