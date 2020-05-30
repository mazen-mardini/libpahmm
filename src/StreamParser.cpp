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


#include "StreamParser.hpp"
#include <iostream>
#include <algorithm>


namespace EBC
{

StreamParser::StreamParser(istream &in)
    : IParser()
{
	//FIXME - check the file structure to avoid errors from phylip, nexus and other non-fasta files

	sequences = new vector<string>();
	names = new vector<string>();

    if(!in)
	{
        throw HmmException(string("Could not read from stream."));
	}
	string tmp;
	string seq = "";
	bool justStarted = true;

    while (std::getline(in, tmp,'\n'))
	{
		//wait for the first sequence description
		if(justStarted)
		{
			if (!isDefinitionLine(tmp))
				continue;
			else
			{
				names->push_back(getSequenceName(tmp));
				justStarted = false;
			}
		}
		else
		{
			if (!isDefinitionLine(tmp))
			{
				seq += tmp;
			}
			else
			{
				names->push_back(getSequenceName(tmp));
				trimWsChars(seq);
				sequences->push_back(seq);
				seq = "";
			}
		}
	}

	if(seq != "")
	{
		trimWsChars(seq);
		sequences->push_back(seq);
		seq = "";
	}

	it=sequences->begin();
	itN=names->begin();

    for(size_t i = 0; i < sequences->size(); i++){
		this->mappedSeqs.insert(make_pair(names->at(i),sequences->at(i)));
	}

	sequences->clear();
	names->clear();

	unsigned int nid = 0;

	for(auto itmap : mappedSeqs){
		DEBUG("Found sequence named " << itmap.first << "\t\twith an index of " << nid );
		nid ++;
		names->push_back(itmap.first);
		sequences->push_back(itmap.second);
	}

}

bool StreamParser::isDefinitionLine(string& s)
{
	std::size_t found = s.find(">");
	return (found !=std::string::npos);
}

string StreamParser::getSequenceName(string& s)
{
	//string name(s);
	s.erase(std::remove_if( s.begin(), s.end(), [](char c){ return (c =='>' || c =='\r' || c =='\t' || c == ' ' || c == '\n');}), s.end() );
	return s;//name;
}

void StreamParser::trimWsChars(string& s)
{
	//Trimming stars (stop codon)
	s.erase(std::remove_if( s.begin(), s.end(), [](char c){ return (c =='\r' || c =='\t' || c == ' ' || c == '\n' || c == '*');}), s.end() );
}

string StreamParser::getNextSequence()
{
	//DEBUG("Get next sequence it: " << *it);
	return *it++;
}

string StreamParser::getNextName()
{
	//DEBUG("Get next sequence it: " << *it);
	return *itN++;
}

unsigned int StreamParser::getSequenceCount()
{
    return static_cast<unsigned int>(sequences->size());
}

string StreamParser::getSequenceAt(unsigned int position)
{
	return sequences->at(position);
}

string StreamParser::getSequenceNameAt(unsigned int position)
{
	return names->at(position);
}

StreamParser::~StreamParser()
{
    delete sequences;
    delete names;
}

} /* namespace EBC */


