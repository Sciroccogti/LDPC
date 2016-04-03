#include "NBLdpcBuilder.h"
#include "GF.h"
#include <random>
#include <fstream>
#include <iostream>
#include <ctime>
#include <algorithm>
#include <vector>
#include <numeric>
#include <set>

using namespace std;

void NBLdpcBuilder::BuildRegular(unsigned Length, unsigned Redundancy, unsigned Extension, unsigned NumOfRowElements, unsigned NumOfColumnElements, std::string specFileName)
{
	if(NumOfRowElements * Redundancy != NumOfColumnElements * Length)
	{
		cerr << "Inconsistent number of row and column nonzero elements\n";
		return;
	}
	unsigned NumOfRowSubs = Length / NumOfRowElements;
	unsigned NumOfColSubs = Redundancy / NumOfColumnElements;
	if(NumOfRowSubs != NumOfColSubs)
	{
		cerr << "Block check matrices cannot be square\n";
		return;		
	}
	GFSymbol *pCheckMatrix = new GFSymbol[Length * Redundancy]();
	vector<unsigned> pPermutation(NumOfRowSubs);
	mt19937 gen(time(nullptr));
	uniform_int_distribution<int> distr(1, (1 << Extension) - 1);
	
	/*iota(pPermutation.begin(), pPermutation.end(), 0);
	vector<vector<unsigned>> pColumnPositions(Length);
	for (unsigned i = 0; i < Redundancy; ++i)
	{
		shuffle(pPermutation.begin(), pPermutation.end(), gen);
		for (unsigned j = 0; j < NumOfRowElements; ++j)
		{
			unsigned ind = pPermutation[j];
			pCheckMatrix[i * Length + ind] = distr(gen);
			pColumnPositions[ind].push_back(i);
		}
	}
	for (unsigned j = 0; j < Length; ++j)
	{
		if(pColumnPositions[j].size() < NumOfColumnElements)
		{
			for (unsigned k = 0; k < Length && pColumnPositions[j].size() < NumOfColumnElements; ++k)
			{
				if (k == j) continue;
				if(pColumnPositions[k].size() > NumOfColumnElements)
				{
					unsigned copyLen = min(NumOfColumnElements - pColumnPositions[j].size(), pColumnPositions[k].size() - NumOfColumnElements);
					unsigned cnt = 0;
					for (auto && e = pColumnPositions[k].begin(); e != pColumnPositions[k].end() && cnt < copyLen;)
					{
						if (pCheckMatrix[(*e) * Length + j] == 0)
						{
							swap(pCheckMatrix[(*e) * Length + j], pCheckMatrix[(*e) * Length + k]);
							pColumnPositions[j].push_back(*e);
							e = pColumnPositions[k].erase(e);
							cnt++;
						}
						else
							++e;
					}					
				}
			}
		}
	}*/

	for (unsigned i = 0; i < Redundancy; i += NumOfColSubs)
	{
		for (unsigned j = 0; j < Length; j += NumOfRowSubs)
		{
			iota(pPermutation.begin(), pPermutation.end(), j);
			shuffle(pPermutation.begin(), pPermutation.end(), gen);
			for (unsigned t = i; t < i + NumOfColSubs; ++t)
				pCheckMatrix[t * Length + pPermutation[t - i]] = distr(gen);
		}
	}
	ofstream ofs(specFileName);
	ofs << Length << ' ' << Redundancy << ' ' << Extension <<  endl;
	ofs << "R" << ' ' << NumOfRowElements << ' ' << NumOfColumnElements << endl;
	for (unsigned i = 0; i < Redundancy; ++i)
	{
		for (unsigned j = 0; j < Length; ++j)
			if(pCheckMatrix[i * Length + j])
				ofs << j << ' ' << (int)pCheckMatrix[i * Length + j] << ' ';
		ofs << endl;
	}
	ofs.close();
	delete[] pCheckMatrix;
}
