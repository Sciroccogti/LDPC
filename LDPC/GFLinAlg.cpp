#include "GFLinAlg.h"

//construct a nullspace basis for a given matrix. 
GFSymbol* Nullspace(GFSymbol*pMatrix  ///the matrix to be considered
	, unsigned Length ///number of columns in the matrix
	, unsigned& Dimension ///number of rows in the input matrix. On return, the dimension of the nullspace
	, unsigned* pPermutation ///preferred column permutation
	, GaloisField &gf
	)
{
	unsigned* pPerm;
	if (pPermutation)
		pPerm = pPermutation;
	else
	{
		pPerm = new unsigned[Length];
		for (unsigned i = 0; i<Length; i++)
			pPerm[i] = i;
	}
	Gauss(pMatrix, Dimension, Length, true, pPerm, gf);
	unsigned NumOfChecks = Length - Dimension;
	GFSymbol* pNullspace = new GFSymbol[NumOfChecks * Length];
	memset(pNullspace, 0, sizeof(GFSymbol) * Length * NumOfChecks);
	for (unsigned i = 0; i < NumOfChecks; i++)
	{
		pNullspace[i * Length + pPerm[Dimension + i]] = 1;
		for (unsigned j = 0; j < Dimension; j++)
			pNullspace[i * Length + pPerm[j]] = pMatrix[j * Length + pPerm[Dimension + i]];
	}
	Dimension = Length - Dimension;
	if (!pPermutation)
		delete[]pPerm;
	return pNullspace;
}

void Gauss(GFSymbol* pMatrix
	, unsigned& NumOfRows
	, unsigned NumOfColumns
	, bool ReversePass
	, unsigned *pPermutation
	, GaloisField& gf)
{
	int gfDeg;
	for (unsigned i = 0; i < NumOfRows; i++)
	{
		//identify the leading column
		unsigned c = i;
		bool Success = false;
		for (c; c < NumOfColumns; c++)
		{
			unsigned C = pPermutation[c];
			for (unsigned j = i; j < NumOfRows; j++)
			{
				if (pMatrix[j * NumOfColumns + C] != 0)
				{
					Success = true;
					gfDeg = gf.pLogTable[pMatrix[i*NumOfColumns + C]];
					for (unsigned k = 0; k < NumOfColumns; ++k)
						pMatrix[j*NumOfColumns + k] = gf.divideConst(pMatrix[j*NumOfColumns + k], gfDeg);
					if (j > i)
					{
						for (unsigned k = 0; k < NumOfColumns; ++k)
							pMatrix[i*NumOfColumns + k] ^= pMatrix[j*NumOfColumns + k];
					}
					break;
				}
			}
			if (Success)
			{
				if (c != i)
					std::swap(pPermutation[c], pPermutation[i]);
				break;
			}
		}
		if (!Success)
		{
			NumOfRows = i;
			break;
		}
		unsigned LoopStart = (ReversePass) ? 0 : (i + 1);
		unsigned C = pPermutation[i];
		for (unsigned j = LoopStart; j < NumOfRows; ++j)
		{
			if (j == i)
				continue;
			if (pMatrix[j*NumOfColumns + C] != 0)
			{
				gfDeg = gf.pLogTable[pMatrix[j*NumOfColumns + C]];
				for (unsigned k = 0; k < NumOfColumns; ++k)
					pMatrix[j*NumOfColumns + k] ^= gf.multiplyConst(pMatrix[i*NumOfColumns + k], gfDeg);
			}
		}
	}
}
