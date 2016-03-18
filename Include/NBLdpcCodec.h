#pragma once

//#include "TannerGraph.h"

#include "GFLinAlg.h"

class NBLdpcCodec
{
	unsigned m_Length;
	unsigned m_Dimension;
	unsigned m_NumOfChecks;
	GFSymbol *m_pGenMatrix;
	GFPair **m_ppCheckConstraints;
	GaloisField m_GF;

public:
	NBLdpcCodec(unsigned Length, unsigned NumOfChecks, unsigned Extension, GFSymbol *pCheckMatrix)
		: m_Length(Length), m_NumOfChecks(NumOfChecks)
	{
		m_Dimension = m_NumOfChecks;
		m_GF.Init(Extension);
		m_ppCheckConstraints = new GFPair*[m_NumOfChecks];
		GFPair *pCurConstraints;
		unsigned cnt;
		for (unsigned i = 0; i < m_NumOfChecks; ++i)
		{
			pCurConstraints = new GFPair[m_Length + 1];
			cnt = 0;
			for (unsigned j = 0; j < m_Length; j++)
			{
				if (pCheckMatrix[i * m_Length + j])
					pCurConstraints[cnt++] = std::make_pair(j, pCheckMatrix[i * m_Length + j]);
			}
			pCurConstraints[cnt++] = std::make_pair(0, 0);
			pCurConstraints = (GFPair *)realloc(pCurConstraints, sizeof(GFPair) * cnt);
			m_ppCheckConstraints[i] = pCurConstraints;
		}
		m_pGenMatrix = Nullspace(pCheckMatrix, Length, m_Dimension, nullptr, m_GF);
	}
	void Encode(GFSymbol *pInfSymbols, GFSymbol *pCodeword)
	{
		GFSymbol *pRow = m_pGenMatrix;
		memset(pCodeword, 0, sizeof(GFSymbol) * m_Length);
		for (unsigned j = 0; j < m_Dimension; ++j)
		{	
			for (unsigned i = 0; i < m_Length; ++i)
				pCodeword[i] ^= m_GF.multiply(pInfSymbols[j], pRow[i]);
			pRow += m_Length;
		}
	}
	void Decode(float *pReceivedSymbols, GFSymbol *pInfSymbols, GFSymbol *pCodeword)
	{

	}

	void Modulate()
	{

	}
	bool VerifyCodeword(GFSymbol *pCodeword)
	{
		bool f = true;
		for (unsigned i = 0; i < m_NumOfChecks; ++i)
		{
			GFSymbol C = 0;
			GFPair *pCurConstraints = m_ppCheckConstraints[i];
			while (pCurConstraints->second)
			{
				C ^= m_GF.multiply(pCodeword[pCurConstraints->first], pCurConstraints->second);
				++pCurConstraints;
			}
			f &= (C == 0);
		}
		return f;
	}
	~NBLdpcCodec()
	{
		delete[] m_pGenMatrix; 
		for (unsigned i = 0; i < m_NumOfChecks; ++i)
			delete[] m_ppCheckConstraints[i];
		delete[] m_ppCheckConstraints;
	}

	unsigned GetDimension() const
	{
		return m_Dimension;
	}
};
