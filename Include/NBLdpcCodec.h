#pragma once

#include "GFLinAlg.h"
#include <fstream>

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
			pCurConstraints = static_cast<GFPair *>(realloc(pCurConstraints, sizeof(GFPair) * cnt));
			m_ppCheckConstraints[i] = pCurConstraints;
		}
		m_pGenMatrix = Nullspace(pCheckMatrix, Length, m_Dimension, nullptr, m_GF);
	}

	explicit NBLdpcCodec(std::string specFile)
	{
		std::ifstream ifs(specFile);
		unsigned Extension;
		ifs >> m_Length >> m_NumOfChecks >> Extension;
		m_GF.Init(Extension);
		m_Dimension = m_NumOfChecks;
		m_ppCheckConstraints = new GFPair*[m_NumOfChecks];
		GFSymbol *pCheckMatrix = new GFSymbol[m_Length * m_NumOfChecks]();
		char type;
		ifs >> type;
		if(type == 'R')
		{
			unsigned NumOfRowElements, NumOfColumnElements;
			ifs >> NumOfRowElements >> NumOfColumnElements;
			for (unsigned i = 0; i < m_NumOfChecks; ++i)
			{
				m_ppCheckConstraints[i] = new GFPair[NumOfRowElements + 1];
				for (unsigned j = 0; j < NumOfRowElements; ++j)
				{
					int tmpSym;
					ifs >> m_ppCheckConstraints[i][j].first >> tmpSym;
					pCheckMatrix[i * m_Length + m_ppCheckConstraints[i][j].first]
						= m_ppCheckConstraints[i][j].second = tmpSym;
				}
				m_ppCheckConstraints[i][NumOfRowElements].first = m_ppCheckConstraints[i][NumOfRowElements].second = 0;
			}
		}
		ifs.close();
		m_pGenMatrix = Nullspace(pCheckMatrix, m_Length, m_Dimension, nullptr, m_GF);
		for (unsigned j = 0; j < m_Dimension; ++j)
		{
			for (unsigned i = 0; i < m_NumOfChecks; ++i)
			{
				GFSymbol C = 0;
				GFPair *pCurConstraints = m_ppCheckConstraints[i];
				while (pCurConstraints->second)
				{
					C ^= m_GF.multiply(m_pGenMatrix[j * m_Length + pCurConstraints->first], pCurConstraints->second);
					++pCurConstraints;
				}
				_ASSERT(C == 0);
			}
		}
		delete[] pCheckMatrix;
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
	unsigned GetLength() const
	{
		return m_Length;
	}
	unsigned GetMaxFieldElement() const
	{
		return m_GF.FieldSize_1;
	}
};
