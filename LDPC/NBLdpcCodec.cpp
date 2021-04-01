#include "NBLdpcCodec.h"
#include <iostream>
#include <cassert>

NBLdpcCodec::NBLdpcCodec(unsigned Length, unsigned NumOfChecks, unsigned Extension, GFSymbol* pCheckMatrix): 
	m_Length(Length), m_NumOfChecks(NumOfChecks), m_Type(L_OTHER)
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
	InitModem();
}

NBLdpcCodec::NBLdpcCodec(std::string specFile): m_Type(L_OTHER)
{
	std::ifstream ifs(specFile);
	if(!ifs)
	{
		std::cerr << "No such file!\n";
		throw std::exception();
	}
	unsigned Extension;
	ifs >> m_Length >> m_NumOfChecks >> Extension;
	m_GF.Init(Extension);
	m_Dimension = m_NumOfChecks;
	m_ppCheckConstraints = new GFPair*[m_NumOfChecks];
	GFSymbol *m_pCheckMatrix = new GFSymbol[m_Length * m_NumOfChecks]();
	char type;
	ifs >> type;
	if(type == 'R')
	{
		m_Type = L_REGULAR;
		m_pAdditionalParams = new regular_params();
		regular_params *pParams = reinterpret_cast<regular_params*>(m_pAdditionalParams);
		unsigned &NumOfRowElements = pParams->m_NumOfRowElements;
		unsigned &NumOfColumnElements = pParams->m_NumOfColumnElements;
		ifs >> NumOfRowElements >> NumOfColumnElements;
		for (unsigned i = 0; i < m_NumOfChecks; ++i)
		{
			m_ppCheckConstraints[i] = new GFPair[NumOfRowElements + 1];
			for (unsigned j = 0; j < NumOfRowElements; ++j)
			{
				int tmpSym;
				ifs >> m_ppCheckConstraints[i][j].first >> tmpSym;
				m_pCheckMatrix[i * m_Length + m_ppCheckConstraints[i][j].first]
					= m_ppCheckConstraints[i][j].second = tmpSym;
			}
			m_ppCheckConstraints[i][NumOfRowElements].first = m_ppCheckConstraints[i][NumOfRowElements].second = 0;
		}
	}
	ifs.close();
	m_pGenMatrix = Nullspace(m_pCheckMatrix, m_Length, m_Dimension, nullptr, m_GF);
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
			assert(C == 0);
		}
	}
	InitModem();
}

bool NBLdpcCodec::VerifyCodeword(GFSymbol* pCodeword)
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
