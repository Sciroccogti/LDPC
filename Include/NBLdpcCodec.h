#pragma once

#include "GFLinAlg.h"
#include <fstream>

struct regular_params
{
	unsigned m_NumOfRowElements;
	unsigned m_NumOfColumnElements;
};

class NBLdpcCodec
{
protected:
	enum ldpc_type_t {L_REGULAR, L_OTHER};
	
	GaloisField m_GF;
	//Code params
	ldpc_type_t m_Type;
	unsigned m_Length;
	unsigned m_Dimension;
	unsigned m_NumOfChecks;
	GFSymbol *m_pGenMatrix;
	GFPair **m_ppCheckConstraints;

	void *m_pAdditionalParams;

	//For modulation
	float *m_pConstellation;
	float m_ScaleFactor;

	void InitModem()
	{
		m_ScaleFactor = sqrtf(3.f / ((1u << (2 * m_GF.Extension)) - 1));
		m_pConstellation = new float[m_GF.FieldSize_1 + 1];
		for (unsigned i = 0; i <= m_GF.FieldSize_1; ++i)
			m_pConstellation[i] = (2.f * i - m_GF.FieldSize_1) * m_ScaleFactor;
	}

public:
	NBLdpcCodec(unsigned Length, unsigned NumOfChecks, unsigned Extension, GFSymbol* pCheckMatrix);

	explicit NBLdpcCodec(std::string specFile);

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

	float Modulate(const GFSymbol& Symbol) const
	{
		return m_pConstellation[Symbol];
	}

	bool VerifyCodeword(GFSymbol* pCodeword);

	virtual ~NBLdpcCodec()
	{
		delete[] m_pGenMatrix;
		for (unsigned i = 0; i < m_NumOfChecks; ++i)
			delete[] m_ppCheckConstraints[i];
		delete[] m_ppCheckConstraints;
		delete[] m_pConstellation;
		delete[] m_pAdditionalParams;
	}

	unsigned GetDimension() const
	{
		return m_Dimension;
	}
	unsigned GetFieldOrder() const
	{
		return m_GF.Extension;
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
