#include "NBLdpcDecoder.h"
#include <algorithm>
#include <functional>

NBLdpcDecoder::NBLdpcDecoder(std::string specFile, unsigned NumOfRemainingLLRs)
	: NBLdpcCodec(specFile), m_NumOfComponents(NumOfRemainingLLRs)
{
	m_ppVar2Checks = new tIndPair*[m_Length];
	m_pCheckDegrees = new unsigned[m_NumOfChecks];
	m_pVarDegrees = new unsigned[m_Length];
	if (m_Type == L_REGULAR)
	{
		regular_params *pParams = reinterpret_cast<regular_params*>(m_pAdditionalParams);
		std::fill(m_pCheckDegrees, m_pCheckDegrees + m_NumOfChecks, pParams->m_NumOfRowElements);
		std::fill(m_pVarDegrees, m_pVarDegrees + m_Length, pParams->m_NumOfColumnElements);

	}
	m_ppVar2Checks = new tIndPair*[m_Length];
	for (unsigned i = 0; i < m_Length; ++i)
		m_ppVar2Checks[i] = new tIndPair[m_pVarDegrees[i]];
	m_ppChecks2Var = new tIndPair*[m_NumOfChecks];
	for (unsigned i = 0; i < m_NumOfChecks; ++i)
		m_ppChecks2Var[i] = new tIndPair[m_pCheckDegrees[i]];
	unsigned *pCnt = new unsigned[m_Length]();
	for (unsigned i = 0; i < m_NumOfChecks; ++i)
	{
		for (unsigned j = 0; j < m_pCheckDegrees[i]; ++j)
		{
			auto a = m_ppCheckConstraints[i][j];
			m_ppChecks2Var[i][j] = std::make_pair(a.first, pCnt[a.first]);
			m_ppVar2Checks[a.first][pCnt[a.first]++] = std::make_pair(i, j);
		}
	}
	delete[] pCnt;

	m_ppCheckInputs = new tLLRPair **[m_NumOfChecks];
	m_ppCheckOutputs = new tLLRPair **[m_NumOfChecks];
	m_ppVarInputs = new tLLRPair **[m_Length];
	m_ppVarOutputs = new tLLRPair **[m_Length];
	for (unsigned i = 0; i < m_Length; ++i)
	{
		m_ppVarInputs[i] = new tLLRPair*[m_pVarDegrees[i]];
		m_ppVarOutputs[i] = new tLLRPair*[m_pVarDegrees[i]];
		for(unsigned j = 0; j < m_pVarDegrees[i]; ++j)
		{
			m_ppVarInputs[i][j] = new tLLRPair[m_NumOfComponents];
			m_ppVarOutputs[i][j] = new tLLRPair[m_NumOfComponents];
		}
	}
	for (unsigned i = 0; i < m_NumOfChecks; ++i)
	{
		m_ppCheckInputs[i] = new tLLRPair*[m_pCheckDegrees[i]];
		m_ppCheckOutputs[i] = new tLLRPair*[m_pCheckDegrees[i]];
		for (unsigned j = 0; j < m_pCheckDegrees[i]; ++j)
		{
			m_ppCheckInputs[i][j] = new tLLRPair[m_NumOfComponents];
			m_ppCheckOutputs[i][j] = new tLLRPair[m_NumOfComponents];
		}
	}
	m_ppInputLLRs = new tLLRPair*[m_Length];
	for (unsigned i = 0; i < m_Length; ++i)
		m_ppInputLLRs[i] = new tLLRPair[m_GF.FieldSize_1 + 1];
}

//compute log(1+exp(-x)), x>=0;
inline double JacobiLog(double x)
{
	if (x<5.5)
		return .131478060787483 + (3.43245248951091 - 1.74351875423031*x) / (6.11620503468962 + (2.21612896897222 + x)*x);
	else
		return 0;
}

inline double LogSum(const double& LogA, const double& LogB)
{
	return (LogA > LogB) ? LogA + JacobiLog(LogA - LogB) : LogB + JacobiLog(LogB - LogA);
}

void BoxPlus(tLLRPair* L1, GFSymbol a1, tLLRPair* L2, GFSymbol a2, GaloisField& gf, tLLRPair* Res)
{
	GFSymbol invA1 = gf.inverseDeg(a1), invA2 = gf.inverseDeg(a2);
	double tmp = LogSum(0, L1[1] + L2[gf.multiplyConst(a1, gf.inverseDeg(a2))]);
	for (unsigned i = 2; i <= gf.FieldSize_1; ++i)
		tmp = LogSum(tmp, L1[i] + L2[gf.multiplyConst(gf.multiplyConst(a1, invA2), gf.pLogTable[i])]);

	Res[0] = 0;
	for(unsigned i = 1; i <= gf.FieldSize_1; ++i)
	{
		Res[i] = LogSum(L2[gf.multiplyConst(i, invA2)], L1[1] + L2[gf.multiplyConst(i ^ a1, invA2)]);
		for (unsigned j = 2; j <= gf.FieldSize_1; ++j)
			Res[i] = LogSum(Res[i], L1[j] + L2[gf.multiplyConst(i ^ gf.multiplyConst(a1, gf.pLogTable[i]), invA2)]);
		Res[i] -= tmp;
	}
}

bool NBLdpcDecoder::Decode(float* pNoisyData, GFSymbol* pCodeword, unsigned NumOfIterations, float m_NoiseVariance)
{
	for (unsigned i = 0; i < m_Length; ++i)
	{
		for (unsigned j = 0; j <= m_GF.FieldSize_1; ++j)
			m_ppInputLLRs[i][j] = (m_pConstellation[0] - m_pConstellation[j])*(m_pConstellation[j] + m_pConstellation[0] - 2 * pNoisyData[i]) / (2 * m_NoiseVariance);//std::make_pair((m_pConstellation[0] - m_pConstellation[j])*(m_pConstellation[i] + m_pConstellation[0] - 2 * pNoisyData[i]) / (2 * m_NoiseVariance), j);
		//std::sort(m_ppInputLLRs[i], m_ppInputLLRs[i] + m_GF.FieldSize_1 + 1, std::greater<tLLRPair>());
		for (unsigned j = 0; j < m_pVarDegrees[i]; ++j)
			memcpy(m_ppVarInputs[i][j], m_ppInputLLRs[i], sizeof(tLLRPair) * m_NumOfComponents);		
	}
	for (unsigned i = 0; i < m_NumOfChecks; ++i)
		for (unsigned j = 0; j < m_pCheckDegrees[i]; ++j)
			memset(m_ppCheckOutputs[i][j], 0, sizeof(tLLRPair) * m_NumOfComponents);

	for(unsigned i = 0; i < NumOfIterations; ++i)
	{
		// 1. Variable node update
		for(unsigned j = 0; j < m_Length; ++j)
		{
			memcpy(m_ppVarOutputs[j][0], m_ppInputLLRs[j], sizeof(tLLRPair) * m_NumOfComponents);
			for (unsigned k = 0; k < m_pVarDegrees[j]; ++k)
			{
				auto&& ind = m_ppVar2Checks[j][k];
				for (unsigned t = 0; t <= m_GF.FieldSize_1; ++t)
					m_ppVarOutputs[j][0][t] += m_ppCheckOutputs[ind.first][ind.second][t];
			}
			pCodeword[j] = std::distance(m_ppVarOutputs[j][0], std::max_element(m_ppVarOutputs[j][0], m_ppVarOutputs[j][0] + m_GF.FieldSize_1 + 1));
		}
		if (VerifyCodeword(pCodeword))
			return true;

		for (unsigned j = 0; j < m_Length; ++j)
		{
			for (unsigned k = 0; k < m_pVarDegrees[j]; ++k)
			{
				memcpy(m_ppVarInputs[j][k], m_ppVarOutputs[j][0], sizeof(tLLRPair) * m_NumOfComponents);
				auto&& ind = m_ppVar2Checks[j][k];
				for (unsigned t = 0; t <= m_GF.FieldSize_1; ++t)
					m_ppVarInputs[j][k][t] -= m_ppCheckOutputs[ind.first][ind.second][t];
			}
		}
		// 2. Permutation step
		/*for (unsigned j = 0; j < m_NumOfChecks; ++j)
		{
			for (unsigned k = 0; k < m_pCheckDegrees[j]; ++k)
			{
				auto ind = m_ppChecks2Var[j][k];
				GFSymbol curConstraint = m_GF.pLogTable[m_ppCheckConstraints[j][k].second];
				for (unsigned t = 0; t < m_NumOfComponents; ++t)
					m_ppCheckInputs[j][k][t].second = m_GF.multiplyConst(m_ppVarOutputs[ind.first][ind.second][t].second, curConstraint);
			}
		}*/

		// 3. Check node update

		for(unsigned j = 0; j < m_NumOfChecks; ++j)
		{
			//Calculate \sigma and \rho distributions (see H. Wymeersch, H. Steendam and M. Moeneclaey "Log-domain decoding of LDPC codes over GF(q)")
			for(unsigned k = 0; k <= m_GF.FieldSize_1; ++k)
			{
				auto&& ind0 = m_ppChecks2Var[j][0];
				auto&& indL = m_ppChecks2Var[j][m_pCheckDegrees[j] - 1];
				m_ppCheckInputs[j][0][m_GF.multiplyConst(k, m_ppCheckConstraints[j][0].second)] = m_ppVarInputs[ind0.first][ind0.second][k];
				m_ppCheckInputs[j][m_pCheckDegrees[j] - 1][m_GF.multiplyConst(k, m_ppCheckConstraints[j][m_pCheckDegrees[j] - 1].second)] = m_ppVarInputs[indL.first][indL.second][k];
			}
			for(unsigned k = 1; k < m_pCheckDegrees[j]; ++k)
			{
				auto&& ind = m_ppChecks2Var[j][k];
				BoxPlus(m_ppCheckInputs[j][k - 1], 1, m_ppVarInputs[ind.first][ind.second], m_GF.pGF[m_ppCheckConstraints[j][k].second], m_GF, m_ppCheckInputs[j][k]);
				BoxPlus(m_ppCheckOutputs[j][m_pCheckDegrees[j] - k - 1], 1, m_ppVarInputs[ind.first][ind.second], m_ppCheckConstraints[j][k].second, m_GF, m_ppCheckOutputs[j][m_pCheckDegrees[j] - k]);
			}

			for (unsigned k = 0; k <= m_GF.FieldSize_1; ++k)
			{
				m_ppCheckOutputs[j][0][m_GF.multiplyConst(k, m_GF.inverseDeg(m_GF.pGF[m_ppCheckConstraints[j][0].second]))] 
					= m_ppCheckOutputs[j][1][k];
				m_ppCheckOutputs[j][m_pCheckDegrees[j] - 1][m_GF.multiplyConst(k, m_GF.inverseDeg(m_GF.pGF[m_ppCheckConstraints[j][m_pCheckDegrees[j] - 1].second]))]
					= m_ppCheckInputs[j][m_pCheckDegrees[j] - 2][k];
			}
			for(unsigned k = 1; k < m_pCheckDegrees[j] - 1; ++k)
			{
				GFSymbol tmpSym = m_GF.inverse(m_GF.pGF[m_ppCheckConstraints[j][k].second]);
				BoxPlus(m_ppCheckOutputs[j][k + 1], tmpSym, m_ppCheckInputs[j][k - 1], tmpSym, m_GF, m_ppCheckOutputs[j][k]);
			}
		}

		// 4. Reverse permutation step
		/*for (unsigned j = 0; j < m_Length; ++j)
		{
			for (unsigned k = 0; k < m_pVarDegrees[j]; ++k)
			{
				auto ind = m_ppVar2Checks[j][k];
				GFSymbol curConstraint = m_GF.inverseDeg(m_ppCheckConstraints[ind.first][ind.second].second);
				for (unsigned t = 0; t < m_NumOfComponents; ++t)
					m_ppVarInputs[j][k][t].second = m_GF.multiplyConst(m_ppCheckOutputs[ind.first][ind.second][t].second, curConstraint);
			}
		}*/
	}
	return false;
}

NBLdpcDecoder::~NBLdpcDecoder()
{
	for (unsigned i = 0; i < m_Length; ++i)
		delete[] m_ppVar2Checks[i];
	delete[] m_ppVar2Checks;
}
