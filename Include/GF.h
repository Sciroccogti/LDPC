#pragma once
#include <utility>

typedef unsigned char tBit;
typedef unsigned char GFSymbol;
typedef std::pair<unsigned, GFSymbol> GFPair;

void BuildGF(unsigned Extension, GFSymbol* pGF, int* pLogTable);

struct GaloisField
{
	unsigned Extension;
	unsigned FieldSize_1;
	GFSymbol* pGF;
	int* pLogTable;

	GaloisField(unsigned extension)
	{
		Init(extension);
	}
	GaloisField() : pGF(nullptr), pLogTable(nullptr)
	{

	}
	void Init(unsigned extension)
	{
		Extension = extension;
		FieldSize_1 = (1 << extension) - 1;
		pGF = new GFSymbol[(FieldSize_1) * 2];
		pLogTable = new int[FieldSize_1 + 1];
		BuildGF(extension, pGF, pLogTable);
	}
	// a * b in GF(2^m)
	inline GFSymbol multiply(GFSymbol a, GFSymbol b) const
	{
		return (a && b) ? pGF[pLogTable[a] + pLogTable[b]] : 0;
	}
	// a * \alpha^x in GF(2^m)
	inline GFSymbol multiplyConst(GFSymbol a, int x) const
	{
		return (a) ? pGF[pLogTable[a] + x] : 0;
	}
	// return y: \alpha^y = a / \alpha^x in GF(2^m)
	inline int divideConstDeg(GFSymbol a, int x) const
	{
		if (a)
		{
			int var = pLogTable[a] - x;
			return (var > 0) ? var : FieldSize_1 + var;
		}
		else
			return -1;
	}
	// a / \alpha^x in GF(2^m)
	inline GFSymbol divideConst(GFSymbol a, int x) const
	{
		if (a)
		{
			int var = pLogTable[a] - x;
			return pGF[(var >= 0) ? var : FieldSize_1 + var];
		}
		else
			return 0;
	}
	// for nonzero a return b: a * b = 1 mod 2^m
	inline GFSymbol inverse(GFSymbol a) const
	{
		return (a != 1) ? pGF[FieldSize_1 - pLogTable[a]] : 1;
	}
	// for nonzero a return x: a * \alpha^x = 1 mod 2^m
	inline GFSymbol inverseDeg(GFSymbol a) const
	{
		return (a != 1) ? FieldSize_1 - pLogTable[a] : 0;
	}
	~GaloisField()
	{
		delete[] pGF;
		delete[] pLogTable;
	}
};