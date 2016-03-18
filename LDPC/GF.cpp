#include "GF.h"

const GFSymbol GFGenerators[] = {
	0, 0, 7, 0xB, 0x13, 0x25, 0x43, 0x83, 0x11D, 0x211, 02011, 04005, 010123, 020033, 21993, 39065, 85245 };

void BuildGF(unsigned Extension, GFSymbol* pGF, int* pLogTable)
{
	GFSymbol GFGenPoly = GFGenerators[Extension];
	unsigned FieldSize = 1 << Extension;

	pGF[0] = 1; pLogTable[0] = -1;
	pLogTable[1] = 0;
	unsigned mask = 1 << Extension;
	for (unsigned i = 1; i < FieldSize - 1; ++i)
	{
		pGF[i] = pGF[i - 1] << 1;
		if (pGF[i] & mask)
			pGF[i] ^= GFGenPoly;
		pLogTable[pGF[i]] = i;
	}
	for (unsigned i = FieldSize - 1; i < 2 * (FieldSize - 1); ++i)
		pGF[i] = pGF[i - FieldSize + 1];
}
