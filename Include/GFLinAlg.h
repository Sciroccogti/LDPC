#pragma once

#include "GF.h"

GFSymbol* Nullspace(GFSymbol*pMatrix, unsigned Length, unsigned& Dimension, unsigned* pPermutation, GaloisField &gf);
void Gauss(GFSymbol* pMatrix, unsigned& NumOfRows, unsigned NumOfColumns, bool ReversePass, unsigned *pPermutation, GaloisField& gf);
