#include <iostream>
#include "NBLdpcCodec.h"

using namespace std;

int main()
{
	const unsigned Length = 3, NumOfChecks = 2;
	GFSymbol pCheckMatrix[] = {
		1, 1, 0,
		3, 0, 2
	};
	NBLdpcCodec codec(Length, NumOfChecks, 2, pCheckMatrix);

	unsigned Dimension = codec.GetDimension();
	GFSymbol pMsg[] = { 3, 2, 1 };
	GFSymbol pCW[Length];

	codec.Encode(pMsg, pCW);
	cout << codec.VerifyCodeword(pCW);
	return 0;
}