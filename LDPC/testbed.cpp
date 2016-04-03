#include <iostream>
#include "NBLdpcCodec.h"
#include "NBLdpcBuilder.h"
#include <algorithm>
#include <random>
#include <ctime>

using namespace std;

char *NextArgument()
{
	static int curArg = 1;
	return __argv[curArg++];
}

void Usage()
{
	cout << "Usage:\n"
		<< "B Length Redundancy Extension [R RowVals ColVals] specFileName\n"
		<< "T specFileName NumOfIterations\n";
}

int main()
{
	if(__argc != 9 && __argc != 4)
	{
		Usage();
		return 0;
	}
	char mode = NextArgument()[0];
	unsigned Length, Extension;
	string specFileName;
	if(mode == 'B')
	{
		Length = stoi(NextArgument());
		unsigned Redundancy = stoi(NextArgument());
		Extension = stoi(NextArgument());
		char type = NextArgument()[0];
		if(type == 'R')
		{
			unsigned RowVals = stoi(NextArgument()),
				ColVals = stoi(NextArgument());
			specFileName = NextArgument();
			NBLdpcBuilder::BuildRegular(Length, Redundancy, Extension, RowVals, ColVals, specFileName);
		}
	}
	else if(mode == 'T')
	{
		specFileName = NextArgument();
		unsigned NumOfIterations = stoi(NextArgument());
		NBLdpcCodec codec(specFileName);

		GFSymbol *pData = new GFSymbol[codec.GetDimension()];
		GFSymbol *pEncoded = new GFSymbol[codec.GetLength()];
		mt19937 gen(time(nullptr));
		uniform_int_distribution<int> distr(1, codec.GetMaxFieldElement());
		int errNum = 0;
		for (unsigned i = 0; i < NumOfIterations; ++i)
		{
			transform(pData, pData + codec.GetDimension(), pData, [&](GFSymbol a) { return (GFSymbol)distr(gen); });
			codec.Encode(pData, pEncoded);
			bool result = codec.VerifyCodeword(pEncoded);
			if (!result)
			{
				cerr << "Verification fail at iteration #" << i << endl;
				errNum++;
			}			
		}
		if (!errNum)
			cout << "Verification successful" << endl;
		else
			cout << "Verification failed" << endl;

		delete[] pData;
		delete[] pEncoded;
	}
	return 0;
}