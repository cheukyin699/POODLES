#include "stdafx.h"

/**********************************
         SVMsmplStatus.cpp
Copyright(C) 2006, Kana Shimizu
**********************************/

void SVMSmplStatus::outToFile(char *filename){
	ofstream outFile(filename);
	if(!outFile){
		cerr<<filename<<"\n";
		exit(1);
	}

	for(int i=0;i<20;i++){
		outFile<<AAlabel[i]<<"\n";
	}

	outFile<<complexity<<"\n";
	outFile<<pairpotential<<"\n";
	outFile<<hydropathy<<"\n";
	outFile<<windowSize<<"\n";
	outFile<<slideSize<<"\n";
	outFile<<startPos<<"\n";
	outFile<<stopPos<<"\n";
	outFile<<minDisLen<<"\n";
	outFile<<maxDisLen<<"\n";
	for(int i=0;i<NUM_OF_FEATURE-2;i++){
		outFile<<ftr.ftr[i]<<"\n";	
	}

	outFile.close();
}

void SVMSmplStatus::InputFromFile(char *filename){
	ifstream inFile(filename);
	if(!inFile){
		cerr<<filename<<"\n";
		exit(1);
	}

	for(int i=0;i<20;i++){
		inFile>>AAlabel[i];
	}

	inFile>>complexity;
	inFile>>pairpotential;
	inFile>>hydropathy;
	inFile>>windowSize;
	inFile>>slideSize;
	inFile>>startPos;
	inFile>>stopPos;
	inFile>>minDisLen;
	inFile>>maxDisLen;
	for(int i=0;i<NUM_OF_FEATURE-2;i++){
		inFile>>ftr.ftr[i];	
	}

	inFile.close();

}
