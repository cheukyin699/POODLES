#include "stdafx.h"

/**********************************
         BioTables.cpp
 Copyright(C) 2006, Kana Shimizu
**********************************/

void BioTables::ReadPairEnergy(string pefile){

    char buf[256];
    int pos=0;
    ifstream inFile(pefile.c_str());

    if(!inFile){
        cout<<pefile+" :can not open!\n";
    }

    for(int cntF=0;cntF<20;cntF++){
        inFile.getline(buf,256,'\n');
        string tmp(buf);
        for(int cntB=0;cntB<20;cntB++){
            pos=tmp.find(',');
            pairenergy[cntF][cntB] = atof(tmp.substr(0,pos).c_str());
            tmp = tmp.substr(pos+1);
        }
    }
}

void BioTables::ReadFeature(string ftrfile){

    char buf[256];
    int pos=0;
    ifstream inFile(ftrfile.c_str());

    if(!inFile){
        cout<<ftrfile+" :can not open!\n";
        exit(1);
    }

    //weight ‚ÌŽæ“¾
    inFile.getline(buf,256,'\n');
    string tmp(buf);
    for(int wCnt=0;wCnt<NUM_OF_FEATURE-1;wCnt++){
        pos = tmp.find('\t');
        weight[wCnt] = atof( tmp.substr(0,pos).c_str() );
        tmp = tmp.substr(pos+1);
    }
    weight[NUM_OF_FEATURE-1] =  atof(tmp.c_str());

    //“Á’¥ ‚ÌŽæ“¾
    int aCnt=0;
    while(!inFile.eof()){
        inFile.getline(buf,256,'\n');
        string tmp(buf);
        for(int pCnt=0;pCnt<NUM_OF_FEATURE-2;pCnt++){
            pos = tmp.find('\t');
            feature[aCnt][pCnt] = atoi( tmp.substr(0,pos).c_str() );
            tmp = tmp.substr(pos+1);
        }
        pos = tmp.find('\t');
        flex[aCnt] = atof( tmp.substr(0,pos).c_str() );
        tmp = tmp.substr(pos+1);

        hydro[aCnt] =  atof(tmp.c_str());
        aCnt++;
    }

    inFile.close();
}

char BioTables::ToAcidChar(int input){

    char ret;

    switch(input){
        case Ala:
            ret ='A';
            break;
        case Cys:
            ret ='C';
            break;
        case Asp:
            ret = 'D';
            break;
        case Glu:
            ret ='E';
            break;
        case Gly:
            ret ='G';
            break;
        case Phe:
            ret ='F';
            break;
        case His:
            ret ='H';
            break;
        case Ile:
            ret ='I';
            break;
        case Lys:
            ret ='K';
            break;
        case Leu:
            ret ='L';
            break;
        case Met:
            ret ='M';
            break;
        case Asn:
            ret ='N';
            break;
        case Pro:
            ret ='P';
            break;
        case Gln:
            ret ='Q';
            break;
        case Arg:
            ret ='R';
            break;
        case Ser:
            ret ='S';
            break;
        case Thr:
            ret = 'T';
            break;
        case Val:
            ret ='V';
            break;
        case Trp:
            ret ='W';
            break;
        case Tyr:
            ret ='Y';
            break;
        default:
            ret = 'Z';
            break;
    }
    return(ret);
}

int BioTables::ToAcidCode(char input){

    int ret;

    switch(input){
        case 'A':
            ret = Ala;
            break;
        case 'C':
            ret =Cys;
            break;
        case 'D':
            ret = Asp;
            break;
        case 'E':
            ret =Glu;
            break;
        case 'G':
            ret =Gly;
            break;
        case 'F':
            ret =Phe;
            break;
        case 'H':
            ret =His;
            break;
        case 'I':
            ret =Ile;
            break;
        case 'K':
            ret =Lys;
            break;
        case 'L':
            ret =Leu;
            break;
        case 'M':
            ret =Met;
            break;
        case 'N':
            ret =Asn;
            break;
        case 'P':
            ret =Pro;
            break;
        case 'Q':
            ret =Gln;
            break;
        case 'R':
            ret =Arg;
            break;
        case 'S':
            ret =Ser;
            break;
        case 'T':
            ret = Thr;
            break;
        case 'V':
            ret =Val;
            break;
        case 'W':
            ret = Trp;
            break;
        case 'Y':
            ret = Tyr;
            break;
        default:
            ret = UNKNOWN;
            break;
    }
    return(ret);
}
