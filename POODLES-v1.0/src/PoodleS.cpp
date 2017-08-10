// PoodleS.cpp : コンソール アプリケーションのエントリ ポイントを定義します。
//

#include "stdafx.h"

/**********************************
         PoodleS.cpp
Copyright(C) 2006, Kana Shimizu
**********************************/

char *POODLES_DIR;
char *SVM_DIR;


int main(int argc, char** argv)
{
    if(argc < 2){
        exit(1);
    }

    if(argv[1][0]=='-'){
        string feature("feature.txt");
        string instdir;

        SVMData sd;
        switch(argv[1][1]){
            case 'c':
                cerr<<"disorder2.exe -4 filename\n";
                instdir=argv[2];
                feature  = instdir +"/POODLES-v1.0/"  + feature;
                sd.ftrfile = feature;

                POODLES_DIR= (char *)malloc(sizeof(char)*(instdir.size()+16));
                SVM_DIR= (char *)malloc(sizeof(char)*(instdir.size()+10));
                sprintf(POODLES_DIR,"%s%s",instdir.c_str(),"/POODLES-v1.0");
                sprintf(SVM_DIR,"%s%s",instdir.c_str(),"/libsvm");
                sd.MKCasp7Test(argv[3]);
                break;
            case 'e':
                instdir=argv[2];
                feature  = instdir +"/POODLES-v1.0/"  + feature;
                sd.ftrfile = feature;
                POODLES_DIR= (char *)malloc(sizeof(char)*(instdir.size()+16));
                SVM_DIR= (char *)malloc(sizeof(char)*(instdir.size()+10));
                sprintf(POODLES_DIR,"%s%s",instdir.c_str(),"/POODLES-v1.0");
                sprintf(SVM_DIR,"%s%s",instdir.c_str(),"/libsvm");
                sd.ExCasp7Test(argv[3],argv[4]);
                break;
            case 'm':
                instdir=argv[2];
                feature  = instdir +"/POODLES-v1.0/"  + feature;
                sd.ftrfile = feature;
                sd.MrgCasp7Test(argv[3]);
                break;
        }
    }
    return 0;
}

