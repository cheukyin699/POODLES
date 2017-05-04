#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<algorithm>
#include<math.h>

using namespace std;

/**********************************
         param.h
Copyright(C) 2006, Kana Shimizu
**********************************/

#define SCLDDATAPATH ".\\scaled_data\\"
#define SCLDDATATESTPATH ".\\scaled_dataTest\\"
#define BLNCDATAPATH ".\\blcData\\"
#define SVMDATAPATH ".\\svmData\\"
#define SVMTRAINPATH ".\\"
#define SVMFTRDATAPATH ".\\svmFtr\\"
#define TMPAREAPATH ".\\tmp\\"
#define TESTDATAPATH ".\\CaspDataTest\\"
#define TESTFTRDATAPATH ".\\CaspDataTestFtr\\"
#define INTEGRATEDPATH ".\\integrated\\"
#define ULBLDDATAPATH ".\\unprtData\\"
#define SCULBLDDATAPATH ".\\scaled_unprtData\\"

#define ESTIMATEPATH ".\\estimated\\"
#define DISOPREDPATH ".\\disopred\\"

#define CIBCB_CROSS_VALID_NUM 5
#define CIBCB_POS_CLASS 7
#define CIBCB_POS_OS_CLASS 3

#define IO_ERR -1
#define IO_NORM 0

#define GLAY 2
#define ORDER 1
#define DISORDER 0

#define INVALID -100

#define CTERM 1
#define NTERM 0

#define TRUE 1
#define FALSE 0

#define DP_TRAN_PENALTY -0.15
#define DISORDER_AVE 0.011
#define ORDER_AVE 0.014

#define MIN_LENGTH_ARRAY 30

//#define DP_TRAN_PENALTY -0.15
//#define DISORDER_AVE 0.009
//#define ORDER_AVE 0.012

#define LENGTHCLASS 7
#define POSITIONCLASS 6

#define NUM_OF_ACIDS 21
#define NUN_OF_PAIRS 3
#define CHOSEN_F 20
#define SGTTHRESHOLD 0
#define SGTTHRESHOLD_ORD 0.5
#define SGTTHRESHOLD_DIS -0.5

#define Ala 0
#define Cys 1
#define Asp 2
#define Glu 3
#define Phe 4
#define Gly 5
#define His 6
#define Ile 7
#define Lys 8
#define Leu 9
#define Met 10
#define Asn 11
#define Pro 12
#define Gln 13
#define Arg 14
#define Ser 15
#define Thr 16
#define Val 17
#define Trp 18
#define Tyr 19
#define UNKNOWN 20

#define HELIX 0
#define SHEET 1
#define COIL 2
#define ANY 3
#define FORTEST 5
#define NUM_OF_SND_STR 4

#define HYDROPHOBIC 0
#define HYDROPHILIC 1
#define CHARGED 2
#define POSITIVE 3
#define NEGATIVE 4
#define AROMATIC 5
#define ALIPHATIC 6
#define TINY 7
#define SMALL 8
#define POLAR 9


#define AAWEIGHT 0.001

#define NUM_OF_FEATURE 12
//#define PROFILE_PARAM 0.2

#define N_WINDOW_SIZE 7
#define M_WINDOW_SIZE 15
#define C_WINDOW_SIZE 7

#define NBAND 10
#define CBAND 10

#define PROFILE(seqPos,ftrPos) profile[seqPos*(NUM_OF_FEATURE-1) + ftrPos ]
#define RAW(seqPos,acid) raw[ seqPos*(NUM_OF_ACIDS) + acid ]
