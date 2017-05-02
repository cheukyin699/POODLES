#include "param.h"

/**********
　構造体
 ***********/

//matrix の一行分の生のデータ
class Elements{
public:
	int value; //ORDER 1, DISORDER 0
	int acid;
	int str;
	double aa[NUM_OF_ACIDS];
	Elements(){
		value=-1;
		acid=-1;
		str=COIL;
	}
};

struct ALLDATA{
	double profile[70];
	int value;
};



/**********
　クラス
 ***********/

class FSdata{
public:
	int acidL[NUM_OF_ACIDS];
	int kpairL[NUN_OF_PAIRS*NUM_OF_ACIDS*NUM_OF_ACIDS];
	int fpairL[NUN_OF_PAIRS*10*10];
	FSdata(){
		memset(acidL,0x00,NUM_OF_ACIDS*sizeof(int));
		memset(kpairL,0x00,NUN_OF_PAIRS*NUM_OF_ACIDS*NUM_OF_ACIDS*sizeof(int));
		memset(fpairL,0x00,NUN_OF_PAIRS*10*10*sizeof(int));
	}
};

class FSdataTmp{
public:
	double acid[NUM_OF_ACIDS];
	double kpair[NUN_OF_PAIRS*NUM_OF_ACIDS*NUM_OF_ACIDS];
	double fpair[NUN_OF_PAIRS*10*10];
	FSdataTmp(){
		memset(acid,0x00,NUM_OF_ACIDS*sizeof(double));
		memset(kpair,0x00,NUN_OF_PAIRS*NUM_OF_ACIDS*NUM_OF_ACIDS*sizeof(double));
		memset(fpair,0x00,NUN_OF_PAIRS*10*10*sizeof(double));
	}
};

class Result{
public:
	double st;
	double sp;
	double cc;
	double fr;

	double g;
	double c;

	double casCc;
	double casSp;
	double casSt;
	double casFr;
	
	Result(){
		st=0.0;sp=0.0;cc=0.0,fr=0.0;
		g=0.0;c=0.0;
		casCc=0.0;casSp=0.0;casSt=0.0;casFr=0.0;
	}
};

class DResult{
public:
	int tp;
	int tn;
	int fp;
	int fn;

	double st;
	double sp;
	double cc;
	double fr;
	DResult(){
		st=0.0;sp=0.0;cc=0.0,fr=0.0;
		tp=0;tn=0;fp=0;fn=0;
	}
	void CalcResult();
};

class FtrData{
public:
	int ftr[NUM_OF_FEATURE-2];
	FtrData(){
		memset(ftr,0x00,sizeof(int)*(NUM_OF_FEATURE-2));
	}
};

class BioTables{
public:
	double weight[NUM_OF_FEATURE];
	int feature[NUM_OF_ACIDS][NUM_OF_FEATURE-2];
	double flex[NUM_OF_ACIDS];
	double hydro[NUM_OF_ACIDS];
	double pairenergy[NUM_OF_ACIDS][NUM_OF_ACIDS];

	char ToAcidChar(int input);
	void ReadFeature(string ftrfile);
	void ReadPairEnergy(string pairenergyfile);
	int ToAcidCode(char input);
};

class Ranges{

public:
	int start[POSITIONCLASS];
	int stop[POSITIONCLASS];
	int min[LENGTHCLASS];
	int max[LENGTHCLASS];

	Ranges(){
		memset(start,0x00,sizeof(int)*POSITIONCLASS);
		memset(stop,0x00,sizeof(int)*POSITIONCLASS);
		memset(min,0x00,sizeof(int)*LENGTHCLASS);
		memset(max,0x00,sizeof(int)*LENGTHCLASS);
	}

	void ReadLength(char *filename);
	void ReadPosition(char *filenaem);
};

class Statistics{
public:
	string ftrfile;
	string listfile;
	void Composition20(int minlength,int maxlength, vector<double> *outData);
	void Composition20(int minlength,int maxlength, vector<int> *outData);
	void Composition20(int start,int stop,int Type);
	void Composition20(int start,int stop,int Type,vector<double> *outData);
	void Composition20();
	void Composition20(string pdbfile, string disprotfile, int start,int stop,int Type,int strID,vector<double> *outData);
	void DrLenDist();
	void DrLenDist(int lenID,int posID);
	void CompositionType(int length);
	void CalcQuaiSquare(int threshold);
	void EnergyDist(int frgLen,int skpLen);
	void FtrDist(int frgLen,int skpLen);
	void CalcAve();
	void CalcAllChiSqFromFile(char *filename);
	void FeatureSelection(int start,int stop);
};

class SeqList{
public:
	vector<string> list;

	void GetList(string listfile);
	void AppList(string listfile);
};

class PBMatrix{
public:
	void ReadMtxFile(string ftrfile,vector<Elements> &seqRawData);
};

class FormatFile{
public:
	void SetFmtFile(string ftrfile,vector<Elements> &seqRawData);
	void ReadFmtFile(string ftrfile,vector<Elements> &seqRawData);
};

class FastaFile{
public:
	void ReadFstFile(string ftrfile,vector<Elements> &seqRawData);
};

class STRFile{
public:
	int ReadStrFile(string ftrfile,vector<Elements> &seqRawData);
};

class PairMatrix{
public:
	string ftrfile;
	string listfile;
	string pairEfile;
	void CalcEnergyMatrix(void);
	void CalcChargeMatrix(void);
	void CalcHydroMatrix(void);
	void CalcHydroPMatrix(void);
	void CalcEnergyEachAcid(void);
};

class PairMtxData{
	double *data;
	int mtxSize;
	int mlFlg;
public:
	PairMtxData(){mlFlg=0;}
	~PairMtxData(){if(mlFlg==1){free(data);}}
	int GetSize();
	double GetValue(int i,int j);
	double GetLineSum(int i,int j);
	void MalocData();
	void ReadMtxData(string filename);
	double GetLectSum(int i, int j);
	double GetLectSabun(int i, int j,double ave);
	void ChangeTabComma(string filename);
};


class SeqWindow{
public:

	double *profile;
	double *raw;
	double flex;
	int value; // order or disorder
	int goodSmple;

	int windowSize;
	int frontWidth; // ----O++++ -の長さ

	SeqWindow(){flex=0.0;goodSmple=FALSE;}
	~SeqWindow(){delete profile,delete raw;}

	void GetProfile(int middlePos,BioTables tbl,vector<Elements> seqRawData);
	void GetRaw(int middlePos,BioTables tbl,vector<Elements> seqRawData);

	void PrintProfile();
	void PrintProfile(string input);
	void PrintWindowProfile();
	void PrintRaw();
	void PrintNetCharge(int middlePos,BioTables tbl,vector<Elements> seqRawData);
};

class NtermWindow: public SeqWindow{
public:
	NtermWindow(){
		windowSize=7;frontWidth=3;
		profile = new double[(N_WINDOW_SIZE)*(NUM_OF_FEATURE-1)];
		raw = new double[(N_WINDOW_SIZE)*(NUM_OF_ACIDS)];

		memset(profile,0x00,sizeof(double)*(N_WINDOW_SIZE)*(NUM_OF_FEATURE-1) );
		memset(raw,0x00,sizeof(double)*(N_WINDOW_SIZE)*(NUM_OF_ACIDS) );
	}
};

class MtermWindow: public SeqWindow{
public:
	MtermWindow(){
		windowSize=M_WINDOW_SIZE; frontWidth=M_WINDOW_SIZE/2;
		profile = new double[(M_WINDOW_SIZE)*(NUM_OF_FEATURE-1)];
		raw = new double[(M_WINDOW_SIZE)*(NUM_OF_ACIDS)];

		memset(profile,0x00,sizeof(double)*(M_WINDOW_SIZE)*(NUM_OF_FEATURE-1) );
		memset(raw,0x00,sizeof(double)*(M_WINDOW_SIZE)*(NUM_OF_ACIDS) );
	}
};

class CtermWindow: public SeqWindow{
public:
	CtermWindow(){
		windowSize=7;frontWidth=3;
		profile = new double[(C_WINDOW_SIZE)*(NUM_OF_FEATURE-1)];
		raw = new double[(C_WINDOW_SIZE)*(NUM_OF_ACIDS)];

		memset(profile,0x00,sizeof(double)*(C_WINDOW_SIZE)*(NUM_OF_FEATURE-1) );
		memset(raw,0x00,sizeof(double)*(C_WINDOW_SIZE)*(NUM_OF_ACIDS) );
	}
};


class SVMSmplStatus{
public:
	int AAlabel[20];
	int complexity;
	int pairpotential;
	int hydropathy;
	int windowSize;
	int slideSize;
	int startPos;
	int stopPos;
	int minDisLen;
	int maxDisLen;
	FtrData ftr;

	SVMSmplStatus(){
		memset(AAlabel,0x00,sizeof(int)*20);
		complexity=0;
		pairpotential=1;
		windowSize=0;
		slideSize=0;
		startPos=0;
		stopPos=0;
		hydropathy=0;
		minDisLen=0;
		maxDisLen=0;
	}
	void clear();
	void outToFile(char *filename);
	void InputFromFile(char *filename);
};


class SVMData{
	BioTables btl;
	int disSttl[POSITIONCLASS*LENGTHCLASS];
	int ordSttl[POSITIONCLASS*LENGTHCLASS];
	int curPos;
	int curLen;
	int disSmpl[CIBCB_POS_CLASS*NUM_OF_SND_STR];
	int ordSmpl[CIBCB_POS_CLASS*NUM_OF_SND_STR];

public:
	string ftrfile;
	string listfile;
	string listfile2;
	string totallistfile;

	//for CASP7
	string pdbfile;
	string disprotfile;

	int INITIALFLAG;
	SVMData(){
		memset(disSttl,0x00,sizeof(int)*POSITIONCLASS*LENGTHCLASS);
		memset(ordSttl,0x00,sizeof(int)*POSITIONCLASS*LENGTHCLASS);
		memset(disSmpl,0x00,sizeof(int)*CIBCB_POS_CLASS*NUM_OF_SND_STR);
		memset(ordSmpl,0x00,sizeof(int)*CIBCB_POS_CLASS*NUM_OF_SND_STR);

	}
	void BalanceSvmSmpl();
	void BalanceSvmSmpl(int red,char *input,char *output);

	void printSVMSmpl(vector<Elements> ele,vector<int> disSize,SVMSmplStatus *sss,PairMtxData *pd,ofstream* outFile);
	void MakeBatFile();
	void MakeGSBatFile();
	Result CVGS(int lenID, int posID, int cvrd);
	void allCVGS();
	void MakeCaspChkBatFile();
	void calcResLibSvm(char* answer, char* libsvm);
	void calcResLibSvm(char *answer, char* libSvm,Result &res);

	void ExecAllCycle();
	void ExecOneCycle(char* filename);
	void MakeOneSvmData(vector<Elements> ele, PairMtxData &pd);
	void MakeOneSvmData(vector<Elements> ele, PairMtxData &pd, char *posfile, char *lenfile,int slideSize);
	void printSVMSmpl(vector<Elements> ele,SVMSmplStatus *sss,PairMtxData *pd,ofstream* outFile);
	void ExecPredict(string option);
	void MergeResult(char* filename, vector<Elements> ele);

	void MakeFtrSmpl();
	void prinSVMFtrSmpl(vector<Elements> ele,SVMSmplStatus *sss,PairMtxData *pd,ofstream* outFile);
	void prinSVMFtrSmpl(vector<Elements> ele,vector<int> disSize,SVMSmplStatus *sss,PairMtxData *pd,ofstream* outFile);

	void MakeFtrTestSmpl();
	void MakeOneFtrData(vector<Elements> ele, PairMtxData &pd);
	void MakeTestSvmData();

	void EstimatePerformance(int posClass);

	void MoveTargetModels(char* directory);
	void FeatureSelection(int start,int stop,SVMSmplStatus &sss);
	void FeatureSelection(int start,int stop,int min, int max,SVMSmplStatus &sss);

	void MakeCIBCBData(char* allistname);
	void CalcEachFSPSCAV(int fs, int ps, Result &res);

	void CalcAllPSAV(int fs);

	void MakeNSData(char* allistname);
	void MakeNSCaspData(char* allistname);
	void CalcNSEachPSAV(int ps, Result &res);
	void CalcOSEachPSAV(int ps, vector<Result> &res);
	void CalcAllOSPAV();


	 //for casp7
	 void CalcDRRateOfSubData();
	 void CalcAllTgtRes(char *listfile,char *dir);
	 int CalcEachTgtRes(char *dir,char *tgtName, char *ansName,DResult &resAll);
	 void CalcAllTgtRes(char *listfile);
	 void CalcEachTgtRes(char *tgtName, char *ansName,DResult &resAll);
	 void MkCasp7Learn();
	 void MKCasp7Test(char *targetFile);
	 void ExCasp7Test(char *targetFile,char *dir);
	 void MrgCasp7Test(char *targetFile);
	 void FeatureSelection(string pdb, string disprot,int start,int stop,int strId,SVMSmplStatus &sss);
	 void printCasp7SVMSmpl(vector<Elements> ele,int strID,SVMSmplStatus *sss,PairMtxData *pd,ofstream *outFile);
	 void TrainCasp7SVM(char *dir);
	 void ScalesCasp7LearnData();
	 void BlncCasp7SVMSmpl();
	 int Saiki(int *seq,int jyouken,int cnt, int start,int stop);
	 void Casp7CV();
};
