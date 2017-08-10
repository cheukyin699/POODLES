#include "stdafx.h"

#define MAX_LINE 10000

/**********************************
Files.cpp
Copyright(C) 2006, Kana Shimizu
**********************************/

//　listに　全てのファイル名をpushする。
void SeqList::GetList(string listfile){
	ifstream inFile(listfile.c_str());
	if(!inFile){
		cerr<<"cannot find"<<listfile<<"\n";
		exit(1);
	}
	list.clear();
	string tmp="0";
	string cmp="1";
	while(!inFile.eof()){
		inFile >> tmp;
		if(cmp!=tmp)
			list.push_back(tmp);
		cmp=tmp;
	}
	inFile.close();
}

void SeqList::AppList(string listfile){
	ifstream inFile(listfile.c_str());
	if(!inFile){
		cerr<<"cannot find"<<listfile<<"\n";
		exit(1);
	}
	string tmp="0";
	string cmp="1";
	while(!inFile.eof()){
		inFile >> tmp;
		if(cmp!=tmp)
			list.push_back(tmp);
		cmp=tmp;
	}
	inFile.close();
}

void PBMatrix::ReadMtxFile(string matrixfile, vector<Elements> &seqRawData){

	int bPos;
	int aCnt=0;
	int mCnt=0;
	char buf[50000];

	matrixfile +=".mtx";
	//	matrixfile = MTXFILEPATH + matrixfile;

	ifstream inFile(matrixfile.c_str());
	if(!inFile){
		cerr<<"cannot find"<<matrixfile<<"\n";
		exit(1);
	}
	string str;
	//読み飛ばし
	Elements ele;
	int lCnt=0;
	while(!inFile.eof()){
		inFile.getline(buf,50000,'\n');
		string tmp(buf);
		if(lCnt>13 && tmp.size()>1){
			aCnt=0;
			for(mCnt=0;mCnt<24;mCnt++){
				if(mCnt!=0 && mCnt!=2 && mCnt!=21 && aCnt<NUM_OF_ACIDS){
					str = tmp.substr(0, tmp.find(" "));
					//		    cerr<<str<<" ";
					ele.aa[aCnt] = atof(str.c_str());
					aCnt++;
				}
				tmp = tmp.substr(tmp.find(" ")+2);
			}
			seqRawData.push_back(ele);
		}
		lCnt++;
	}
	inFile.close();
}

void FormatFile::SetFmtFile(string formatfile, vector<Elements> &seqRawData){

	char buf[MAX_LINE];
	int cnt=0;

	formatfile += ".FORMAT";
	//	formatfile = FMTFILEPATH + formatfile;

	ifstream inFile(formatfile.c_str());
	if(!inFile){
		cerr<<"cannot find"<<formatfile<<"\n";
		exit(1);
	}
	Elements ele;
	while(!inFile.eof()){
		inFile.getline(buf,MAX_LINE,'\n');
		if(buf[0]==NULL){}
		else if(buf[0]!='>'){
			seqRawData.push_back(ele);
		}
	}
	inFile.close();
}


/*format file はmatrixの後に読むこと！*/
/*もしくはformat file はsetFmtFileのあとに読むこと！*/
void FormatFile::ReadFmtFile(string formatfile, vector<Elements> &seqRawData){

	BioTables tbl;
	char buf[MAX_LINE];
	int cnt=0;
	char acid;
	int odf;//order or disorder

	int length=seqRawData.size();

	formatfile += ".FORMAT";
	//	formatfile = FMTFILEPATH + formatfile;

	ifstream inFile(formatfile.c_str());
	if(!inFile){
		cerr<<"cannot find"<<formatfile<<"\n";
		exit(1);
	}
	while(!inFile.eof()){
		inFile.getline(buf,MAX_LINE,'\n');
		if(buf[0]==NULL){}
		else if(buf[0]!='>'){
			acid = buf[0];
			odf = atoi(&buf[2]);
			seqRawData[cnt].value = odf;
			seqRawData[cnt].acid = tbl.ToAcidCode(acid);
			cnt++;
		}
	}
	if(cnt<2){
		cerr<<"irregular file (empty)"<<formatfile<<"\n";
		exit(1);
	}
	inFile.close();
	for(;cnt<length;cnt++){
		seqRawData.pop_back();
	}
}

int STRFile::ReadStrFile(string strfile, vector<Elements> &seqRawData){

	BioTables tbl;
	char buf[MAX_LINE];
	int cnt=0;
	char acid;
	int strID;//order or disorder

	int length=seqRawData.size();

	strfile += ".STR";
	//	strfile = STRFILEPATH + strfile;

	ifstream inFile(strfile.c_str());
	if(!inFile){
		cerr<<"cannot find"<<strfile<<"\n";
		return(0);
	}
	while(!inFile.eof()){
		inFile.getline(buf,MAX_LINE,'\n');
		if(buf[0]==NULL){}
		else if(buf[0]!='>'){
			acid = buf[0];
			switch(buf[2]){
				case 'H':
					strID =HELIX;
					break;

				case 'E':
					strID =SHEET;
					break;
				default:
					strID = COIL;
					break;
			}
			seqRawData[cnt].str = strID;
			cnt++;
		}
	}
	if(cnt<2){
		cerr<<"irregular file (empty)"<<strfile<<"\n";
		exit(1);
	}
	inFile.close();
	for(;cnt<length;cnt++){
		seqRawData.pop_back();
	}
}

/*fst file はmatrixの後に読むこと！*/
void FastaFile::ReadFstFile(string fstfile, vector<Elements> &seqRawData){

	BioTables tbl;
	char buf[MAX_LINE];
	int cnt=0;
	char acid;
	int odf;//order or disorder

	int length=seqRawData.size();

	//	fstfile = FSTFILEPATH + fstfile;
	fstfile += ".fst";

	ifstream inFile(fstfile.c_str());
	if(!inFile){
		cerr<<"cannot find"<<fstfile<<"\n";
		exit(1);
	}
	while(!inFile.eof()){
		inFile.getline(buf,MAX_LINE);
		string tmp=buf;
		if(tmp.find(">")==0){}else{
			for(int i=0;i<tmp.size();i++){
				acid = (tmp.substr(i,1).c_str())[0];
				seqRawData[cnt].acid = tbl.ToAcidCode(acid);
				cnt++;
			}
		}
		tmp.clear();
	}
	if(cnt<2){
		cerr<<"irregular file (empty)"<<fstfile<<"\n";
		exit(1);
	}
	inFile.close();
	for(;cnt<length;cnt++){
		seqRawData.pop_back();
	}
}
