#include "stdafx.h"

/**********************************
         server.cpp
Copyright(C) 2006, Kana Shimizu
**********************************/

extern char *POODLES_DIR;
extern char *SVM_DIR;


void SVMData::MrgCasp7Test(char *targetFile){
	char buf[256];
	int startA[CIBCB_POS_CLASS];
	int stopA[CIBCB_POS_CLASS];
	double *outVal;

	vector<Elements> ele;
	PBMatrix pm;
	FastaFile fs;

	pm.ReadMtxFile(targetFile,ele);
	fs.ReadFstFile(targetFile,ele);
	for(int i=0;i<ele.size();i++){
		cerr<<btl.ToAcidChar(ele[i].acid)<<",";	
	}

	string outfilename="./mgr/";
	outfilename += targetFile;
	outfilename +=".txt";
	ofstream outFile(outfilename.c_str());

	outfilename = "./caspFmt/";
	outfilename += targetFile;
	outfilename +=".txt";
	ofstream caspFile(outfilename.c_str());

	outVal = (double *)malloc(sizeof(double)*CIBCB_POS_CLASS*NUM_OF_SND_STR*ele.size());
	memset(outVal,0x00,sizeof(double)*CIBCB_POS_CLASS*NUM_OF_SND_STR*ele.size());

	startA[0]=2;stopA[0]=10;
	startA[1]=5;stopA[1]=15;
	startA[2]=10;stopA[2]=25;
	startA[3]=20;stopA[3]=-35;
	startA[4]=-15;stopA[4]=-40;
	startA[5]=-5;stopA[5]=-20;
	startA[6]=-2;stopA[6]=-10;

	string tmp;
	for(int ps=0;ps<CIBCB_POS_CLASS;ps++){
		if(startA[ps]<0){
			startA[ps] =ele.size()+stopA[ps];
		}

		for(int strID=0;strID<NUM_OF_SND_STR;strID++){
			sprintf(buf,"./rs_targets/%s%d-%d.txt",targetFile,ps,strID);
			ifstream inFile(buf);
			if(!inFile){
				cerr<<buf<<"\n";
				exit(1);
			}
			int cnt=0;
			curPos=0;
			while(!inFile.eof()){
				inFile>>tmp;
				if(cnt%3 == 2 && cnt!=2){
					outVal[ ps*(NUM_OF_SND_STR*ele.size())+strID*ele.size()+ curPos +startA[ps] ] = atof(tmp.c_str());
					curPos++;
				}
				cnt++;
			}
			inFile.close();
		}
	}

	for(int strID=0;strID<NUM_OF_SND_STR;strID++){
		outVal[ 0*(NUM_OF_SND_STR*ele.size())+strID*ele.size()+ 0] = outVal[ 0*(NUM_OF_SND_STR*ele.size())+strID*ele.size()+ 2];
		outVal[ 0*(NUM_OF_SND_STR*ele.size())+strID*ele.size()+ 1] = outVal[ 0*(NUM_OF_SND_STR*ele.size())+strID*ele.size()+ 2];
	}

	for(int strID=0;strID<NUM_OF_SND_STR;strID++){
		outVal[ 6*(NUM_OF_SND_STR*ele.size())+strID*ele.size()+ ele.size()-1] = outVal[ 6*(NUM_OF_SND_STR*ele.size())+strID*ele.size()+ ele.size()-4];
		outVal[ 6*(NUM_OF_SND_STR*ele.size())+strID*ele.size()+ ele.size()-2] = outVal[ 6*(NUM_OF_SND_STR*ele.size())+strID*ele.size()+ ele.size()-4];
		outVal[ 6*(NUM_OF_SND_STR*ele.size())+strID*ele.size()+ ele.size()-3] = outVal[ 6*(NUM_OF_SND_STR*ele.size())+strID*ele.size()+ ele.size()-4];
	}


	double *cspOut;
	cspOut = (double *)malloc(sizeof(double)*NUM_OF_SND_STR*ele.size());
	memset(cspOut,0x00,sizeof(double)*NUM_OF_SND_STR*ele.size());
	double tmpD;
	double tmpCnt;
	for(int i=0;i<ele.size();i++){
		for(int strID=0;strID<NUM_OF_SND_STR;strID++){
			tmpD=0;
			tmpCnt=0;
			for(int ps=0;ps<CIBCB_POS_CLASS;ps++){
				if(outVal[ ps*(NUM_OF_SND_STR*ele.size())+strID*ele.size()+ i]==0){
				}else{
					tmpD += outVal[ ps*(NUM_OF_SND_STR*ele.size())+strID*ele.size()+ i];
					tmpCnt++;
				}
			}
			if(tmpCnt==1 || tmpCnt==2){
				cspOut[strID*ele.size() + i] = tmpD/tmpCnt;
			}else{
				cerr<<"not 1 or 2\n";
			}
		}
	}
	
	//smoothing
	vector<double> cspFnl;
	for(int i=0;i<ele.size();i++){
		tmpD=0;
		tmpCnt=0;
		for(int j=i-2;j<i+3;j++){
			if(j>=0 && j<ele.size()){
				tmpCnt++;
				tmpD +=cspOut[3*ele.size() + j];
			}
		}
		cspFnl.push_back(tmpD/tmpCnt);
	}

	//MRG out
	for(int i=0;i<ele.size();i++){
		outFile<<btl.ToAcidChar(ele[i].acid)<<",";	
		for(int ps=0;ps<CIBCB_POS_CLASS;ps++){
			for(int strID=0;strID<NUM_OF_SND_STR;strID++){
				outFile.precision(2);
				outFile.width(3);
				outFile<<outVal[ ps*(NUM_OF_SND_STR*ele.size())+strID*ele.size()+ i]<<",";
			}
		}
		for(int strID=0;strID<NUM_OF_SND_STR;strID++){
			outFile<<cspOut[strID*ele.size() + i]<<",";
		}
		outFile<<cspFnl[i]<<"\n";
	}
	outFile.close();

	free(cspOut);


	caspFile<<"PFRMAT DR\n";
	caspFile<<"REMARK K. Shimizu, Y. Muraoka, S. Hirose, and T. Noguchi\n";
	caspFile<<"REMARK \"Feature Selection Based on Physicochemical Properties of\n";

	caspFile<<"REMARK Redefined N-term Region and C-term Regions for Predicting Disorder\"\n";
	caspFile<<"REMARK Proc. of IEEE CIBCB 2005, pp262-267. \n";
	caspFile<<"METHOD Prediction for short disorder using modified PSSM\n";
	caspFile<<"METHOD ------\n";

	for(int i=0;i<ele.size();i++){
		caspFile<<btl.ToAcidChar(ele[i].acid)<<" ";
		if(cspFnl[i]>0.5){
			caspFile<<"D ";
		}else{
			caspFile<<"O ";
		}
		caspFile.precision(3);
		caspFile<<cspFnl[i]<<"\n";	
	}

	caspFile<<"END\n";

	caspFile<<"\n\n";
	caspFile<<"Copyright:\nComputational Biology Research Center (CBRC), AIST\n";
	caspFile<<"Pharma Design, Inc.\n\n";
	
	caspFile.close();
}

void SVMData::ExCasp7Test(char *targetFile,char *dir){
	char buf[1024];
	ofstream scaleFile("oneScale.sh");
	for(int ps=0;ps<CIBCB_POS_CLASS;ps++){
		for(int strID=0;strID<NUM_OF_SND_STR;strID++){
			sprintf(buf,"%s/svm-scale -r %s/scales_param/sc_param_%d-%d.txt ./targets/%s%d-%d.txt >  ./sc_targets/%s%d-%d.txt \n",SVM_DIR,POODLES_DIR,ps,strID,targetFile,ps,strID,targetFile,ps,strID);
			string tmp(buf);
			scaleFile<<tmp<<"\n";
		}
	}
	scaleFile.close();


	ofstream predFile("predict.sh");
	for(int ps=0;ps<CIBCB_POS_CLASS;ps++){
		for(int strID=0;strID<NUM_OF_SND_STR;strID++){
			sprintf(buf,"%s/svm-predict -b 1 ./sc_targets/%s%d-%d.txt %s/models/%s/model%d-%d.txt ./rs_targets/%s%d-%d.txt \n",SVM_DIR,targetFile,ps,strID,POODLES_DIR,dir,ps,strID,targetFile,ps,strID);
			string tmp(buf);
			predFile<<tmp<<"\n";
		}
	}
	predFile.close();
}


void SVMData::MKCasp7Test(char *targetFile){
	char buf[50];
	btl.ReadFeature(ftrfile);

	SVMSmplStatus sss[CIBCB_POS_CLASS*NUM_OF_SND_STR];


	for(int ps=0;ps<CIBCB_POS_CLASS;ps++){
		for(int strID=0;strID<NUM_OF_SND_STR;strID++){
			sprintf(buf,"%s/statusFile/%d-%d.txt\0",POODLES_DIR,ps,strID);
			sss[NUM_OF_SND_STR*ps+strID].InputFromFile(buf);
		}
	}

	cerr<<"CREATE casp FILE.\n";
	ofstream caspFile[CIBCB_POS_CLASS*NUM_OF_SND_STR];

	for(int ps=0;ps<CIBCB_POS_CLASS;ps++){
		for(int strID=0;strID<NUM_OF_SND_STR;strID++){
			sprintf(buf,"./targets/%s%d-%d.txt",targetFile,ps,strID);
			caspFile[NUM_OF_SND_STR*ps+strID].open(buf);
			if(!caspFile[NUM_OF_SND_STR*ps+strID]){
			    cerr<<"can not open "<<buf<<"\n";
			    exit(1);
			}
		}
	}

	int w[CIBCB_POS_CLASS];
	int startA[CIBCB_POS_CLASS];
	int stopA[CIBCB_POS_CLASS];

	for(int ps=0;ps<CIBCB_POS_CLASS;ps++){
		for(int strID=0;strID<NUM_OF_SND_STR;strID++){
			sss[NUM_OF_SND_STR*ps+strID].pairpotential=0;
			sss[NUM_OF_SND_STR*ps+strID].slideSize =1;
			sss[NUM_OF_SND_STR*ps+strID].minDisLen=0;
			sss[NUM_OF_SND_STR*ps+strID].maxDisLen=100000;
			sss[NUM_OF_SND_STR*ps+strID].windowSize=5;
		}
	}

	cerr<<"sat param\n";

	w[0]=5;	w[1]=5; w[2]=10;	w[3]=15;	w[4]=10;	w[5]=10;	w[6]=5;	

	startA[0]=0;stopA[0]=12; //0-10+2
	startA[1]=3;stopA[1]=17; //5-15 -2+2
	startA[2]=5;stopA[2]=29; //-5 +4
	startA[3]=13;stopA[3]=-28;//-7 +7
	startA[4]=-11;stopA[4]=-45;
	startA[5]=-1;stopA[5]=-25;
	startA[6]=0;stopA[6]=-12;

	//dummy
	PairMtxData pd;
	//dummy
	curLen=0;curPos=0;

	string target(targetFile);

	vector<Elements> ele;
	PBMatrix pm;
	FastaFile ff;
	STRFile st;

	cerr<<"read mtx ";
	pm.ReadMtxFile(target,ele);
	cerr<<"... done\nread fasta ";
	ff.ReadFstFile(target,ele);
	cerr<<"... done\n";
	st.ReadStrFile(target,ele);
	cout<<target<<" "<<ele.size()<<"\n";
	cout.flush();

	//createdata
	for(int ps=0;ps<CIBCB_POS_CLASS;ps++){
		for(int strID=0;strID<NUM_OF_SND_STR;strID++){
			sss[NUM_OF_SND_STR*ps+strID].windowSize = w[ps];
			if(stopA[ps]>0){
				sss[NUM_OF_SND_STR*ps+strID].startPos =startA[ps];
				sss[NUM_OF_SND_STR*ps+strID].stopPos =stopA[ps];
			}else if(startA[ps]>0){
				sss[NUM_OF_SND_STR*ps+strID].startPos =startA[ps];
				sss[NUM_OF_SND_STR*ps+strID].stopPos =ele.size()+stopA[ps];
			}else{
				sss[NUM_OF_SND_STR*ps+strID].startPos =ele.size()+stopA[ps];
				sss[NUM_OF_SND_STR*ps+strID].stopPos =ele.size()+startA[ps];
			}
			if(sss[NUM_OF_SND_STR*ps+strID].stopPos<0){
				sss[NUM_OF_SND_STR*ps+strID].stopPos=0;
			}
			printCasp7SVMSmpl(ele,FORTEST,&(sss[NUM_OF_SND_STR*ps+strID]),&pd,&caspFile[NUM_OF_SND_STR*ps+strID]);
		}
	}

	for(int ps=0;ps<CIBCB_POS_CLASS;ps++){
		for(int strID=0;strID<NUM_OF_SND_STR;strID++){
			caspFile[NUM_OF_SND_STR*ps+strID].close();
		}
	}
}

void SVMData::MkCasp7Learn(){
	
	char buf[256];
	btl.ReadFeature(ftrfile);

	SVMSmplStatus sss[CIBCB_POS_CLASS*NUM_OF_SND_STR];

	for(int strID=0;strID<NUM_OF_SND_STR;strID++){
		FeatureSelection(pdbfile,disprotfile,0,10,strID,sss[NUM_OF_SND_STR*0+strID]);
		FeatureSelection(pdbfile,disprotfile,5,15,strID,sss[NUM_OF_SND_STR*1+strID]);
		FeatureSelection(pdbfile,disprotfile,10,25,strID,sss[NUM_OF_SND_STR*2+strID]);
		FeatureSelection(pdbfile,disprotfile,20,-35,strID,sss[NUM_OF_SND_STR*3+strID]);
		FeatureSelection(pdbfile,disprotfile,-15,-40,strID,sss[NUM_OF_SND_STR*4+strID]);
		FeatureSelection(pdbfile,disprotfile,-5,-20,strID,sss[NUM_OF_SND_STR*5+strID]);
		FeatureSelection(pdbfile,disprotfile,0,-10,strID,sss[NUM_OF_SND_STR*6+strID]);
	}

	for(int ps=0;ps<CIBCB_POS_CLASS;ps++){
		for(int strID=0;strID<NUM_OF_SND_STR;strID++){
			if(strID==COIL || strID==ANY){
				sss[NUM_OF_SND_STR*ps+strID].hydropathy=1;
				sss[NUM_OF_SND_STR*ps+strID].complexity=1;
			}
			sprintf(buf,".\\CASP7\\statusFile\\%d-%d.txt\0",ps,strID);
			sss[NUM_OF_SND_STR*ps+strID].outToFile(buf);
		}
	}

	cerr<<"CREATE LEARN FILE\n";
	ofstream learnFile[CIBCB_POS_CLASS*NUM_OF_SND_STR];

	for(int ps=0;ps<CIBCB_POS_CLASS;ps++){
		for(int strID=0;strID<NUM_OF_SND_STR;strID++){
			sprintf(buf,".\\CASP7\\LearnOrg%d-%d.txt",ps,strID);
			learnFile[NUM_OF_SND_STR*ps+strID].open(buf);
		}
	}

	int w[CIBCB_POS_CLASS];
	int startA[CIBCB_POS_CLASS];
	int stopA[CIBCB_POS_CLASS];

	for(int ps=0;ps<CIBCB_POS_CLASS;ps++){
		for(int strID=0;strID<NUM_OF_SND_STR;strID++){
			sss[NUM_OF_SND_STR*ps+strID].pairpotential=0;
			sss[NUM_OF_SND_STR*ps+strID].slideSize =1;
			sss[NUM_OF_SND_STR*ps+strID].minDisLen=0;
			sss[NUM_OF_SND_STR*ps+strID].maxDisLen=100000;
			sss[NUM_OF_SND_STR*ps+strID].windowSize=5;
		}
	}

	w[0]=5;	w[1]=5; w[2]=10;	w[3]=15;	w[4]=10;	w[5]=10;	w[6]=5;	

	startA[0]=0;stopA[0]=10;
	startA[1]=5;stopA[1]=15;
	startA[2]=10;stopA[2]=25;
	startA[3]=20;stopA[3]=-35;
	startA[4]=-15;stopA[4]=-40;
	startA[5]=-5;stopA[5]=-20;
	startA[6]=0;stopA[6]=-10;

	//dummy
	PairMtxData pd;
	//dummy
	curLen=0;curPos=0;

	SeqList list;
	list.GetList(pdbfile);
	list.AppList(disprotfile);

	for(int listCnt=0;listCnt<list.list.size();listCnt++){

		vector<Elements> ele;
		PBMatrix pm;
		FormatFile ff;
		STRFile st;

		pm.ReadMtxFile(list.list[listCnt],ele);
		ff.ReadFmtFile(list.list[listCnt],ele);
		st.ReadStrFile(list.list[listCnt],ele);

		cout<<listCnt<<" "<<list.list[listCnt]<<" "<<ele.size()<<"\n";
		cout.flush();

		//createdata
		for(int ps=0;ps<CIBCB_POS_CLASS;ps++){
			for(int strID=0;strID<NUM_OF_SND_STR;strID++){
				sss[NUM_OF_SND_STR*ps+strID].windowSize = w[ps];
				if(ps==3 && strID==ANY){
//					sss[NUM_OF_SND_STR*ps+strID].slideSize =7;		
				}
				if(stopA[ps]>0){
					sss[NUM_OF_SND_STR*ps+strID].startPos =startA[ps];
					sss[NUM_OF_SND_STR*ps+strID].stopPos =stopA[ps];
				}else if(startA[ps]>0){
					sss[NUM_OF_SND_STR*ps+strID].startPos =startA[ps];
					sss[NUM_OF_SND_STR*ps+strID].stopPos =ele.size()+stopA[ps];
				}else{
					sss[NUM_OF_SND_STR*ps+strID].startPos =ele.size()+stopA[ps];
					sss[NUM_OF_SND_STR*ps+strID].stopPos =ele.size()+startA[ps];
				}
				if(sss[NUM_OF_SND_STR*ps+strID].stopPos<0){
					sss[NUM_OF_SND_STR*ps+strID].stopPos=0;
				}
				curPos=ps;
				printCasp7SVMSmpl(ele,strID,&(sss[NUM_OF_SND_STR*ps+strID]),&pd,&learnFile[NUM_OF_SND_STR*ps+strID]);
			}
		}
	}

	for(int ps=0;ps<CIBCB_POS_CLASS;ps++){
		for(int strID=0;strID<NUM_OF_SND_STR;strID++){
			cerr<<ps<<"-"<<strID<<":"<<disSmpl[NUM_OF_SND_STR*ps+strID]<<"/"<<ordSmpl[NUM_OF_SND_STR*ps+strID]<<"\n";
			learnFile[NUM_OF_SND_STR*ps+strID].close();
		}
	}

}

void SVMData::FeatureSelection(string pdb, string disprot,int start,int stop,int strId,SVMSmplStatus &sss){


	vector<double> dis;vector<double> ord;
	Statistics st;

	st.ftrfile = ftrfile;
	btl.ReadFeature(ftrfile);

	st.Composition20(pdb,disprot,start,stop,DISORDER,COIL,&dis);
	st.Composition20(pdb,disprot,start,stop,ORDER,strId,&ord);

	for(int i=0;i<20;i++){
	//	cout<<dis[i]/ord[i]<<",";
	}
//	cout<<"\n";

	double total=0,sabun=0;
	int fs[NUM_OF_FEATURE-2];
	int AA[20];
	memset(AA,0x00,sizeof(int)*20);

	for(int j=0;j<NUM_OF_FEATURE-2;j++){
		for(int i=0;i<20;i++){
			sabun += btl.feature[i][j]*(ord[i]-dis[i]);
			total += btl.feature[i][j]*fabs(ord[i]-dis[i]);
		}
		cout<<sabun/total<<",";
		if((fabs(sabun)/total)<1.0){
			fs[j]=0;
		}else{
			fs[j]=1;
		}
		total=0;sabun=0;
	}
	cout<<"\n";
	/**
	for(int j=0;j<NUM_OF_FEATURE-2;j++){
		sss.ftr.ftr[j] = fs[j];
		if(fs[j]==1){
			for(int i=0;i<20;i++){
				AA[i] += btl.feature[i][j];
			}
		}
	}
	int rdCnt=0;
	for(int i=0;i<20;i++){
		if(AA[i]==0){
			sss.AAlabel[i]=1;
			rdCnt++;
		}
	}
	cout<<"\n"<<20-rdCnt<<"\n";
	**/
}

void Statistics::Composition20(string pdbfile, string disprotfile, int start,int stop,int Type,int strID,vector<double> *outData){
	BioTables btl;
	SeqList list;
	double alpha=0;

	list.GetList(pdbfile);
	list.AppList(disprotfile);
	
	int totalAcids[NUM_OF_ACIDS];

	memset(totalAcids,0x00,sizeof(int)*NUM_OF_ACIDS);

	int startPos,stopPos;
	int incriment=1;
	int band=0;

	startPos = start;
	stopPos = stop;

	for(int i=0;i<list.list.size();i++){

		vector<Elements> ele;
		STRFile st;
		FormatFile ff;

		ff.SetFmtFile(list.list[i],ele);
		ff.ReadFmtFile(list.list[i],ele);
		st.ReadStrFile(list.list[i],ele);

		if(stop>0){
			startPos = start;
			stopPos = stop;
		}else if(start>0){
			startPos = start;
			stopPos = ele.size()+stop;			
		}else{
			startPos = ele.size()+stop;
			stopPos = ele.size()+start;
		}

		if(stop<0){
			startPos = 0;
		}

		for(int aCnt=startPos;aCnt<stopPos;aCnt++){
			if(ele[aCnt].acid>=NUM_OF_ACIDS){
			}
			else if(strID==ANY){
				if(ele[aCnt].value==Type){
					totalAcids[ele[aCnt].acid]++;				
				}
			}
			else if(ele[aCnt].value==Type &&ele[aCnt].str==strID){
				totalAcids[ele[aCnt].acid]++;			
			}
		}
	}

	int totalSize=0;

	for(int aCnt=0;aCnt<NUM_OF_ACIDS;aCnt++){
		totalSize+=totalAcids[aCnt];
	}

	cerr<<"totalSizeof"<<strID<<": "<<totalSize<<"\n";
	for(int aCnt=0;aCnt<NUM_OF_ACIDS;aCnt++){
		outData->push_back((double)totalAcids[aCnt]/totalSize);
		cerr<<(double)totalAcids[aCnt]/totalSize<<",";
	}
	cerr<<"\n";
}

void SVMData::printCasp7SVMSmpl(vector<Elements> ele,int strID,SVMSmplStatus *sss,PairMtxData *pd,ofstream *outFile){

	int startPos,stopPos;
	double ftrVl[NUM_OF_FEATURE-2];
	int dis,ord,outLabel;
	double hydropathy;
	int complex[21];

	int sCnt=0,aCnt=0;
	for(sCnt=0;;sCnt++){
		startPos = sss->startPos+sCnt*sss->slideSize;
		stopPos = startPos + sss->windowSize;

		if(startPos<0){
			break;
		}
		if(stopPos>=ele.size()){
			break;
		}
		if(stopPos>sss->stopPos){
			break;
		}

		dis=0;ord=0;
		for(aCnt=startPos;aCnt<stopPos;aCnt++){
			if(ele[aCnt].value==DISORDER){
				dis++;
			}else if(strID==ANY && ele[aCnt].value==ORDER ){
				ord++;
			}
			else if(ele[aCnt].value==ORDER && ele[aCnt].str==strID){
				ord++;
			}
		}

		outLabel=INVALID;
		if(strID==FORTEST){
			outLabel=ORDER;
		}else if(dis==sss->windowSize){
			outLabel=DISORDER;
			disSmpl[NUM_OF_SND_STR*curPos + strID]++;
		}else if(ord==sss->windowSize){
			outLabel=ORDER;
			ordSmpl[NUM_OF_SND_STR*curPos + strID]++;
		}else if(sss->windowSize==10){
			if(ord>6){
				outLabel=ORDER;
				ordSmpl[NUM_OF_SND_STR*curPos + strID]++;
			}
		}else if(sss->windowSize==15){
			if(ord>8){
				outLabel=ORDER;
				ordSmpl[NUM_OF_SND_STR*curPos + strID]++;
			}
		}

		if(outLabel!=INVALID){
			(*outFile)<<outLabel<<" ";
			for(aCnt=startPos;aCnt<stopPos;aCnt++){
				//calc sample
				memset(ftrVl,0x00,sizeof(double)*(NUM_OF_FEATURE-2));
				for(int aa=0;aa<20;aa++){
					for(int ff=0;ff<NUM_OF_FEATURE-2;ff++){
						ftrVl[ff] += btl.feature[aa][ff]*ele[aCnt].aa[aa];
					}
				}
				//print sample
				//ƒAƒ~ƒmŽ_
				for(int aa=0;aa<20;aa++){
					if(sss->AAlabel[aa]==1){
						(*outFile)<<10+(aCnt-startPos)*100+aa<<":"<<ele[aCnt].aa[aa]<<" ";

					}
				}
				//1000 -- 
				for(int ff=0;ff<NUM_OF_FEATURE-2;ff++){
					if(sss->ftr.ftr[ff]==1){
						(*outFile)<<50+(aCnt-startPos)*100+ff<<":";
						(*outFile)<<ftrVl[ff]<<" ";
					}
				}
			}
			if(sss->complexity==1){
				memset(complex,0x00,sizeof(int)*21);
				int complexity=0;
				for(aCnt=startPos;aCnt<stopPos;aCnt++){
					complex[ele[aCnt].acid]++;
				}
				for(aCnt=0;aCnt<20;aCnt++){
					if(complex[aCnt]>0){
						complexity++;
					}
				}
				(*outFile)<<"1900:"<<complexity<<" ";
			}
			if(sss->hydropathy==1){
				hydropathy=0;
				for(aCnt=startPos;aCnt<stopPos;aCnt++){
					hydropathy += btl.hydro[ele[aCnt].acid];
				}
				(*outFile)<<"2000:"<<hydropathy<<" ";
			}
#ifdef READ_MTXDATA
			if(sss->pairpotential==1){
				(*outFile)<<"2100:"<<pd->GetLectSum(startPos,stopPos-1)<<" ";
			}
#endif
			(*outFile)<<"\n";
		}
	}
}

void SVMData::BlncCasp7SVMSmpl(){
	//File Open
	ifstream inFile[CIBCB_POS_CLASS*NUM_OF_SND_STR];
	ofstream outFile[CIBCB_POS_CLASS*NUM_OF_SND_STR];

	char buf[256];

	for(int ps=0;ps<CIBCB_POS_CLASS;ps++){
		for(int strID=0;strID<NUM_OF_SND_STR;strID++){
			sprintf(buf,".\\CASP7\\LearnOrg%d-%d.txt\0",ps,strID);
			string iFileName=buf;
			inFile[NUM_OF_SND_STR*ps+strID].open(iFileName.c_str());
			if(!inFile[NUM_OF_SND_STR*ps+strID]){
				exit(1);
			}
			sprintf(buf,".\\CASP7\\Learn%d-%d.txt\0",ps,strID);
			string oFileName=buf;
			outFile[NUM_OF_SND_STR*ps+strID].open(oFileName.c_str());
			if(!outFile[NUM_OF_SND_STR*ps+strID]){
				exit(1);
			}
		}
	}

	char tmp[20000];
	int ord,dis,rCnt,lCnt;
	double times;
	for(int ps=0;ps<CIBCB_POS_CLASS;ps++){
		cout<<"---"<<ps<<"---"<<"\n";
		for(int strID=0;strID<NUM_OF_SND_STR;strID++){
			switch(ps){
			case 0:
				times = 2;
				break;
			case 1:
				times = 2;
				break;
			case 2:
				times = 2;
				break;
			case 3:
				times = 3;
				break;
			case 4:
				times = 2;
				break;
			case 5:
				times = 2;
				break;
			case 6:
				times = 2;
				break;
			}

			ord=0;dis=0;
			while(!inFile[NUM_OF_SND_STR*ps+strID].eof()){
				inFile[NUM_OF_SND_STR*ps+strID].getline(tmp,20000,'\n');
				if(tmp[0]=='0'){
					dis++;
				}else if(tmp[0]=='1'){
					ord++;
				}
			}
			inFile[NUM_OF_SND_STR*ps+strID].clear();
			inFile[NUM_OF_SND_STR*ps+strID].seekg(0,ios::beg);

			double bairitsu = (double)ord/times/dis;
			cout<<dis<<"/"<<ord<<":"<<bairitsu<<"\n";

			int ordFrom=10;
			int ordSlct=10;
			
			int disFrom=10;
			int disSlct=10;

			if(bairitsu>1){
				ordFrom = bairitsu*10;
				ordSlct = 10;
			}else{
				disFrom = 10.0/bairitsu;
				disSlct = 10;
			}
//			cerr<<disFrom<<","<<disSlct<<"\n";
			while(1){
				if((ordSlct*ord/ordFrom + disSlct*dis/disFrom)>30000){
					ordFrom += 1;
					disFrom += 1;
				}else{
					break;
				}
			}

			int *ordLbl;
			int *disLbl;

			ordLbl = (int *)malloc(sizeof(int)*ordFrom);
			disLbl = (int *)malloc(sizeof(int)*disFrom);

			int tCnt=0;	
			Saiki(ordLbl,100,tCnt,0,ordFrom);
			tCnt=0;			
			Saiki(disLbl,100,tCnt,0,disFrom);

			int pCnt=0;
			int tol=0;
			while(tol<disSlct){
				tol += pow((double)2,(double)pCnt);
				pCnt++;
			}
			tol -= pow((double)2,(double)(pCnt-1));
			pCnt--;
			for(int test=0;test<disFrom;test++){
				if(disLbl[test]<pCnt){
					disLbl[test]=1;
				}else if(disLbl[test]==pCnt && disSlct-tol>0){
					disLbl[test]=1;		
					tol++;
				}else{
					disLbl[test]=0;
				}

			}

			pCnt=0;
			tol=0;
			while(tol<ordSlct){
				tol += pow((double)2,(double)pCnt);
				pCnt++;
			}
			tol -= pow((double)2,(double)(pCnt-1));
			pCnt--;
			for(int test=0;test<ordFrom;test++){
				if(ordLbl[test]<pCnt){
					ordLbl[test]=1;
				}else if(ordLbl[test]==pCnt && ordSlct-tol>0){
					ordLbl[test]=1;		
					tol++;
				}else{
					ordLbl[test]=0;
				}

			}

			rCnt=0;	ord=0;dis=0; lCnt=0;
			while(!inFile[NUM_OF_SND_STR*ps+strID].eof()){
				inFile[NUM_OF_SND_STR*ps+strID].getline(tmp,20000,'\n');
				if(tmp[0]=='0'){
					if(disLbl[lCnt%disFrom]==1){
						outFile[NUM_OF_SND_STR*ps+strID]<<tmp<<"\n";
						dis++;
					}
					lCnt++;
				}else if(tmp[0]=='1'){
					if(ordLbl[rCnt%ordFrom]==1 ){
						outFile[NUM_OF_SND_STR*ps+strID]<<tmp<<"\n";
						ord++;
					}
					rCnt++;	
				}
			}
			cout<<dis<<"/"<<ord<<"\n";
END_FILE_WRITE:;
		}
	}

	for(int ps=0;ps<CIBCB_POS_CLASS;ps++){
		for(int strID=0;strID<NUM_OF_SND_STR;strID++){
			inFile[NUM_OF_SND_STR*ps+strID].close();
			outFile[NUM_OF_SND_STR*ps+strID].close();
		}
	}
}

int SVMData::Saiki(int *seq,int jyouken,int cnt, int start,int stop){
	if(start==stop || start>stop){
		return(0);
	}
	seq[start+(stop-start)/2]=cnt;
//	cerr<<start<<"-"<<stop<<","<<start+(stop-start)/2<<"|";
	cnt++;
	if(cnt>jyouken){
		return(0);
	}
	Saiki(seq,jyouken,cnt,start,start+(stop-start)/2);
	Saiki(seq,jyouken,cnt,start+(stop-start)/2+1,stop);
}
