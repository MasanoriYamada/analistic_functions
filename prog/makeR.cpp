/*2点相関関数(confを抜いたやつ)とBS波動関数(confを抜いたやつ)を読み込んでR相関関数を作る*/

#include<complex>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "../include/analys.h"

using namespace std;
typedef std::complex<double> COMPLEX;
#define I       std::complex<double>(0.0,1.0)

static const int datasize=XYZnodeSites;



int call_data(int,int,COMPLEX[],COMPLEX[]);
int call_prop(char[],COMPLEX[]);
int call_omegaprop(char[],COMPLEX[]);
void make_R(COMPLEX[],COMPLEX[],COMPLEX[]);
void out_data(int,int,COMPLEX[]);
int out_R(char[],COMPLEX[]);



main(int argc, char* argv[]){

	dir_path=argv[1];
	cout <<"Directory path  ::"<<dir_path<<endl;
        in_dir_path = dir_path;
	out_dir_path = dir_path;
	
	root_mkdir(out_dir_path.c_str());
        out_dir_path = out_dir_path + "/Rcor";
	root_mkdir(out_dir_path.c_str());
        out_dir_path = out_dir_path + "/binR";
	root_mkdir(out_dir_path.c_str());
        out_dir_path = out_dir_path + "/xyz";
	root_mkdir(out_dir_path.c_str());
	
	
	for ( int j=0; j<binnumber; j++) {
	  cout <<"conf"<<j<<"を計算中・・・"<<endl;
	  for (int it=T_in; it<T_fi+1; it++) {
	    COMPLEX* omega_sub_prop = new COMPLEX[1];
	    COMPLEX* ave_sub_sub = new COMPLEX[datasize];
		call_data(it,j,&(ave_sub_sub[0]),&(omega_sub_prop[0]));
		COMPLEX* R_sub=new COMPLEX[datasize];
		make_R(&(ave_sub_sub[0]),&(omega_sub_prop[0]),&(R_sub[0]));
			out_data(it,j,&(R_sub[0]));
		delete []ave_sub_sub;
		delete []omega_sub_prop;		
		delete []R_sub;
	}}
		cout << "finished"<<endl;
	return 0;
}

//メインの計算の所、時間と空間をいれてomegaの相関関数を計算する
/************************************************************************************************************************************************/
//call_q_propを使いquarkプロパゲータを呼び出す。call_dataでcall_q_propの適用先をconf i番目にする no (fanameをいじると読み込むファイル,データを入れる配列[Maxsize])  //
/************************************************************************************************************************************************/
int call_data(int it,int b,COMPLEX local[datasize],COMPLEX local1[1]){
	char fnameBS[200]={0};
	char fname2point[200]={0};
	sprintf(fnameBS,"%s/Projwave/binProjwave/xyz/binBS.%s.%06d-%06d.it%03d",in_dir_path.c_str(),base,binnumber,b,it);
	sprintf(fname2point,"%s/tcor2pt/bintcor2pt/%02d/tcor2pt.%s.%06d-%06d.it%02d",in_dir_path.c_str(),it,base,binnumber,b,it);
	call_prop(&(fnameBS[0]),&(local[0]));
	call_omegaprop(&(fname2point[0]),&(local1[0]));
	return 0;
}
/******************************************************************************/
//ファイルからぬいた4点データを呼び出す	no (読み込むファイルパス,読みだしたout用の配列(double)	  //
/******************************************************************************/
int call_prop(char fname[200],COMPLEX data[datasize]){
	int switch1=0;
	if (switch1==0) {
		fstream infile;
		
		infile.open(fname,ios::in|ios::binary);
		if (!infile.is_open()) {
            cout << "ERROR file can't open (no exist)	::"<<fname<<endl;
            exit(1);
            return EXIT_FAILURE;
        }
        if (infile.fail()){
            cout << "ERROR file size is 0 (can open)	::"<<fname<<endl;
            exit(1);
            return EXIT_FAILURE;
        }
		int id=0;
		while(!infile.eof()){
			infile.read( ( char * ) &data[id], sizeof( COMPLEX ) );
			id=id+1; 
			//cout << id<<endl;
		}
		static int tmp=0;
		if (tmp==0) {
			cout <<"reading data size is	;;"<<id<<endl;
			tmp=tmp+1;
		}
		//endian_convert((double*)data,datasize*2);
		for (int point=0; point<id; point++) {
			//cout << data[point]<<endl;
		}}
if (switch1==1) {
	std::ifstream ifs( fname , ios::in );			/* テキストモードで */
		if(! ifs  )
		{
			cout << fname<<"ファイルが開きません"<<endl;
			return 1;
		}
		//		getline( ifs,ss );
		int tmpx[datasize];
		int tmpy[datasize];
		int tmpz[datasize];		
		for(int	id=0; id<datasize; ++id)  {			
			ifs>>std::setw(3)>>tmpx[id]>>std::setw(3)>>tmpy[id]>>std::setw(3)>>tmpz[id]>>setprecision(15)>> data[id].real()>>setprecision(15)>>data[id].imag();
			//cout<<std::setw(3)<<tmpx[id]<<std::setw(3)<<tmpy[id]<<std::setw(3)<<tmpz[id]<<setprecision(15)<< local[id].real()<<setprecision(15)<<local[id].imag()<<endl;
		}
}
	return 0;
}
/******************************************************************************/
//ファイルから抜いた２点データを呼び出す	no (読み込むファイルパス,読みだしたout用の配列(double)	  //
/******************************************************************************/
int call_omegaprop(char fname[200],COMPLEX data[1]){
	int switch1=0;
	if (switch1==0) {
		fstream infile;
		
		infile.open(fname,ios::in|ios::binary);
		if (!infile.is_open()) {
		  cout << "ERROR file can't open (no exist)   ::"<<fname<<endl;
			return EXIT_FAILURE;
		}		
		int id=0;
		while(!infile.eof()){
			infile.read( ( char * ) &data[id], sizeof( COMPLEX ) );
			id=id+1; 
			//cout << id<<endl;
		}
		static int tmp=0;
		if (tmp==0) {
			cout <<"reading data size is	;;"<<id<<endl;
			tmp=tmp+1;
		}
		//endian_convert((double*)data,datasize*2);
		for (int point=0; point<id; point++) {
			//cout << data[point]<<endl;
		}
		infile.close();
	}
if (switch1==1) {
	std::ifstream ifs( fname , ios::in );			/* テキストモードで */
	if(! ifs  )
	{
		cout << fname<<"ファイルが開きません"<<endl;
		return 1;
	}
	//		getline( ifs,ss );
	int tmpx[1];
	int tmpy[1];
	int tmpz[1];		
	for(int	id=0; id<1; ++id)  {
		ifs>>std::setw(15)>> data[id].real()>>std::setw(15)>> data[id].imag();
		//cout << std::setw(15)<<local[id].real()<<std::setw(15)<<local[id].img()<<endl;
		//ifs>>std::setw(3)>>tmpx[id]>>std::setw(3)>>tmpy[id]>>std::setw(3)>>tmpz[id]>>setprecision(15)>> local[id].real()>>setprecision(15)>>local[id].imag();
		//cout<<std::setw(3)<<tmpx[id]<<std::setw(3)<<tmpy[id]<<std::setw(3)<<tmpz[id]<<setprecision(15)<< local[id].real()<<setprecision(15)<<local[id].imag()<<endl;
	}}
	return 0;
}
/************************************************************************************************************************************************/
//out_BSwaveを使いBS波動関数を書き出す。out_dataでout_BSwaveの適用先をconf j番目にする no (fanameをいじると出力ファイル,入力データの配列[Maxsize])  //
/************************************************************************************************************************************************/
void out_data(int it,int j,COMPLEX local[datasize]){
  char fname[300]={0};
  sprintf(fname,"%s/binRwave.%s.%06d-%06d.it%03d",out_dir_path.c_str(),base,binnumber,j,it);
		out_R(&(fname[0]),&(local[0]));
}

/**************************************************************************/
//結果を書き出す	no (書きだすファイルパス,BS波動関数)	  //
/**************************************************************************/
int out_R(char fname[300],COMPLEX  R_sub[datasize]){
	ofstream ofs(fname,ios::binary|ios::trunc);
	if (!ofs.is_open()) {
		cout << "ERROR output file can't open (can't create)"<<endl;
		return EXIT_FAILURE;
	}
	for (int z=0; z<ZnodeSites; z++) {
		for (int y=0; y<YnodeSites; y++) {
			for (int x=0; x<XnodeSites; x++) {
				ofs.write((const char*) &(R_sub[(x) +XnodeSites*((y) + YnodeSites*((z)))]),sizeof( COMPLEX ));
				//cout << x<<"	"<<y<<"	"<<z<<"	"<<R_sub[(x) +XnodeSites*((y) + YnodeSites*((z)))]<<endl;
			}
		}
	}
	ofs.close();
	
	/*
	std::ofstream ofs( &(fname[0]));
	for (int z=0; z<ZnodeSites; z++) {
	for (int y=0; y<YnodeSites; y++) {
	for (int x=0; x<XnodeSites; x++) {
				ofs<<std::setw(3)<<x<<std::setw(3)<<y<<std::setw(3)<<z<<setprecision(15)<< R_sub[(x) +XnodeSites*((y) + YnodeSites*((z)))].real()<<setprecision(15)<< R_sub[(x) +XnodeSites*((y) + YnodeSites*((z)))].imag()<<endl;
	}}}*/
		return 0;
	
}
/******************************************************************************/
//R_correlatorを作る。 BS/(2point*2point)	  //
/******************************************************************************/
void make_R(COMPLEX jack_ave_sub_sub[datasize],COMPLEX omega_sub_prop[1],COMPLEX R_sub[datasize]){
for (int id=0; id<datasize; id++) {
	R_sub[id]=jack_ave_sub_sub[id]/((omega_sub_prop[0])*(omega_sub_prop[0]));
	//cout << "jack_ave"<<jack_ave_sub_sub[id]<<endl;
	//cout << "omegaprop"<<omega_sub_prop[0]<<endl;
	//cout << "R_sub"<<R_sub[id]<<endl;
}
}










