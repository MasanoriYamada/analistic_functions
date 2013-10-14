#include "../include/analys.h"


static const int datasize=TnodeSites;									//max size of it direction array

#define I       std::complex<double>(0.0,1.0)
using namespace std;
typedef std::complex<double> COMPLEX;


int call_data(int,int,COMPLEX[]);
int call_jack(char[],COMPLEX[]);
void out_file(int,char[],COMPLEX[],COMPLEX[]);
void out_data(int,COMPLEX[],COMPLEX[]);
COMPLEX jack_ave_calc(COMPLEX[]);
COMPLEX jack_err_calc(COMPLEX[]);



int main(int argc,char** argv){

  //make out put dir
  dir_path=argv[1];
  cout <<"Directory path  ::"<<dir_path<<endl;
  in_dir_path = dir_path;
  out_dir_path = dir_path;
  root_mkdir(out_dir_path.c_str());
  out_dir_path = out_dir_path + "/tcor2pt";
  root_mkdir(out_dir_path.c_str());
  out_dir_path = out_dir_path + "/jack_error";
  root_mkdir(out_dir_path.c_str());


  COMPLEX* jack_ave_sub = new COMPLEX[datasize*binnumber];
  for (int it=T_in; it < (T_fi +1); it++) {
    cout << "時間"<<it<<"を計算中・・・"<<endl;
    for (int b=0; b<binnumber; b++) {
      call_data(b,it,&(jack_ave_sub[0]));
    }}
  COMPLEX* avesub = new COMPLEX[datasize];
  COMPLEX* errsub = new COMPLEX[datasize];
  COMPLEX jack_ave_sub_sub[binnumber];
  for (int it=T_in; it < (T_fi +1); it++) {
    COMPLEX ave=0.0;
    COMPLEX err=0.0;
    for (int b=0; b<binnumber; b++) {				
      jack_ave_sub_sub[b]=jack_ave_sub[it + datasize*b];
    }
    ave=jack_ave_calc(&(jack_ave_sub_sub[0]));
    err=jack_err_calc(&(jack_ave_sub_sub[0]));
    cout << ave<<endl;
    cout << err<<endl;
    avesub[it]=ave;
    errsub[it]=err;
  }
  for (int it=T_in; it < (T_fi +1); it++) {
    out_data(it,&(avesub[0]),&(errsub[0]));
  }
  delete[] avesub;
  delete[] errsub;
  delete[] jack_ave_sub;
  cout << "finished"<<endl;
  return 0;
	
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/************************************************************************************************************************************************/
//call_q_propを使いquarkプロパゲータを呼び出す。call_dataでcall_q_propの適用先をconf i番目にする no (fanameをいじると読み込むファイル,データを入れる配列[Maxsize])  //
/************************************************************************************************************************************************/
int call_data(int j,int it,COMPLEX local[datasize*binnumber]){
	char fname[300]={0};
	COMPLEX data[1];
	sprintf(fname,"%s/tcor2pt/bintcor2pt/%02d/tcor2pt.%s.%06d-%06d.it%02d",in_dir_path.c_str(),it,base,binnumber,j,it);
	call_jack(fname,data);
	//cout << fname<<endl;
	//cout << data[0]<<endl;
	local[it+ datasize*j]=data[0];
	//cout << it<<j<<local[it+datasize*j]<<endl;
	
	return 0;
}
/******************************************************************************/
//ファイルからデータを呼び出す	no (読み込むファイルパス,読みだしたout用の配列(double)	  //
/******************************************************************************/
int call_jack(char fname[200],COMPLEX data[1]){

	fstream infile;
	//COMPLEX data[datasize];
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
		infile.read( ( char * ) &(data[id]), sizeof( COMPLEX ) );
		id=id+1; 
		//cout << id<<endl;
	}
	static int tmp=0;
	if (tmp==0) {
		cout <<"reading data size is	;;"<<id<<endl;
		tmp=tmp+1;
	}
	//endian_convert((double*)data,datasize*2);
	//cout<<data[0]<<endl;
	infile.close();
	return 0;
}	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*******************************************************************************************************************************************/
//call_fileを利用して、複数のファイルから呼び出したデータを配列としてまとめる (読みだしたout用の配列(double)	この中のfanameをいじると入力ファイルを変えられる。  //
/******************************************************************************************************************************************/
void out_data(int it,COMPLEX avesub[datasize],COMPLEX errsub[datasize]){
	
	char fname[200]={0};
	sprintf(fname,"%s/tcor2pt.%s.%06d",out_dir_path.c_str(),base,binnumber);
	out_file(it,fname,&(avesub[0]),&(errsub[0]));
}

/**************************************************************************/
//結果を書き出す	no (書きだすファイルパス,BS波動関数)	  //(rについて書くファイル名,xについて書くファイル名,書きだす内容)
/**************************************************************************/
void out_file(int it,char fname[200],COMPLEX avesub[datasize], COMPLEX errsub[datasize]){
	std::ofstream ofs;
	
	ofs.open(&(fname[0]));
	for (int it=0; it<datasize; it++) {
	  ofs<<it<<" "<<avesub[it].real()<<"	"<<errsub[it].real()<<endl;
	}
	ofs.close();
	}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/********************************************************************************************************/
//jackナイフのためにconf[j]を取り除き平均をとったomega_prop[j]をだす no(入力のデータ,出力のデータ)
/*******************************************************************************************************/
COMPLEX jack_ave_calc(COMPLEX jack_ave_sub[binnumber]){
	COMPLEX buffer=0.0;
	COMPLEX ave=0.0;
	for (int b=0; b<binnumber; b++) {
		buffer=buffer+jack_ave_sub[b];
	}
	ave=buffer/(double)binnumber;
	return ave;
}

/*******************************************************************************************************/
//jackナイフのためにconf[j]を取り除き平均をとったomega_prop[j]をだす no(入力のデータ,出力のデータ)
/*******************************************************************************************************/
COMPLEX jack_err_calc(COMPLEX jack_ave_sub[binnumber]){
	COMPLEX ave1=0;
	COMPLEX ave2=0;
	COMPLEX err=0;
	COMPLEX double_jack_ave_sub[binnumber];
	for (int b=0; b<binnumber; b++) {
		double_jack_ave_sub[b]=jack_ave_sub[b]*jack_ave_sub[b];
	}
	ave1=jack_ave_calc(&(jack_ave_sub[0]));
	ave2=jack_ave_calc(&(double_jack_ave_sub[0]));
	err= sqrt(((double)binnumber -1.0)*((ave2)-(ave1)*(ave1)));
	return err;
}

