
#include "../include/analys.h"


static const int datasize=TnodeSites;									//number of data in a file
static const int switch1=1;                          //0;bainary	1;text

#define data(id,j) data[j+(Confsize*id)]					//id:number ofdata		
#define I       std::complex<double>(0.0,1.0)
using namespace std;
typedef std::complex<double> COMPLEX;



int call_file(char[],COMPLEX[]);
void call_data(COMPLEX[]);
int out_file(int,int,char[],COMPLEX[]);
void out_data(int,COMPLEX[]);
void jack_avesub_calc(int,COMPLEX[],COMPLEX[]);



int main(int argc,char** argv){

  //make out put dir
	dir_path=argv[1];
	cout <<"Directory path  ::"<<dir_path<<endl;
        in_dir_path = dir_path;
	out_dir_path = dir_path;
	root_mkdir(out_dir_path.c_str());
        out_dir_path = out_dir_path + "/tcor2pt";
	root_mkdir(out_dir_path.c_str());
        out_dir_path = out_dir_path + "/bintcor2pt";
	root_mkdir(out_dir_path.c_str());

	COMPLEX *data= new COMPLEX[datasize*Confsize];
	call_data(&(data[0]));
		COMPLEX *jack_ave_sub_sub= new COMPLEX[datasize*binnumber];
		for (int id=0; id < datasize; id++) {
			COMPLEX jack_ave_sub[binnumber];
			jack_avesub_calc(id,&(data[0]),&(jack_ave_sub[0]));
			for (int b=0; b<binnumber; b++) {
				jack_ave_sub_sub[id+datasize*(b)]=jack_ave_sub[b];
			}
			out_data(id,&(jack_ave_sub_sub[0]));
			
		}
		cout << "finished"<<endl;
	return 0;
	
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*******************************************************************************************************************************************/
//call_fileを利用して、複数のファイルから呼び出したデータを配列としてまとめる (読みだしたout用の配列(double)	この中のfanameをいじると入力ファイルを変えられる。  //
/******************************************************************************************************************************************/
void call_data(COMPLEX data[Confsize*datasize]){
  in_dir_path = in_dir_path + "/tcor2pt";
  COMPLEX *buffer = new COMPLEX[Confsize*datasize];
  for (int j=0; j<Confsize; j++) {
    for (int tmp=0; tmp<Tshiftsize; tmp++) {

      char fname[200]={0};
      COMPLEX local[datasize]={0};
      /* set1
	 if(0<=j && j<500){
	 sprintf(fname,"%s/ts%02d/omega_correlator.+%03d.%s-1-%05d0",in_dir_path.c_str(),tshift[tmp],tshift[tmp],base,j+41);//hachi
	 }
	 else if(500<=j && j<600){
	 sprintf(fname,"%s/ts%02d/omega_correlator.+%03d.%s-2-%05d0",in_dir_path.c_str(),tshift[tmp],tshift[tmp],base,j+51-500);//hachi
	 }
	 else if(600<=j && j<700){
	 sprintf(fname,"%s/ts%02d/omega_correlator.+%03d.%s-3-%05d0",in_dir_path.c_str(),tshift[tmp],tshift[tmp],base,j+41-600);//hachi
	 }
	 else {
	 cout << "ERR open file name is worng"<<endl;
	 }
	 set2 */
      if(0<=j && j<550){
	sprintf(fname,"%s/raw/omega_correlator.+%03d.%s-1-%05d0",in_dir_path.c_str(),tshift[tmp],base,j+41);//hachi
      }
      else if(550<=j && j<800){
	sprintf(fname,"%s/raw/omega_correlator.+%03d.%s-2-%05d0",in_dir_path.c_str(),tshift[tmp],base,j-499);//hachi
      }
      else {
	cout << "ERR open file name is worng"<<endl;
      }

            
      call_file(&(fname[0]),&(local[0]));
      for (int id=0; id < datasize; id++) {
	buffer[Confsize*id+j] +=local[id];
      }}
    for (int id=0; id < datasize; id++) {
      data(id,j)=buffer[Confsize*id+j]/((double)Tshiftsize);
    }}
  //debug
  cout<<"data[0]"<<data(0,0)<<"buffer[0]"<<buffer[0]<<" Tshift  "<<Tshiftsize<<endl;
  //debug end

}
/*******************************************************************************************************/
//ファイルからデータを呼び出す	no (読み込むファイルパス,読みだしたout用の配列(double)	  //
/******************************************************************************************************/
int call_file(char fname[200],COMPLEX proj_wave[datasize]){

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
		//while(!infile.eof()){
		for (int id=0; id < datasize; id++) {
			infile.read( ( char * ) &proj_wave[id], sizeof( COMPLEX ) );
			//cout << id<<endl;
		}
		/*
		static int tmp=0;
		if (tmp==0) {
			cout <<"reading data size is	;;"<<id<<endl;
			tmp=tmp+1;
		}
		 */
		//endian_convert((double*)proj_wave,datasize*2);
		//for (int point=0; point<id; point++) {
			//cout << proj_wave[point]<<endl;
		}//}
	//-----------------------------------------------------------------------------------------------------------------------------------------------------------
	if (switch1==1) {std::ifstream ifs( fname , ios::in );			/* テキストモードで */
		if(! ifs  )
		{
			cout << fname<<"ファイルが開きません"<<endl;
			return 1;
		}
		//		getline( ifs,ss );
		int tmptmp=0;
		for(int	id=0; id<datasize; ++id)  {
			ifs>>tmptmp>>proj_wave[id].real()>>proj_wave[id].imag();
			//cout<< proj_wave[id].real()<<"	"<<proj_wave[id].imag()<<endl;
		}}
	return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*******************************************************************************************************************************************/
//call_fileを利用して、複数のファイルから呼び出したデータを配列としてまとめる (読みだしたout用の配列(double)	この中のfanameをいじると入力ファイルを変えられる。  //
/******************************************************************************************************************************************/
void out_data(int it,COMPLEX jack_ave_sub_sub[datasize*binnumber]){
    char fname[200]={0};

    sprintf(fname,"%s/%02d/",out_dir_path.c_str(),it);//hachi
    root_mkdir(fname);

	for (int b=0; b<binnumber; b++) {
		char fname1[200]={0};
		sprintf(fname1,"%s/%02d/tcor2pt.%s.%06d-%06d.it%02d",out_dir_path.c_str(),it,base,binnumber,b,it);
		out_file(it,b,fname1,&(jack_ave_sub_sub[0]));
	}
}
/*******************************************************************************************************/
//結果を書き出す	no (書きだすファイルパス,平均(データの大きさ),誤差(データの大きさ)	  //
/******************************************************************************************************/
int out_file(int it,int b, char fname[200],COMPLEX jack_ave_sub_sub[binnumber*datasize]){
	ofstream ofs;
	//ofs( &(fname[0]),std::ios::out | std::ios::trunc);	
	//ofs<< jack_ave_sub_sub[it].real()<<setprecision(15)<<jack_ave_sub_sub[it].imag()<<endl;
	//cout <<std::setw(3)<<x<<std::setw(3)<<y<<std::setw(3)<<z<<setprecision(15)<< jack_ave_sub_sub[(x) +XnodeSites*((y) + YnodeSites*((z)))].real()<<setprecision(15)<<jack_ave_sub_sub[(x) +XnodeSites*((y) + YnodeSites*((z)))].imag()<<endl;

	
	ofs.open(fname,ios::out|ios::binary|ios::trunc);
	if (!ofs.is_open()) {
		cout << "ERROR output file can't open (no exist)"<<endl;
		return EXIT_FAILURE;
	}
	ofs.write((const char*) &(jack_ave_sub_sub[it+datasize*(b)]),sizeof( COMPLEX ));
	//cout << it <<"	"<<b<<"	"<<jack_ave_sub_sub[it+datasize*(b)]<<endl;
	ofs.close();
	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*******************************************************************************************************/
//jackナイフのためにconf[j]を取り除き平均をとったomega_prop[j]をだす no(入力のデータ,出力のデータ)
/*******************************************************************************************************/
void jack_avesub_calc(int id,COMPLEX data[datasize*Confsize],COMPLEX jack_ave_sub[binnumber]){
	
	COMPLEX full=0.0;
	for (int j=0; j<Confsize; j++) {
		full=full+data(id,j);
	}
	for (int b=0; b<binnumber; b++) {
		COMPLEX subfull[Confsize];
		COMPLEX localfull=0;
		for (int j=b*binsize; j<(b+1)*binsize; j++) {
			localfull=localfull+data(id,j);
		}
		subfull[b]=full-localfull;
		jack_ave_sub[b]=subfull[b]/((double)Confsize-(double)binsize);
	}
}

