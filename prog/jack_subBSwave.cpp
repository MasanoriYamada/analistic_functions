#include "../include/analys.h"

static const int datasize=XnodeSites*YnodeSites*ZnodeSites;
static const int switch1=0;                          //0:bainary read　1:text read
#define data(id,j) data[Confsize*id+j]
#define I       std::complex<double>(0.0,1.0)
using namespace std;
typedef std::complex<double> COMPLEX;



int call_file(char[],COMPLEX[]);
void call_data(int,COMPLEX[]);
int out_file(int,char[],COMPLEX[]);
void out_data(int,COMPLEX[]);
void jack_avesub_calc(int,COMPLEX[],COMPLEX[]);
static void endian_convert(double* buf,int len);




int main(int argc , char** argv){
  //make out put dir
	dir_path=argv[1];
	cout <<"Directory path  ::"<<dir_path<<endl;
	in_dir_path = dir_path;
	out_dir_path = dir_path;
	root_mkdir(dir_path.c_str());
    out_dir_path=out_dir_path + "/Projwave";
	root_mkdir(out_dir_path.c_str());
    out_dir_path = out_dir_path + "/binProjwave";
	root_mkdir(out_dir_path.c_str());
    out_dir_path = out_dir_path + "/xyz";
	root_mkdir(out_dir_path.c_str());

	for (int it=T_in; it < (T_fi +1); it++) {
		cout << "time"<<it<<endl;
		COMPLEX *data= new COMPLEX[datasize*Confsize];
		COMPLEX *jack_ave_sub_sub= new COMPLEX[datasize*binnumber];
		call_data(it,&(data[0]));
	for (int id=0; id < datasize; id++) {
			COMPLEX jack_ave_sub[binnumber];
			jack_avesub_calc(id,&(data[0]),&(jack_ave_sub[0]));
			for (int b=0; b<binnumber; b++) {
				jack_ave_sub_sub[id+datasize*(b)]=jack_ave_sub[b];
			}}
		out_data(it,&(jack_ave_sub_sub[0]));
	}	
	cout << "finish"<<endl;
	return 0;
	
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*******************************************************************************************************************************************/
//call_fileを利用して、複数のファイルから呼び出したデータを配列としてまとめる (読みだしたout用の配列(double)	この中のfanameをいじると入力ファイルを変えられる。  //
/******************************************************************************************************************************************/
void call_data(int it,COMPLEX data[Confsize*datasize]){
	for (int j=0; j<Confsize; j++) {
 	    char fname[200]={0};
	    COMPLEX local[datasize]={0};
	    sprintf(fname,"%s/impBSwave/aveTshift/OmgOmgwave_PH1.+%03d.%s-%06d",in_dir_path.c_str(),it,base,j);
		cout << fname<<"reading now"<<endl;
		call_file(&(fname[0]),&(local[0]));
		for (int id=0; id < datasize; id++) {
			data(id,j)=local[id];
		}
	}
}
/*******************************************************************************************************/
//ファイルからデータを呼び出す	no (読み込むファイルパス,読みだしたout用の配列(double)	  //
/******************************************************************************************************/
int call_file(char fname[300],COMPLEX data[datasize]){
	//string ss;
	int tmpx[datasize];
	int tmpy[datasize];
	int tmpz[datasize];
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
	//-----------------------------------------------------------------------------------------------------------------------------------------------------------
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
			ifs>>std::setw(3)>>tmpx[id]>>std::setw(3)>>tmpy[id]>>std::setw(3)>>tmpz[id]>>std::setw(15)>> data[id].real()>>std::setw(15)>>data[id].imag();
		 //	cout<<std::setw(3)<<tmpx[id]<<std::setw(3)<<tmpy[id]<<std::setw(3)<<tmpz[id]<<std::setw(15)<< data[id].real()<<std::setw(15)<<data[id].imag()<<endl;
		}}
	return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*******************************************************************************************************************************************/
//call_fileを利用して、複数のファイルから呼び出したデータを配列としてまとめる (読みだしたout用の配列(double)	この中のfanameをいじると入力ファイルを変えられる。  //
/******************************************************************************************************************************************/
void out_data(int it,COMPLEX jack_ave_sub_sub[datasize*binnumber]){
	for (int b=0; b<binnumber; b++) {
		char fname[200]={0};
	sprintf(fname,"%s/binBS.%s.%06d-%06d.it%03d",out_dir_path.c_str(),base,binnumber,b,it);
			out_file(b,fname,&(jack_ave_sub_sub[0]));
	}
}
/*******************************************************************************************************/
//結果を書き出す	no (書きだすファイルパス,平均(データの大きさ),誤差(データの大きさ)	  //
/******************************************************************************************************/
int out_file(int b, char fname[200],COMPLEX jack_ave_sub_sub[binnumber*datasize]){
	COMPLEX outdata[datasize];
	if (switch1==0) {
		ofstream ofs_xyz;		
		ofs_xyz.open(fname,ios::out|ios::binary|ios::trunc);
		if (!ofs_xyz.is_open()) {
			cout << "ERROR output file can't open (no exist)"<<endl;
            exit(1);
			return EXIT_FAILURE;
		}
		for (int z=0; z<ZnodeSites; z++) {
		for (int y=0; y<YnodeSites; y++) {
		for (int x=0; x<XnodeSites; x++) {
			int	id=(x) +XnodeSites*((y) + YnodeSites*((z)));
					outdata[id]=jack_ave_sub_sub[id+datasize*(b)];
					ofs_xyz.write((const char*) &(outdata[id]),sizeof( COMPLEX ));
		}}}
	}
	if (switch1==1) {
	std::ofstream ofs( &(fname[0]),std::ios::out | std::ios::trunc);
	for (int z=0; z<ZnodeSites; z++) {
	for (int y=0; y<YnodeSites; y++) {
	for (int x=0; x<XnodeSites; x++) {
	int	id=(x) +XnodeSites*((y) + YnodeSites*((z)));
		outdata[id]=jack_ave_sub_sub[id+datasize*(b)];
		ofs<<std::setw(3)<<x<<std::setw(3)<<y<<std::setw(3)<<z<<std::setw(21)<<setprecision(15)<< outdata[id].real()<<std::setw(21)<<setprecision(15)<<outdata[id].imag()<<endl;
		//cout<<std::setw(3)<<x<<std::setw(3)<<y<<std::setw(3)<<z<<std::setw(21) <<setprecision(15)<< outdata[id].real()<<std::setw(21) <<setprecision(15)<<outdata[id].imag()<<endl;
	}}}}
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
/*******************************************************************/
//endian convert
/*******************************************************************/

static void endian_convert(double* buf,int len)
{
	for(int i=0; i<len; i++){
		char tmp[8];
		((double*)tmp)[0] = *buf;
		for(int j=0; j<8; j++){
			((char*)buf)[j] = ((char*)tmp)[7-j];
		}
		buf++;
	}
}


