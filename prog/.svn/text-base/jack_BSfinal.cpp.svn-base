#include "../include/analys.h"

static const int datasize=XYZnodeSites;									//１つのファイルにあるデータの数
#define data(id,j) data[binnumber*id+j]					//id:dataの数		
#define I       std::complex<double>(0.0,1.0)
using namespace std;
typedef std::complex<double> COMPLEX;


int call_file(int,int,double[]);
int call_data(int,char[],double[]);
void out_file(char[],char[],double[],double[]);
void out_data(int,double[],double[]);
double jack_ave_calc(int,double[]);
double jack_err_calc(int,double[]);
static void endian_convert(double* buf,int len);



int main(int argc,char** argv){

  //make out put dir
	dir_path=argv[1];
	cout <<"Directory path  ::"<<dir_path<<endl;
        in_dir_path = dir_path;
	out_dir_path = dir_path;
	root_mkdir(out_dir_path.c_str());
        out_dir_path = out_dir_path + "/Projwave";
	root_mkdir(out_dir_path.c_str());
        out_dir_path = out_dir_path + "/jack_error";
	root_mkdir(out_dir_path.c_str());
        xyz_out_dir_path = out_dir_path + "/xyz";
	root_mkdir(xyz_out_dir_path.c_str());
        r_out_dir_path = out_dir_path + "/r";
	root_mkdir(r_out_dir_path.c_str());


	for (int it=T_in; it < (T_fi +1); it++) {
		cout << "時間"<<it<<"を読み込んでいます・・・"<<endl;
		double* jack_ave_sub=new double[datasize*binnumber];
	for (int b=0; b<binnumber; b++) {
	call_file(b,it,&(jack_ave_sub[0]));
	}
		double avesub[datasize]={0};
		double errsub[datasize]={0};
	for (int id=0; id < datasize; id++) {
		double ave=0;
		double err=0;
	ave=jack_ave_calc(id,&(jack_ave_sub[0]));
	err=jack_err_calc(id,&(jack_ave_sub[0]));
		avesub[id]=ave;
		errsub[id]=err;
	}
	out_data(it,&(avesub[0]),&(errsub[0]));
		delete []jack_ave_sub;
	}
	cout << "finished"<<endl;
	return 0;

	}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/************************************************************************************************************************************************/
//call_q_propを使いquarkプロパゲータを呼び出す。call_dataでcall_q_propの適用先をconf i番目にする no (fanameをいじると読み込むファイル,データを入れる配列[Maxsize])  //
/************************************************************************************************************************************************/
int call_file(int b,int it,double local[binnumber*datasize]){
	char fname[200]={0};
	sprintf(fname,"%s/Projwave/binProjwave/xyz/binBS.%s.%06d-%06d.it%03d",in_dir_path.c_str(),base,binnumber,b,it);

	call_data(b,&(fname[0]),&(local[0]));
	return 0;
}
/******************************************************************************/
//ファイルからデータを呼び出す	no (読み込むファイルパス,読みだしたout用の配列(double)	  //
/******************************************************************************/
int call_data(int b,char fname[200], double jack_ave_sub[binnumber*datasize]){
	fstream infile;
	COMPLEX data[datasize];
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
	for (int z=0; z<ZnodeSites; z++) {
	for (int y=0; y<YnodeSites; y++) {
	for (int x=0; x<XnodeSites; x++) {
	int	id=(x) +XnodeSites*((y) + YnodeSites*((z)));
		//cout << data[point]<<endl;
		jack_ave_sub[id+datasize*(b)]=data[id].real();
		//cout << x<<y<<z<<data[id]<<endl;
	}}}
	/*
	
	
	COMPLEX proj_BSwave[datasize];
	std::ifstream ifs( fname , ios::in );			// テキストモードで 
	if(! ifs  )
	{
		cout << fname<<"ファイルが開きません"<<endl;
		return 1;
	}
	int tmpx[XYZnodeSites];
	int tmpy[XYZnodeSites];
	int tmpz[XYZnodeSites];
	for (int z=0; z<ZnodeSites; z++) {
	for (int y=0; y<YnodeSites; y++) {
	for (int x=0; x<XnodeSites; x++) {
			ifs>>std::setw(3)>>tmpx[(x) +XnodeSites*((y) + YnodeSites*((z)))]>>std::setw(3)>>tmpy[(x) +XnodeSites*((y) + YnodeSites*((z)))]>>std::setw(3)>>tmpz[(x) +XnodeSites*((y) + YnodeSites*((z)))]>>std::setw(21) >>setprecision(15)>> proj_BSwave[(x) +XnodeSites*((y) + YnodeSites*((z)))].real()>>std::setw(21) >>setprecision(15)>>proj_BSwave[(x) +XnodeSites*((y) + YnodeSites*((z)))].imag();
							//	cout << proj_BSwave[(x) +XnodeSites*((y) + YnodeSites*((z)))].real()<<"///"<<proj_BSwave[(x) +XnodeSites*((y) + YnodeSites*((z)))].imag()<<endl;
		jack_ave_sub[(x) +XnodeSites*((y) + YnodeSites*((z)))+datasize*(b)]=proj_BSwave[(x) +XnodeSites*((y) + YnodeSites*((z)))].real();
	}
	}
	}
	*/
		return 0;
	}	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*******************************************************************************************************************************************/
//call_fileを利用して、複数のファイルから呼び出したデータを配列としてまとめる (読みだしたout用の配列(double)	この中のfanameをいじると入力ファイルを変えられる。  //
/******************************************************************************************************************************************/
void out_data(int it,double avesub[datasize],double errsub[datasize]){
	
	char fnamer[200]={0};
	char fnamexyz[200]={0};
		sprintf(fnamer,"%s/ProjBS.%06d.%s.it%02d",r_out_dir_path.c_str(),binnumber,base,it);
	sprintf(fnamexyz,"%s/ProjBS.%06d.%s.it%02d",xyz_out_dir_path.c_str(),binnumber,base,it);
	out_file(fnamer,fnamexyz,&(avesub[0]),&(errsub[0]));
}

/**************************************************************************/
//結果を書き出す	no (書きだすファイルパス,BS波動関数)	  //(rについて書くファイル名,xについて書くファイル名,書きだす内容)
/**************************************************************************/
void out_file(char fname[200],char fname1[200],double avesub[datasize], double errsub[datasize]){
#define radius_sq(x,y,z) ((x)*(x) + (y)*(y) + (z)*(z))
#define min(a,b) (((a) < (b)) ? (a) : (b))
		int r_sq=0;
		int r_sq_count[r_sq_max]={0};
	double ave[r_sq_max];
	double err[r_sq_max];
	std::ofstream ofsxyz;
	std::ofstream ofsr;
		
	ofsxyz.open(&(fname1[0]));
	ofsr.open( &(fname[0]));

	ofsxyz.setf(ios::scientific);
	ofsr.setf(ios::scientific);

	
		//	ofs<<"#"<<std::setw(7)<<"r"<<std::setw(21) <<setprecision(15)<<"BSwave(r)real part"<<std::setw(21) <<setprecision(15)<<"BSwave(r)imaginary part"<<endl;
		for (int z=0; z<ZnodeSites; z++) {
		for (int y=0; y<YnodeSites; y++) {
		for (int x=0; x<XnodeSites; x++) {
		ofsxyz<<x<<"	"<<y<<"	"<<z<<"	"<< avesub[(x) +XnodeSites*((y) + YnodeSites*((z)))]<<"	"<<errsub[(x) +XnodeSites*((y) + YnodeSites*((z)))]<<endl;
					int r_sq = radius_sq( min(x,XnodeSites-x), min(y,YnodeSites-y), min(z,ZnodeSites-z) );
					r_sq_count[r_sq]= r_sq_count[r_sq] +1;
			ave[r_sq] =ave[r_sq] + avesub[(x) +XnodeSites*((y) + YnodeSites*((z)))];
			err[r_sq] =err[r_sq] + errsub[(x) +XnodeSites*((y) + YnodeSites*((z)))];
			//float rad = sqrt((float)r_sq);
			//ofsr<<rad<<"	"<< avesub[(x) +XnodeSites*((y) + YnodeSites*((z)))]<<"	"<<errsub[(x) +XnodeSites*((y) + YnodeSites*((z)))]<<endl;
			
				}
			}
		}
		ofsxyz.close();
		for (int r_sq=0; r_sq<r_sq_max; r_sq++) {
			
			if ( r_sq_count[r_sq] == 0 ) continue;
			
			ave[r_sq] =ave[r_sq]/ ((double) r_sq_count[r_sq]);
			err[r_sq] =err[r_sq]/ ((double) r_sq_count[r_sq]);
			
			float rad = sqrt((float)r_sq);
			ofsr<<rad<<"	"<< ave[r_sq]<<"	"<<err[r_sq]<<endl;
			//ofsr<<std::setw(7)<<rad<<std::setw(21)<<setprecision(15)<< ave[r_sq]<<std::setw(21)<<setprecision(15)<<err[r_sq]<<endl;
		}
		ofsr.close();
	}
	

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/********************************************************************************************************/
//jackナイフのためにconf[j]を取り除き平均をとったomega_prop[j]をだす no(入力のデータ,出力のデータ)
/*******************************************************************************************************/
double jack_ave_calc(int id,double jack_ave_sub[binnumber*datasize]){
	double buffer=0.0;
	double ave=0.0;
	for (int b=0; b<binnumber; b++) {
			buffer=buffer+jack_ave_sub[id+datasize*(b)];
		}
		ave=buffer/(double)binnumber;
	return ave;
}

/*******************************************************************************************************/
//jackナイフのためにconf[j]を取り除き平均をとったomega_prop[j]をだす no(入力のデータ,出力のデータ)
/*******************************************************************************************************/
double jack_err_calc(int id,double jack_ave_sub[binnumber*datasize]){
	double ave1=0;
	double ave2=0;
	double err=0;
	double* double_jack_ave_sub=new double[binnumber*datasize];
		for (int b=0; b<binnumber; b++) {
		double_jack_ave_sub[id+datasize*(b)]=jack_ave_sub[id+datasize*(b)]*jack_ave_sub[id+datasize*(b)];
		}
		ave1=jack_ave_calc(id,&(jack_ave_sub[0]));
		ave2=jack_ave_calc(id,&(double_jack_ave_sub[0]));
		err= sqrt(((double)binnumber -1.0)*((ave2)-(ave1)*(ave1)));
	delete []double_jack_ave_sub;
	return err;
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



