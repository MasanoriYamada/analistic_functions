/***********************************************/
/*-1*時間1階微分を計算するプログラム*/
/***********************************************/
/*information
 16^3*32
 conf100
 t=0~14
 in ./results/proj.omega-omega1.3.3/xyz/
 out ./results/proj.omega-omega1.3.3/laplacian/xyz/
 out ./results/proj.omega-omega1.3.3/laplacian/r/
 上書き
 version info 0.2 normalizationを取る
 
 読み込みは成功
 差分演算は成功
 最後まで確認した
 */


#include "../include/analys.h"


#define proj_BSwave(x,y,z) proj_BSwave[(x) +XnodeSites*((y) + YnodeSites*((z)))]
#define	sub_BSwave(x,y,z,it)   sub_BSwave[(x) +XnodeSites*((y) + YnodeSites*((z)+ZnodeSites*(it)))]
#define	sub_sub_BSwave(x,y,z,it)   sub_sub_BSwave[(x) +XnodeSites*((y) + YnodeSites*((z)+ZnodeSites*(it)))]
#define	Tdep_BSwave(x,y,z,it) Tdep_BSwave[(x) +XnodeSites*((y) + YnodeSites*((z)+ZnodeSites*(it)))]


#define I       std::complex<double>(0.0,1.0)
using namespace std;
typedef std::complex<double> COMPLEX;


int call_data(int,int,COMPLEX[]);
int call_BSwave(char[],COMPLEX[]);
void out_data(int,int,COMPLEX[]);
int out_TdepBSwave(int,char[],COMPLEX[]);
void calc_Tdep(COMPLEX[],COMPLEX[]);



int main(int argc,char** argv){

	dir_path=argv[1];
	cout <<"Directory path  ::"<<dir_path<<endl;
        in_dir_path = dir_path;
	out_dir_path = dir_path;

	root_mkdir(out_dir_path.c_str());
        out_dir_path = out_dir_path + "/Tder";
	root_mkdir(out_dir_path.c_str());
        out_dir_path = out_dir_path + "/binTder";
	root_mkdir(out_dir_path.c_str());
        out_dir_path = out_dir_path + "/xyz";
	root_mkdir(out_dir_path.c_str());
	in_dir_path = in_dir_path + "/Rcor/binR/xyz";

		
	for ( int j=0; j<binnumber; j++) {
		COMPLEX* sub_BSwave= new COMPLEX[XYZnodeSites*TnodeSites];
		for (int it=T_in; it<T_fi+1; it++) {
				COMPLEX* proj_BSwave = new COMPLEX[XYZnodeSites];
			//cout<<"今conf"<<j<<"、時間"<<it<<"データの読み出し中"<<endl;
			call_data(j,it,&(proj_BSwave[0]));
			//cout<<"今conf"<<j<<"、時間"<<it<<"データの読み出し完了"<<endl;
			for (int z=0; z<ZnodeSites; z++) {
			for (int y=0; y<YnodeSites; y++) {
			for (int x=0; x<XnodeSites; x++) {
			sub_BSwave(x,y,z,it)=proj_BSwave(x,y,z);
				//cout <<"	"<<it<<"	"<<sub_BSwave(x,y,z,it)<<endl;
			}}}
				delete []proj_BSwave;
		}
			COMPLEX Tdep_BSwave[XYZnodeSites*TnodeSites];
			calc_Tdep(&(sub_BSwave[0]),&(Tdep_BSwave[0]));
		for (int it=T_in; it<T_fi+1; it++) {
			out_data(j,it,&(Tdep_BSwave[0]));
		}
		delete []sub_BSwave;	
		cout<<"conf"<<j<<"	finish"<<endl;
	}
	cout << "finished"<<endl;
	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/************************************************************************************************************************************************/
//call_q_propを使いquarkプロパゲータを呼び出す。call_dataでcall_q_propの適用先をconf i番目にする no (fanameをいじると読み込むファイル,データを入れる配列[Maxsize])  //
/************************************************************************************************************************************************/
int call_data(int b,int it,COMPLEX local[XYZnodeSites]){
        char fname[200];
	sprintf(fname,"%s/binRwave.%s.%06d-%06d.it%03d",in_dir_path.c_str(),base,binnumber,b,it);
	call_BSwave(&(fname[0]),&(local[0]));
	return 0;
}
/******************************************************************************/
//ファイルからデータを呼び出す	no (読み込むファイルパス,読みだしたout用の配列(double)	  //
/******************************************************************************/
int call_BSwave(char fname[200],COMPLEX proj_BSwave[XYZnodeSites]){
	
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
			infile.read( ( char * ) &proj_BSwave[id], sizeof( COMPLEX ) );
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
					//cout <<"	"<<x<<y<<z<<"	"<<proj_BSwave[(x) +XnodeSites*((y) + YnodeSites*((z)))]<<endl;
				}}}
	}
if (switch1==1) {
	std::ifstream ifs( fname , ios::in );			/* テキストモードで */
	if(! ifs  )
	{
		cout << fname<<"ファイルが開きません";
		return 1;
	}
	int tmpx[XYZnodeSites];
	int tmpy[XYZnodeSites];
	int tmpz[XYZnodeSites];
	for (int z=0; z<ZnodeSites; z++) {
	for (int y=0; y<YnodeSites; y++) {
	for (int x=0; x<XnodeSites; x++) {
				ifs>>std::setw(3)>>tmpx[(x) +XnodeSites*((y) + YnodeSites*((z)))]>>std::setw(3)>>tmpy[(x) +XnodeSites*((y) + YnodeSites*((z)))]>>std::setw(3)>>tmpz[(x) +XnodeSites*((y) + YnodeSites*((z)))]>>setprecision(15)>> proj_BSwave(x,y,z).real()>>setprecision(15)>>proj_BSwave(x,y,z).imag();
							//cout << proj_BSwave(x,y,z).real()<<"///"<<proj_BSwave(x,y,z).imag()<<endl;
	}
	}
	}
}
	return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/************************************************************************************************************************************************/
//out_TdepBSwaveを使い1回微分したBS波動関数を書き出す。out_dataでTdep_BSwaveの適用先をconf j番目にする no (fanameをいじると出力ファイル,入力データの配列[Maxsize])  //
/************************************************************************************************************************************************/
void out_data(int j,int it,COMPLEX local[XYZnodeSites*TnodeSites]){
	char fname[200]={0};
	sprintf(fname,"%s/binTder.%06d-%06d.%s.it%02d",out_dir_path.c_str(),binnumber,j,base,it);
	out_TdepBSwave(it,&(fname[0]),&(local[0]));
}
/**************************************************************************/
//結果を書き出す	no (書きだすファイルパス,BS波動関数)	  //
/**************************************************************************/
int out_TdepBSwave(int it,char fname[200],COMPLEX Tdep_BSwave[XYZnodeSites*TnodeSites]){
	int switch1=0;
	if (switch1==0) {
	ofstream ofs_xyz;		
	ofs_xyz.open(fname,ios::out|ios::binary|ios::trunc);
	if (!ofs_xyz.is_open()) {
		cout << "ERROR output file can't open (can't create)"<<endl;
		return EXIT_FAILURE;
	}
	for (int z=0; z<ZnodeSites; z++) {
	for (int y=0; y<YnodeSites; y++) {
	for (int x=0; x<XnodeSites; x++) {
				ofs_xyz.write((const char*) &(Tdep_BSwave(x,y,z,it)),sizeof( COMPLEX ));
					//cout << x<<"	"<<y<<"	"<<z<<"	"<<it<<"	"<<Tdep_BSwave(x,y,z,it)<<endl;
			}
		}
	}
	ofs_xyz.close();
	}
	
if (switch1==1) {
	std::ofstream ofs;
	ofs.open(&(fname[0]));
	//	ofs<<"#"<<std::setw(7)<<"r"<<setprecision(15)<<"BSwave(r)real part"<<setprecision(15)<<"BSwave(r)imaginary part"<<endl;
	for (int z=0; z<ZnodeSites; z++) {
		for (int y=0; y<YnodeSites; y++) {
			for (int x=0; x<XnodeSites; x++) {
				ofs<<std::setw(3)<<x<<std::setw(3)<<y<<std::setw(3)<<z<<setprecision(15)<< Tdep_BSwave(x,y,z,it).real()<<setprecision(15)<<Tdep_BSwave(x,y,z,it).imag()<<endl;
			}
		}
	}
	ofs.close();}
	return 0;
	}

/*******************************************************************/
//1階の時間差分を計算する -{f(x+1)-f(x-1)}/2を用いた
/*******************************************************************/
void calc_Tdep(COMPLEX sub_BSwave[XYZnodeSites*TnodeSites],COMPLEX Tdep_BSwave[XYZnodeSites*TnodeSites]){
	COMPLEX* sub_sub_BSwave=new COMPLEX[XYZnodeSites*TnodeSites];
	for (int z=0; z<ZnodeSites; z++) {
	for (int y=0; y<YnodeSites; y++) {
	for (int x=0; x<XnodeSites; x++) {
	for (int it=T_in; it<T_fi+1;it++ ) {
		sub_sub_BSwave(x,y,z,it)=log(sub_BSwave(x,y,z,it));
	}
	for (int it=T_in; it<T_fi+1;it++ ) {
		Tdep_BSwave(x,y,z,it)=(sub_sub_BSwave(x,y,z,it-1)-sub_sub_BSwave(x,y,z,it+1))/2.0;
	}
	}}}
	delete []sub_sub_BSwave;
}


