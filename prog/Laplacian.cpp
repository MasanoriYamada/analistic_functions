/***********************************************/
/*Laplacian part(1/mR*(Δψ))*/
/***********************************************/
#include "../include/analys.h"


/**************************************************************************************/
#define proj_BSwave(x,y,z) proj_BSwave[(x) +XnodeSites*((y) + YnodeSites*((z)))]
#define Laplacian_BSwave(x,y,z) Laplacian_BSwave[(x) +XnodeSites*((y) + YnodeSites*((z)))]											//ただし最初と最後には値が入ってない。	
/***********************************synkでprojectionを取る場合に必要*********************************/

#define I       std::complex<double>(0.0,1.0)
using namespace std;
typedef std::complex<double> COMPLEX;


int call_file(int,int,COMPLEX[]);
int call_data(char[],COMPLEX[]);
void out_data(int,int,COMPLEX[]);
int out_LaplacianBSwave(char[],COMPLEX[]);
void calc_Laplacian(COMPLEX[],COMPLEX[]);
void divide(COMPLEX[],COMPLEX[]);



int main(int argc , char** argv){
	
	dir_path=argv[1];
	cout <<"Directory path  ::"<<dir_path<<endl;
        in_dir_path = dir_path;
	out_dir_path = dir_path;

	root_mkdir(out_dir_path.c_str());
        out_dir_path=out_dir_path + "/Laplacian";
	root_mkdir(out_dir_path.c_str());
        out_dir_path = out_dir_path + "/binLaplacian";
	root_mkdir(out_dir_path.c_str());
        out_dir_path = out_dir_path + "/xyz";
	root_mkdir(out_dir_path.c_str());
	in_dir_path = in_dir_path + "/Rcor/binR/xyz";


	for ( int j=0; j<binnumber; j++) {
		for (int it=T_in; it<T_fi+1; it++) {
			COMPLEX* proj_BSwave = new COMPLEX[XYZnodeSites];
			COMPLEX* Laplacian_BSwave=new COMPLEX[XYZnodeSites];
			//cout<<"今conf"<<j<<"、時間"<<it<<"データの読み出し中"<<endl;
			call_file(j,it,&(proj_BSwave[0]));
			//cout<<"今conf"<<j<<"、時間"<<it<<"データの読み出し完了"<<endl;
			//cout<<"今conf"<<j<<"、時間"<<it<<"の計算を実行中"<<endl;
			calc_Laplacian(&(proj_BSwave[0]),&(Laplacian_BSwave[0]));
			divide(&(Laplacian_BSwave[0]),&(proj_BSwave[0]));
			delete []proj_BSwave;
			out_data(j,it,&(Laplacian_BSwave[0]));
			delete []Laplacian_BSwave;
		}
	}
	cout << "finished"<<endl;
	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/************************************************************************************************************************************************/
//call_q_propを使いquarkプロパゲータを呼び出す。call_dataでcall_q_propの適用先をconf i番目にする no (fanameをいじると読み込むファイル,データを入れる配列[Maxsize])  //
/************************************************************************************************************************************************/
int call_file(int b,int it,COMPLEX local[XYZnodeSites]){
	char fname[300]={0};
	sprintf(fname,"%s/binRwave.%s.%06d-%06d.it%03d",in_dir_path.c_str(),base,binnumber,b,it);
	call_data(&(fname[0]),&(local[0]));
	return 0;
}
/******************************************************************************/
//ファイルからデータを呼び出す	no (読み込むファイルパス,読みだしたout用の配列(double)	  //
/******************************************************************************/
int call_data(char fname[300], COMPLEX data[XYZnodeSites]){
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
	infile.close();
		return 0;
}	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/************************************************************************************************************************************************/
//out_LaplacianBSwaveを使い２回微分したBS波動関数を書き出す。out_dataでout_LaplacianBSwaveの適用先をconf j番目にする no (fanameをいじると出力ファイル,入力データの配列[Maxsize])  //
/************************************************************************************************************************************************/
void out_data(int j,int it,COMPLEX local[XYZnodeSites]){
	char fname[200]={0};
	sprintf(fname,"%s/Lap.%06d-%06d.%s.it%02d",out_dir_path.c_str(),binnumber,j,base,it);
out_LaplacianBSwave(&(fname[0]),&(local[0]));
}
/**************************************************************************/
//結果を書き出す	no (書きだすファイルパス,BS波動関数)	  //
/**************************************************************************/
int out_LaplacianBSwave(char fname[200],COMPLEX Laplacian_BSwave[XYZnodeSites]){
	ofstream ofs_xyz;		
	ofs_xyz.open(fname,ios::out|ios::binary|ios::trunc);
	if (!ofs_xyz.is_open()) {
		cout << "ERROR output file can't open (can't create)"<<endl;
		return EXIT_FAILURE;
	}
	//	ofs<<"#"<<std::setw(7)<<"r"<<setprecision(15)<<"BSwave(r)real part"<<setprecision(15)<<"BSwave(r)imaginary part"<<endl;
	for (int z=0; z<ZnodeSites; z++) {
	for (int y=0; y<YnodeSites; y++) {
	for (int x=0; x<XnodeSites; x++) {
		ofs_xyz.write((const char*) &(Laplacian_BSwave(x,y,z)),sizeof( COMPLEX ));
		//cout << x<<"	"<<y<<"	"<<z<<"	"<<Laplacian_BSwave(x,y,z)<<endl;
			}
		}
	}
	ofs_xyz.close();
	return 0;
	}

/*******************************************************************/
//2階差分を計算する  f(x-1)+f(x+1)-2f(x)なので配列は Xnodesites→Xnodesites-2の長さになる
/*******************************************************************/
void calc_Laplacian(COMPLEX proj_BSwave[XYZnodeSites],COMPLEX Laplacian_BSwave[XYZnodeSites]){
	for (int z=0; z<ZnodeSites; z++) {
		for (int y=0; y<YnodeSites; y++) {
			for (int x=0; x<XnodeSites; x++) {
				//cout <<"projBS"<<x<<y<<z<<proj_BSwave(x,y,z)<<endl;
				/*if (x==0) {
					if (y==0) {
						if (z==0) {Laplacian_BSwave(x,y,z)=proj_BSwave(x+1,y,z)+proj_BSwave(XnodeSites-1,y,z)+proj_BSwave(x,y+1,z)+proj_BSwave(x,YnodeSites-1,z)+proj_BSwave(x,y,z+1)+proj_BSwave(x,y,ZnodeSites-1)-6.0*proj_BSwave(x,y,z);}
						if (z==ZnodeSites-1) {Laplacian_BSwave(x,y,z)=proj_BSwave(x+1,y,z)+proj_BSwave(XnodeSites-1,y,z)+proj_BSwave(x,y+1,z)+proj_BSwave(x,YnodeSites-1,z)+proj_BSwave(x,y,0)+proj_BSwave(x,y,z-1)-6.0*proj_BSwave(x,y,z);}
						if (z!=0 && z!=ZnodeSites-1){Laplacian_BSwave(x,y,z)=proj_BSwave(x+1,y,z)+proj_BSwave(XnodeSites-1,y,z)+proj_BSwave(x,y+1,z)+proj_BSwave(x,YnodeSites-1,z)+proj_BSwave(x,y,z+1)+proj_BSwave(x,y,z-1)-6.0*proj_BSwave(x,y,z);}
					}
					if (y==YnodeSites-1) {
						if (z==0) {Laplacian_BSwave(x,y,z)=proj_BSwave(x+1,y,z)+proj_BSwave(XnodeSites-1,y,z)+proj_BSwave(x,0,z)+proj_BSwave(x,y-1,z)+proj_BSwave(x,y,z+1)+proj_BSwave(x,y,ZnodeSites-1)-6.0*proj_BSwave(x,y,z);}
						if (z==ZnodeSites-1) {Laplacian_BSwave(x,y,z)=proj_BSwave(x+1,y,z)+proj_BSwave(XnodeSites-1,y,z)+proj_BSwave(x,0,z)+proj_BSwave(x,y-1,z)+proj_BSwave(x,y,0)+proj_BSwave(x,y,z-1)-6.0*proj_BSwave(x,y,z);}
						if (z!=0 && z!=ZnodeSites-1){Laplacian_BSwave(x,y,z)=proj_BSwave(x+1,y,z)+proj_BSwave(XnodeSites-1,y,z)+proj_BSwave(x,0,z)+proj_BSwave(x,y-1,z)+proj_BSwave(x,y,z+1)+proj_BSwave(x,y,z-1)-6.0*proj_BSwave(x,y,z);}
					}
					if (y!=0 && y!=YnodeSites-1){
						if (z==0) {Laplacian_BSwave(x,y,z)=proj_BSwave(x+1,y,z)+proj_BSwave(XnodeSites-1,y,z)+proj_BSwave(x,y+1,z)+proj_BSwave(x,y-1,z)+proj_BSwave(x,y,z+1)+proj_BSwave(x,y,ZnodeSites-1)-6.0*proj_BSwave(x,y,z);}
						if (z==ZnodeSites-1) {Laplacian_BSwave(x,y,z)=proj_BSwave(x+1,y,z)+proj_BSwave(XnodeSites-1,y,z)+proj_BSwave(x,y+1,z)+proj_BSwave(x,y-1,z)+proj_BSwave(x,y,0)+proj_BSwave(x,y,z-1)-6.0*proj_BSwave(x,y,z);}
						if (z!=0 && z!=ZnodeSites-1){Laplacian_BSwave(x,y,z)=proj_BSwave(x+1,y,z)+proj_BSwave(XnodeSites-1,y,z)+proj_BSwave(x,y+1,z)+proj_BSwave(x,y-1,z)+proj_BSwave(x,y,z+1)+proj_BSwave(x,y,z-1)-6.0*proj_BSwave(x,y,z);}
					}
				}
				if (x==XnodeSites-1){
					if (y==0) {
						if (z==0) {Laplacian_BSwave(x,y,z)=proj_BSwave(0,y,z)+proj_BSwave(x-1,y,z)+proj_BSwave(x,y+1,z)+proj_BSwave(x,YnodeSites-1,z)+proj_BSwave(x,y,z+1)+proj_BSwave(x,y,ZnodeSites-1)-6.0*proj_BSwave(x,y,z);}
						if (z==ZnodeSites-1) {Laplacian_BSwave(x,y,z)=proj_BSwave(0,y,z)+proj_BSwave(x-1,y,z)+proj_BSwave(x,y+1,z)+proj_BSwave(x,YnodeSites-1,z)+proj_BSwave(x,y,0)+proj_BSwave(x,y,z-1)-6.0*proj_BSwave(x,y,z);}
						if (z!=0 && z!=ZnodeSites-1){Laplacian_BSwave(x,y,z)=proj_BSwave(0,y,z)+proj_BSwave(x-1,y,z)+proj_BSwave(x,y+1,z)+proj_BSwave(x,YnodeSites-1,z)+proj_BSwave(x,y,z+1)+proj_BSwave(x,y,z-1)-6.0*proj_BSwave(x,y,z);}
					}
					if (y==YnodeSites-1) {
						if (z==0) {Laplacian_BSwave(x,y,z)=proj_BSwave(0,y,z)+proj_BSwave(x-1,y,z)+proj_BSwave(x,0,z)+proj_BSwave(x,y-1,z)+proj_BSwave(x,y,z+1)+proj_BSwave(x,y,ZnodeSites-1)-6.0*proj_BSwave(x,y,z);}
						if (z==ZnodeSites-1) {Laplacian_BSwave(x,y,z)=proj_BSwave(0,y,z)+proj_BSwave(x-1,y,z)+proj_BSwave(x,0,z)+proj_BSwave(x,y-1,z)+proj_BSwave(x,y,0)+proj_BSwave(x,y,z-1)-6.0*proj_BSwave(x,y,z);}
						if (z!=0 && z!=ZnodeSites-1){Laplacian_BSwave(x,y,z)=proj_BSwave(0,y,z)+proj_BSwave(x-1,y,z)+proj_BSwave(x,0,z)+proj_BSwave(x,y-1,z)+proj_BSwave(x,y,z+1)+proj_BSwave(x,y,z-1)-6.0*proj_BSwave(x,y,z);}
					}
					if (y!=0 && y!=YnodeSites-1){
						if (z==0) {Laplacian_BSwave(x,y,z)=proj_BSwave(0,y,z)+proj_BSwave(x-1,y,z)+proj_BSwave(x,y+1,z)+proj_BSwave(x,y-1,z)+proj_BSwave(x,y,z+1)+proj_BSwave(x,y,ZnodeSites-1)-6.0*proj_BSwave(x,y,z);}
						if (z==ZnodeSites-1) {Laplacian_BSwave(x,y,z)=proj_BSwave(0,y,z)+proj_BSwave(x-1,y,z)+proj_BSwave(x,y+1,z)+proj_BSwave(x,y-1,z)+proj_BSwave(x,y,0)+proj_BSwave(x,y,z-1)-6.0*proj_BSwave(x,y,z);}
						if (z!=0 && z!=ZnodeSites-1){Laplacian_BSwave(x,y,z)=proj_BSwave(0,y,z)+proj_BSwave(x-1,y,z)+proj_BSwave(x,y+1,z)+proj_BSwave(x,y-1,z)+proj_BSwave(x,y,z+1)+proj_BSwave(x,y,z-1)-6.0*proj_BSwave(x,y,z);}
					}
				}
				if(x!=0 && x!=XnodeSites-1){
					if (y==0) {
						if (z==0) {Laplacian_BSwave(x,y,z)=proj_BSwave(x+1,y,z)+proj_BSwave(x-1,y,z)+proj_BSwave(x,y+1,z)+proj_BSwave(x,YnodeSites-1,z)+proj_BSwave(x,y,z+1)+proj_BSwave(x,y,ZnodeSites-1)-6.0*proj_BSwave(x,y,z);}
						if (z==ZnodeSites-1) {Laplacian_BSwave(x,y,z)=proj_BSwave(x+1,y,z)+proj_BSwave(x-1,y,z)+proj_BSwave(x,y+1,z)+proj_BSwave(x,YnodeSites-1,z)+proj_BSwave(x,y,0)+proj_BSwave(x,y,z-1)-6.0*proj_BSwave(x,y,z);}
						if (z!=0 && z!=ZnodeSites-1){Laplacian_BSwave(x,y,z)=proj_BSwave(x+1,y,z)+proj_BSwave(x-1,y,z)+proj_BSwave(x,y+1,z)+proj_BSwave(x,YnodeSites-1,z)+proj_BSwave(x,y,z+1)+proj_BSwave(x,y,z-1)-6.0*proj_BSwave(x,y,z);}
					}
					if (y==YnodeSites-1) {
						if (z==0) {Laplacian_BSwave(x,y,z)=proj_BSwave(x+1,y,z)+proj_BSwave(x-1,y,z)+proj_BSwave(x,0,z)+proj_BSwave(x,y-1,z)+proj_BSwave(x,y,z+1)+proj_BSwave(x,y,ZnodeSites-1)-6.0*proj_BSwave(x,y,z);}
						if (z==ZnodeSites-1) {Laplacian_BSwave(x,y,z)=proj_BSwave(x+1,y,z)+proj_BSwave(x-1,y,z)+proj_BSwave(x,0,z)+proj_BSwave(x,y-1,z)+proj_BSwave(x,y,0)+proj_BSwave(x,y,z-1)-6.0*proj_BSwave(x,y,z);}
						if (z!=0 && z!=ZnodeSites-1){Laplacian_BSwave(x,y,z)=proj_BSwave(x+1,y,z)+proj_BSwave(x-1,y,z)+proj_BSwave(x,0,z)+proj_BSwave(x,y-1,z)+proj_BSwave(x,y,z+1)+proj_BSwave(x,y,z-1)-6.0*proj_BSwave(x,y,z);}
					}
					if (y!=0 && y!=YnodeSites-1){
						if (z==0) {Laplacian_BSwave(x,y,z)=proj_BSwave(x+1,y,z)+proj_BSwave(x-1,y,z)+proj_BSwave(x,y+1,z)+proj_BSwave(x,y-1,z)+proj_BSwave(x,y,z+1)+proj_BSwave(x,y,ZnodeSites-1)-6.0*proj_BSwave(x,y,z);}
						if (z==ZnodeSites-1) {Laplacian_BSwave(x,y,z)=proj_BSwave(x+1,y,z)+proj_BSwave(x-1,y,z)+proj_BSwave(x,y+1,z)+proj_BSwave(x,y-1,z)+proj_BSwave(x,y,0)+proj_BSwave(x,y,z-1)-6.0*proj_BSwave(x,y,z);}
						if (z!=0 && z!=ZnodeSites-1){Laplacian_BSwave(x,y,z)=proj_BSwave(x+1,y,z)+proj_BSwave(x-1,y,z)+proj_BSwave(x,y+1,z)+proj_BSwave(x,y-1,z)+proj_BSwave(x,y,z+1)+proj_BSwave(x,y,z-1)-6.0*proj_BSwave(x,y,z);}
					}
				}
				 */
				Laplacian_BSwave(x,y,z)=proj_BSwave((x+1)%XnodeSites,y,z)+proj_BSwave((x+XnodeSites-1)%XnodeSites,y,z)+proj_BSwave(x,(y+1)%YnodeSites,z)+proj_BSwave(x,(y+YnodeSites-1)%YnodeSites,z)+proj_BSwave(x,y,(z+1)%ZnodeSites)+proj_BSwave(x,y,(z-1+ZnodeSites)%ZnodeSites)-6.0*proj_BSwave(x,y,z);
				//cout <<"LaplacianBS"<<x<<y<<z<<Laplacian_BSwave(x,y,z)<<endl;
			}}}			
}
/*******************************************************************/
//BS波動関数で割る　1/ψ*∇ψを行う
/*******************************************************************/
void divide(COMPLEX Laplacian_BSwave[XYZnodeSites],COMPLEX proj_BSwave[XYZnodeSites]){
	for (int id=0; id<XYZnodeSites; id++) {
		Laplacian_BSwave[id]=Laplacian_BSwave[id]/(effectivemass*proj_BSwave[id]);
	}
}


