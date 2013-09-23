#include "../include/analys.h"


#define proj_BSwave(x,y,z) proj_BSwave[(x) +XnodeSites*((y) + YnodeSites*((z)))]
#define Laplacian_BSwave(x,y,z) Laplacian_BSwave[(x) +XnodeSites*((y) + YnodeSites*((z)))]			   

#define I       std::complex<double>(0.0,1.0)
using namespace std;
typedef std::complex<double> COMPLEX;


int call_data(int,int,COMPLEX[],COMPLEX[]);
int call_BSwave(char[],COMPLEX[]);
void out_data(int,int,COMPLEX[]);
int out_LaplacianBSwave(char[],char[],COMPLEX[]);
void sum(COMPLEX[],COMPLEX[],COMPLEX[]);
void divide(COMPLEX[],COMPLEX[]);



int main(int argc,char** argv){

	dir_path=argv[1];
	cout <<"Directory path  ::"<<dir_path<<endl;
        in_dir_path = dir_path;
	out_dir_path = dir_path;

	root_mkdir(out_dir_path.c_str());
        out_dir_path = out_dir_path + "/Potential";
	root_mkdir(out_dir_path.c_str());
        out_dir_path = out_dir_path + "/binPotential";
	root_mkdir(out_dir_path.c_str());
    xyz_out_dir_path = out_dir_path + "/xyz";
	root_mkdir(xyz_out_dir_path.c_str());
    r_out_dir_path = out_dir_path + "/r";
	root_mkdir(r_out_dir_path.c_str());

	
	for ( int j=0; j<binnumber; j++) {
		for (int it=T_in; it<T_fi+1; it++) {
			COMPLEX* Tdep_BSwave= new COMPLEX[XYZnodeSites];
			COMPLEX* Laplacian_BSwave=new COMPLEX[XYZnodeSites];
			//cout<<"今conf"<<j<<"、時間"<<it<<"データの読み出し中"<<endl;
			call_data(j,it,&(Laplacian_BSwave[0]),&(Tdep_BSwave[0]));
			//cout<<"今conf"<<j<<"、時間"<<it<<"データの読み出し完了"<<endl;
			cout<<"今conf"<<j<<"、時間"<<it<<"の計算を実行中"<<endl;
			COMPLEX* potencial=new COMPLEX[XYZnodeSites];
			sum(&(Laplacian_BSwave[0]),&(Tdep_BSwave[0]),&(potencial[0]));
			delete []Laplacian_BSwave;
			delete []Tdep_BSwave;
			out_data(j,it,&(potencial[0]));
			delete []potencial;
		}
	}
	cout << "finished"<<endl;
	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/************************************************************************************************************************************************/
//call_q_propを使いquarkプロパゲータを呼び出す。call_dataでcall_q_propの適用先をconf i番目にする no (fanameをいじると読み込むファイル,データを入れる配列[Maxsize])  //
/************************************************************************************************************************************************/
int call_data(int j,int it,COMPLEX local[XYZnodeSites],COMPLEX local1[XYZnodeSites]){
	char fnamelap[200]={0};
	char fnameTdep[200]={0};
	
	sprintf(fnamelap,"%s/Laplacian/binLaplacian/xyz/Lap.%06d-%06d.%s.it%02d",in_dir_path.c_str(),binnumber,j,base,it);
	sprintf(fnameTdep,"%s/Tder/binTder/xyz/binTder.%06d-%06d.%s.it%02d",in_dir_path.c_str(),binnumber,j,base,it);
	
	call_BSwave(&(fnamelap[0]),&(local[0]));
	call_BSwave(&(fnameTdep[0]),&(local1[0]));

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
		  cout << "ERROR file can't open (no exist)   ::"<<fname<<endl;
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
			//cout << proj_BSwave[id]<<endl;
			id=id+1; 
			//cout << id<<endl;
		}
		static int tmp=0;
		if (tmp==0) {
			cout <<"reading data size is	;;"<<id<<endl;
			tmp=tmp+1;
		}
		//endian_convert((double*)data,datasize*2);
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
				//			cout << proj_BSwave(x,y,z).real()<<"///"<<proj_BSwave(x,y,z).imag()<<endl;
			}
		}
	}
}
	return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/************************************************************************************************************************************************/
//out_LaplacianBSwaveを使い２回微分したBS波動関数を書き出す。out_dataでout_LaplacianBSwaveの適用先をconf j番目にする no (fanameをいじると出力ファイル,入力データの配列[Maxsize])  //
/************************************************************************************************************************************************/
void out_data(int j,int it,COMPLEX local[XYZnodeSites]){
    char fnamer[200]={0};
	char fnamexyz[200]={0};

	sprintf(fnamexyz,"%s/Potential.%06d-%06d.%s.it%02d",xyz_out_dir_path.c_str(),binnumber,j,base,it);
	sprintf(fnamer,"%s/Potential.%06d-%06d.%s.it%02d",r_out_dir_path.c_str(),binnumber,j,base,it);
	out_LaplacianBSwave(&(fnamer[0]),&(fnamexyz[0]),&(local[0]));
}
/**************************************************************************/
//結果を書き出す	no (書きだすファイルパス,BS波動関数)	  //
/**************************************************************************/
int out_LaplacianBSwave(char fnamer[200],char fnamexyz[200],COMPLEX Laplacian_BSwave[XYZnodeSites]){
#define radius_sq(x,y,z) ((x)*(x) + (y)*(y) + (z)*(z))
#define min(a,b) (((a) < (b)) ? (a) : (b))
    int r_sq=0;
    int r_sq_count[r_sq_max]={0};
    COMPLEX Laplacian_BSwave_ave[r_sq_max];
    
    
	int switch1=0;
	std::ofstream ofsxyz;
	if (switch1==0) {
	ofsxyz.open(fnamexyz,ios::out|ios::binary|ios::trunc);
	if (!ofsxyz.is_open()) {
		cout << "ERROR output file can't open (no exist)"<<endl;
        exit(1);
		return EXIT_FAILURE;
	}
	for (int z=0; z<ZnodeSites; z++) {
	for (int y=0; y<YnodeSites; y++) {
	for (int x=0; x<XnodeSites; x++) {
				ofsxyz.write((const char*) &(Laplacian_BSwave(x,y,z)),sizeof( COMPLEX ));
		//cout <<x<<"	"<<y<<"	"<<z<<"	"<< Laplacian_BSwave(x,y,z).real()<<"	"<<Laplacian_BSwave(x,y,z).imag()<<endl;
        int r_sq = radius_sq( min(x,XnodeSites-x), min(y,YnodeSites-y), min(z,ZnodeSites-z) );
        r_sq_count[r_sq]= r_sq_count[r_sq] +1;
        Laplacian_BSwave_ave[r_sq] =Laplacian_BSwave_ave[r_sq] + Laplacian_BSwave[(x) +XnodeSites*((y) + YnodeSites*((z)))];
    }}}
				ofsxyz.close();
        std::ofstream ofsr;
        ofsr.open( &(fnamer[0]));
        ofsr.setf(ios::scientific);
        
        for (int r_sq=0; r_sq<r_sq_max; r_sq++) {
            
            if ( r_sq_count[r_sq] == 0 ) continue;
            
            Laplacian_BSwave_ave[r_sq] =Laplacian_BSwave_ave[r_sq]/ ((double) r_sq_count[r_sq]);
            
            float rad = sqrt((float)r_sq);
            ofsr<<rad<<"	"<< Laplacian_BSwave_ave[r_sq].real()<<"	"<<Laplacian_BSwave_ave[r_sq].imag()<<endl;
        }
        ofsr.close();
        
        return 0;

    }
	if (switch1==1) {
	ofsxyz.open(&(fnamexyz[0]));
	//	ofs<<"#"<<std::setw(7)<<"r"<<setprecision(15)<<"BSwave(r)real part"<<setprecision(15)<<"BSwave(r)imaginary part"<<endl;
	for (int z=0; z<ZnodeSites; z++) {
	for (int y=0; y<YnodeSites; y++) {
	for (int x=0; x<XnodeSites; x++) {
				ofsxyz<<std::setw(3)<<x<<std::setw(3)<<y<<std::setw(3)<<z<<setprecision(15)<< Laplacian_BSwave(x,y,z).real()<<setprecision(15)<<Laplacian_BSwave(x,y,z).imag()<<endl;
        int r_sq = radius_sq( min(x,XnodeSites-x), min(y,YnodeSites-y), min(z,ZnodeSites-z) );
        r_sq_count[r_sq]= r_sq_count[r_sq] +1;
        Laplacian_BSwave_ave[r_sq] =Laplacian_BSwave_ave[r_sq] + Laplacian_BSwave[(x) +XnodeSites*((y) + YnodeSites*((z)))];
		
            }
		}
	}
		ofsxyz.close();
        std::ofstream ofsr;
        ofsr.open( &(fnamer[0]));
        ofsr.setf(ios::scientific);
        
        for (int r_sq=0; r_sq<r_sq_max; r_sq++) {
            
            if ( r_sq_count[r_sq] == 0 ) continue;
            
            Laplacian_BSwave_ave[r_sq] =Laplacian_BSwave_ave[r_sq]/ ((double) r_sq_count[r_sq]);
            
            float rad = sqrt((float)r_sq);
            ofsr<<rad<<"	"<< Laplacian_BSwave_ave[r_sq].real()<<"	"<<Laplacian_BSwave_ave[r_sq].imag()<<endl;
        }
        ofsr.close();
        
        return 0;

    }
    }

/*******************************************************************/
//2階差分を計算する  f(x-1)+f(x+1)-2f(x)なので配列は Xnodesites→Xnodesites-2の長さになる
/*******************************************************************/
void sum(COMPLEX Laplacian_BSwave[XYZnodeSites],COMPLEX Tdep_BSwave[XYZnodeSites],COMPLEX potencial[XYZnodeSites]){
	for (int z=0; z<ZnodeSites; z++) {
	for (int y=0; y<YnodeSites; y++) {
	for (int x=0; x<XnodeSites; x++) {
		potencial[(x) +XnodeSites*((y) + YnodeSites*((z)))]=Laplacian_BSwave[(x) +XnodeSites*((y) + YnodeSites*((z)))]+Tdep_BSwave[(x) +XnodeSites*((y) + YnodeSites*((z)))];		
			}}}}
