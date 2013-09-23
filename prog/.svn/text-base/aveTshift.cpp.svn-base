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
    out_dir_path=out_dir_path + "/impBSwave";
	root_mkdir(out_dir_path.c_str());
    out_dir_path = out_dir_path + "/aveTshift";
	root_mkdir(out_dir_path.c_str());


	for (int it=T_in; it < (T_fi +1); it++) {
		cout << "time"<<it<<endl;
		COMPLEX *data= new COMPLEX[datasize*Confsize];
		call_data(it,&(data[0]));
		out_data(it,data);
	}	
	cout << "finish"<<endl;
	return 0;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*******************************************************************************************************************************************/
//call_fileを利用して、複数のファイルから呼び出したデータを配列としてまとめる (読みだしたout用の配列(double)	この中のfanameをいじると入力ファイルを変えられる。  //
/******************************************************************************************************************************************/
void call_data(int it,COMPLEX data[Confsize*datasize]){
        COMPLEX *buffer = new COMPLEX[Confsize*datasize];
	for (int j=0; j<Confsize; j++) {
    for (int tmp=0; tmp<Tshiftsize; tmp++) {
	    char fname[200]={0};
	    COMPLEX local[datasize]={0};
	    sprintf(fname,"%s/impBSwave/noise_reduction/OmgOmgwave_PH1.+%03d+%03d.%s-%06d",in_dir_path.c_str(),it,tshift[tmp],base,j);//hachi
		cout << fname<<"reading now"<<endl;
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
void out_data(int it,COMPLEX data[datasize*Confsize]){
	for (int j=0; j<Confsize; j++) {
		char fname[200]={0};
	sprintf(fname,"%s/OmgOmgwave_PH1.+%03d.%s-%06d",out_dir_path.c_str(),it,base,j);
			out_file(j,fname,data);
	}
}
/*******************************************************************************************************/
//結果を書き出す	no (書きだすファイルパス,平均(データの大きさ),誤差(データの大きさ)	  //
/******************************************************************************************************/
int out_file(int j, char fname[200],COMPLEX data[Confsize*datasize]){
	COMPLEX outdata[datasize];
	if (switch1==0) {
		ofstream ofs_xyz;		
		ofs_xyz.open(fname,ios::out|ios::binary|ios::trunc);
		if (!ofs_xyz.is_open()) {
			cout << "ERROR output file can't open (no exist)"<<endl;
			return EXIT_FAILURE;
		}
		for (int z=0; z<ZnodeSites; z++) {
		for (int y=0; y<YnodeSites; y++) {
		for (int x=0; x<XnodeSites; x++) {
			int	id=(x) +XnodeSites*((y) + YnodeSites*((z)));
					outdata[id]=data(id,j);
					ofs_xyz.write((const char*) &(outdata[id]),sizeof( COMPLEX ));
		}}}
	}
	if (switch1==1) {
	std::ofstream ofs( &(fname[0]),std::ios::out | std::ios::trunc);
	for (int z=0; z<ZnodeSites; z++) {
	for (int y=0; y<YnodeSites; y++) {
	for (int x=0; x<XnodeSites; x++) {
	int	id=(x) +XnodeSites*((y) + YnodeSites*((z)));
		outdata[id]=data(id,j);
		ofs<<std::setw(3)<<x<<std::setw(3)<<y<<std::setw(3)<<z<<std::setw(21)<<setprecision(15)<< outdata[id].real()<<std::setw(21)<<setprecision(15)<<outdata[id].imag()<<endl;
		//cout<<std::setw(3)<<x<<std::setw(3)<<y<<std::setw(3)<<z<<std::setw(21) <<setprecision(15)<< outdata[id].real()<<std::setw(21) <<setprecision(15)<<outdata[id].imag()<<endl;
	}}}}
	return 0;
}
//*******************************************************************/
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


