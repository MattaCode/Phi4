/*
***********
Phi4Base.h
***********

Tartalmazza a sz�ks�ges oszt�lydefin�ci�kat

*/

#ifndef PHI4BASE_H
#define PHI4BASE_H

#include<iostream>
#include<fstream>
#include<vector>
#include<string.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include <cstdlib>
#include <ctime> 
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include "Array.h"
#include "fftw3.h"
#include "fftw++.h"

using namespace std;
using namespace Array;
using namespace fftwpp;

/*Seg�doszt�ly: 3d array*/

namespace fftwpp {

std::ifstream fftw::ifWisdom;
std::ofstream fftw::ofWisdom;
bool fftw::Wise=false;
//Original
//const double fftw::twopi=2.0*acos(-1.0);
//Changed:
const double fftw::twopi=2.0*M_PI;



// User settings:
unsigned int fftw::effort=FFTW_MEASURE;
const char *fftw::WisdomName="wisdom3.txt";
}


//random gener�tor
boost::mt19937 gen(time(NULL));

//egyenletes vals�ggel sorsol
//kell neki a tartom�nyhat�r
double GetRandom(double x){
	double randnum;
	boost::random::uniform_real_distribution<> dist(-x, x);


	return randnum=dist(gen);
}


/*
Field oszt�ly:
a t�r, amin a Fi4 modellt vizsg�ljuk
*/

class Field{
	
public: 
	//adattagok
	double a; //r�cs�lland�
	bool isempty; //�res-e a t�r
	bool isemptyFFT; //�res-e az FFT t�r

	//a t�r kiterjed�se
	int xmax;
	int ymax;
	int zmax;

	//az FFT miatt sz�ks�ges v�ltoz�k
	unsigned int nzp;
	size_t align1;
	size_t align2;
	size_t align3;
	array3<double> phifield; //3D t�r Phi �rt�keit t�rolja
	array3<double> pifield;  //3D t�r Pi �rt�keit t�rolja
	array3<Complex> phiFFTfield; //3D t�r, Fourier-transzform�lt Phi �rt�kei, nem szimmetrikus
	array3<Complex> piFFTfield;  //3D t�r, Fourier-transzform�lt Pi �rt�kei, nem szimmetrikus
	array3<Complex> dataFFTfield;//3D t�r, Fourier-transzform�lt data �rt�kei, nem szimmetrikus
	
	
	//friend class EvolveModel;

//public f�ggv�nyek:

	//Konstruktorok
	//3D m�retek �s r�cs�lland� megad�s�val, �resen l�trehozva
	Field(  int xm,   int ym,   int zm, double grid):a(grid),xmax(xm),ymax(ym),zmax(zm),
		nzp(zmax/2+1),align1(sizeof(Complex)),align2(sizeof(Complex)),align3(sizeof(Complex)),phifield(xmax,ymax,zmax,align1),pifield(xmax,ymax,zmax,align2),phiFFTfield(xmax,ymax,nzp,align1),
		piFFTfield(xmax,ymax,nzp,align2),dataFFTfield(xmax,ymax,nzp,align3){
			//a seg�dterekr�l destruktorban is gondoskodni kell
			//Forward3Phi=new rcfft3d(zmax,phifield,phiFFTfield);
			//Forward3Pi=new rcfft3d(zmax,pifield,piFFTfield);

			//a terek kezdetben �resek
			isempty=true;
			isemptyFFT=true;
			//cout<<"FIELD KONSTR"<<endl;

	}
	//Default konstruktor
	Field():a(1),xmax(50),ymax(50),zmax(50),nzp(26),
		align1(sizeof(Complex)),align2(sizeof(Complex)),align3(sizeof(Complex)),phifield(xmax,ymax,zmax,align1),pifield(xmax,ymax,zmax,align2),phiFFTfield(xmax,ymax,nzp,align1),
		piFFTfield(xmax,ymax,nzp,align2),dataFFTfield(xmax,ymax,nzp,align3){
			//a seg�dterekr�l destruktorban is gondoskodni kell
			//Forward3Phi=new rcfft3d(zmax,phifield,phiFFTfield);
			//Forward3Pi=new rcfft3d(zmax,pifield,piFFTfield);

			//a terek kezdetben �resek
			isempty=true;
			isemptyFFT=true;
			//cout<<"FIELD KONSTR NO PARA"<<endl;
	}

	//felt�lt�tt t�rrel f�jlb�l l�trehozva
	Field(  int xm,   int ym,   int zm, double grid,const char* filename):a(grid),xmax(xm),ymax(ym),zmax(zm),
		nzp(zmax/2+1),align1(sizeof(Complex)),align2(sizeof(Complex)),align3(sizeof(Complex)),phifield(xmax,ymax,zmax,align1),pifield(xmax,ymax,zmax,align2),phiFFTfield(xmax,ymax,nzp,align1),
		piFFTfield(xmax,ymax,nzp,align2),dataFFTfield(xmax,ymax,nzp,align3){
			int dimensioncheck;

			//Forward3Phi=new rcfft3d(zmax,phifield,phiFFTfield);
			//Forward3Pi=new rcfft3d(zmax,pifield,piFFTfield);

			dimensioncheck = this->FillField(filename);
			//ha nem siker�lt -teljesen-felt�lteni a teret
			if(dimensioncheck==0 || dimensioncheck!=xmax*ymax*zmax+1){
				std::cout<<"Warning!!! Field still EMPTY or file dimension error!"<<std::endl;
				isempty=true;
			}
			//ha a teret felt�lt�tt�k, m�r nem �res
			else isempty=false;
			//de a Fourier-transzform�lt m�g nincs el��ll�tva
			isemptyFFT=true;
			//cout<<"FIELD KONSTR FILE"<<endl;
	}

		//a t�r lem�sol�sa
	Field& operator=(const Field& field){
		if(this!=&field){
			xmax=field.xmax;
			ymax=field.ymax;
			zmax=field.zmax;
			align1=field.align1;
			align2=field.align2;
			align3=field.align3;
			nzp=field.nzp;
			isempty=field.isempty;
			isemptyFFT=field.isemptyFFT;
			a=field.a;
			phifield=field.phifield;
			pifield=field.pifield;
			phiFFTfield=field.phiFFTfield;
			piFFTfield=field.piFFTfield;
			dataFFTfield=field.dataFFTfield;
		}
		return *this;
	}

	
	//Copykonstruktor �rt�kad�s felhaszn�l�s�val
	Field(const Field& field){
	//cout<<"FIELD COPY"<<endl;
			*this=field;
	}

	//direkt fieldek felt�lt�se f�jlb�l - konstruktorban is ezt haszn�ljuk
	//visszat�r a beolvasott sorok sz�m�val
	int FillField(const char* filename){
		std::ifstream ifs;
		int linenumber=0; //sz�moljuk a beolvasott sorokat
		//f�jl nyit�sa
		ifs.open(filename,std::ios::in);

		//ha nem siker�lt nyitni
		if(!(ifs.is_open()))
		{
			throw "Fajl nyitas nem sikerult";
		}

		//siker�lt a f�jl nyit�s
		if(ifs.is_open())
		{
			int tempx,tempy,tempz; //koordin�t�k beolvas�s�hoz
			double tempq;      //phi �rt�k beolvas�s�hoz
			double tempp;      //pi �rt�k beolvas�s�hoz

			//f�jl v�g�ig olvas be
			while(!(ifs.eof()))
			{

				ifs>>tempx;
				ifs>>tempy;
				ifs>>tempz;
				ifs>>tempq;
				ifs>>tempp;
				(phifield(tempx,tempy,tempz))=tempq; //phi t�r megfelel� hely�re betessz�k a beolvasott adatot
				(pifield(tempx,tempy,tempz))=tempp;  //pi t�r megfelel� hely�re betessz�k a beolvasott adatot
				ifs.ignore(1); //sorv�ge karakter nem kell
				linenumber++;  //beolvastunk egy sort

			}

			ifs.close();
		}
		return linenumber;
	}
	
	//direkt field adatok ki�r�sa tetsz�leges kimenetre
	void writeField(std::ostream &os)const
	{
		os.precision(6);
		//amennyiben t�nylegesen adatokat t�rolunk
		if(!isempty){
			for(  int i=0;i<zmax;i++){ //for Z koordin�t�kra
				for(  int j=0;j<ymax;j++){//for Y koordin�t�kra
					for(  int k=0;k<xmax;k++){ //for X koordin�t�kra
						os<<std::fixed<<k<<'\t'<<j<<'\t'<<i<<'\t'<<phifield(k,j,i)<<'\t'<<pifield(k,j,i)<<std::endl;
					}
				}
			}
		}

	}

	//direkt fieldek f�jlba �r�sa
	void ToFileField(const char* filename)const{
		std::ofstream myfile;

		myfile.open (filename,std::ios::out);
		this->writeField(myfile);
		myfile.close();
	}

	//FFT fieldek - adatok ki�r�sa os-ra
	//nem szimmetrikus terek
	void writeFFTField(std::ostream &os)const{
		os.precision(6);
		if(!isemptyFFT){

			for(int i=0;i<(zmax/2+1);i++){ //for Z koordin�t�kra

				for(int j=0;j<ymax;j++){//for Y koordin�t�kra

					for(int k=0;k<xmax;k++){ //for X koordin�t�kra
						os<<std::fixed<<k<<'\t'<<j<<'\t'<<i<<'\t'<<phiFFTfield(k,j,i)<<'\t'<<piFFTfield(k,j,i)<<std::endl;
					}
				}
			}
		}
	}

	//FFT fieldek - adatok ki�r�sa f�jlba
	void ToFileFFTField(const char* filename)const{
		std::ofstream myfile;

		myfile.open (filename,std::ios::out);
		this->writeFFTField(myfile);
		myfile.close();
	}

	//FFT h�v�s
	void RunFFT(){
		
		//az FFT miatt sz�ks�ges tov�bbi v�ltoz�k
		rcfft3d Forward3Phi(this->zmax,this->phifield,this->phiFFTfield);
		rcfft3d Forward3Pi(this->zmax,this->pifield,this->piFFTfield);
		Forward3Phi.fft(phifield,phiFFTfield);
		Forward3Pi.fft(pifield,piFFTfield);
		//ezent�l m�r adatok vannak az FFT terekben is
		isemptyFFT=false;
	}

	double GetA(){
		return this->a;
	}

	int GetX(){
		return this->xmax;
	}
	int GetY(){
		return ymax;
	}

	int GetZ(){
		return zmax;
	}

	

	array3<Complex>& GetPhiFTfield(){
		return phiFFTfield;
	}

	array3<double>& GetPifield(){
		return pifield;
	}

	size_t Getalign1(){
		return align1;
	}

	size_t Getalign3(){
		return align3;
	}

	array3<Complex>& GetDataFFTfield(){
		return dataFFTfield;
	}

	//destruktor
	~Field(){
		/*if(!isemptyFFT){
			cout<<"FFT DELETE"<<endl;
		free(Forward3Phi);
		free(Forward3Pi);
             }*/
	}

};

/**********Field v�ge***************************/

/*
Model oszt�ly:
a t�ren v�grehajta az id�fejl�d�st
*/


class Model{
	/*adattagok*/
	double dt;         //id�l�p�s egys�ge
	int    dr;         //szakaszelem
	double maxtime;    //eddig megy az id�l�ptet�s
	double m2;          //t�megn�gyzet
	double lambda;     //lambda
	int    measuretime; //ennyi dt-nk�nt m�r�nk
	double g;	//phi6 csatolasi allando

	bool measenergy;   //m�rj�nk-e �sszenergi�t
	bool measpress;    //m�rj�nk-e nyom�st

	double randphi;	   //random intervallum
	double randpi;

	double xiint;      //xi random intervallum
	double gamma;

	

public:
	Field discrete; //a ter�nk, kezdeti felt�telekkel felt�ltve
	//friend class EvolveModel;
	/*f�ggv�nyei*/


	//konstruktor
	//id�l�p�s, m�retek �s r�cs�lland� megad�s�val
	Model(double t,   int xm,   int ym,   int zm, double grid,double rndphi, double rndpi):dt(t),randphi(rndphi),randpi(rndpi),discrete(xm,ym,zm,grid){
	  lambda=0;
	g=0;
	  m2=1;
	  maxtime=dt*100;
	  measuretime=1;
	  measenergy=false;
	  measpress=false;
	  xiint=this->randpi;
	  gamma=randpi/500;
	  //cout<<"MODEL KONSTR"<<endl;
	}
	
	//Default konstruktor
	Model():dt(0.1),randphi(0),randpi(1),discrete(50,50,50,1){
	  lambda=0;
	g=0;
	  m2=1;
	  maxtime=dt*100;
	  measuretime=1;
	  measenergy=false;
	  measpress=false;
	   xiint=this->randpi;
	  gamma=randpi/500;
//	  cout<<"MODEL DEF KONSTR"<<endl;
	}

	//konstruktor f�jlb�l val� l�trehoz�shoz
	Model(double t,   int xm,   int ym,   int zm, double grid, double rndphi, double rndpi, const char* filename):dt(t),randphi(rndphi),randpi(rndpi),discrete(xm,ym,zm,grid,filename){
	  lambda=0;
	  g=0;
	m2=1;
	  maxtime=dt*100;
	  measuretime=1;
	   measenergy=false;
	  measpress=false;
	   xiint=this->randpi;
	  gamma=randpi/500;
//	  cout<<"MODEL KONSTR FILE"<<endl;
	}
	
	//Copykonstruktor �rt�kad�s felhaszn�l�s�val
	Model(const Model& model){
		*this=model;
//		cout<<"MODEL COPY"<<endl;
	}

	//a t�r lem�sol�sa
	Model& operator=(const Model& model){
		if(this!=&model){

			dt=model.dt;         //id�l�p�s egys�ge
			dr=model.dr;         //szakaszelem
			maxtime=model.maxtime;    //eddig megy az id�l�ptet�s
			m2=model.m2;          //t�megn�gyzet
			lambda=model.lambda;     //lambda
			g=model.g;
			measuretime=model.measuretime; //ennyi dt-nk�nt m�r�nk

			measenergy=model.measenergy;   //m�rj�nk-e �sszenergi�t
			measpress=model.measpress;    //m�rj�nk-e nyom�st

			randphi=model.randphi;	   //random intervallum
			randpi=model.randpi;
			discrete=model.discrete;
			xiint=model.xiint;
			gamma=model.gamma;

		}
		return *this;
	}

	//Laplace Phi sz�m�t�sa
	//visszat�r Laplace(Phi) �rt�k�vel
	double LaplacePhi(int x, int y, int z){
		double value;
		  int xindex=(x+discrete.xmax)%discrete.xmax;
		  int yindex=(y+discrete.ymax)%discrete.ymax;
		  int zindex=(z+discrete.zmax)%discrete.zmax;
		  int x_plus=(x+1+discrete.xmax)%discrete.xmax;
		  int x_minus=(x-1+discrete.xmax)%discrete.xmax;
		  int y_plus=(y+1+discrete.ymax)%discrete.ymax;
		  int y_minus=(y-1+discrete.ymax)%discrete.ymax;
		  int z_plus=(z+1+discrete.zmax)%discrete.zmax;
		  int z_minus=(z-1+discrete.zmax)%discrete.zmax;
		
			value=(discrete.phifield(x_plus,yindex,zindex)-2*discrete.phifield(xindex,yindex,zindex)+discrete.phifield(x_minus,yindex,zindex))/(discrete.a*discrete.a);
			value+=(discrete.phifield(xindex,y_plus,zindex)-2*discrete.phifield(xindex,yindex,zindex)+discrete.phifield(xindex,y_minus,zindex))/(discrete.a*discrete.a);
			value+=(discrete.phifield(xindex,yindex,z_plus)-2*discrete.phifield(xindex,yindex,zindex)+discrete.phifield(xindex,yindex,z_minus))/(discrete.a*discrete.a);
		
		return value;
	}

	//T12 sz�mol�sa
	double T12(int x, int y, int z){
		double value;
		double tempdiffx; //Phi(x,y,z) x ir�ny� deriv�ltja
		double tempdiffy; //Phi(x,y,z) y ir�ny� deriv�ltja
		  int xindex=(x+discrete.xmax)%discrete.xmax;
		  int yindex=(y+discrete.ymax)%discrete.ymax;
		  int zindex=(z+discrete.zmax)%discrete.zmax;
		  int x_plus=(x+1+discrete.xmax)%discrete.xmax;
		  int x_minus=(x-1+discrete.xmax)%discrete.xmax;
		  int y_plus=(y+1+discrete.ymax)%discrete.ymax;
		  int y_minus=(y-1+discrete.ymax)%discrete.ymax;
		  int z_plus=(z+1+discrete.zmax)%discrete.zmax;
		  int z_minus=(z-1+discrete.zmax)%discrete.zmax;

		  tempdiffx=1./(2*discrete.a)*(discrete.phifield(x_plus,yindex,zindex)-discrete.phifield(x_minus,yindex,zindex));
		  tempdiffy=1./(2*discrete.a)*(discrete.phifield(xindex,y_plus,zindex)-discrete.phifield(xindex,y_minus,zindex));
		  value=tempdiffx*tempdiffy;
		  return value;
	}


	//t�r friss�t�se - egy id�l�ptet�s elv�gz�se
	void RefreshField(){
		
		//els� l�p�s: Phi(to+dt) sz�mol�sa
		for(  int i=0;i<discrete.xmax;i++) //els� for
		{
			for(  int j=0;j<discrete.ymax;j++) //m�sodik for
			{
				for(  int k=0;k<discrete.zmax;k++) //harmadik for
				{
					/*Fi az adott helyen*/ discrete.phifield(i,j,k)=discrete.phifield(i,j,k)+dt*discrete.pifield(i,j,k);
			    }
			}

		}//els� for v�ge

		//ezek ut�n F(xi) sz�mol�sa �s p m�dos�t�sa
		for(  int i=0;i<discrete.xmax;i++) //els� for
		{
			for(  int j=0;j<discrete.ymax;j++) //m�sodik for
			{
				for(  int k=0;k<discrete.zmax;k++) //harmadik for
				{
					discrete.pifield(i,j,k)=discrete.pifield(i,j,k)+(LaplacePhi(i,j,k)-m2*discrete.phifield(i,j,k)-lambda*discrete.phifield(i,j,k)*discrete.phifield(i,j,k)*discrete.phifield(i,j,k)*1./6-g*discrete.phifield(i,j,k)*discrete.phifield(i,j,k)*discrete.phifield(i,j,k)*discrete.phifield(i,j,k)*discrete.phifield(i,j,k)*1./120)*dt;
			    }
			}

		}//els� for v�ge

	}

	//t�r friss�t�se zajjal terhelve
	void RefreshNoiseField(){
		double xi;
		//els� l�p�s: Phi(to+dt) sz�mol�sa
		for(  int i=0;i<discrete.xmax;i++) //els� for
		{
			for(  int j=0;j<discrete.ymax;j++) //m�sodik for
			{
				for(  int k=0;k<discrete.zmax;k++) //harmadik for
				{
					/*Fi az adott helyen*/ discrete.phifield(i,j,k)=discrete.phifield(i,j,k)+dt*discrete.pifield(i,j,k);
			    }
			}

		}//els� for v�ge

		//ezek ut�n F(xi) sz�mol�sa �s p m�dos�t�sa
		for(  int i=0;i<discrete.xmax;i++) //els� for
		{
			for(  int j=0;j<discrete.ymax;j++) //m�sodik for
			{
				for(  int k=0;k<discrete.zmax;k++) //harmadik for
				{
					xi=GetRandom(xiint);
					
					discrete.pifield(i,j,k)=discrete.pifield(i,j,k)+(LaplacePhi(i,j,k)-m2*discrete.phifield(i,j,k)-lambda*discrete.phifield(i,j,k)*discrete.phifield(i,j,k)*discrete.phifield(i,j,k)*1./6+xi-gamma*discrete.pifield(i,j,k)-g*discrete.phifield(i,j,k)*discrete.phifield(i,j,k)*discrete.phifield(i,j,k)*discrete.phifield(i,j,k)*discrete.phifield(i,j,k)*1./120)*dt;
			    }
			}

		}//els� for v�ge
	}


	//id�l�ptet�s v�ltoztat�sa
	void Setdt(double t){
		dt=t;
	}
	
	//szakaszelem v�ltoztat�sa
	void Setdr(int r){
		dr=r;
	}

	//r�cs�lland� v�ltoztat�sa
	void Setgrid(double a){
		discrete.a=a;
	}

	//t�meg v�ltoztat�sa
	void SetM2(double tomeg2){
		m2=tomeg2;
	}

	//lambda v�ltoztat�sa
	void Setlambda(double l){
		lambda=l;
	}

	//g valtoztatasa
	void Setg(double csatg){
		g=csatg;
	}

	void Setmaxtime(double max){
		maxtime=max;
	}

	//xiint be�ll�t�sa
	void SetXiint(double Xiint){
		this->xiint=Xiint;
	}

	//gamma be�ll�t�sa
	void SetGamma(double Gamma){
		this->gamma=Gamma;
	}

	//id�l�ptet�s lek�rdez�se
	double Getdt(){
		return dt;
	}

	//szakaszelem lek�r�se
	int Getdr(){
		return dr;
	}

	//xi random intervall lek�rdez�se
	double Getxiint(){
		return xiint;
	}

	//gamma lek�rdez�se
	double GetGamma(){
		return gamma;
	}

	//R�cs�lland�
	double Geta(){
		return this->discrete.GetA();
	}

	int GetXmax(){
		return discrete.GetX();
	}
	int GetYmax(){
		return discrete.GetY();
	}

	int GetZmax(){
		return discrete.GetZ();
	}

	

	array3<Complex>& GetphiFTfield(){
		return discrete.GetPhiFTfield();
	}

	array3<double>& GetPIField(){
		return discrete.GetPifield();
	}

	array3<Complex>& GetdataFFTfield(){
		return discrete.GetDataFFTfield();
	}

	Complex GetDataelement(unsigned x,unsigned y, unsigned z){
		return discrete.dataFFTfield(x,y,z);
	}

	double GetPhiElement(unsigned x,unsigned y, unsigned z){
		return discrete.phifield(x,y,z);
	}
	double GetPiElement(unsigned x,unsigned y, unsigned z){
		return discrete.pifield(x,y,z);
	}

	Complex GetphiFFTelement(unsigned x, unsigned y, unsigned z){
		return discrete.phiFFTfield(x,y,z);
	}
	Complex GetpiFFTelement(unsigned x,unsigned y, unsigned z){
		return discrete.piFFTfield(x,y,z);
	}

	size_t GeTAlign1(){
		return discrete.Getalign1();
	}

	size_t GeTAlign3(){
		return discrete.Getalign3();
	}

	//max id� lek�rdez�se
	double GetMaxtime(){
		return maxtime;
	}

	//m�r�sgyakoris�g lek�rdez�se
	int GetMeasTime(){
		return measuretime;
	}
	
	//r�cs�lland� lek�rdez�se
	double GetGrid(){
		return discrete.a;
	}
	
	//t�meg �s lambda lek�rdez�sei
	double GetM2(){
		return m2;
	}
	double GetL(){
		return lambda;
	}
	double GetG(){
		return g;
	}
	double Getrandphi(){
		return randphi;
	}
	double Getrandpi(){
		return randpi;
	}

	//direkt terek f�jlba �r�sa
	void ToFile(const char* filename)const{
		cout<<"Valos ter fajlba irasa..."<<endl;
		discrete.ToFileField(filename);
	}

	//FFT terek f�jlba �r�sa
	void ToFileFFT(const char* filename)const{
		cout<<"FFT ter fajlba irasa"<<endl;
		discrete.ToFileFFTField(filename);
	}

	//GradPhi^2 - seg�df�ggv�ny a sz�mol�sokhoz
	//visszat�r gradPhi abs.n�gyzet �rt�k�vel
	double GradPhi_sq(  int x,   int y,   int z){
		double xvalue;
		double yvalue;
		double zvalue;
		  int xindex=(x+discrete.xmax)%discrete.xmax;
		  int yindex=(y+discrete.ymax)%discrete.ymax;
		  int zindex=(z+discrete.zmax)%discrete.zmax;
		  int x_plus=(x+1+discrete.xmax)%discrete.xmax;
		  int y_plus=(y+1+discrete.ymax)%discrete.ymax;
		  int z_plus=(z+1+discrete.zmax)%discrete.zmax;

		//Forward differencia
		xvalue=1./(discrete.a)*(discrete.phifield(x_plus,yindex,zindex)-discrete.phifield(xindex,yindex,zindex));
		yvalue=1./(discrete.a)*(discrete.phifield(xindex,y_plus,zindex)-discrete.phifield(xindex,yindex,zindex));
		zvalue=1./(discrete.a)*(discrete.phifield(xindex,yindex,z_plus)-discrete.phifield(xindex,yindex,zindex));

		//gradphi^2
		double value=xvalue*xvalue+yvalue*yvalue+zvalue*zvalue;
        return value;

	}

	//Kisz�moljuk az �sszenergi�t
	double CountSumEnergy(){
		double energy; //ez a v�ltoz� fogja tartalmazni az �sszenergi�t
		double a=discrete.a;
		
		energy=0;

		for(int i=0;i<discrete.xmax;i++){
			for(int j=0;j<discrete.ymax;j++){
				for(int k=0;k<discrete.zmax;k++){
					//seg�df�ggv�ny h�v�sa
					energy+=(this->CountEnergy(i,j,k));
				}
			}
		}//for x v�ge
		
		return energy;
	}
	

	//Energia-tagok osszeget szamoljuk
	void CountSumEnergies(double (&value)[5]){
	 value[0]=0;
	 value[1]=0;
	 value[2]=0;
	 value[3]=0;
	 value[4]=0;
	 double a=discrete.a;
	 for(int i=0;i<discrete.xmax;i++){
                        for(int j=0;j<discrete.ymax;j++){
                                for(int k=0;k<discrete.zmax;k++){
                                        //seg�df�ggv�ny h�v�sa
                                        value[0]+=(this->CountPiEnergy(i,j,k)); //pi
					value[1]+=(this->CountGradEnergy(i,j,k)); //gradphi**2
					value[2]+=(this->CountSquareEnergy(i,j,k)); //phi**2
					value[3]+=(this->CountQuadrEnergy(i,j,k)); //phi**4
					value[4]+=(this->CountSixEnergy(i,j,k));  //phi**6
                                }
                        }
                }//for x v�ge


	}

	//Kisz�moljuk a teljes nyom�st
	double CountSumPressure(){
		double press; //ez a v�ltoz� fogja tartalmazni az �sszenergi�t
		double a=discrete.a;
		
		press=0;

		for(int i=0;i<discrete.xmax;i++){
			for(int j=0;j<discrete.ymax;j++){
				for(int k=0;k<discrete.zmax;k++){
					//seg�df�ggv�ny h�v�sa
					press+=(this->CountPressure(i,j,k));
				}
			}
		}//for x v�ge
		
		return press;
	}
	//T�r inicializ�l�s - random
	void RandomFill(){
		 srand((int)time(0)); 
		 for(int i=0;i<discrete.xmax;i++){
			 for(int j=0;j<discrete.ymax;j++){
				 for(int k=0;k<discrete.zmax;k++){

					 if((randphi!=0)&&(randpi!=0)){

						 discrete.phifield(i,j,k)=((double)rand()-(RAND_MAX/2))/(double)RAND_MAX*2*randphi; //random phi
						 discrete.pifield(i,j,k)=((double)rand()-(RAND_MAX/2))/(double)RAND_MAX*2*randpi;  //impulzus
					 }
					 else{
						 if(randpi==0){
							 discrete.phifield(i,j,k)=((double)rand()-(RAND_MAX/2))/(double)RAND_MAX*2*randphi; //random phi
							 discrete.pifield(i,j,k)=0;  //impulzus
						 }
						 if(randphi==0){
							 discrete.phifield(i,j,k)=0; //random phi
							 discrete.pifield(i,j,k)=((double)rand()-(RAND_MAX/2))/(double)RAND_MAX*2*randpi;  //impulzus
						 }
					 }
					  				 }
			 }
		 }//x for v�ge
		 //m�r nem �res a t�r
		discrete.isempty=false;
	}


	//T�r inicializ�l�s - random - Lorenz a=1
        void RandomLorenzFill(){
                 srand((int)time(0)); 
                 for(int i=0;i<discrete.xmax;i++){
                         for(int j=0;j<discrete.ymax;j++){
                                 for(int k=0;k<discrete.zmax;k++){

                                         if((randphi!=0)&&(randpi!=0)){
						double pilorenz;
						double philorenz;
                                                 philorenz=((double)rand()-(RAND_MAX/2))/(double)RAND_MAX*2*randphi; //impulzus
                                                 pilorenz=((double)rand()-(RAND_MAX/2))/(double)RAND_MAX*2*randpi;  //random phi
						philorenz=tan(M_PI*philorenz);
						pilorenz=tan(M_PI*pilorenz);
						discrete.phifield(i,j,k)=philorenz*exp(-abs(philorenz));
						discrete.pifield(i,j,k)=pilorenz*exp(-abs(pilorenz));
                                         }
                                         else{
                                                 if(randpi==0){
							 double philorenz;
							 philorenz=((double)rand()-(RAND_MAX/2))/(double)RAND_MAX*2*randphi;  //random phi
							 philorenz=tan(M_PI*philorenz);
							discrete.phifield(i,j,k)=philorenz*exp(-abs(philorenz));
                                                         discrete.pifield(i,j,k)=0;  //impulzus
                                                 }
                                                 if(randphi==0){
							double pilorenz;
							
							pilorenz=((double)rand())/(double)RAND_MAX; //impulzus
							 pilorenz=1/M_PI*log(tan(M_PI*pilorenz/2));
							
							
							discrete.pifield(i,j,k)=pilorenz;
							discrete.phifield(i,j,k)=0; //random phi
							
                                                 }
                                         }
                                                                         }
                         }
                 }//x for v�ge
                 //m�r nem �res a t�r
                discrete.isempty=false;
        }


	//T�r inicializ�l�s - sin �s cos fgv. alapj�n
	void SinCosFill(){
		
		for(  int i=0;i<discrete.xmax;i++){
			for(  int j=0;j<discrete.ymax;j++){
				for(  int k=0;k<discrete.zmax;k++){
					
					discrete.phifield(i,j,k)=(sin(2*M_PI*i/discrete.xmax));
					discrete.pifield(i,j,k)=(cos(2*M_PI*i/discrete.xmax));
				}
			}
		}//x for v�ge
		discrete.isempty=false;
	}


	//Param�terbe�ll�t�sok: r�cs�lland�, id�l�p�s, max id�, t�meg �s lambda, m�r�sgyakoris�g
	void ModelInit(){
		char swal;
		std::cout<<"Kerem a racsallandot, idolepest, max idot, tomegnegyzetet, lambdat es meresgyakorisagot es gt vesszovel elvalasztva:"<<std::endl;
		std::cin>>discrete.a>>swal>>dt>>swal>>maxtime>>swal>>m2>>swal>>lambda>>swal>>measuretime>>swal>>g;
	}

	//Mit m�rj�nk?
	void MeasInit(){
		char swal;
		std::cout<<"Merjek energiat es nyomast is?<y/n,y/n>"<<endl;
		std::cin>>swal;
		if(swal=='y'){ measenergy=true; }
		std::cin>>swal;
		std::cin>>swal;
		if(swal=='y'){ measpress=true; }
	}

	//Energia kisz�m�t�sa egy adott indexen
	double CountEnergy(int i, int j, int k){
		double a=discrete.a;
		double phithis=discrete.phifield(i,j,k);
		double energy=
		(1./2*discrete.pifield(i,j,k)*discrete.pifield(i,j,k)+1./2*this->GradPhi_sq(i,j,k)+m2*1./2*phithis*phithis+lambda*1./24*phithis*phithis*phithis*phithis+g/720*phithis*phithis*phithis*phithis*phithis*phithis)*a*a*a;
		return energy;
	}

	//energiatag: pi**2
	double CountPiEnergy(int i, int j, int k){
                double a=discrete.a;
                double energy=
                (1./2*discrete.pifield(i,j,k)*discrete.pifield(i,j,k))*a*a*a;
                return energy;
        }

	//energiatag: gradphi**2
	double CountGradEnergy(int i, int j, int k){
                double a=discrete.a;
                double phithis=discrete.phifield(i,j,k);
                double energy=
                (1./2*this->GradPhi_sq(i,j,k))*a*a*a;
                return energy;
        }

	//energiatag: phi**2
	double CountSquareEnergy(int i, int j, int k){
                double a=discrete.a;
                double phithis=discrete.phifield(i,j,k);
                double energy=
                (m2*1./2*phithis*phithis)*a*a*a;
                return energy;
        }

	//energiatag: phi**4
	double CountQuadrEnergy(int i, int j, int k){
                double a=discrete.a;
                double phithis=discrete.phifield(i,j,k);
                double energy=
                (lambda*1./24*phithis*phithis*phithis*phithis)*a*a*a;
                return energy;
        }

	//energiatag: phi**6
	double CountSixEnergy(int i, int j, int k){
                double a=discrete.a;
                double phithis=discrete.phifield(i,j,k);
                double energy=
                (g/720*phithis*phithis*phithis*phithis*phithis*phithis)*a*a*a;
                return energy;
        }
	//Nyom�s sz�m�t�sa
	double CountPressure(int i, int j, int k){
		double a=discrete.a;
		double phithis=discrete.phifield(i,j,k);

		double pressure=
			(1./2*discrete.pifield(i,j,k)*discrete.pifield(i,j,k)-1./6*this->GradPhi_sq(i,j,k)-m2*1./2*phithis*phithis-lambda*1./24*phithis*phithis*phithis*phithis)*a*a*a;
		return pressure;
	}

	//phik^2 sz�m�t�sa
	double CountPhik2(int i, int j, int k){
		double value=0;
		if(!discrete.isemptyFFT){
			const Complex phithisk=discrete.phiFFTfield(i,j,k);
			value=real(phithisk)*real(phithisk)+imag(phithisk)*imag(phithisk);	
		}
		return value;
	}


	//pik^2 sz�m�t�sa
	double CountPik2(int i, int j, int k){
		double value=0;
		if(!discrete.isemptyFFT){
			 const Complex pithisk=discrete.piFFTfield(i,j,k);
			value=(real(pithisk)*real(pithisk))+(imag(pithisk)*imag(pithisk));	
		}
		return value;
	}


	
	//Hisztogram m�r�sek
	//visszat�r�si �rt�k: 
	//vektor, melyben az egyes cell�kon m�rt �rt�k(ek) vannak - minden eleme vektor, mivel egy cell�nak t�bb jellemz�je is van
	vector<vector<double> > HistoMeas(){
		vector<vector<double> > value; //visszat�r�si �rt�k
		double dr3=dr*dr*dr;          //elemi cell�ban lev� r�cspontok
		
		//nagy ciklus dr egys�gekkel - cell�kon haladunk
		for(int X=0;X<=discrete.xmax-dr;X+=dr){
			for(int Y=0;Y<=discrete.ymax-dr;Y+=dr){
				for(int Z=0;Z<=discrete.zmax-dr;Z+=dr){

					double cellenergy=0; //cellaenergia
					double cellpressure=0; //cellanyom�s
				
					//bels� for a cell�kra
					for(int x=0;x<dr;x++){
						for(int y=0;y<dr;y++){
							for(int z=0;z<dr;z++){

								//ha m�r�nk cellaenergi�t
								if(this->measenergy){
									cellenergy+=this->CountEnergy((X+x),(Y+y),(Z+z)); 
								}

								//ha m�r�nk nyom�st
								if(this->measpress){
									cellpressure+=this->CountPressure((X+x),(Y+y),(Z+z));
								}

							}//z for
						}//y for
					}//x for

					vector<double> temp;
					temp.push_back(cellenergy);
					temp.push_back(cellpressure);
					temp.resize(2);
					value.push_back(temp);
				}//Z for
			}//Y for
		}//X for

		return value;
	}

	//phi v�rhat� �rt�k �s sz�r�sn�gyzet
	vector<double> MomentMeas(){

		double meanphi=0;
		double sigma2=0;
		vector<double> value;

		for(int i=0;i<discrete.xmax;i++){
			for(int j=0;j<discrete.ymax;j++){
				for(int k=0;k<discrete.zmax;k++){

					meanphi+=discrete.phifield(i,j,k);

				}//z for
			}//y for
		}//x for

		meanphi=meanphi/(discrete.xmax*discrete.ymax*discrete.zmax);

		for(int i=0;i<discrete.xmax;i++){
			for(int j=0;j<discrete.ymax;j++){
				for(int k=0;k<discrete.zmax;k++){

					sigma2+=(discrete.phifield(i,j,k)-meanphi)*(discrete.phifield(i,j,k)-meanphi);

				}//z for
			}//y for
		}//x for

		sigma2=sigma2/(discrete.xmax*discrete.ymax*discrete.zmax);

		value.push_back(meanphi);
		value.push_back(sigma2);
		return value;

	}

    //cellaenergia hisztogram dV f�gg�s
	//visszat�r�si �rt�k:
	//vektor: egyes cellafelbont�sra m�rt �rt�kek
	//de ezek is vektorok, mivel minden cell�hoz tartozik egy �rt�k
	vector<vector<double> > CellEnergyDV(){

	vector<vector<double> > value;

		//k�ls� ciklus cellam�ret n�vel�s�re
		for(int d=1;d<discrete.xmax;d+=1){
			dr=d;
			vector<double> temp;

		//nagy ciklus dr egys�gekkel - cell�kon haladunk
		for(int X=0;X<discrete.xmax-dr;X+=dr){
			for(int Y=0;Y<discrete.ymax-dr;Y+=dr){
				for(int Z=0;Z<discrete.zmax-dr;Z+=dr){
					
					double cellenergy=0;
					//bels� for a cella�tlagol�sra
					for(int x=0;x<dr;x++){
						for(int y=0;y<dr;y++){
							for(int z=0;z<dr;z++){

                                cellenergy+=this->CountEnergy((X+x),(Y+y),(Z+z));
						
							}//z for
						}//y for
					}//x for

					temp.push_back(cellenergy); //a vektor v�g�re tessz�k
				}//Z for
			}//Y for
		}//X for

		value.push_back(temp); //a vektor k�vetkez� eleme

		}

		return value;

	}

	//M�r�sek a k t�rben
	//Param�ter a 2D t�r referenci�ja, melyet kit�lt�nk �rt�kekkel
	//Ezt Inicializ�ltan kapja!!!
	//phikabs^2  �s pikabs^2 v�rhat� �rt�k k^2 f�ggv�ny�ben
    //�j javaslat phikabs^2 v�rhat� �rt�k eo-ra
	void phikfunctionmod(array2<double>& value, int N){
		
		int *kx=new int[discrete.xmax];
		int *ky=new int[discrete.ymax];
		double k2;
		int idx;
		double kmax=sqrt((double)discrete.xmax/2*discrete.xmax/2+(discrete.ymax/2*discrete.ymax/2)+(discrete.zmax/2*discrete.zmax/2));
		double phik2;
		double pik2;
		int xdim=discrete.xmax*discrete.ymax*(discrete.zmax/2+1);
	
		//indext�mb�k felt�lt�se
		int i=0;
		for(i;i<(discrete.xmax/2)+1;i++){
			kx[i]=i;
		}
		for(i;i<discrete.xmax;i++){
			kx[i]=(i-discrete.xmax);
		}

		i=0;
		for(i;i<(discrete.xmax/2)+1;i++){
			ky[i]=i;
		}
		for(i;i<discrete.xmax;i++){
			ky[i]=(i-discrete.xmax);
		}

		//kit�ltj�k az els� sort a k (hat�r)�rt�kekkel
		for(int n=0;n<xdim;n++){
			value(n,0)=kmax/N*(n+1);
		}

		//bej�rjuk a k teret
		
		for(i=0;i<discrete.xmax;i++){
			for(int j=0;j<discrete.ymax;j++){
				for(int kz=0;kz<(discrete.zmax/2)+1;kz++){

					k2=kx[i]*kx[i]+ky[j]*ky[j]+kz*kz;
					phik2=CountPhik2(i,j,kz);
					pik2=CountPik2(i,j,kz);

					if(sqrt(k2)<(kmax/(sqrt((double)2)))){
					//a megfelel� index
					idx=N*sqrt(k2)/kmax; 

				value(idx,1)+=phik2;
				value(idx,2)+=pik2;
				value(idx,3)++;
					}
				}//kzfor
			
			}//kyfor
	
		}//kxfor
		cout<<"ciklus kesz"<<endl;
		
		delete[] ky;
		delete[] kx;
		//return value;
	}

	//phikfunctionmod de nem k f�ggv�ny�ben, hanem sum_i(4sin^2(Pi k_i)/N) fgv�ben (itt N=xmax)
	//int N: s�vsz�m
	void phikfunctionsin(array2<double>& value, int N){
		
		int *kx=new int[discrete.xmax];
		int *ky=new int[discrete.ymax];
		double site;
		int idx;
		double sitemax=sin(M_PI/2)*sin(M_PI/2)+sin(M_PI/2)*sin(M_PI/2)+sin(M_PI/2)*sin(M_PI/2);
		double phik2;
		double pik2;
		int xdim=discrete.xmax*discrete.ymax*(discrete.zmax/2+1);
	
		//indext�mb�k felt�lt�se
		int i=0;
		for(i;i<(discrete.xmax/2)+1;i++){
			kx[i]=i;
		}
		for(i;i<discrete.xmax;i++){
			kx[i]=(i-discrete.xmax);
		}

		i=0;
		for(i;i<(discrete.xmax/2)+1;i++){
			ky[i]=i;
		}
		for(i;i<discrete.xmax;i++){
			ky[i]=(i-discrete.xmax);
		}

		//kit�ltj�k az els� sort a site (hat�r)�rt�kekkel
		for(int n=0;n<xdim;n++){
			value(n,0)=sitemax/N*(n+1);
		}

		//bej�rjuk a k teret
		
		for(i=0;i<discrete.xmax;i++){
			for(int j=0;j<discrete.ymax;j++){
				for(int kz=0;kz<(discrete.zmax/2)+1;kz++){

					site=(sin(2*M_PI*kx[i]/discrete.xmax)*sin(2*M_PI*kx[i]/discrete.xmax))+(sin(2*M_PI*ky[j]/discrete.xmax)*sin(2*M_PI*ky[j]/discrete.xmax))+(sin(2*M_PI*kz/discrete.xmax)*sin(2*M_PI*kz/discrete.xmax));
					phik2=CountPhik2(i,j,kz);
					pik2=CountPik2(i,j,kz);

					//if(sqrt(k2)<(kmax/(sqrt((double)2)))){
					//a megfelel� index
					idx=N*site/sitemax; 

				value(idx,1)+=phik2;
				value(idx,2)+=pik2;
				value(idx,3)++;
					//}
				}//kzfor
			
			}//kyfor
	
		}//kxfor
		cout<<"ciklus kesz"<<endl;
		
		delete[] ky;
		delete[] kx;
		//return value;
	}

	//H�m�rs�klet meghat�roz�sa
	double CountTemperature(){
	
		double value=0;
		for(int i=0;i<discrete.xmax;i++){
			for(int j=0;j<discrete.ymax;j++){
				for(int kz=0;kz<discrete.zmax/2+1;kz++){
				
				value+=this->CountPik2(i,j,kz);
				
				}//for k
			}//for j
		}//for i

		//�tlagol�s
		value=value/(discrete.xmax*discrete.ymax*(discrete.zmax/2+1));
		return value;
	}

	//Fourier transzform�ci�
	void FFTModel(){
	//	cout<<"FFT szamolasa..."<<endl;
		discrete.RunFFT();
	}

	void SetPiFieldElement(unsigned i, unsigned j, unsigned k, double element){
		discrete.pifield(i,j,k)=element;
	}

	

	//Destruktor
	~Model(){//cout<<"MODEL DESTR"<<endl;
	}

};


class EvolveModel{
	//adattagjai
	//array4<double> phifield;
	//array4<double> pifield;
	double dV;
	double dr;
//	vector<Model> Modelset1;
//	vector<Model> Modelset2;
	vector<Model> *Modelset1;
	
	int pointnumber;
	int* meastime2;
	double zerotime;
	int measincrease;
	int measmin;

	//f�ggv�nyei
public:
	//konstruktor
	/*EvolveModel(Model& model, int timedim, int xmax, int ymax, int zmax):phifield(xmax,ymax,zmax,timedim),pifield(xmax,ymax,zmax,timedim){
		for(unsigned timestep=0;timestep<phifield.N4();timestep++){

			//adott id�re a 3D alt�r felt�lt�se
			for(unsigned i=0;i<phifield.Nx();i++){
				for(unsigned j=0;j<phifield.Ny();j++){
					for(unsigned k=0;k<phifield.Nz();k++){
						phifield(i,j,k,timestep)=model.discrete.phifield(i,j,k);
						pifield(i,j,k,timestep)=model.discrete.pifield(i,j,k);
					}//z for
				}//y for
			}//x for

			//model id�l�ptet�se
			model.RefreshField();

		}//model id�l�ptet�s v�ge

		dV=model.discrete.a*model.discrete.a*model.discrete.a*model.Getdt();
		dr=model.discrete.a;
	}*/

	
	//timedim: ennyi dt-nk�nti als�k van
	//Konstruktor
	EvolveModel(Model& model, double zerotime, int pointnumber,int measincrease,int measmin,bool noise):zerotime(zerotime),pointnumber(pointnumber),measincrease(measincrease),measmin(measmin){
		
		//pointerfoglal�sok ah�ny t-re +1 (t'-nek)
		meastime2= new int[pointnumber+1];
		Modelset1= new vector<Model>[pointnumber+1];

		//m�r�sgyakoris�gok �s bel�l�k r�gt�n counter
		meastime2[0]=0; //counter, a t�bbi counteri
		for(int i=0;i<(pointnumber);i++){
			meastime2[i+1]=0-(measmin+(i*measincrease));
		}//for m�r�sid�k felt�lt�se


		//id�fejl�d�s �s felt�lt�s
		//int counter=0;
		//int counter2=counter-meastime2;
		double dt=model.Getdt();
		double maxtime=model.GetMaxtime();
		int modmeastime=model.GetMeasTime();

	//DEBUG	std::ofstream myfile;
	//DEBUG	myfile.open("D:/doc/phi4modelinfo.dat",std::ios::out);
		for(double timer=zerotime;timer<maxtime;timer+=dt){

			//csak m�r�sgyakoris�gonk�nt t�lt�nk fel
			for(int i=0;i<(pointnumber+1);i++){
				//ellen�rizni kell, hogy a counter m�r el�rte-e a null�t
				if(meastime2[i]>=0){
					if((meastime2[i]%modmeastime)==0){
						
				Modelset1[i].push_back(model);
				cout<<"modell import "<<i<<": "<<Modelset1[i].size()<<"timer: "<<timer<<endl;
		//DEBUG		myfile<<"modell import "<<i<<": "<<Modelset1[i].size()<<"timer: "<<timer<<endl;
					}//ha counteri m�r�s van
				}//ha counteri>=0
			}//for valyon counteri >=0?
			cout<<"Frissites ciklus: "<<timer<<endl;
			//counterek n�vel�se
			for(int i=0;i<(pointnumber+1);i++){
				meastime2[i]+=1;
			}//counterek n�vel�se
			
			//fejleszt�nk
		if(noise){
			model.RefreshNoiseField();
		}
		else{
			model.RefreshField();	
		}
			
		}
		cout<<"EvolveModel feltoltve"<<endl;
		//DEBUG myfile.close();
		//m�retez�s
		//for(int i=0;i<(pointnumber+1);i++){
		//	Modelset1[i].shrink_to_fit();
		//}//m�retez�s
		
		//ezt sz�mol�s hely�n fogjuk ellen�rizni
		/*if(Modelset1.size()>Modelset2.size()){

			Modelset1.resize(Modelset2.size());
		}
		else if(Modelset2.size()>Modelset1.size()){
			Modelset2.resize(Modelset1.size());
		}*/
		
		dr=model.Geta();
		dV=dr*dr*dr*model.Getdt();
	}

/*	double Viscosity(){
		double value=0;
		double tempval=0;
		double pisquare=0;
		//k�ls� for 4D
		for(int Time=0; Time<timedim;Time++){
			for(int I=0; I<Modelset[Time].discrete.xmax;I++){
				for(int J=0;J<Modelset[Time].discrete.ymax;J++){
					for(int K=0;K<Modelset[Time].discrete.zmax;K++){

						//bels� for
						double Ttx=this->CountTtx(I,J,K,Time);
						double Tcoma;
						for(int t=0;t<timedim;t++){
							for(int i=0;i<Modelset[Time].discrete.xmax;i++){
								for(int j=0;j<Modelset[Time].discrete.ymax;j++){
									for(int k=0;k<Modelset[Time].discrete.zmax;k++){
										Tcoma=this->CountTtx(i,j,k,t);
										tempval+=Ttx*Tcoma*dV;

									}//bels� z for
								}//bels� y for
							}//bels� x for
						}//bels� t for

						value+=tempval*dV;
						pisquare+=Modelset[Time].discrete.pifield(I,J,K)*Modelset[Time].discrete.pifield(I,J,K)*dV;

						cout<<"Belso ciklus kesz."<<std::endl;
						

					}//k�ls� z for
					cout<<"Belso z-re kesz"<<std::endl;
				}//K�ls� y for
				cout<<"Belso y-re kesz"<<std::endl;
			}//k�ls� x for
			cout<<Time<<"*dt "<<" kesz."<<std::endl;
		}//k�ls� t for

		value=value/pisquare;

		return value;
	}*/

	//visszat�r�s: els� sor: phik2k t=0
	//m�sodik sor: Gtk val�s r�sz
	//harmadik sor: gtk im. r�sz
	//negyedik sor: 4sumsin^2
	//value legyen null�ra inicializ�lva
	void GTK(array3<double>& value){
		int counter=0;//phik2k �tlagol�shoz
		int xmax=Modelset1[0][0].GetXmax();
		int ymax=Modelset1[0][0].GetYmax();
		int zmax=Modelset1[0][0].GetZmax();
		int *kx=new int[xmax];
		int *ky=new int[ymax];
		int i=0;
		int idx=0;
		Complex temp;
		int limit=Modelset1[0].size();//minden ablakra �ll�tani fogjuk
		int meas1=Modelset1[0][0].GetMeasTime();
		double dt=Modelset1[0][0].Getdt();

		for(i;i<xmax/2+1;i++){
			kx[i]=i;
		}
		for(i;i<xmax;i++){
			kx[i]=(i-xmax);
		}

		i=0;
		for(i;i<(ymax/2)+1;i++){
			ky[i]=i;
		}
		for(i;i<ymax;i++){
			ky[i]=(i-ymax);
		}
		
		i=0;

		for(int timecounter=0;timecounter<limit;timecounter++){
			//els�sz�r fourier teret �ll�tunk el�
			Modelset1[0][timecounter].FFTModel();

		}//timefor v�ge
		
		//adatment�shez tempfile-t hozunk l�tre
		std::ofstream tempfile;
		tempfile.precision(6);
		tempfile.open("gtk_prop_temp.dat");
		cout<<"Letrehozzuk az ideiglenes fajlt, adatmentes celjabol."<<endl;

		//itt kezd�dik amit minden t el kell v�gezni
		for(int l=1;l<(pointnumber+1);l++){

			//limit helyes be�ll�t�sa: amelyik a kisebb
			if(Modelset1[l].size()<Modelset1[0].size()){
				limit=Modelset1[l].size();
			}
			else{
				limit=Modelset1[0].size();
			}

		for(int timecounter=0;timecounter<limit;timecounter++){
			//els�sz�r fourier teret �ll�tunk el�
			Modelset1[l][timecounter].FFTModel();
		}//timefor v�ge

		cout<<"FFT terek eloalltak, Ablak: "<<l<<endl;

		
		//bej�rjuk a k teret �s sz�moljuk sum:phi(t+t',k)*phi(*t',k)-t
		//valamint phik2k-t
		i=0;
		for(i;i<xmax;i++){
			for(int j=0;j<ymax;j++){
				for(int kz=0;kz<(zmax/2)+1;kz++){
				
					// ha l�tez� index
					if(idx<(int)value.Ny()){
					//integr�lunk t' szerint
					for(int timecounter=0;timecounter<limit;timecounter++){
						
						temp=Modelset1[l][timecounter].GetphiFFTelement(i,j,kz)*conj(Modelset1[0][timecounter].GetphiFFTelement(i,j,kz))*(meas1*dt);
						value(l-1,idx,1)+=real(temp);
						value(l-1,idx,2)+=imag(temp);
						//id��tlagolva
						value(l-1,idx,0)+=(Modelset1[l][timecounter].CountPhik2(i,j,kz));
						counter+=1;
						
					}
					value(l-1,idx,0)/=counter;
					value(l-1,idx,1)/=(counter*dt*meas1);
					value(l-1,idx,2)/=(counter*dt*meas1);
					counter=0;
					
					

					value(l-1,idx,3)+=(4*sin(M_PI*kx[i]/xmax)*sin(M_PI*kx[i]/xmax))+(4*sin(M_PI*ky[j]/xmax)*sin(M_PI*ky[j]/xmax))+(4*sin(M_PI*kz/xmax)*sin(M_PI*kz/xmax));
					tempfile<<i<<'\t'<<j<<'\t'<<kz<<'\t'<<value(l-1,idx,0)<<'\t'<<value(l-1,idx,1)<<'\t'<<value(l-1,idx,2)<<'\t'<<value(l-1,idx,3)<<(l-1)<<'\t'<<endl;
					
					cout<<"t' integralas volt adott k helyen. Index: "<<idx<<"Ablak: "<<l<<endl;
					idx++;
					}
					else{
						cout<<"HIBAS INDEX"<<endl;
					}

					}//kzfor
			
			}//kyfor
	
		}//kxfor
		cout<<l<<"-1 ablakra kesz."<<endl;
		idx=0;
	}//minden t ablakra
		tempfile.close();
	}

	void phi4(array3<double>& value){
		int counter=0;//phik2k �tlagol�shoz
		int xmax=Modelset1[0][0].GetXmax();
		int ymax=Modelset1[0][0].GetYmax();
		int zmax=Modelset1[0][0].GetZmax();
		int *kx=new int[xmax];
		int *ky=new int[ymax];
		int i=0;
		int idx=0;
		Complex temp;
		int limit=Modelset1[0].size();//minden ablakra �ll�tani fogjuk
		int meas1=Modelset1[0][0].GetMeasTime();
		double dt=Modelset1[0][0].Getdt();

		for(i;i<(xmax/2)+1;i++){
			kx[i]=i;
		}
		for(i;i<xmax;i++){
			kx[i]=(i-xmax);
		}

		i=0;
		for(i;i<(ymax/2)+1;i++){
			ky[i]=i;
		}
		for(i;i<ymax;i++){
			ky[i]=(i-ymax);
		}
		
		i=0;

		
		//0. index (t' Modelt�mb) phi2 �rt�kei
		array3<double> phi2array0(xmax,ymax,zmax,Modelset1[0][0].GeTAlign1());

		//0. index modelljeire el��ll�tjuk phi^2 FFT tereket
		for(int timecounter=0;timecounter<limit;timecounter++){
			
			//phi2 el��ll�t�sa
			for(int I=0;I<Modelset1[0][0].discrete.xmax;I++){
				for (int J=0;J<Modelset1[0][0].discrete.ymax;J++){
					for (int K=0;K<Modelset1[0][0].discrete.zmax;K++){

						phi2array0[I][J][K]=Modelset1[0][timecounter].GetPhiElement(I,J,K)*Modelset1[0][timecounter].GetPhiElement(I,J,K);

					}//Kfor
				}//Jfor
			}//Ifor

			//Fourier transzform�ci�
			// �s igen, belepiszk�lunk, mert phiFFTfield-et haszn�ljuk Fourier t�rnek
			rcfft3d ForwardPhi(Modelset1[0][timecounter].GetZmax(),phi2array0,Modelset1[0][timecounter].GetphiFTfield());
			ForwardPhi.fft(phi2array0,Modelset1[0][timecounter].GetphiFTfield());
			//((Modelset1[0][timecounter].discrete.Forward3Phi)).fft(phi2array0,Modelset1[0][timecounter].discrete.phiFFTfield);
				
		}//timefor v�ge
		
		//adatment�shez tempfile-t hozunk l�tre
		std::ofstream tempfile;
		tempfile.precision(6);
		tempfile.open("gtk_prop_temp2.dat");
		cout<<"Letrehozzuk az ideiglenes fajlt, adatmentes celjabol."<<endl;

		//itt kezd�dik amit minden t el kell v�gezni
		for(int l=1;l<(pointnumber+1);l++){

			//limit helyes be�ll�t�sa: amelyik a kisebb
			if(Modelset1[l].size()<Modelset1[0].size()){
				limit=Modelset1[l].size();
			}
			else{
				limit=Modelset1[0].size();
			}

		for(int timecounter=0;timecounter<limit;timecounter++){
			//els�sz�r fourier teret �ll�tunk el�
			//phi2 el��ll�t�sa
			for(int I=0;I<xmax;I++){
				for (int J=0;J<ymax;J++){
					for (int K=0;K<zmax;K++){
						
						phi2array0[I][J][K]=Modelset1[l][timecounter].GetPhiElement(I,J,K)*Modelset1[l][timecounter].GetPhiElement(I,J,K);

					}//Kfor
				}//Jfor
			}//Ifor

			//Fourier transzform�ci�
			// �s igen, belepiszk�lunk, mert phiFFTfield-et haszn�ljuk Fourier t�rnek
			rcfft3d ForwardPhi(Modelset1[l][timecounter].GetZmax(),phi2array0,Modelset1[l][timecounter].GetphiFTfield());
			ForwardPhi.fft(phi2array0,Modelset1[l][timecounter].GetphiFTfield());
		
			//((Modelset1[l][timecounter].discrete.Forward3Phi)).fft(phi2array0,Modelset1[l][timecounter].discrete.phiFFTfield);
		}//timefor v�ge

		cout<<"FFT terek eloalltak, Ablak: "<<l<<endl;

		
		//bej�rjuk a k teret �s sz�moljuk sum:phi(t+t',k)*phi(*t',k)-t
		//valamint phik2k-t
		i=0;
		for(i;i<xmax;i++){
			for(int j=0;j<ymax;j++){
				for(int kz=0;kz<(zmax/2)+1;kz++){
				
					// ha l�tez� index
					if(idx<(int)value.Ny()){
					//integr�lunk t' szerint
					for(int timecounter=0;timecounter<limit;timecounter++){
						
						temp=Modelset1[l][timecounter].GetphiFFTelement(i,j,kz)*conj(Modelset1[0][timecounter].GetphiFFTelement(i,j,kz))*(meas1*dt);
						value(l-1,idx,1)+=real(temp);
						value(l-1,idx,2)+=imag(temp);
						//id��tlagolva
						value(l-1,idx,0)+=(Modelset1[l][timecounter].CountPhik2(i,j,kz));
						counter+=1;
						
					}
					value(l-1,idx,0)/=counter;
					value(l-1,idx,1)/=(counter*dt*meas1);
					value(l-1,idx,2)/=(counter*dt*meas1);
					counter=0;
					
					value(l-1,idx,3)+=(4*sin(M_PI*kx[i]/xmax)*sin(M_PI*kx[i]/xmax))+(4*sin(M_PI*ky[j]/ymax)*sin(M_PI*ky[j]/ymax))+(4*sin(M_PI*kz/zmax)*sin(M_PI*kz/zmax));
					

					tempfile<<i<<'\t'<<j<<'\t'<<kz<<'\t'<<value(l-1,idx,0)<<'\t'<<value(l-1,idx,1)<<'\t'<<value(l-1,idx,2)<<'\t'<<value(l-1,idx,3)<<(l-1)<<'\t'<<endl;
					
					cout<<"t' integralas volt adott k helyen. Index: "<<idx<<"Ablak: "<<l<<endl;
					idx++;
					}
					else{
						cout<<"HIBAS INDEX"<<endl;
					}

					}//kzfor
			
			}//kyfor
	
		}//kxfor
		cout<<l<<"-1 ablakra kesz."<<endl;
		idx=0;
	}//minden t ablakra

		tempfile.close();

	}

	void visco(array3<double>& value){
		int counter=0;//phik2k �tlagol�shoz
		int xmax=Modelset1[0][0].GetXmax();
		int ymax=Modelset1[0][0].GetYmax();
		int zmax=Modelset1[0][0].GetZmax();
		int *kx=new int[xmax];
		int *ky=new int[ymax];
		int i=0;
		int idx=0;
		Complex temp;
		int limit=Modelset1[0].size();//minden ablakra �ll�tani fogjuk
		int meas1=Modelset1[0][0].GetMeasTime();
		double dt=Modelset1[0][0].Getdt();

		for(i;i<(xmax/2)+1;i++){
			kx[i]=i;
		}
		for(i;i<xmax;i++){
			kx[i]=(i-xmax);
		}

		i=0;
		for(i;i<(ymax/2)+1;i++){
			ky[i]=i;
		}
		for(i;i<ymax;i++){
			ky[i]=(i-ymax);
		}
		
		i=0;

		
		//0. index (t' Modelt�mb) T12 �rt�kei
		array3<double> T12array0(xmax,ymax,zmax,Modelset1[0][0].GeTAlign3());

		//0. index modelljeire el��ll�tjuk T12 FFT tereket
		for(int timecounter=0;timecounter<limit;timecounter++){
			
			//T12 t�r el��ll�t�sa
			for(int I=0;I<xmax;I++){
				for (int J=0;J<ymax;J++){
					for (int K=0;K<zmax;K++){
						
						T12array0[I][J][K]=Modelset1[0][timecounter].T12(I,J,K);

					}//Kfor
				}//Jfor
			}//Ifor

			//Fourier transzform�ci�
			// �s igen MOD:NEM!!!(, belepiszk�lunk, mert phiFFTfield-et haszn�ljuk Fourier t�rnek)
			rcfft3d ForwardPhi(Modelset1[0][timecounter].GetZmax(),T12array0,Modelset1[0][timecounter].GetdataFFTfield());
			ForwardPhi.fft(T12array0,Modelset1[0][timecounter].GetdataFFTfield());
			
			//((Modelset1[0][timecounter].discrete.Forward3Phi)).fft(T12array0,Modelset1[0][timecounter].discrete.phiFFTfield);
				
		}//timefor v�ge
		
		//adatment�shez tempfile-t hozunk l�tre
		std::ofstream tempfile;
		tempfile.precision(6);
		tempfile.open("gtk_prop_tempT12.dat");
		cout<<"Letrehozzuk az ideiglenes fajlt, adatmentes celjabol."<<endl;

		//itt kezd�dik amit minden t el kell v�gezni
		for(int l=1;l<(pointnumber+1);l++){

			//limit helyes be�ll�t�sa: amelyik a kisebb
			if(Modelset1[l].size()<Modelset1[0].size()){
				limit=Modelset1[l].size();
			}
			else{
				limit=Modelset1[0].size();
			}

		for(int timecounter=0;timecounter<limit;timecounter++){
			//els�sz�r fourier teret �ll�tunk el�
			//T12 el��ll�t�sa
			for(int I=0;I<xmax;I++){
				for (int J=0;J<ymax;J++){
					for (int K=0;K<zmax;K++){

						T12array0[I][J][K]=Modelset1[l][timecounter].T12(I,J,K);

					}//Kfor
				}//Jfor
			}//Ifor

			//Fourier transzform�ci�
			// �s igen MOD: NEM!!!(, belepiszk�lunk, mert phiFFTfield-et haszn�ljuk Fourier t�rnek)
			rcfft3d ForwardPhi(Modelset1[l][timecounter].GetZmax(),T12array0,Modelset1[l][timecounter].GetdataFFTfield());
			ForwardPhi.fft(T12array0,Modelset1[l][timecounter].GetdataFFTfield());
			
			//((Modelset1[l][timecounter].discrete.Forward3Phi)).fft(T12array0,Modelset1[l][timecounter].discrete.phiFFTfield);
		}//timefor v�ge

		cout<<"FFT terek eloalltak, Ablak: "<<l<<endl;

		
		//bej�rjuk a k teret �s sz�moljuk sum:T12(t+t',k)*T12(*t',k)-t
		//valamint phik2k-t
		i=0;
		for(i;i<xmax;i++){
			for(int j=0;j<ymax;j++){
				for(int kz=0;kz<(zmax/2)+1;kz++){
				
					// ha l�tez� index
					if(idx<(int)value.Ny()){
					//integr�lunk t' szerint
					for(int timecounter=0;timecounter<limit;timecounter++){

						temp=Modelset1[l][timecounter].GetDataelement(i,j,kz)*conj(Modelset1[0][timecounter].GetDataelement(i,j,kz))*(meas1*dt);
						value(l-1,idx,1)+=real(temp);
						value(l-1,idx,2)+=imag(temp);
						//id��tlagolva
					//Mostm�r nem lesz j�, dataFFT elemeinek absz�rtn�gyzete kell	value(l-1,idx,0)+=(Modelset1[l][timecounter].CountPhik2(i,j,kz));
						//Ink�bb null�zzuk
						value(l-1,idx,0)+=0;
						counter+=1;
						
					}
					value(l-1,idx,0)/=counter;
					value(l-1,idx,1)/=(counter*dt*meas1);
					value(l-1,idx,2)/=(counter*dt*meas1);
					counter=0;
					
					value(l-1,idx,3)+=(4*sin(M_PI*kx[i]/xmax)*sin(M_PI*kx[i]/xmax))+(4*sin(M_PI*ky[j]/ymax)*sin(M_PI*ky[j]/ymax))+(4*sin(M_PI*kz/zmax)*sin(M_PI*kz/zmax));
					

					tempfile<<i<<'\t'<<j<<'\t'<<kz<<'\t'<<value(l-1,idx,0)<<'\t'<<value(l-1,idx,1)<<'\t'<<value(l-1,idx,2)<<'\t'<<value(l-1,idx,3)<<(l-1)<<'\t'<<endl;
					
					cout<<"t' integralas volt adott k helyen. Index: "<<idx<<"Ablak: "<<l<<endl;
					idx++;
					}
					else{
						cout<<"HIBAS INDEX"<<endl;
					}

					}//kzfor
			
			}//kyfor
	
		}//kxfor
		cout<<l<<"-1 ablakra kesz."<<endl;
		idx=0;
	}//minden t ablakra

		tempfile.close();
	}


	void evolveCombo(array3<double>& value){
		int counter=0;//phik2k �tlagol�shoz
		int xmax=Modelset1[0][0].GetXmax();
		int ymax=Modelset1[0][0].GetYmax();
		int zmax=Modelset1[0][0].GetZmax();
		int *kx=new int[xmax];
		int *ky=new int[ymax];
		int i=0;
		int idx=0;
		Complex temp;
		Complex temp2;
		int limit=Modelset1[0].size();//minden ablakra �ll�tani fogjuk
		int meas1=Modelset1[0][0].GetMeasTime();
		double dt=Modelset1[0][0].Getdt();

		for(i;i<(xmax/2)+1;i++){
			kx[i]=i;
		}
		for(i;i<xmax;i++){
			kx[i]=(i-xmax);
		}

		i=0;
		for(i;i<(ymax/2)+1;i++){
			ky[i]=i;
		}
		for(i;i<ymax;i++){
			ky[i]=(i-ymax);
		}
		
		i=0;

		
		//0. index (t' Modelt�mb) T12 �rt�kei
		array3<double> T12array0(xmax,ymax,zmax,Modelset1[0][0].GeTAlign3());

		//0. index modelljeire el��ll�tjuk T12 FFT tereket
		for(int timecounter=0;timecounter<limit;timecounter++){
			
			//T12 t�r el��ll�t�sa
			for(int I=0;I<xmax;I++){
				for (int J=0;J<ymax;J++){
					for (int K=0;K<zmax;K++){
						
						T12array0[I][J][K]=Modelset1[0][timecounter].T12(I,J,K);

					}//Kfor
				}//Jfor
			}//Ifor

			//Fourier transzform�ci�
			// �s igen MOD:NEM!!!(, belepiszk�lunk, mert phiFFTfield-et haszn�ljuk Fourier t�rnek)
			rcfft3d ForwardPhi(Modelset1[0][timecounter].GetZmax(),T12array0,Modelset1[0][timecounter].GetdataFFTfield());
			ForwardPhi.fft(T12array0,Modelset1[0][timecounter].GetdataFFTfield());
			
			//els�sz�r fourier teret �ll�tunk el�
			Modelset1[0][timecounter].FFTModel();




			//((Modelset1[0][timecounter].discrete.Forward3Phi)).fft(T12array0,Modelset1[0][timecounter].discrete.phiFFTfield);
				
		}//timefor v�ge
		
		//adatment�shez tempfile-t hozunk l�tre
		/*std::ofstream tempfile;
		tempfile.precision(6);
		tempfile.open("gtk_prop_tempT12.dat");
		cout<<"Letrehozzuk az ideiglenes fajlt, adatmentes celjabol."<<endl;*/

		//itt kezd�dik amit minden t el kell v�gezni
		for(int l=1;l<(pointnumber+1);l++){

			//limit helyes be�ll�t�sa: amelyik a kisebb
			if(Modelset1[l].size()<Modelset1[0].size()){
				limit=Modelset1[l].size();
			}
			else{
				limit=Modelset1[0].size();
			}

		for(int timecounter=0;timecounter<limit;timecounter++){
			//els�sz�r fourier teret �ll�tunk el�
			//T12 el��ll�t�sa
			for(int I=0;I<xmax;I++){
				for (int J=0;J<ymax;J++){
					for (int K=0;K<zmax;K++){

						T12array0[I][J][K]=Modelset1[l][timecounter].T12(I,J,K);

					}//Kfor
				}//Jfor
			}//Ifor

			//Fourier transzform�ci�
			// �s igen MOD: NEM!!!(, belepiszk�lunk, mert phiFFTfield-et haszn�ljuk Fourier t�rnek)
			rcfft3d ForwardPhi(Modelset1[l][timecounter].GetZmax(),T12array0,Modelset1[l][timecounter].GetdataFFTfield());
			ForwardPhi.fft(T12array0,Modelset1[l][timecounter].GetdataFFTfield());
			
			Modelset1[l][timecounter].FFTModel();


			//((Modelset1[l][timecounter].discrete.Forward3Phi)).fft(T12array0,Modelset1[l][timecounter].discrete.phiFFTfield);
		}//timefor v�ge

		cout<<"FFT terek eloalltak, Ablak: "<<l<<endl;

		
		//bej�rjuk a k teret �s sz�moljuk sum:T12(t+t',k)*T12(*t',k)-t
		//valamint phik2k-t
		i=0;
		for(i;i<xmax;i++){
			for(int j=0;j<ymax;j++){
				for(int kz=0;kz<(zmax/2)+1;kz++){
				
					// ha l�tez� index
					if(idx<(int)value.Ny()){
					//integr�lunk t' szerint
					for(int timecounter=0;timecounter<limit;timecounter++){

						temp=Modelset1[l][timecounter].GetDataelement(i,j,kz)*conj(Modelset1[0][timecounter].GetDataelement(i,j,kz))*(meas1*dt);
						value(l-1,idx,1)+=real(temp);
						value(l-1,idx,2)+=imag(temp);

						temp2=Modelset1[l][timecounter].GetphiFFTelement(i,j,kz)*conj(Modelset1[0][timecounter].GetphiFFTelement(i,j,kz))*(meas1*dt);
						value(l-1,idx,4)+=real(temp2);
						value(l-1,idx,5)+=imag(temp2);
						//id��tlagolva
						value(l-1,idx,0)+=(Modelset1[l][timecounter].CountPhik2(i,j,kz)); //Phikabs2

						//id��tlagolva
					//Mostm�r nem lesz j�, dataFFT elemeinek absz�rtn�gyzete kell	value(l-1,idx,0)+=(Modelset1[l][timecounter].CountPhik2(i,j,kz));
						
						counter+=1;
						
					}
					value(l-1,idx,0)/=counter;
					value(l-1,idx,1)/=(counter*dt*meas1);
					value(l-1,idx,2)/=(counter*dt*meas1);
					value(l-1,idx,4)/=(counter*dt*meas1);
					value(l-1,idx,5)/=(counter*dt*meas1);
					counter=0;
					
					value(l-1,idx,3)+=(4*sin(M_PI*kx[i]/xmax)*sin(M_PI*kx[i]/xmax))+(4*sin(M_PI*ky[j]/ymax)*sin(M_PI*ky[j]/ymax))+(4*sin(M_PI*kz/zmax)*sin(M_PI*kz/zmax));
					

					//tempfile<<i<<'\t'<<j<<'\t'<<kz<<'\t'<<value(l-1,idx,0)<<'\t'<<value(l-1,idx,1)<<'\t'<<value(l-1,idx,2)<<'\t'<<value(l-1,idx,3)<<(l-1)<<'\t'<<endl;
					
					cout<<"t' integralas volt adott k helyen. Index: "<<idx<<"Ablak: "<<l<<endl;
					idx++;
					}
					else{
						cout<<"HIBAS INDEX"<<endl;
					}

					}//kzfor
			
			}//kyfor
	
		}//kxfor
		cout<<l<<"-1 ablakra kesz."<<endl;
		idx=0;
	}//minden t ablakra

		//tempfile.close();
	}

/*	double CountTtx(int x, int y, int z, int time){
		double Ttx=0;
		double xvalue=0;
		double yvalue=0;
		Model& tempmodel=Modelset[time];

		int xindex=(x+tempmodel.discrete.xmax)%tempmodel.discrete.xmax;
		  int yindex=(y+tempmodel.discrete.ymax)%tempmodel.discrete.ymax;
		  int zindex=(z+tempmodel.discrete.zmax)%tempmodel.discrete.zmax;
		  int x_plus=(x+1+tempmodel.discrete.xmax)%tempmodel.discrete.xmax;
		  int x_minus=(x-1+tempmodel.discrete.xmax)%tempmodel.discrete.xmax;
		  int y_plus=(y+1+tempmodel.discrete.ymax)%tempmodel.discrete.ymax;
		  int y_minus=(y-1+tempmodel.discrete.ymax)%tempmodel.discrete.ymax;
		  int z_plus=(z+1+tempmodel.discrete.zmax)%tempmodel.discrete.zmax;
		  int z_minus=(z-1+tempmodel.discrete.zmax)%tempmodel.discrete.zmax;

		  //K�zponti differencia
		  xvalue=1./(2*dr)*(tempmodel.discrete.phifield(x_plus,yindex,zindex)-tempmodel.discrete.phifield(x_minus,yindex,zindex));
		yvalue=1./(2*dr)*(tempmodel.discrete.phifield(xindex,y_plus,zindex)-tempmodel.discrete.phifield(xindex,y_minus,zindex));
		return xvalue*yvalue;

	}*/
	~EvolveModel(){
		delete[] meastime2;
		delete[] Modelset1;
	}

};

	
#endif
