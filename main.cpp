//Memory check
#ifdef _MSC_VER
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#endif

#include<iostream>
#include<fstream>
#include<string.h>
#include<vector>
#include<math.h>
#include <cstdlib>
#include <ctime> 
#include <string.h>
#include <iomanip>

#include"Phi4Base.h"
#include"Array.h"
#include"fftw++.h"


using namespace std;
using namespace Array;
using namespace fftwpp;

//SVN Próba

//Tesztesetek


//FFT tesztelése
void test_1(){
	string dir;
	Model phi1(0.1,10,10,10,1,1,0);

	//célkönyvtár megadása
	cout<<"Kerem a celkonyvtarat (win: vegen perjellel): "<<endl;
	cin>>dir;

	phi1.SinCosFill();
	phi1.ToFile((dir+"mivanbennem_init.dat").c_str());

	phi1.FFTModel();
	phi1.ToFileFFT((dir+"FFT.dat").c_str());
}

//Random kezdõállapot tesztelése
void test_2(){
	string dir;

	cout<<"Kerem a celkonyvtarat (win: vegen perjellel): "<<endl;
	cin>>dir;

	Model phi1(0.1,50,50,50,1,1,0);
	phi1.RandomFill();

	phi1.ToFile((dir+"mivanbennem_init.dat").c_str());
	phi1.FFTModel();
	phi1.ToFileFFT((dir+"FFT.dat").c_str());
}

//Fájlból feltöltés tesztje
void test_3(){
	
	const char* readin="D:\\doc\\8félév\\Szakdoga\\Mérések\\TesztFájl\\beolvasasutanirtam.dat";

	string dir;
	cout<<"Kerem a celkonyvtarat (win: vegen perjellel): "<<endl;
	cin>>dir;

	//egy teret random feltöltünk
	cout<<"Model letrehozasa es random feltoltese."<<endl;
	Model phi1(0.1,50,50,50,1,1,0);
	phi1.RandomFill();
	phi1.ToFile((dir+"kiirtam.dat").c_str());
	double sum0=phi1.CountSumEnergy();
	cout<<"Az osszenergia:\n"<<sum0<<endl;
	//ebbe fogjuk beolvasni
	cout<<"Masodik modell letrehozasa az elozobol kiirt adatokbol."<<endl;
	Model phi2(0.1,50,50,50,1,1,0,(dir+"kiirtam.dat").c_str());
	phi2.ToFile((dir+"beolvasasutanirtam.dat").c_str());
	double sum=phi2.CountSumEnergy();
	cout<<"Az osszenergia a masodik modellbol: \n"<<sum<<endl;
	cout<<"varok egy entert."<<endl;
	getchar();
}

//dV lefutás, két fájl: kezdeti idõpontra és idõfejlõdés végén
void test_4(){

	//elõzõekhez hasonlóan
	double duration;
	double timestep;
	int meastime;
	int counter=0;
	string dir;

	cout<<"Kerem a celkonyvtarat (win: vegen perjellel): "<<endl;
	cin>>dir;

	
	Model phi1(0.1,50,50,50,1,1,0);
	const char* outfilename=(dir+"SumEnergy.dat").c_str();
	const char* outdV=(dir+"initialdV.dat").c_str();
	const char* outdVfinal=(dir+"finaldV.dat").c_str();
	const char* outinit=(dir+"mivanbennem_init.dat").c_str();
	const char* outinitFFT=(dir+"mivanbennemFFT_init.dat").c_str();
	const char* outfinal=(dir+"mivanbennem_evolved.dat").c_str();
	const char* outfinalFFT=(dir+"mivanbennemFFT_evolved.dat").c_str();
	const char* infofile=(dir+"info.dat").c_str();

	//inicializálás
	phi1.RandomFill();
	cout<<"Kezdeti ter fajlba irasa..."<<endl;
	phi1.ToFile((dir+"mivanbennem_init.dat").c_str());
	cout<<"...kesz."<<endl;
	phi1.ModelInit();
	
	phi1.FFTModel();
	cout<<"Kezdeti FFT ter fajlba irasa..."<<endl;
	phi1.ToFileFFT((dir+"mivanbennemFFT_init.dat").c_str());
	cout<<"...kesz."<<endl;
	duration=phi1.GetMaxtime();
	timestep=phi1.Getdt();
	meastime=phi1.GetMeasTime();


	//szükséges mennyiségek
	double sumenergy;
	vector<vector<double> > rundV; //cellaméret függvényében energiahisztogram

	//a fejlõdés ciklus elején és végén kiírjuk az összenergiát
	std::ofstream myfile1;
	myfile1.precision(6);
	myfile1.open((dir+"SumEnergy.dat").c_str(),std::ios::out);
	sumenergy=phi1.CountSumEnergy();
	myfile1<<"0"<<"\t"<<sumenergy<<std::endl;
	
    //a fejlõdés elején kiírjuk az energiahisztogram dV lefutását
	std::ofstream myfile2;
	myfile2.precision(6);
	myfile2.open((dir+"initialdV.dat").c_str(),std::ios::out);
	
	//itt számoljuk a lényeget - még az idõfejlõdés elõtt
	cout<<"Idofejlodes elotti meres..."<<endl;
	rundV=phi1.CellEnergyDV();

	//fájlba írás
	cout<<"Eredmenyek fajlba irasa..."<<endl;
	for (vector<vector<double> >::iterator it = rundV.begin(); it!=rundV.end(); it++) {
			vector<double> temp=*it;
			for(vector<double>::iterator it2 =temp.begin(); it2!=temp.end();it2++){
				double v=*it2;
				myfile2<<fixed<<v<<"\t";
			}
			myfile2<<endl;
			
		}//iterátor vége
	

	//ide jön az idõfejlõdés és kiírás
	for(double i=0;i<duration;i+=timestep,counter++){
		
		cout<<i<<endl;
		
		cout<<"Refreshing..."<<endl;
		//fejlesztünk
		phi1.RefreshField();
		
		
    }
	//itt véget ér az idõfejlõdés és kiírás

	//végeredmények kiírása
	phi1.ToFile((dir+"mivanbennem_evolved.dat").c_str());
	phi1.FFTModel();
	phi1.ToFileFFT((dir+"mivanbennemFFT_evolved.dat").c_str());
	sumenergy=phi1.CountSumEnergy();
	myfile1<<duration<<"\t"<<sumenergy<<std::endl;

	std::ofstream myfile3;
	myfile3.precision(6);
	myfile3.open((dir+"finaldV.dat").c_str(),std::ios::out);

	//itt a lényeg számolása - idõfejlõdés végén
	cout<<"Meres az idofejlodes vegen..."<<endl;
	rundV=phi1.CellEnergyDV();
	
	cout<<"Eredmenyek fajlba irasa..."<<endl;
	for (vector<vector<double> >::iterator it = rundV.begin(); it!=rundV.end(); it++) {
			vector<double> temp=*it;
			for(vector<double>::iterator it2 =temp.begin(); it2!=temp.end();it2++){
				double v=*it2;
				myfile3<<fixed<<v<<"\t";
			}
			myfile3<<endl;
			
		}//iterátor vége

	//eddig nyitott fájlok zárása
	myfile1.close();
	myfile2.close();
	myfile3.close();
	

	//info fájl elkészítése
	std::ofstream myfile;
	myfile.open ((dir+"info.dat").c_str(),std::ios::out);
	 myfile<<"Racsallando: "<<phi1.GetGrid()<<"\n"<<"Idolepes: "<<timestep<<"\n"<<"Max ido: "<<duration<<"\n"<<"Tomegnegyzet: "<<phi1.GetM2()<<"\n"
		 <<"Lambda: "<<phi1.GetL()<<"\n"<<"Randphi: "<<phi1.Getrandphi()<<"\n"<<"Randpi: "<<phi1.Getrandpi()<<endl;
	 myfile.close();
}


//szimpla mérések - fájlból
void simpletest(string opendir,string dir,double randphi,double randpi){

	//hogy ne kelljen mindig lekérni a modellbõl
	
	int dr;          //szakaszelem
	int band=1;
	
	//Modell létrehozása
	Model phi1(0.1,50,50,50,1,randphi,randpi,opendir.c_str());
	//kimeneti fájlnevek
	
	char swal;
	double grid,tomeg2,lambda,g;
	char sumen,histo,kpill,momentum,temperature,kteruj,pres;

	cout<<"Kerem a racsallandot, tomegnegyzetet es lambdat es csat gt vesszovel elvalasztva: "<<endl;
	cin>>grid>>swal>>tomeg2>>swal>>lambda>>g;

	phi1.Setgrid(grid);
	phi1.SetM2(tomeg2);
	phi1.Setlambda(lambda);
	phi1.Setg(g);

	cout<<"Mit merjek? Osszenergia, hisztogram, k terben, momentumok, homerseklet, k ter uj,nyomas: <y/n,...,y/n>"<<endl;
	cin>>sumen>>swal>>histo>>swal>>kpill>>swal>>momentum>>swal>>temperature>>swal>>kteruj>>swal>>pres;

	if(histo=='y'){
		cout<<"Kerem a szakaszelemet: "<<endl;
		cin>>dr;
		phi1.Setdr(dr);
		phi1.MeasInit();
	}
	if((kpill=='y')||(kteruj=='y')){
		
        cout<<"Hany sav legyen?"<<endl;
		cin>>band;
	}

	//kezdeti terek kiírása
	cout<<"A kezdeti teret nem irom ismet fajlba."<<endl;
	phi1.FFTModel();
	phi1.ToFileFFT((dir+"mivanbennemFFT_init.dat").c_str());
	cout<<"...kesz."<<endl;

	cout<<"A meresinicializalas kesz. Kezdem a munkat."<<endl;

	//info fájl nyitása
	std::ofstream myfile;
	myfile.open((dir+"info.dat").c_str(),std::ios::out);
	myfile<<"Racsallando: "<<grid<<endl;
	myfile<<"Tomegnegyzet: "<<tomeg2<<endl;
	myfile<<"Lambda: "<<lambda<<endl;
	myfile<<"Randphi: "<<phi1.Getrandphi()<<endl;
	myfile<<"Randpi: "<<phi1.Getrandpi()<<endl;
	myfile<<"csat g: "<<phi1.GetG()<<endl;
	//kimenet nyitása
	std::ofstream myfile1;
	std::ofstream myfile2;
	std::ofstream myfile3;
	myfile1.precision(6);
	myfile1.open((dir+"measdata.dat").c_str(),std::ios::out);
	if(sumen=='y'){
		double sumenergy=phi1.CountSumEnergy();
		myfile1<<"Osszenergia: "<<sumenergy<<endl;
		myfile<<"Mertem osszenergiat."<<endl;
	}
	if(pres=='y'){
		double sumpres=phi1.CountSumPressure();
		myfile1<<"Ossznyomás: "<<sumpres<<endl;
		myfile<<"Mertem ossznyomast."<<endl;
	}
	if(momentum=='y'){
		vector<double> value=phi1.MomentMeas();
		myfile1<<"Phi atlag: "<<value[0]<<endl;
		myfile1<<"Phi 2. momentum: "<<value[1]<<endl;
		myfile<<"Mertem phi atlagot es szorast"<<endl;
	}
	//ha mérünk hõmérsékletet
	if(temperature=='y'){
		double temp=phi1.CountTemperature();
		myfile1<<"Homerseklet: "<<temp<<endl;
		myfile<<"Mertem homersekletet."<<endl;
	}
	myfile1.close();

	if(histo=='y'){
		
		phi1.FFTModel();
		vector<vector<double> > histomeas=phi1.HistoMeas();
		
		myfile1.precision(6);
		myfile1.open((dir+"Energy_hist.dat").c_str(),std::ios::out);

		
		myfile2.precision(6);
		myfile2.open((dir+"Pressure_hist.dat").c_str(),std::ios::out);

		cout<<"Hisztogram adatok fajlba irasa..."<<endl;
		for (vector<vector<double> >::iterator it = histomeas.begin(); it!=histomeas.end(); it++) {
			vector<double> temp=*it;
			myfile1<<fixed<<temp[0]<<"\t";
			myfile2<<fixed<<temp[1]<<"\t";


		}//iterátor vége
			myfile1<<endl;
			myfile2<<endl;
			myfile1.close();
			myfile2.close();
			myfile<<"Mertem hisztogramot."<<endl;
			myfile<<"Szakaszelem: "<<dr<<endl;
	}//hiszto vége
	
	//ha mérünk k térben
	if(kpill=='y'){

		int xdim=phi1.GetXmax()*phi1.GetYmax()*(phi1.GetZmax()/2+1);

		array2<double> kmeas(xdim,4,sizeof(Complex));
	//kmeas-t nullára inicializáljuk
		for(int i=0;i<xdim;i++){
			for(int j=0;j<4;j++){
				kmeas(i,j)=0;
			}//yfor
		}//xfor

		phi1.phikfunctionmod(kmeas,band);

		myfile1.precision(6);
		myfile1.open((dir+"phik2k.dat").c_str(),std::ios::out);

		myfile2.precision(6);
		myfile2.open((dir+"pik2k.dat").c_str(),std::ios::out);

		for(int i=0;i<xdim;i++){
		if(kmeas(i,3)!=0){
			kmeas(i,1)=kmeas(i,1)/kmeas(i,3);
			kmeas(i,2)=kmeas(i,2)/kmeas(i,3);
			myfile1<<kmeas(i,0)<<"\t"<<kmeas(i,1)<<endl;
			myfile2<<kmeas(i,0)<<"\t"<<kmeas(i,2)<<endl;
		}
	}//kiírás vége

		myfile1.close();
		myfile2.close();

		myfile<<"Mertem a k terben."<<endl;
		myfile<<"k savok szama: "<<band<<endl;
	}//mérés ktérben

	//ha mérünk k térben új függvénnyel
	if(kteruj=='y'){
		int xdim=phi1.GetXmax()*phi1.GetYmax()*(phi1.GetZmax()/2+1);

		array2<double> kmeas(xdim,4,sizeof(Complex));
	//kmeas-t nullára inicializáljuk
		for(int i=0;i<xdim;i++){
			for(int j=0;j<4;j++){
				kmeas(i,j)=0;
			}//yfor
		}//xfor

		phi1.phikfunctionsin(kmeas,band);

		myfile1.precision(6);
		myfile1.open((dir+"phik2ku.dat").c_str(),std::ios::out);

		myfile2.precision(6);
		myfile2.open((dir+"pik2ku.dat").c_str(),std::ios::out);

		for(int i=0;i<xdim;i++){
		if(kmeas(i,3)!=0){
			kmeas(i,1)=kmeas(i,1)/kmeas(i,3);
			kmeas(i,2)=kmeas(i,2)/kmeas(i,3);
			myfile1<<kmeas(i,0)<<"\t"<<kmeas(i,1)<<endl;
			myfile2<<kmeas(i,0)<<"\t"<<kmeas(i,2)<<endl;
		}
	}//kiírás vége

		myfile1.close();
		myfile2.close();

		myfile<<"Mertem a k terben sin fuggvenyben."<<endl;
		myfile<<"savok szama: "<<band<<endl;
	}//mérés ktérben

	myfile.close();
	

}

//szimpla mérés fájlból vagy random
void test_5(){
	
	string opendir;
	string dir;
	double randphi;
	double randpi;

	cout<<"Udvozol a Szimpla Meres Program!\n"<<endl;
	cout<<"Ha random kezdeti feltetelt szeretnel, ird: RANDOM"<<endl;
	cout<<"Kulonben ird be a beolvasando fajlnevet eleresi uttal egyutt."<<endl;
	cout<<"Ide varom a valaszt: ";
	cin>>opendir;
	cout<<"Kerem a celkonyvtarat (win: vegen perjellel): "<<endl;
	cin>>dir;

	//ha random kezdeti feltételt szeretnénk
	if(!(opendir.compare("RANDOM"))){
		cout<<"Mi legyen randphi?"<<endl;
		cin>>randphi;
		cout<<"Mi legyen randpi?"<<endl;
		cin>>randpi;
		Model phirandom(0.1,50,50,50,1,randphi,randpi);
		phirandom.RandomFill();
		cout<<"Fajlba irom a generalt teret..."<<endl;
		phirandom.ToFile((dir+"randomgeninit.dat").c_str());
		cout<<"...fajlba iras kesz."<<endl;
		opendir=(dir+"randomgeninit.dat");
	}

	//mostantól jön a mérés.
	simpletest(opendir,dir,randphi,randpi);


}

//idõfejlõdéses mérés - fájlból
void evoltest(string opendir, string dir, double randphi, double randpi,bool noise){
	int dr;          //szakaszelem
	int band=1;
	
	//Modell létrehozása
	Model phi1(0.1,50,50,50,1,randphi,randpi,opendir.c_str());
	const char* infofile=(dir+"info.dat").c_str();
	
	char swal;
	double dt,maxtime,avgtime,starttime,zerotime;
	int meastime,counter=0;
	char sumen,histo,kpill,momentum,temperature,histoavg,kteruj,pres,iswrite;
	phi1.ModelInit();

	dt=phi1.Getdt();
	maxtime=phi1.GetMaxtime();
	meastime=phi1.GetMeasTime();

	cout<<"Mit merjek? Osszenergia,hisztogram,k terben, momentumok, homerseklet, atlagolt hisztogram, kterben sinfgvben,nyomas: <y/n,...,y/n>"<<endl;
	cin>>sumen>>swal>>histo>>swal>>kpill>>swal>>momentum>>swal>>temperature>>swal>>histoavg>>swal>>kteruj>>swal>>pres;

	cout<<"Kiirjam-e a tereket meresidonkent?: <y/n>"<<endl;
	cin>>iswrite;

	cout<<"Ha korabbi merest folytatunk, mi a kezdeti idopont? (Random esetben 0)"<<endl;
	cin>>zerotime;

	cout<<"Mikor kezdodjon a meres? "<<endl;
	cin>>starttime;

	

	//info fájl nyitása
	std::ofstream myfile;
	myfile.open((dir+"info.dat").c_str(),std::ios::out);
	myfile<<"Racsallando: "<<phi1.GetGrid()<<endl;
	myfile<<"Tomegnegyzet: "<<phi1.GetM2()<<endl;
	myfile<<"Lambda: "<<phi1.GetL()<<endl;
	myfile<<"Idolepes: "<<phi1.Getdt()<<endl;
	myfile<<"Meresgyakorisag: "<<phi1.GetMeasTime();
	myfile<<"Meres vege: "<<phi1.GetMaxtime()<<endl;
	myfile<<"Randphi: "<<phi1.Getrandphi()<<endl;
	myfile<<"Randpi: "<<phi1.Getrandpi()<<endl;
	myfile<<"csat g: "<<phi1.GetG()<<endl;

	if(iswrite=='y'){
	myfile<<"kiirom a tereket is!"<<endl;
	}


	if(noise){
		double xi,gamma;
	cout<<"Mennyi legyen xi es gamma?"<<endl;
	cin>>xi>>swal>>gamma;
	phi1.SetGamma(gamma);
	phi1.SetXiint(xi);
	myfile<<"Zajos meres!!!"<<endl;
	myfile<<"xi: "<<phi1.Getxiint()<<endl;
	myfile<<"gamma: "<<phi1.GetGamma()<<endl;
	}


	if(histoavg=='y'){
	cout<<"Mikor kezdjunk atlagolni? (Mereskezdes vagy nagyobb)"<<endl;
	cin>>avgtime;
	myfile<<"Atlagolas kezdete: "<<avgtime<<endl;
	}
	else{
		avgtime=maxtime;
	}

	if(histo=='y' || histoavg=='y'){
	cout<<"Kerem a szakaszelemet: "<<endl;
				cin>>dr;
				phi1.Setdr(dr);
				phi1.MeasInit();
				myfile<<"Szakaszelem: "<<dr<<endl;
	}
	if(kpill=='y'|| kteruj=='y'){
	cout<<"Hany sav legyen?"<<endl;
	cin>>band;
	myfile<<"savok szama: "<<band<<endl;
	}

	cout<<"A meresinicializalas kesz, kezdem a munkat."<<endl;

	//Fájlba írjuk miket mérünk
	if(sumen=='y'){
		myfile<<"Merem az osszenergia idofuggeset."<<endl;
	}
	if(histo=='y'){
		myfile<<"Merek hisztogramot minden meresidoben."<<endl;
	}
	if(kpill=='y'){
		myfile<<"Merem a k terbeli fuggest minden meresidoben."<<endl;
	}
	if(momentum=='y'){
		myfile<<"Merem phi varhato erteket es szorasat az ido fuggvenyeben."<<endl;
	}
	if(temperature=='y'){
		myfile<<"Merem a homerseklet idofuggeset."<<endl;
	}
	if(histoavg=='y'){
		myfile<<"Merek idoatlagolt hisztogramot."<<endl;
	}
	
	if(kteruj=='y'){
		myfile<<"Merek k terben sin^2 fuggvenyeben."<<endl;
	}
	if(pres=='y'){
		myfile<<"Merek ossznyomast az ido fuggvenyeben."<<endl;
	}


	//kezdeti terek kiírása
	cout<<"A kezdeti teret nem irom ismet fajlba."<<endl;
	phi1.FFTModel();
	cout<<"Kezdeti FFT ter fajlba irasa..."<<endl;
	phi1.ToFileFFT((dir+"mivanbennemFFT_init.dat").c_str());
	cout<<"...kesz."<<endl;

	//kimeneti fájlok nyitása

	//idõ-energia
	std::ofstream myfile1;
	myfile1.precision(6);
	myfile1.open((dir+"Time_SumEnergy.dat").c_str());

	//idõ-nyomás
	std::ofstream myfile1B;
	myfile1B.precision(6);
	myfile1B.open((dir+"Time_SumPressure.dat").c_str());

	//idõ-hõmérséklet
	std::ofstream myfile2;
	myfile2.precision(6);
	myfile2.open((dir+"Ido_Homerseklet.dat").c_str());

	//idõ-momentumok
	std::ofstream myfile3;
	myfile3.precision(6);
	myfile3.open((dir+"time_phimean.dat").c_str());

	std::ofstream myfile4;
	myfile4.precision(6);
	myfile4.open((dir+"time_phi2nd.dat").c_str());
	
	//idõ-hisztogram
	std::ofstream myfile5;
	myfile5.precision(6);
	myfile5.open((dir+"Time_Energy_hist.dat").c_str());

	std::ofstream myfile6;
	myfile6.precision(6);
	myfile6.open((dir+"Time_Pressure_hist.dat").c_str());

	//idõ - ktérben
	std::ofstream myfile7;
	myfile7.precision(6);
	myfile7.open((dir+"time_phik2k.dat").c_str());

	std::ofstream myfile8;
	myfile8.precision(6);
	myfile8.open((dir+"time_pik2k.dat").c_str());

	//idõátlagolt mérésekhez
	std::ofstream myfile10;
	myfile10.precision(6);
	myfile10.open((dir+"TimeAVG_Energy_hist.dat").c_str());

	std::ofstream myfile11;
	myfile11.precision(6);
	myfile11.open((dir+"TimeAVG_Pressure_hist.dat").c_str());

	std::ofstream myfile15;
	myfile15.precision(6);
	myfile15.open((dir+"time_phik2ku.dat").c_str());

	std::ofstream myfile16;
	myfile16.precision(6);
	myfile16.open((dir+"time_pik2ku.dat").c_str());

	std::ofstream myfile21;
	myfile21.precision(6);
	myfile21.open((dir+"kteruj_variance.dat").c_str());


	//a mérésekhez szükséges:
	double sumenergy,temper,sumpres;
	double energies[5];
	int xdim=phi1.GetXmax()*phi1.GetYmax()*(phi1.GetZmax()/2+1);
	array2<double> kmeas(xdim,4,sizeof(Complex));
	//kmeas-t nullára inicializáljuk
		for(int i=0;i<xdim;i++){
			for(int j=0;j<4;j++){
				kmeas(i,j)=0;
			}//yfor
		}//xfor

	array2<double> kmeasuj(xdim,4,sizeof(Complex));
	//kmeas-t nullára inicializáljuk
		for(int i=0;i<xdim;i++){
			for(int j=0;j<4;j++){
				kmeasuj(i,j)=0;
			}//yfor
		}//xfor


	//idõbeli szórás számolásához
	array2<double> kmeasujvar(xdim,2,sizeof(Complex));
	//kmeas-t nullára inicializáljuk
		for(int i=0;i<xdim;i++){
			for(int j=0;j<2;j++){
				kmeasujvar(i,j)=0;
			}//yfor
		}//xfor
	int count=0;

    vector<vector<double> > histomeas;
	vector<double> mom;
	//mérések és idõfejlõdés


	for(double i=zerotime;i<maxtime;i+=dt){
		
		cout<<i<<endl;
		
		//ha mérésidõ van
		if(i>=starttime){
		if(counter%meastime==0){
			cout<<"Meres!"<<endl;
			phi1.FFTModel();

			if(iswrite=='y'){
			stringstream ss;
			ss<< counter;
			string outnum=ss.str();
			cout<<"terek kiirasa..."<<outnum<<endl;
			 phi1.ToFile((dir+"valoster_evolved"+outnum+".dat").c_str());
			cout<<"...sikeres"<<endl;

			}
			if(sumen=='y'){
			sumenergy=phi1.CountSumEnergy();
			phi1.CountSumEnergies(energies);
			myfile1<<i<<"\t"<<sumenergy<<"\t"<<energies[0]<<"\t"<<energies[1]<<"\t"<<energies[2]<<"\t"<<energies[3]<<"\t"<<energies[4]<<endl;
			}
			if(pres=='y'){
				sumpres=phi1.CountSumPressure();
				myfile1B<<i<<"\t"<<sumpres<<endl;
			}
			if(temperature=='y'){
			temper=phi1.CountTemperature();
			myfile2<<i<<"\t"<<temper<<endl;
			}
			if(momentum=='y'){
				mom=phi1.MomentMeas();
				myfile3<<i<<"\t"<<mom[0]<<endl;
				myfile4<<i<<"\t"<<mom[1]<<endl;
			}
			if(kpill=='y'){
				phi1.phikfunctionmod(kmeas,band);

				cout<<"k terbeli adatok fajlba irasa..."<<endl;

				for(int j=0;j<xdim;j++){
					if(kmeas(j,3)!=0){
						kmeas(j,1)=kmeas(j,1)/kmeas(j,3);
						kmeas(j,2)=kmeas(j,2)/kmeas(j,3);
						myfile7<<kmeas(j,0)<<"\t"<<kmeas(j,1)<<";";
						myfile8<<kmeas(j,0)<<"\t"<<kmeas(j,2)<<";";
					}
				}//forj
				myfile7<<endl;
				myfile8<<endl;
				
			
			//kmeas-t újra inicializáljuk
			for(int l=0;l<xdim;l++){
				for(int z=0;z<4;z++){
					kmeas(l,z)=0;
				}//yfor
			}//xfor
		}//kpill

		if(kteruj=='y'){
				phi1.phikfunctionsin(kmeasuj,band);

				cout<<"k terbeli adatok fajlba irasa..."<<endl;

				for(int j=0;j<xdim;j++){
					if(kmeasuj(j,3)!=0){
						kmeasuj(j,1)=kmeasuj(j,1)/kmeasuj(j,3);
						kmeasuj(j,2)=kmeasuj(j,2)/kmeasuj(j,3);
						myfile15<<kmeasuj(j,0)<<"\t"<<kmeasuj(j,1)<<";";
						myfile16<<kmeasuj(j,0)<<"\t"<<kmeasuj(j,2)<<";";

						//szórásszámoláshoz
						kmeasujvar(j,0)+=kmeasuj(j,1);
						kmeasujvar(j,1)+=kmeasuj(j,1)*kmeasuj(j,1);
						
						count++;

					}
				}//forj
				myfile15<<endl;
				myfile16<<endl;
				
			//kmeas-t újra inicializáljuk
			for(int l=0;l<xdim;l++){
				for(int z=0;z<4;z++){
					kmeasuj(l,z)=0;
				}//yfor
			}//xfor
		}//kpill

			if(histo=='y'){
				
				histomeas=phi1.HistoMeas();

				cout<<"Hisztogram adatok fajlba irasa..."<<endl;
				for (vector<vector<double> >::iterator it = histomeas.begin(); it!=histomeas.end(); it++) {
					vector<double> temp=*it;
					myfile5<<fixed<<temp[0]<<"\t";
					myfile6<<fixed<<temp[1]<<"\t";
					if(histoavg=='y' && i>avgtime){
						myfile10<<fixed<<temp[0]<<"\t";
						myfile11<<fixed<<temp[1]<<"\t";
					}
			}//iterátor vége
			myfile5<<i<<endl;
			myfile6<<i<<endl;
			}//histo
			
			//ha csak idõátlagoltak vannak
			if(i>avgtime){
				if(histo!='y' && histoavg=='y'){
					histomeas=phi1.HistoMeas();

					cout<<"Hisztogram adatok fajlba irasa..."<<endl;
				for (vector<vector<double> >::iterator it = histomeas.begin(); it!=histomeas.end(); it++) {
					vector<double> temp=*it;
					myfile10<<fixed<<temp[0]<<"\t";
					myfile11<<fixed<<temp[1]<<"\t";
				}//iterátor vége
				}
				
			}

			
		}//ha mérésidõ vége
		counter++;
	}
		cout<<"Refreshing..."<<endl;
		
		
		//fejlesztünk
		if(noise){
			phi1.RefreshNoiseField();
		}
		else{
			phi1.RefreshField();	
		}
    }
	myfile10<<endl;
	myfile11<<endl;

	//kteruj szórásának fájlbaírása
	if(kteruj=='y'){
		double temp;
		for(int j=0;j<xdim;j++){
			kmeasujvar(j,0)/=count;
			kmeasujvar(j,1)/=count;
			temp=kmeasujvar(j,1)-kmeasujvar(j,0)*kmeasujvar(j,0);
			myfile21<<kmeasujvar(j,0)<<"\t"<<kmeasujvar(j,1)<<"\t"<<temp<<endl;
						
		}//forj
	}

	phi1.ToFile((dir+"valoster_evolved.dat").c_str());
	phi1.ToFileFFT((dir+"reciprokter_evolved.dat").c_str());


	myfile.close();
	myfile1.close();
	myfile1B.close();
	myfile2.close();
	myfile3.close();
	myfile4.close();
	myfile5.close();
	myfile6.close();
	myfile7.close();
	myfile8.close();
	myfile10.close();
	myfile11.close();
	myfile15.close();
	myfile16.close();
	myfile21.close();
}

//idõfejlõdéses mérések
void test_6(){
	string opendir;
	string dir;
	double randphi;
	double randpi;
	bool noise=false;
	char noi;
	char noisetype;

	cout<<"Udvozol az Idofejlodeses Meres Program!\n"<<endl;
	cout<<"Ha random LORENZ kezdeti feltetelt szeretnel, ird: RANDOM"<<endl;
	cout<<"Kulonben ird be a beolvasando fajlnevet eleresi uttal egyutt."<<endl;
	cout<<"Ide varom a valaszt: ";
	cin>>opendir;
	cout<<"Kerem a celkonyvtarat (win: vegen perjellel): "<<endl;
	cin>>dir;

	cout<<"Zajos meres? <y,n>"<<endl;
	cin>>noi;

	if(noi=='y'){
		noise=true;
	}

	cout<<"Mi legyen randphi?"<<endl;
	cin>>randphi;
	cout<<"Mi legyen randpi?"<<endl;
	cin>>randpi;

	//ha random kezdeti feltételt szeretnénk
	if(!(opendir.compare("RANDOM"))){
		cout<<"Zajtipus hiperbol vagy uniform? <h,u>"<<endl;
                cin>>noisetype;
		Model phirandom(0.1,50,50,50,1,randphi,randpi);
		if(noisetype=='h'){
                phirandom.RandomLorenzFill();
                }
                else if(noisetype=='u'){
                phirandom.RandomFill();
                }
		cout<<"Fajlba irom a generalt teret..."<<endl;
		phirandom.ToFile((dir+"randomgeninit.dat").c_str());
		cout<<"...fajlba iras kesz."<<endl;
		opendir=(dir+"randomgeninit.dat");
	}

	//mostantól jön a mérés.
	evoltest(opendir,dir,randphi,randpi,noise);
}


//randomint hõmérséklet összefüggés mérése
void test_7(){
	string opendir;
	string dir;
	double randphi=0;
	

	cout<<"Udvozol az Idofejlodeses Meres Program!\n"<<endl;
	

	cout<<"Randphi erteke nulla, randpi 4-tol indul 20-ig, felenkent."<<endl;
	
	for(double randpi=0.2;randpi<=20;randpi+=0.2){

		cout<<"Kerem a celkonyvtarat (win: vegen perjellel): "<<endl;
		cin>>dir;
		Model phirandom(0.1,50,50,50,1,randphi,randpi);
		phirandom.RandomFill();
		cout<<"Fajlba irom a generalt teret..."<<endl;
		phirandom.ToFile((dir+"randomgeninit.dat").c_str());
		cout<<"...fajlba iras kesz."<<endl;
		opendir=(dir+"randomgeninit.dat");
	

	//mostantól jön a mérés.
	evoltest(opendir,dir,randphi,randpi,false);

	}
		
		
}


//G(t,k_) meghatározása
void GTKclever(string opendir, string dir, double randphi, double randpi,bool noise){

	//Modell létrehozása
	Model phi1(0.1,50,50,50,1,randphi,randpi,opendir.c_str());
	const char* infofile=(dir+"info.dat").c_str();
	
	double dt,maxtime,zerotime;
	int meastime,measmin,timepoint,tinc,counter=0;
	int bins;

	cout<<"Maxtime legyen oszthato dt-vel!"<<endl;
	phi1.ModelInit();

	dt=phi1.Getdt();
	maxtime=phi1.GetMaxtime();
	meastime=phi1.GetMeasTime();

	cout<<"Ha korabbi merest folytatunk, mi a kezdeti idopont? (Random esetben 0)"<<endl;
	cin>>zerotime;

	cout<<"Mekkora legyen a minimális ablakmeret? (dt egysegben) "<<endl;
	cin>>measmin;

	cout<<"Mennyivel noveljuk az ablakmeretet (dt egysegben)? "<<endl;
	cin>>tinc;

	cout<<"Hany kulonbozo ablakot vegyunk fel? "<<endl;
	cin>>timepoint;

	cout<<"Bin-ek száma? "<<endl;
	cin>>bins;

	//info fájl nyitása
	std::ofstream myfile;
	myfile.open((dir+"info.dat").c_str(),std::ios::out);
	myfile<<"Racsallando: "<<phi1.GetGrid()<<endl;
	myfile<<"Tomegnegyzet: "<<phi1.GetM2()<<endl;
	myfile<<"Lambda: "<<phi1.GetL()<<endl;
	myfile<<"Idolepes: "<<phi1.Getdt()<<endl;
	myfile<<"Meresgyakorisag: "<<phi1.GetMeasTime()<<endl;
	myfile<<"Meres vege: "<<phi1.GetMaxtime()<<endl;
	myfile<<"Randphi: "<<phi1.Getrandphi()<<endl;
	myfile<<"Randpi: "<<phi1.Getrandpi()<<endl;
	myfile<<"Min. ablakmeret: "<<measmin<<endl;
	myfile<<"Ablaknoveles: "<<tinc<<endl;
	myfile<<"Ablakszam: "<<timepoint<<endl;

	myfile.close();
	cout<<"A meresinicializalas kesz, kezdem a munkat."<<endl;

	//kezdeti terek kiírása
	cout<<"A kezdeti teret nem irom ismet fajlba."<<endl;
	phi1.FFTModel();
	cout<<"Kezdeti FFT ter fajlba irasa..."<<endl;
	phi1.ToFileFFT((dir+"mivanbennemFFT_init.dat").c_str());
	cout<<"...kesz."<<endl;

	
	//a mérésekhez szükséges:
	int xdim=phi1.GetXmax()*phi1.GetYmax()*(phi1.GetZmax()/2+1);
	array3<double> gtk(1,xdim,4,sizeof(Complex));
	//kmeas-t nullára inicializáljuk
	for(int l=0;l<1;l++){	
	for(int i=0;i<xdim;i++){
			for(int j=0;j<4;j++){
				gtk(l,i,j)=0;
			}//yfor
		}//xfor
	}//lfor
		array3<double> gtk2(1,xdim,6,sizeof(Complex));
	//kmeas-t nullára inicializáljuk
		for(int l=0;l<1;l++){
		for(int i=0;i<xdim;i++){
			for(int j=0;j<6;j++){
				gtk2(l,i,j)=0;
			}//yfor
		}//xfor
		}//lfor
	
		Model temp(1,50,50,50,1,1,0);

			temp=phi1;
		//DEBUG	temp.ToFile((dir+"temperedeti.dat").c_str());
		std::ofstream myfile1;
		myfile1.precision(6);
		myfile1.open((dir+"Propagator_starttime_k.dat").c_str());

		std::ofstream myfile2;
		myfile2.precision(6);
		myfile2.open((dir+"Propagator_Bin.dat").c_str());
		
	//ismétléssel
	for(int count=0;count<timepoint;count++){
	
	//létrehozás és idõfejlõdés
	phi1=temp;
	//DEBUG phi1.ToFile((dir+"vajoneredeti").c_str());
	EvolveModel myevol(phi1,zerotime,1,tinc,measmin+count*tinc,noise);

	cout<<"Idofejlodes kesz es memoriaban van."<<endl;

	
	cout<<"Most kezdjuk a propagatort szamolni!"<<endl;
	myevol.GTK(gtk);
	cout<<"Propagatort kiszamoltuk. Fajlba iras kovetkezik!"<<endl;

	//oszlopokba írjuk ki, ne sorokba
	for(int l=0;l<1;l++){
	for(int i=0;i<xdim;i++){
		myfile1<<measmin+count*tinc<<"\t";
		for(int j=0;j<4;j++){
			myfile1<<gtk(l,i,j)<<"\t";
		}//jfor
		myfile1<<endl;
	}//ifor
	}//lfor
	
	cout<<"Meres vegeredmenye fajlba irva!"<<endl;
	cout<<"Az ideiglenes fajlt toroljuk!"<<endl;
	if( remove( "gtk_prop_temp.dat" ) != 0 )
    perror( "Error deleting file" );
	else
    puts( "File successfully deleted" );
	cout<<"Bin szamolas következik... "<<endl;

		int idx;
		double sitemax=4*sin(M_PI/2)*sin(M_PI/2)+4*sin(M_PI/2)*sin(M_PI/2)+4*sin(M_PI/2)*sin(M_PI/2);
		
		//kitöltjük az elsõ sort a site (határ)értékekkel
		for(int l=0;l<1;l++){
		for(int n=0;n<xdim;n++){
			gtk2(l,n,0)=sitemax/bins*(n+1);
		}
		}//lfor
		for(int l=0;l<1;l++){
		for(int i=0;i<xdim;i++){

			idx=bins*gtk(l,i,3)/sitemax;
			gtk2(l,idx,1)+=gtk(l,i,0);//phi
			gtk2(l,idx,2)+=gtk(l,i,1);//re
			gtk2(l,idx,3)+=gtk(l,i,2);//im
			gtk2(l,idx,4)+=gtk(l,i,3);//sin ?
			gtk2(l,idx,5)++;

		}
		}//lfor
		
	cout<<"...bin-elve rendeztuk. Atlagolas es fajlba iras kovetkezik..."<<endl;
	
	for(int l=0;l<1;l++){
	for(int i=0;i<xdim;i++){
		if(gtk2(l,i,5)!=0){
		gtk2(l,i,1)/=gtk2(l,i,5);
		gtk2(l,i,2)/=gtk2(l,i,5);
		gtk2(l,i,3)/=gtk2(l,i,5);
		gtk2(l,i,4)/=gtk2(l,i,5);
		myfile2<<measmin+count*tinc<<"\t"<<gtk2(l,i,0)<<"\t"<<gtk2(l,i,1)<<"\t"<<gtk2(l,i,2)<<"\t"<<gtk2(l,i,3)<<"\t"<<gtk2(l,i,4)<<endl;
		}//if
		
	}//for
	}//lfor
	
	

	//kmeas-t nullára inicializáljuk
	for(int l=0;l<1;l++){	
	for(int i=0;i<xdim;i++){
			for(int j=0;j<4;j++){
				gtk(l,i,j)=0;
			}//yfor
		}//xfor
	}//lfor

	for(int l=0;l<1;l++){
		for(int i=0;i<xdim;i++){
			for(int j=0;j<6;j++){
				gtk2(l,i,j)=0;
			}//yfor
		}//xfor
		}//lfor

	}//következõ ablak

	myfile1.close();
	myfile2.close();

	cout<<"Minden ablakra kesz!\nFajlba irjuk a vegso valos es reciprokteret!"<<endl;
	phi1.ToFile((dir+"valoster_evolved.dat").c_str());
	phi1.ToFileFFT((dir+"reciprokter_evolved.dat").c_str());

	cout<<"...fajlbairas kesz.!"<<endl;
	
}

void Gphi4kclever(string opendir, string dir, double randphi, double randpi,bool noise){

	//Modell létrehozása
	Model phi1(0.1,50,50,50,1,randphi,randpi,opendir.c_str());
	const char* infofile=(dir+"info.dat").c_str();
	
	double dt,maxtime,zerotime;
	int meastime,measmin,timepoint,tinc,counter=0;
	int bins;

	cout<<"Maxtime legyen oszthato dt-vel!"<<endl;
	phi1.ModelInit();

	dt=phi1.Getdt();
	maxtime=phi1.GetMaxtime();
	meastime=phi1.GetMeasTime();

	cout<<"Ha korabbi merest folytatunk, mi a kezdeti idopont? (Random esetben 0)"<<endl;
	cin>>zerotime;

	cout<<"Mekkora legyen a minimális ablakmeret? (dt egysegben) "<<endl;
	cin>>measmin;

	cout<<"Mennyivel noveljuk az ablakmeretet (dt egysegben)? "<<endl;
	cin>>tinc;

	cout<<"Hany kulonbozo ablakot vegyunk fel? "<<endl;
	cin>>timepoint;

	cout<<"Bin-ek száma? "<<endl;
	cin>>bins;

	//info fájl nyitása
	std::ofstream myfile;
	myfile.open((dir+"info.dat").c_str(),std::ios::out);
	myfile<<"Racsallando: "<<phi1.GetGrid()<<endl;
	myfile<<"Tomegnegyzet: "<<phi1.GetM2()<<endl;
	myfile<<"Lambda: "<<phi1.GetL()<<endl;
	myfile<<"Idolepes: "<<phi1.Getdt()<<endl;
	myfile<<"Meresgyakorisag: "<<phi1.GetMeasTime()<<endl;
	myfile<<"Meres vege: "<<phi1.GetMaxtime()<<endl;
	myfile<<"Randphi: "<<phi1.Getrandphi()<<endl;
	myfile<<"Randpi: "<<phi1.Getrandpi()<<endl;
	myfile<<"Min. ablakmeret: "<<measmin<<endl;
	myfile<<"Ablaknoveles: "<<tinc<<endl;
	myfile<<"Ablakszam: "<<timepoint<<endl;

	myfile.close();
	cout<<"A meresinicializalas kesz, kezdem a munkat."<<endl;

	//kezdeti terek kiírása
	cout<<"A kezdeti teret nem irom ismet fajlba."<<endl;
	phi1.FFTModel();
	cout<<"Kezdeti FFT ter fajlba irasa..."<<endl;
	phi1.ToFileFFT((dir+"mivanbennemFFT_init.dat").c_str());
	cout<<"...kesz."<<endl;

	
	//a mérésekhez szükséges:
	int xdim=phi1.GetXmax()*phi1.GetYmax()*(phi1.GetZmax()/2+1);
	array3<double> gtk(1,xdim,4,sizeof(Complex));
	//kmeas-t nullára inicializáljuk
	for(int l=0;l<1;l++){	
	for(int i=0;i<xdim;i++){
			for(int j=0;j<4;j++){
				gtk(l,i,j)=0;
			}//yfor
		}//xfor
	}//lfor
		array3<double> gtk2(1,xdim,6,sizeof(Complex));
	//kmeas-t nullára inicializáljuk
		for(int l=0;l<1;l++){
		for(int i=0;i<xdim;i++){
			for(int j=0;j<6;j++){
				gtk2(l,i,j)=0;
			}//yfor
		}//xfor
		}//lfor
	
		Model temp=phi1;
		std::ofstream myfile1;
		myfile1.precision(6);
		myfile1.open((dir+"Propagator_starttime_k.dat").c_str());

		std::ofstream myfile2;
		myfile2.precision(6);
		myfile2.open((dir+"Propagator_Bin.dat").c_str());
		
	//ismétléssel
	for(int count=0;count<timepoint;count++){
	
	//létrehozás és idõfejlõdés
		phi1=temp;
	EvolveModel myevol(phi1,zerotime,1,tinc,measmin+count*tinc,noise);

	cout<<"Idofejlodes kesz es memoriaban van."<<endl;

	
	cout<<"Most kezdjuk a propagatort szamolni!"<<endl;
	myevol.phi4(gtk);
	cout<<"Propagatort kiszamoltuk. Fajlba iras kovetkezik!"<<endl;

	//oszlopokba írjuk ki, ne sorokba
	for(int l=0;l<1;l++){
	for(int i=0;i<xdim;i++){
		myfile1<<measmin+count*tinc<<"\t";
		for(int j=0;j<4;j++){
			myfile1<<gtk(l,i,j)<<"\t";
		}//jfor
		myfile1<<endl;
	}//ifor
	}//lfor
	
	cout<<"Meres vegeredmenye fajlba irva!"<<endl;
	cout<<"Az ideiglenes fajlt toroljuk!"<<endl;
	if( remove( "gtk_prop_temp2.dat" ) != 0 )
    perror( "Error deleting file" );
	else
    puts( "File successfully deleted" );
	cout<<"Bin szamolas következik... "<<endl;

		int idx;
		double sitemax=4*sin(M_PI/2)*sin(M_PI/2)+4*sin(M_PI/2)*sin(M_PI/2)+4*sin(M_PI/2)*sin(M_PI/2);
		
		//kitöltjük az elsõ sort a site (határ)értékekkel
		for(int l=0;l<1;l++){
		for(int n=0;n<xdim;n++){
			gtk2(l,n,0)=sitemax/bins*(n+1);
		}
		}//lfor
		for(int l=0;l<1;l++){
		for(int i=0;i<xdim;i++){

			idx=bins*gtk(l,i,3)/sitemax;
			gtk2(l,idx,1)+=gtk(l,i,0);//phi
			gtk2(l,idx,2)+=gtk(l,i,1);//re
			gtk2(l,idx,3)+=gtk(l,i,2);//im
			gtk2(l,idx,4)+=gtk(l,i,3);//sin ?
			gtk2(l,idx,5)++;

		}
		}//lfor
		
	cout<<"...bin-elve rendeztuk. Atlagolas es fajlba iras kovetkezik..."<<endl;
	
	for(int l=0;l<1;l++){
	for(int i=0;i<xdim;i++){
		if(gtk2(l,i,5)!=0){
		gtk2(l,i,1)/=gtk2(l,i,5);
		gtk2(l,i,2)/=gtk2(l,i,5);
		gtk2(l,i,3)/=gtk2(l,i,5);
		gtk2(l,i,4)/=gtk2(l,i,5);
		myfile2<<measmin+count*tinc<<"\t"<<gtk2(l,i,0)<<"\t"<<gtk2(l,i,1)<<"\t"<<gtk2(l,i,2)<<"\t"<<gtk2(l,i,3)<<"\t"<<gtk2(l,i,4)<<endl;
		}//if
		
	}//for
	}//lfor
	
	
	//kmeas-t nullára inicializáljuk
	for(int l=0;l<1;l++){	
	for(int i=0;i<xdim;i++){
			for(int j=0;j<4;j++){
				gtk(l,i,j)=0;
			}//yfor
		}//xfor
	}//lfor

	for(int l=0;l<1;l++){
		for(int i=0;i<xdim;i++){
			for(int j=0;j<6;j++){
				gtk2(l,i,j)=0;
			}//yfor
		}//xfor
		}//lfor

	}//következõ ablak

	myfile1.close();
	myfile2.close();

	cout<<"Minden ablakra kesz!\nFajlba irjuk a vegso valos es reciprokteret!"<<endl;
	phi1.ToFile((dir+"valoster_evolved.dat").c_str());
	phi1.ToFileFFT((dir+"reciprokter_evolved.dat").c_str());

	cout<<"...fajlbairas kesz.!"<<endl;
	
}

//Gtk meghatározása fájlból vagy randomból
void GTK(){
	string opendir;
	string dir;
	double randphi;
	double randpi;
	bool noise=false;
	char noi;

	cout<<"Udvozol a G(t,k_) Meres Program!\n"<<endl;
	cout<<"Ha random kezdeti feltetelt szeretnel, ird: RANDOM"<<endl;
	cout<<"Kulonben ird be a beolvasando fajlnevet eleresi uttal egyutt."<<endl;
	cout<<"Ide varom a valaszt: ";
	cin>>opendir;
	cout<<"Kerem a celkonyvtarat (win: vegen visszaperjellel): "<<endl;
	cin>>dir;

	cout<<"Zajos meres? <y,n>"<<endl;
	cin>>noi;

	if(noi=='y'){
		noise=true;
	}

	cout<<"Mi legyen randphi?"<<endl;
	cin>>randphi;
	cout<<"Mi legyen randpi?"<<endl;
	cin>>randpi;

	//ha random kezdeti feltételt szeretnénk
	if(!(opendir.compare("RANDOM"))){
		
		Model phirandom(0.1,50,50,50,1,randphi,randpi);
		phirandom.RandomFill();
		cout<<"Fajlba irom a generalt teret..."<<endl;
		phirandom.ToFile((dir+"randomgeninit.dat").c_str());
		cout<<"...fajlba iras kesz."<<endl;
		opendir=(dir+"randomgeninit.dat");
	}

	//mostantól jön a mérés.
			GTKclever(opendir,dir,randphi,randpi,noise);
}

void Gphi4k(){
	string opendir;
	string dir;
	double randphi;
	double randpi;
	bool noise=false;
	char noi;

	cout<<"Udvozol a G(t,k_) PHI4 Meres Program!\n"<<endl;
	cout<<"Ha random kezdeti feltetelt szeretnel, ird: RANDOM"<<endl;
	cout<<"Kulonben ird be a beolvasando fajlnevet eleresi uttal egyutt."<<endl;
	cout<<"Ide varom a valaszt: ";
	cin>>opendir;
	cout<<"Kerem a celkonyvtarat (win: vegen visszaperjellel): "<<endl;
	cin>>dir;

	cout<<"Zajos meres? <y,n>"<<endl;
	cin>>noi;

	if(noi=='y'){
		noise=true;
	}

	cout<<"Mi legyen randphi?"<<endl;
	cin>>randphi;
	cout<<"Mi legyen randpi?"<<endl;
	cin>>randpi;

	//ha random kezdeti feltételt szeretnénk
	if(!(opendir.compare("RANDOM"))){
		
		Model phirandom(0.1,50,50,50,1,randphi,randpi);
		phirandom.RandomFill();
		cout<<"Fajlba irom a generalt teret..."<<endl;
		phirandom.ToFile((dir+"randomgeninit.dat").c_str());
		cout<<"...fajlba iras kesz."<<endl;
		opendir=(dir+"randomgeninit.dat");
	}

	//mostantól jön a mérés.
		Gphi4kclever(opendir,dir,randphi,randpi,noise);
}


void T12T12clever(string opendir, string dir, double randphi, double randpi,bool noise){

	//Modell létrehozása
	Model phi1(0.1,50,50,50,1,randphi,randpi,opendir.c_str());
	const char* infofile=(dir+"info.dat").c_str();
	
	double dt,maxtime,zerotime;
	int meastime,measmin,timepoint,tinc,counter=0;
	int bins;

	cout<<"Maxtime legyen oszthato dt-vel!"<<endl;
	phi1.ModelInit();

	dt=phi1.Getdt();
	maxtime=phi1.GetMaxtime();
	meastime=phi1.GetMeasTime();

	cout<<"Ha korabbi merest folytatunk, mi a kezdeti idopont? (Random esetben 0)"<<endl;
	cin>>zerotime;

	cout<<"Mekkora legyen a minimális ablakmeret? (dt egysegben) "<<endl;
	cin>>measmin;

	cout<<"Mennyivel noveljuk az ablakmeretet (dt egysegben)? "<<endl;
	cin>>tinc;

	cout<<"Hany kulonbozo ablakot vegyunk fel? "<<endl;
	cin>>timepoint;

	cout<<"Bin-ek száma? "<<endl;
	cin>>bins;

	//info fájl nyitása
	std::ofstream myfile;
	myfile.open((dir+"info.dat").c_str(),std::ios::out);
	myfile<<"Racsallando: "<<phi1.GetGrid()<<endl;
	myfile<<"Tomegnegyzet: "<<phi1.GetM2()<<endl;
	myfile<<"Lambda: "<<phi1.GetL()<<endl;

	myfile<<"Idolepes: "<<phi1.Getdt()<<endl;
	myfile<<"Meresgyakorisag: "<<phi1.GetMeasTime()<<endl;
	myfile<<"Meres vege: "<<phi1.GetMaxtime()<<endl;
	myfile<<"Randphi: "<<phi1.Getrandphi()<<endl;
	myfile<<"Randpi: "<<phi1.Getrandpi()<<endl;
	myfile<<"Min. ablakmeret: "<<measmin<<endl;
	myfile<<"Ablaknoveles: "<<tinc<<endl;
	myfile<<"Ablakszam: "<<timepoint<<endl;

	myfile.close();
	cout<<"A meresinicializalas kesz, kezdem a munkat."<<endl;

	//kezdeti terek kiírása
	cout<<"A kezdeti teret nem irom ismet fajlba."<<endl;
	phi1.FFTModel();
	cout<<"Kezdeti FFT ter fajlba irasa..."<<endl;
	phi1.ToFileFFT((dir+"mivanbennemFFT_init.dat").c_str());
	cout<<"...kesz."<<endl;

	
	//a mérésekhez szükséges:
	int xdim=phi1.GetXmax()*phi1.GetYmax()*(phi1.GetZmax()/2+1);
	array3<double> gtk(1,xdim,4,sizeof(Complex));
	//kmeas-t nullára inicializáljuk
	for(int l=0;l<1;l++){	
	for(int i=0;i<xdim;i++){
			for(int j=0;j<4;j++){
				gtk(l,i,j)=0;
			}//yfor
		}//xfor
	}//lfor
		array3<double> gtk2(1,xdim,6,sizeof(Complex));
	//kmeas-t nullára inicializáljuk
		for(int l=0;l<1;l++){
		for(int i=0;i<xdim;i++){
			for(int j=0;j<6;j++){
				gtk2(l,i,j)=0;
			}//yfor
		}//xfor
		}//lfor
	
		Model temp=phi1;
		std::ofstream myfile1;
		myfile1.precision(6);
		myfile1.open((dir+"Propagator_starttime_k.dat").c_str());

		std::ofstream myfile2;
		myfile2.precision(6);
		myfile2.open((dir+"Propagator_Bin.dat").c_str());
		
	//ismétléssel
	for(int count=0;count<timepoint;count++){
	
	//létrehozás és idõfejlõdés
		phi1=temp;
	EvolveModel myevol(phi1,zerotime,1,tinc,measmin+count*tinc,noise);

	cout<<"Idofejlodes kesz es memoriaban van."<<endl;

	
	cout<<"Most kezdjuk a propagatort szamolni!"<<endl;
	myevol.visco(gtk);
	cout<<"Propagatort kiszamoltuk. Fajlba iras kovetkezik!"<<endl;

	//oszlopokba írjuk ki, ne sorokba
	for(int l=0;l<1;l++){
	for(int i=0;i<xdim;i++){
		myfile1<<measmin+count*tinc<<"\t";
		for(int j=0;j<4;j++){
			myfile1<<gtk(l,i,j)<<"\t";
		}//jfor
		myfile1<<endl;
	}//ifor
	}//lfor
	
	cout<<"Meres vegeredmenye fajlba irva!"<<endl;
	cout<<"Az ideiglenes fajlt toroljuk!"<<endl;
	if( remove( "gtk_prop_temp2.dat" ) != 0 )
    perror( "Error deleting file" );
	else
    puts( "File successfully deleted" );
	cout<<"Bin szamolas következik... "<<endl;

		int idx;
		double sitemax=4*sin(M_PI/2)*sin(M_PI/2)+4*sin(M_PI/2)*sin(M_PI/2)+4*sin(M_PI/2)*sin(M_PI/2);
		
		//kitöltjük az elsõ sort a site (határ)értékekkel
		for(int l=0;l<1;l++){
		for(int n=0;n<xdim;n++){
			gtk2(l,n,0)=sitemax/bins*(n+1);
		}
		}//lfor
		for(int l=0;l<1;l++){
		for(int i=0;i<xdim;i++){

			idx=bins*gtk(l,i,3)/sitemax;
			gtk2(l,idx,1)+=gtk(l,i,0);//phi
			gtk2(l,idx,2)+=gtk(l,i,1);//re
			gtk2(l,idx,3)+=gtk(l,i,2);//im
			gtk2(l,idx,4)+=gtk(l,i,3);//sin ?
			gtk2(l,idx,5)++;

		}
		}//lfor
		
	cout<<"...bin-elve rendeztuk. Atlagolas es fajlba iras kovetkezik..."<<endl;
	
	for(int l=0;l<1;l++){
	for(int i=0;i<xdim;i++){
		if(gtk2(l,i,5)!=0){
		gtk2(l,i,1)/=gtk2(l,i,5);
		gtk2(l,i,2)/=gtk2(l,i,5);
		gtk2(l,i,3)/=gtk2(l,i,5);
		gtk2(l,i,4)/=gtk2(l,i,5);
		myfile2<<measmin+count*tinc<<"\t"<<gtk2(l,i,0)<<"\t"<<gtk2(l,i,1)<<"\t"<<gtk2(l,i,2)<<"\t"<<gtk2(l,i,3)<<"\t"<<gtk2(l,i,4)<<endl;
		}//if
		
	}//for
	}//lfor
	
	
	//kmeas-t nullára inicializáljuk
	for(int l=0;l<1;l++){	
	for(int i=0;i<xdim;i++){
			for(int j=0;j<4;j++){
				gtk(l,i,j)=0;
			}//yfor
		}//xfor
	}//lfor

	for(int l=0;l<1;l++){
		for(int i=0;i<xdim;i++){
			for(int j=0;j<6;j++){
				gtk2(l,i,j)=0;
			}//yfor
		}//xfor
		}//lfor

	}//következõ ablak

	myfile1.close();
	myfile2.close();

	cout<<"Minden ablakra kesz!\nFajlba irjuk a vegso valos es reciprokteret!"<<endl;
	phi1.ToFile((dir+"valoster_evolved.dat").c_str());
	phi1.ToFileFFT((dir+"reciprokter_evolved.dat").c_str());

	cout<<"...fajlbairas kesz.!"<<endl;
	
}
void T12T12(){
	string opendir;
	string dir;
	double randphi;
	double randpi;
	bool noise=false;
	char noi;

	cout<<"Udvozol a T12T12(t,k_) Viszko Meres Program!\n"<<endl;
	cout<<"Ha random kezdeti feltetelt szeretnel, ird: RANDOM"<<endl;
	cout<<"Kulonben ird be a beolvasando fajlnevet eleresi uttal egyutt."<<endl;
	cout<<"Ide varom a valaszt: ";
	cin>>opendir;
	cout<<"Kerem a celkonyvtarat (win: vegen visszaperjellel): "<<endl;
	cin>>dir;

	cout<<"Zajos meres? <y,n>"<<endl;
	cin>>noi;

	if(noi=='y'){
		noise=true;
	}

	cout<<"Mi legyen randphi?"<<endl;
	cin>>randphi;
	cout<<"Mi legyen randpi?"<<endl;
	cin>>randpi;
	

	//ha random kezdeti feltételt szeretnénk
	if(!(opendir.compare("RANDOM"))){
		
		Model phirandom(0.1,50,50,50,1,randphi,randpi);
		phirandom.RandomFill();
		cout<<"Fajlba irom a generalt teret..."<<endl;
		phirandom.ToFile((dir+"randomgeninit.dat").c_str());
		cout<<"...fajlba iras kesz."<<endl;
		opendir=(dir+"randomgeninit.dat");
	}

	//mostantól jön a mérés.
		T12T12clever(opendir,dir,randphi,randpi,noise);

}
//void ComboClever(string opendir, string dir, double randphi, double randpi){
//
//	//Modell létrehozása
//	Model phi1(0.1,50,50,50,1,randphi,randpi,opendir.c_str());
//	const char* infofile=(dir+"info.dat").c_str();
//	
//	double dt,maxtime,zerotime;
//	int meastime,measmin,timepoint,tinc,counter=0;
//	int bins;
//
//	cout<<"Maxtime legyen oszthato dt-vel!"<<endl;
//	phi1.ModelInit();
//
//	dt=phi1.Getdt();
//	maxtime=phi1.GetMaxtime();
//	meastime=phi1.GetMeasTime();
//
//	cout<<"Ha korabbi merest folytatunk, mi a kezdeti idopont? (Random esetben 0)"<<endl;
//	cin>>zerotime;
//
//	cout<<"Mekkora legyen a minimális ablakmeret? (dt egysegben) "<<endl;
//	cin>>measmin;
//
//	cout<<"Mennyivel noveljuk az ablakmeretet (dt egysegben)? "<<endl;
//	cin>>tinc;
//
//	cout<<"Hany kulonbozo ablakot vegyunk fel? "<<endl;
//	cin>>timepoint;
//
//	cout<<"Bin-ek száma? "<<endl;
//	cin>>bins;
//
//	//info fájl nyitása
//	std::ofstream myfile;
//	myfile.open((dir+"info.dat").c_str(),std::ios::out);
//	myfile<<"Racsallando: "<<phi1.GetGrid()<<endl;
//	myfile<<"Tomegnegyzet: "<<phi1.GetM2()<<endl;
//	myfile<<"Lambda: "<<phi1.GetL()<<endl;
//	myfile<<"Idolepes: "<<phi1.Getdt()<<endl;
//	myfile<<"Meresgyakorisag: "<<phi1.GetMeasTime()<<endl;
//	myfile<<"Meres vege: "<<phi1.GetMaxtime()<<endl;
//	myfile<<"Randphi: "<<phi1.Getrandphi()<<endl;
//	myfile<<"Randpi: "<<phi1.Getrandpi()<<endl;
//	myfile<<"Min. ablakmeret: "<<measmin<<endl;
//	myfile<<"Ablaknoveles: "<<tinc<<endl;
//	myfile<<"Ablakszam: "<<timepoint<<endl;
//
//	myfile.close();
//	cout<<"A meresinicializalas kesz, kezdem a munkat."<<endl;
//
//	//kezdeti terek kiírása
//	cout<<"A kezdeti teret nem irom ismet fajlba."<<endl;
//	phi1.FFTModel();
//	cout<<"Kezdeti FFT ter fajlba irasa..."<<endl;
//	phi1.ToFileFFT((dir+"mivanbennemFFT_init.dat").c_str());
//	cout<<"...kesz."<<endl;
//
//	
//	//a mérésekhez szükséges:
//	int xdim=50*50*(25+1);
//	array3<double> gtk(1,xdim,4,sizeof(Complex));
//	//kmeas-t nullára inicializáljuk
//	for(int l=0;l<1;l++){	
//	for(int i=0;i<xdim;i++){
//			for(int j=0;j<4;j++){
//				gtk(l,i,j)=0;
//			}//yfor
//		}//xfor
//	}//lfor
//		array3<double> gtk2(1,xdim,6,sizeof(Complex));
//	//kmeas-t nullára inicializáljuk
//		for(int l=0;l<1;l++){
//		for(int i=0;i<xdim;i++){
//			for(int j=0;j<6;j++){
//				gtk2(l,i,j)=0;
//			}//yfor
//		}//xfor
//		}//lfor
//	
////
//		array3<double> Ttk(1,xdim,4,sizeof(Complex));
//	//kmeas-t nullára inicializáljuk
//	for(int l=0;l<1;l++){	
//	for(int i=0;i<xdim;i++){
//			for(int j=0;j<4;j++){
//				Ttk(l,i,j)=0;
//			}//yfor
//		}//xfor
//	}//lfor
//		array3<double> Ttk2(1,xdim,6,sizeof(Complex));
//	//kmeas-t nullára inicializáljuk
//		for(int l=0;l<1;l++){
//		for(int i=0;i<xdim;i++){
//			for(int j=0;j<6;j++){
//				Ttk2(l,i,j)=0;
//			}//yfor
//		}//xfor
//		}//lfor
//
//		Model temp=phi1;
//		std::ofstream myfile1;
//		myfile1.precision(6);
//		myfile1.open((dir+"GTKPropagator_starttime_k.dat").c_str());
//
//		std::ofstream myfile2;
//		myfile2.precision(6);
//		myfile2.open((dir+"GTKPropagator_Bin.dat").c_str());
//
//		std::ofstream myfile3;
//		myfile3.precision(6);
//		myfile3.open((dir+"T12Propagator_starttime_k.dat").c_str());
//
//		std::ofstream myfile4;
//		myfile4.precision(6);
//		myfile4.open((dir+"T12Propagator_Bin.dat").c_str());
//		
//	//ismétléssel
//	for(int count=0;count<timepoint;count++){
//	
//	//létrehozás és idõfejlõdés
//		phi1=temp;
//	EvolveModel myevol(phi1,zerotime,1,tinc,measmin+count*tinc);
//
//	cout<<"Idofejlodes kesz es memoriaban van."<<endl;
//
//	
//	cout<<"Most kezdjuk a propagatort szamolni!"<<endl;
//	myevol.GTK(gtk);
//	myevol.visco(Ttk);
//	cout<<"Propagatort kiszamoltuk. Fajlba iras kovetkezik!"<<endl;
//
//	//oszlopokba írjuk ki, ne sorokba
//	for(int l=0;l<1;l++){
//	for(int i=0;i<xdim;i++){
//		myfile1<<measmin+count*tinc<<"\t";
//		for(int j=0;j<4;j++){
//			myfile1<<gtk(l,i,j)<<"\t";
//		}//jfor
//		myfile1<<endl;
//	}//ifor
//	}//lfor
//	
//	cout<<"Meres vegeredmenye fajlba irva!"<<endl;
//	cout<<"Az ideiglenes fajlt toroljuk!"<<endl;
//	if( remove( "gtk_prop_temp2.dat" ) != 0 )
//    perror( "Error deleting file" );
//	else
//    puts( "File successfully deleted" );
//
//	//oszlopokba írjuk ki, ne sorokba
//	for(int l=0;l<1;l++){
//	for(int i=0;i<xdim;i++){
//		myfile3<<measmin+count*tinc<<"\t";
//		for(int j=0;j<4;j++){
//			myfile3<<Ttk(l,i,j)<<"\t";
//		}//jfor
//		myfile3<<endl;
//	}//ifor
//	}//lfor
//
//	cout<<"Meres vegeredmenye fajlba irva!"<<endl;
//	cout<<"Az ideiglenes fajlt toroljuk!"<<endl;
//	if( remove( "gtk_prop_tempT12.dat" ) != 0 )
//    perror( "Error deleting file" );
//	else
//    puts( "File successfully deleted" );
//	cout<<"Bin szamolas következik... "<<endl;
//
//		int idx;
//		double sitemax=4*sin(M_PI/2)*sin(M_PI/2)+4*sin(M_PI/2)*sin(M_PI/2)+4*sin(M_PI/2)*sin(M_PI/2);
//		
//		//kitöltjük az elsõ sort a site (határ)értékekkel
//		for(int l=0;l<1;l++){
//		for(int n=0;n<xdim;n++){
//			gtk2(l,n,0)=sitemax/bins*(n+1);
//		}
//		}//lfor
//
//		//kitöltjük az elsõ sort a site (határ)értékekkel
//		for(int l=0;l<1;l++){
//		for(int n=0;n<xdim;n++){
//			Ttk2(l,n,0)=sitemax/bins*(n+1);
//		}
//		}//lfor
//		for(int l=0;l<1;l++){
//		for(int i=0;i<xdim;i++){
//
//			idx=bins*gtk(l,i,3)/sitemax;
//			gtk2(l,idx,1)+=gtk(l,i,0);//phi
//			gtk2(l,idx,2)+=gtk(l,i,1);//re
//			gtk2(l,idx,3)+=gtk(l,i,2);//im
//			gtk2(l,idx,4)+=gtk(l,i,3);//sin ?
//			gtk2(l,idx,5)++;
//
//		}
//		}//lfor
//		
//		for(int l=0;l<1;l++){
//		for(int i=0;i<xdim;i++){
//
//			idx=bins*Ttk(l,i,3)/sitemax;
//			Ttk2(l,idx,1)+=Ttk(l,i,0);//phi
//			Ttk2(l,idx,2)+=Ttk(l,i,1);//re
//			Ttk2(l,idx,3)+=Ttk(l,i,2);//im
//			Ttk2(l,idx,4)+=Ttk(l,i,3);//sin ?
//			Ttk2(l,idx,5)++;
//
//		}
//		}//lfor
//
//	cout<<"...bin-elve rendeztuk. Atlagolas es fajlba iras kovetkezik..."<<endl;
//	
//	for(int l=0;l<1;l++){
//	for(int i=0;i<xdim;i++){
//		if(gtk2(l,i,5)!=0){
//		gtk2(l,i,1)/=gtk2(l,i,5);
//		gtk2(l,i,2)/=gtk2(l,i,5);
//		gtk2(l,i,3)/=gtk2(l,i,5);
//		gtk2(l,i,4)/=gtk2(l,i,5);
//		myfile2<<measmin+count*tinc<<"\t"<<gtk2(l,i,0)<<"\t"<<gtk2(l,i,1)<<"\t"<<gtk2(l,i,2)<<"\t"<<gtk2(l,i,3)<<"\t"<<gtk2(l,i,4)<<endl;
//		}//if
//		
//	}//for
//	}//lfor
//	
//	
//	for(int l=0;l<1;l++){
//	for(int i=0;i<xdim;i++){
//		if(Ttk2(l,i,5)!=0){
//		Ttk2(l,i,1)/=Ttk2(l,i,5);
//		Ttk2(l,i,2)/=Ttk2(l,i,5);
//		Ttk2(l,i,3)/=Ttk2(l,i,5);
//		Ttk2(l,i,4)/=Ttk2(l,i,5);
//		myfile4<<measmin+count*tinc<<"\t"<<Ttk2(l,i,0)<<"\t"<<Ttk2(l,i,1)<<"\t"<<Ttk2(l,i,2)<<"\t"<<Ttk2(l,i,3)<<"\t"<<Ttk2(l,i,4)<<endl;
//		}//if
//		
//	}//for
//	}//lfor
//
//	//kmeas-t nullára inicializáljuk
//	for(int l=0;l<1;l++){	
//	for(int i=0;i<xdim;i++){
//			for(int j=0;j<4;j++){
//				gtk(l,i,j)=0;
//			}//yfor
//		}//xfor
//	}//lfor
//
//	for(int l=0;l<1;l++){
//		for(int i=0;i<xdim;i++){
//			for(int j=0;j<6;j++){
//				gtk2(l,i,j)=0;
//			}//yfor
//		}//xfor
//		}//lfor
//
//	//kmeas-t nullára inicializáljuk
//	for(int l=0;l<1;l++){	
//	for(int i=0;i<xdim;i++){
//			for(int j=0;j<4;j++){
//				Ttk(l,i,j)=0;
//			}//yfor
//		}//xfor
//	}//lfor
//
//	for(int l=0;l<1;l++){
//		for(int i=0;i<xdim;i++){
//			for(int j=0;j<6;j++){
//				Ttk2(l,i,j)=0;
//			}//yfor
//		}//xfor
//		}//lfor
//
//	}//következõ ablak
//
//	myfile1.close();
//	myfile2.close();
//	myfile3.close();
//	myfile4.close();
//
//	cout<<"Minden ablakra kesz!\nFajlba irjuk a vegso valos es reciprokteret!"<<endl;
//	phi1.ToFile((dir+"valoster_evolved.dat").c_str());
//	phi1.ToFileFFT((dir+"reciprokter_evolved.dat").c_str());
//
//	cout<<"...fajlbairas kesz.!"<<endl;
//	
//}

void ComboExtra(string opendir, string dir, double randphi, double randpi,bool noise){

	//Modell létrehozása
	Model phi1(0.1,50,50,50,1,randphi,randpi,opendir.c_str());
	const char* infofile=(dir+"info.dat").c_str();
	
	double dt,maxtime,zerotime;
	double xi,gamma;
	int meastime,measmin,timepoint,tinc,counter=0;
	int bins;
	char swal;
	cout<<"Maxtime legyen oszthato dt-vel!"<<endl;
	phi1.ModelInit();

	dt=phi1.Getdt();
	maxtime=phi1.GetMaxtime();
	meastime=phi1.GetMeasTime();

	cout<<"Ha korabbi merest folytatunk, mi a kezdeti idopont? (Random esetben 0)"<<endl;
	cin>>zerotime;

	cout<<"Mekkora legyen a minimális ablakmeret? (dt egysegben) "<<endl;
	cin>>measmin;

	cout<<"Mennyivel noveljuk az ablakmeretet (dt egysegben)? "<<endl;
	cin>>tinc;

	cout<<"Hany kulonbozo ablakot vegyunk fel? "<<endl;
	cin>>timepoint;

	cout<<"Bin-ek száma? "<<endl;
	cin>>bins;

	cout<<"Mennyi legyen xi,gamma?"<<endl;
	cin>>xi>>swal>>gamma;
	phi1.SetGamma(gamma);
	phi1.SetXiint(xi);

	//info fájl nyitása
	std::ofstream myfile;
	myfile.open((dir+"info.dat").c_str(),std::ios::out);
	myfile<<"Racsallando: "<<phi1.GetGrid()<<endl;
	myfile<<"Tomegnegyzet: "<<phi1.GetM2()<<endl;
	myfile<<"Lambda: "<<phi1.GetL()<<endl;
	myfile<<"Idolepes: "<<phi1.Getdt()<<endl;
	myfile<<"Meresgyakorisag: "<<phi1.GetMeasTime()<<endl;
	myfile<<"Meres vege: "<<phi1.GetMaxtime()<<endl;
	myfile<<"Randphi: "<<phi1.Getrandphi()<<endl;
	myfile<<"Randpi: "<<phi1.Getrandpi()<<endl;
	myfile<<"Min. ablakmeret: "<<measmin<<endl;
	myfile<<"Ablaknoveles: "<<tinc<<endl;
	myfile<<"Ablakszam: "<<timepoint<<endl;
	myfile<<"Xiint: "<<phi1.Getxiint()<<endl;
	myfile<<"Gamma: "<<phi1.GetGamma()<<endl;

	myfile.close();
	cout<<"A meresinicializalas kesz, kezdem a munkat."<<endl;

	//kezdeti terek kiírása
	cout<<"A kezdeti teret nem irom ismet fajlba."<<endl;
	phi1.FFTModel();
	cout<<"Kezdeti FFT ter fajlba irasa..."<<endl;
	phi1.ToFileFFT((dir+"mivanbennemFFT_init.dat").c_str());
	cout<<"...kesz."<<endl;

	
	//a mérésekhez szükséges:
	int xdim=phi1.GetXmax()*phi1.GetYmax()*(phi1.GetZmax()/2+1);
	array3<double> gtk(1,xdim,6,sizeof(Complex));
	//kmeas-t nullára inicializáljuk
	for(int l=0;l<1;l++){	
	for(int i=0;i<xdim;i++){
			for(int j=0;j<6;j++){
				gtk(l,i,j)=0;
			}//yfor
		}//xfor
	}//lfor
		array3<double> gtk2(1,xdim,8,sizeof(Complex));
	//kmeas-t nullára inicializáljuk
		for(int l=0;l<1;l++){
		for(int i=0;i<xdim;i++){
			for(int j=0;j<8;j++){
				gtk2(l,i,j)=0;
			}//yfor
		}//xfor
		}//lfor
	
//

		Model temp=phi1;
		std::ofstream myfile1;
		myfile1.precision(6);
		myfile1.open((dir+"GTKPropagator_starttime_k.dat").c_str());

		std::ofstream myfile2;
		myfile2.precision(6);
		myfile2.open((dir+"GTKPropagator_Bin.dat").c_str());

		
		
	//ismétléssel
	for(int count=0;count<timepoint;count++){
	
	//létrehozás és idõfejlõdés
		phi1=temp;
	EvolveModel myevol(phi1,zerotime,1,tinc,measmin+count*tinc,noise);

	cout<<"Idofejlodes kesz es memoriaban van."<<endl;

	
	cout<<"Most kezdjuk a propagatort szamolni!"<<endl;
	myevol.evolveCombo(gtk);
	cout<<"Propagatort kiszamoltuk. Fajlba iras kovetkezik!"<<endl;

	//oszlopokba írjuk ki, ne sorokba
	for(int l=0;l<1;l++){
	for(int i=0;i<xdim;i++){
		myfile1<<measmin+count*tinc<<"\t";
		for(int j=0;j<6;j++){
			myfile1<<gtk(l,i,j)<<"\t";
		}//jfor
		myfile1<<endl;
	}//ifor
	}//lfor
	
	cout<<"Meres vegeredmenye fajlba irva!"<<endl;
	/*cout<<"Az ideiglenes fajlt toroljuk!"<<endl;
	if( remove( "gtk_prop_temp2.dat" ) != 0 )
    perror( "Error deleting file" );
	else
    puts( "File successfully deleted" );*/

	

		int idx;
		double sitemax=4*sin(M_PI/2)*sin(M_PI/2)+4*sin(M_PI/2)*sin(M_PI/2)+4*sin(M_PI/2)*sin(M_PI/2);
		
		//kitöltjük az elsõ sort a site (határ)értékekkel
		for(int l=0;l<1;l++){
		for(int n=0;n<xdim;n++){
			gtk2(l,n,0)=sitemax/bins*(n+1);
		}
		}//lfor

		
		for(int l=0;l<1;l++){
		for(int i=0;i<xdim;i++){

			idx=bins*gtk(l,i,3)/sitemax;
			gtk2(l,idx,1)+=gtk(l,i,0);//phi
			gtk2(l,idx,2)+=gtk(l,i,1);//re T12
			gtk2(l,idx,3)+=gtk(l,i,2);//im T12
			gtk2(l,idx,4)+=gtk(l,i,3);//sin ?
			gtk2(l,idx,5)++; //counter
			gtk2(l,idx,6)+=gtk(l,i,4);//re PhikPhik
			gtk2(l,idx,7)+=gtk(l,i,5);//im PhikPhik

		}
		}//lfor
		

	cout<<"...bin-elve rendeztuk. Atlagolas es fajlba iras kovetkezik..."<<endl;
	
	for(int l=0;l<1;l++){
	for(int i=0;i<xdim;i++){
		if(gtk2(l,i,5)!=0){
		gtk2(l,i,1)/=gtk2(l,i,5);
		gtk2(l,i,2)/=gtk2(l,i,5);
		gtk2(l,i,3)/=gtk2(l,i,5);
		gtk2(l,i,4)/=gtk2(l,i,5);
		gtk2(l,i,6)/=gtk2(l,i,5);
		gtk2(l,i,7)/=gtk2(l,i,5);
		myfile2<<measmin+count*tinc<<"\t"<<gtk2(l,i,0)<<"\t"<<gtk2(l,i,1)<<"\t"<<gtk2(l,i,2)<<"\t"<<gtk2(l,i,3)<<"\t"<<gtk2(l,i,4)<<"\t"<<gtk2(l,i,6)<<"\t"<<gtk2(l,i,7)<<endl;
		}//if
		
	}//for
	}//lfor
	
	
	

	//kmeas-t nullára inicializáljuk
	for(int l=0;l<1;l++){	
	for(int i=0;i<xdim;i++){
			for(int j=0;j<6;j++){
				gtk(l,i,j)=0;
			}//yfor
		}//xfor
	}//lfor

	for(int l=0;l<1;l++){
		for(int i=0;i<xdim;i++){
			for(int j=0;j<8;j++){
				gtk2(l,i,j)=0;
			}//yfor
		}//xfor
		}//lfor

	

	}//következõ ablak

	myfile1.close();
	myfile2.close();
	
	cout<<"Minden ablakra kesz!\nFajlba irjuk a vegso valos es reciprokteret!"<<endl;
	phi1.ToFile((dir+"valoster_evolved.dat").c_str());
	phi1.ToFileFFT((dir+"reciprokter_evolved.dat").c_str());

	cout<<"...fajlbairas kesz.!"<<endl;
	
}

//void Combo(){
//	string opendir;
//	string dir;
//	double randphi;
//	double randpi;
//
//
//	cout<<"Udvozol a Combo Meres Program!\n"<<endl;
//	cout<<"Ha random kezdeti feltetelt szeretnel, ird: RANDOM"<<endl;
//	cout<<"Kulonben ird be a beolvasando fajlnevet eleresi uttal egyutt."<<endl;
//	cout<<"Ide varom a valaszt: ";
//	cin>>opendir;
//	cout<<"Kerem a celkonyvtarat (win: vegen visszaperjellel): "<<endl;
//	cin>>dir;
//
//	cout<<"Zajos meres? <y,n>"<<endl;
//	cin>>noi;
//
//
//
//	cout<<"Mi legyen randphi?"<<endl;
//	cin>>randphi;
//	cout<<"Mi legyen randpi?"<<endl;
//	cin>>randpi;
//	
//
//	//ha random kezdeti feltételt szeretnénk
//	if(!(opendir.compare("RANDOM"))){
//		
//		Model phirandom(0.1,50,50,50,1,randphi,randpi);
//		phirandom.RandomFill();
//		cout<<"Fajlba irom a generalt teret..."<<endl;
//		phirandom.ToFile((dir+"randomgeninit.dat").c_str());
//		cout<<"...fajlba iras kesz."<<endl;
//		opendir=(dir+"randomgeninit.dat");
//	}
//
//	//mostantól jön a mérés.
//		ComboClever(opendir,dir,randphi,randpi);
//
//}

void ComboE(){
	string opendir;
	string dir;
	double randphi;
	double randpi;
	bool noise=false;
	char noi;


	cout<<"Udvozol a Combo Meres Program!\n"<<endl;
	cout<<"Ha random kezdeti feltetelt szeretnel, ird: RANDOM"<<endl;
	cout<<"Kulonben ird be a beolvasando fajlnevet eleresi uttal egyutt."<<endl;
	cout<<"Ide varom a valaszt: ";
	cin>>opendir;
	cout<<"Kerem a celkonyvtarat (win: vegen visszaperjellel): "<<endl;
	cin>>dir;

	cout<<"Zajos meres? <y,n>"<<endl;
	cin>>noi;

	if(noi=='y'){
		noise=true;
	}

	cout<<"Mi legyen randphi?"<<endl;
	cin>>randphi;
	cout<<"Mi legyen randpi?"<<endl;
	cin>>randpi;
	

	//ha random kezdeti feltételt szeretnénk
	if(!(opendir.compare("RANDOM"))){
		
		Model phirandom(0.1,50,50,50,1,randphi,randpi);
		phirandom.RandomFill();
		cout<<"Fajlba irom a generalt teret..."<<endl;
		phirandom.ToFile((dir+"randomgeninit.dat").c_str());
		cout<<"...fajlba iras kesz."<<endl;
		opendir=(dir+"randomgeninit.dat");
	}

	//mostantól jön a mérés.
		ComboExtra(opendir,dir,randphi,randpi,noise);

}

void fieldtest(){
	string opendir;
	string dir;
	double randphi;
	double randpi;
	string idoido="n";
	string clever="y";

	cout<<"Udvozol a G(t,k_) Meres Program!\n"<<endl;
	cout<<"Ha random kezdeti feltetelt szeretnel, ird: RANDOM"<<endl;
	cout<<"Kulonben ird be a beolvasando fajlnevet eleresi uttal egyutt."<<endl;
	cout<<"Ide varom a valaszt: ";
	cin>>opendir;
	cout<<"Kerem a celkonyvtarat (win: vegen visszaperjellel): "<<endl;
	cin>>dir;

	Field test(50,50,50,1,opendir.c_str());
	test.RunFFT();
	test.ToFileField((dir+"eredeti.txt").c_str());
	
	Field test2(50,50,50,1);
	
	test2=test;
	test2.ToFileField((dir+"opegyenlo.txt").c_str());
	
}
void modeltest(){
	string opendir;
	string dir;
	double randphi;
	double randpi;
	string idoido="n";
	string clever="y";

	cout<<"Udvozol a G(t,k_) Meres Program!\n"<<endl;
	cout<<"Ha random kezdeti feltetelt szeretnel, ird: RANDOM"<<endl;
	cout<<"Kulonben ird be a beolvasando fajlnevet eleresi uttal egyutt."<<endl;
	cout<<"Ide varom a valaszt: ";
	cin>>opendir;
	cout<<"Kerem a celkonyvtarat (win: vegen visszaperjellel): "<<endl;
	cin>>dir;

	Model test(1,50,50,50,1,0,1,opendir.c_str());
	
	test.ToFile((dir+"eredeti.txt").c_str());
	
	Model test2(1,50,50,50,1,0,1);
	
	test2=test;
	test.RefreshField();
	test.ToFile((dir+"evolved.txt").c_str());
	test2.ToFile((dir+"opegyenlo.txt").c_str());
	test=test2;
	test.ToFile((dir+"opegyenlo2.txt").c_str());
}




/* Tömeg Fixáló programcsokor*/

/*Tömegfix segéd: Energia beálltát ellenõrzi*/
void EFixer(int &a,
	std::ofstream &myfile, std::ofstream &myfile1,std::ofstream &myfile2,std::ofstream &myfile3,std::ofstream &myfile4,
	double &zerotime,const double &maxtime,const double &dt,Model& phi1){
		int counter=0; //idõciklushoz
		int meastime=phi1.GetMeasTime();
		//a mérésekhez szükséges:
		double sumenergy,temper;
		
		vector<double> mom;
		
		//EFixer specifikus
		double data[50];
		double avgdata=0;
		int count=0; //0-tól 50-ig számol
		double maxDiff=0;
		double tempDiff=0;

		myfile<<"Hello: EFixer"<<endl;

		//idõfejlõdés
	for(double i=zerotime;(i<maxtime)&&(a==0);i+=dt){
		cout<<i<<endl;
		
		//ha mérésidõ van
		if((counter%meastime==0)){
			cout<<"Meres!"<<endl;
			phi1.FFTModel();
			
			sumenergy=phi1.CountSumEnergy();
			myfile1<<i<<"\t"<<sumenergy<<endl;
			avgdata+=sumenergy; //átlagoláshoz ha 50 kigyûlt és számoltunk, lenullázzuk

			temper=phi1.CountTemperature();
			myfile2<<i<<"\t"<<temper<<endl;
			
			mom=phi1.MomentMeas();
			myfile3<<i<<"\t"<<mom[0]<<endl;
			myfile4<<i<<"\t"<<mom[1]<<endl;

			//ha kevesebb mint 50 van, adatot gyûjtünk
			if(count<50){
				data[count]=sumenergy;
			}
			else if(count==50){//ha kigyûlt az 50 adat
				avgdata/=count; //átlagolunk
				maxDiff=0;
				tempDiff=0;
				for(int j=0;j<count;j++){
					tempDiff=abs((data[j]-avgdata));
					if(tempDiff>maxDiff){
						maxDiff=tempDiff;
					}
				}//legnagyobb eltérés keresése

				if((maxDiff/avgdata)<0.08){
					a=1;
					zerotime=i;
					myfile<<"a="<<a<<", aktualis ido="<<zerotime<<", Mert energia="<<avgdata<<", Megyek PhiTesterhez!"<<endl;
					return;
				}//ha konfidencián belül vagyunk, visszatérünk, a-t állítjuk

				//ha még nem állt be
				avgdata=0;
				count=-1;


			}//ha kigyûlt 50 adat

			count++;
		}//ha mérésidõ vége
		counter++;
	
		cout<<"Refreshing..."<<endl;
		//fejlesztünk
			phi1.RefreshNoiseField();
		
		
    }//idõfejlõdés

	//ha lejárt az idõ, akkor lépjünk ki
	a=5;
	zerotime=maxtime;
	myfile<<"a="<<a<<", Kifutottunk az idobol!"<<endl;
	
}

/*Tömegfix segéd: ellenõrzi, hogy nem a spontán sértett állapotban vagyunk-e*/
void PhiTester(int &a,
	 std::ofstream &myfile, std::ofstream &myfile1, std::ofstream &myfile2, std::ofstream &myfile3,std::ofstream &myfile4,
	double &zerotime,const double &maxtime,const double &dt,Model& phi1){
		int counter=0; //idõciklushoz
		int meastime=phi1.GetMeasTime();
		//a mérésekhez szükséges:
		double sumenergy,temper;
		vector<double> mom=phi1.MomentMeas();
		
		//EFixer specifikus
		double avgdata=0;
		int count=0; //0-tól 8000-ig számol
		double maxDiff=0;
		
		myfile<<"Hello: PhiTester"<<endl;

		//idõfejlõdés
	for(double i=zerotime;(i<maxtime)&&(a==1);i+=dt){
		cout<<i<<endl;
		
		//ha mérésidõ van
		if((counter%meastime==0)){
			cout<<"Meres!"<<endl;
			phi1.FFTModel();
			
			sumenergy=phi1.CountSumEnergy();
			myfile1<<i<<"\t"<<sumenergy<<endl;
			

			temper=phi1.CountTemperature();
			myfile2<<i<<"\t"<<temper<<endl;
			
			mom=phi1.MomentMeas();
			avgdata+=mom[0];//átlagoláshoz ha 8000 kigyûlt és számoltunk, lenullázzuk
			myfile3<<i<<"\t"<<mom[0]<<endl;
			myfile4<<i<<"\t"<<mom[1]<<endl;

			if(count==8000){//ha kigyûlt az 8000 adat
				avgdata/=count; //átlagolunk
				maxDiff=abs(avgdata);

				if(maxDiff<0.05){
					a=2; //mehet a tömegcalc
					zerotime=i;
					myfile<<"a="<<a<<", aktualis ido="<<zerotime<<", Mert phimean="<<avgdata<<", Megyek MCalchoz!"<<endl;
					return;
				}//ha konfidencián belül vagyunk, visszatérünk, a-t állítjuk

				//ha még nem állt be
				a=10;//elküldjük feedback 2-nek
				zerotime=i;
				myfile<<"a="<<a<<", aktualis ido="<<zerotime<<", Mert phimean="<<avgdata<<", Megyek FeedB2-höz!"<<endl;
				return;

			}//ha kigyûlt 8000 adat

			count++;
		}//ha mérésidõ vége
		counter++;
	
		cout<<"Refreshing..."<<endl;
		//fejlesztünk
			phi1.RefreshNoiseField();
		
		
    }//idõfejlõdés

	//ha lejárt az idõ, akkor lépjünk ki
	a=6;
	zerotime=maxtime;
	myfile<<"a="<<a<<", aktualis ido="<<zerotime<<", Kifutottunk az idobol!"<<endl;
}

/*Tömegfix segéd: zérushelyet keres phimean-ben*/
double FindZH(int &a,
	  std::ofstream &myfile, std::ofstream &myfile1, std::ofstream &myfile2, std::ofstream &myfile3,std::ofstream &myfile4,
	double &zerotime,const double &maxtime,const double &dt,Model& phi1,bool &sign){
		int counter=0; //idõciklushoz
		int meastime=phi1.GetMeasTime();
		//a mérésekhez szükséges:
		double sumenergy,temper;
		vector<double> mom;
		
		//FindZH specifikus
		double tless=zerotime; //elõjelváltás elõtt
		double tmore=maxtime; //elõjelváltás után
		const bool flag=sign; //ezzel ellenõrizzük az elõjelváltást

		myfile<<"Hello: FindZH"<<endl;

		//idõfejlõdés
	for(double i=zerotime;(i<maxtime)&&(a==2);i+=dt){
		cout<<i<<endl;
		
		//ha mérésidõ van
		if((counter%meastime==0)){
			cout<<"Meres!"<<endl;
			phi1.FFTModel();
			
			sumenergy=phi1.CountSumEnergy();
			myfile1<<i<<"\t"<<sumenergy<<endl;
			

			temper=phi1.CountTemperature();
			myfile2<<i<<"\t"<<temper<<endl;
			
			mom=phi1.MomentMeas();
			myfile3<<i<<"\t"<<mom[0]<<endl;
			myfile4<<i<<"\t"<<mom[1]<<endl;
			//ha eredetileg pozitív volt
			if(flag){
				if(mom[0]<0){
					sign=false;//most már negatív
				}
				if(sign){
					tless=i;
				}
				else{//megtaláltuk a zérushelyet
					tmore=i;
					
					double ZH=zerotime+((tmore-tless)/2);
					zerotime=i;
					myfile<<"a="<<a<<", aktualis ido="<<zerotime<<", ZH:"<<ZH<<"-nél"<<endl;
					return ZH;
					
				}
			}
			else{//ha eredetileg negatív volt
				if(mom[0]>0){
					sign=true; //mostmár pozitív
				}
				if((!sign)){
					tless=i;
				}
				else{//megtaláltuk a zérushelyet
					tmore=i;
					
					double ZH=zerotime+((tmore-tless)/2);
					zerotime=i;
					myfile<<"a="<<a<<", aktualis ido="<<zerotime<<", ZH:"<<ZH<<"-nél"<<endl;
					return ZH;
				}
			}
		}//ha mérésidõ vége
		counter++;
	
		cout<<"Refreshing..."<<endl;
		//fejlesztünk
			phi1.RefreshNoiseField();
		
		
    }//idõfejlõdés

	

	//ha lejárt az idõ, akkor lépjünk ki
	a=7;
	zerotime=maxtime;
	myfile<<"a="<<a<<", kifutottunk az idobol!"<<endl;
	return -500;
}

/*Tömegfix segéd: ha spontán sértésben vagyunk, melegít vagy hût*/
void FeedB2(int &a,
	 std::ofstream &myfile, std::ofstream &myfile1, std::ofstream &myfile2, std::ofstream &myfile3,std::ofstream &myfile4,
	 double &zerotime,const double &maxtime,const double &dt,Model& phi1,const bool & sign){
		 myfile<<"Hello: FeedB2"<<endl;
		 double gamma=phi1.GetGamma();
		 double temp2;
		 if(sign){//Fûteni kell
			 phi1.SetGamma(gamma*0.9);
			 myfile<<"Futottem! ";
		 }
		 else{
			 temp2=gamma*1.1;
			 phi1.SetGamma(temp2);
			 myfile<<"Hutottem! ";
		 }
		 a=0;//visszaküldjük energiafixálónak
		 myfile<<"a="<<a<<", aktualis ido="<<zerotime<<", Megyek EFixerhez"<<endl;
}

/*Tömegfix segéd: kiszámolja a tömeg aktuális értékét*/
double MCalc(int &a,
	std::ofstream &myfile, std::ofstream &myfile1,std::ofstream &myfile2,std::ofstream &myfile3,std::ofstream &myfile4,std::ofstream &myfile5,
	double& zerotime,const double &maxtime,const double &dt,Model& phi1,const double & mwanted, const double& errorlimit){
		double tstart; //1. zh idõpillanata
		double tact; //újonnan kapott zh
		double tfin;//n. zh idõpillanata
		int counter=0; //talált zh-ek száma
		double tint1=0; //aktuális két zh közti táv
		double tint100=0; //100 zh közti táv összege
		double actint; //aktuális intervallum
		bool sign; //kezdeti elõjele phimean-nek
		double T; //egy periódus ideje
		double M=-500;
		double data[100]; //intervallum hosszok
		double diff[100]; //eltérések az átlag intervall.hossztól
		double mom=(phi1.MomentMeas())[0];
		if(mom>0){
		sign=true;
		}
		else{
			sign=false;
		}

		myfile<<"Hello: MCalc"<<endl;

		tstart=FindZH(a,myfile,myfile1,myfile2,myfile3,myfile4,zerotime,maxtime,dt,phi1,sign);
		tfin=tstart;

		////counter=0 elõször begyûjtünk 5-öt ->5 intervallum
		//for(counter;(counter<5)&&(a!=6);counter++){
		//	tact=FindZH(a,myfile, myfile1,myfile2,myfile3,myfile4,zerotime,maxtime,dt,phi1,sign);
		//	tint1+=(tact-tfin);
		//	tfin=tact;
		//}//zh gyûjtés
		//if(a==6){
		//	myfile<<"a="<<a<<", kifutottunk az idobol es nem szamoltunk tomeget!"<<endl;
		//	return M;
		//}
		//tint1/=5; //5 intervallumból átlagolva egy tipikus intervall hossz
		//
		//myfile<<"Tipikus intervallum: "<<tint1<<endl;

		////majd begyûjtünk továbbiakat, max 25-öt
		//for(counter;counter<30;counter++){
		//	tact=FindZH(a,myfile,myfile1,myfile2,myfile3,myfile4,zerotime,maxtime,dt,phi1,sign);
		//	actint=(tact-tfin);
		//	if((((actint-tint1)/tint1)<0.05)&&(a!=6)){//ha tipikus intervallum
		//		tfin=tact;
		//	}
		//	else{//ha nem tipikus intervallum
		//		T=((tfin-tstart)*2/(counter+1));
		//		M=((2*M_PI)/T);
		//		if(((M-mwanted)/mwanted)<errorlimit){
		//			//jó a tömeg
		//			a=4;//exit
		//			myfile<<"a="<<a<<", aktualis ido="<<zerotime<<", Nem tipikus intevallum jött! Eddigi alapján JO a tomeg es M="<<M<<", Kesz vagyok!"<<endl;
		//			return M;
		//		}
		//		else{
		//			// nem jó még a tömeg ->feedback
		//			a=3;
		//			myfile<<"a="<<a<<", aktualis ido="<<zerotime<<", Nem tipikus intevallum jött! Eddigi alapján ROSSZ a tomeg es M="<<M<<", Megyek FeedBackhez!"<<endl;
		//			return M;
		//		}
		//	}
		//}//zh gyûjtés


		//100 zérushelyet gyûjtünk be
		//counter=0 elõször begyûjtünk 5-öt ->5 intervallum
		for(counter=0;(counter<100)&&(a!=7);counter++){
			tact=FindZH(a,myfile, myfile1,myfile2,myfile3,myfile4,zerotime,maxtime,dt,phi1,sign);
			tint1=(tact-tfin);
			data[counter]=tint1;
			tint100+=tint1;
			tfin=tact;
		}//zh gyûjtés

		

		if(a==7){
			myfile<<"a="<<a<<", kifutottunk az idobol es nem szamoltunk tomeget!"<<endl;
			return M;
		}

		T=tint100/50; //100 intervallumból átlagolva egy intervall hossz kétszerese

		myfile5<<endl;

		for(counter=0;counter<100;counter++){
			diff[counter]=abs(data[counter]-(tint100/100));
			myfile5<<data[counter]<<"\t"<<diff[counter]<<endl;
		}

		//T=((tfin-tstart)*2/(counter+1));
		M=((2*M_PI)/T);
				if((abs(M-mwanted)/mwanted)<errorlimit){
					//jó a tömeg
					a=4;//exit
					myfile<<"a="<<a<<", aktualis ido="<<zerotime<<", Lemértem! M="<<M<<", Kesz vagyok!"<<endl;
				}
				else{
					// nem jó még a tömeg ->feedback
					a=3;
					myfile<<"a="<<a<<", aktualis ido="<<zerotime<<", Lemértem! M="<<M<<", Megyek FeedBackhez!"<<endl;
				}
return M;
}

/*Tömegfix segéd: tömegértékek összehasonlítása, módosítja a zajt*/
void FeedBack(int &a,
	std::ofstream &myfile, std::ofstream &myfile1,std::ofstream &myfile2,std::ofstream &myfile3,std::ofstream &myfile4,
	double &zerotime,const double &maxtime,const double &dt,Model& phi1, 
	const double & mwanted, const double & mcalcul){
		double gamma=phi1.GetGamma();
		double difference=abs(mcalcul-mwanted)/mwanted;

		myfile<<"Hello: FeedBack"<<endl;

		if(mcalcul>mwanted){//hûteni kell
			difference+=1;
			myfile<<"Hutottem! ";
		}
		else{
			difference=1-difference;
			myfile<<"Futottem! ";
		}
		phi1.SetGamma(gamma*difference);
		a=0;
		myfile<<"a="<<a<<", aktualis ido="<<zerotime<<" Megyek EFixerhez!"<<endl;
}

/*MFixerConstE segéd: a kezdeti energiához próbálja igazítani az energiát*/
void EChanger(int &a,
	std::ofstream &myfile, std::ofstream &myfile1,std::ofstream &myfile2,std::ofstream &myfile3,std::ofstream &myfile4,
	double &zerotime,const double &maxtime,const double &dt,Model& phi1, const double & initE){
		double actualE=phi1.CountSumEnergy();
		double xiint=phi1.Getxiint();

		myfile<<"Hello: EChanger!"<<endl;

		//Fontos: Konfidencia int. szélesebb legyen mint efixernél

		if(actualE<(initE*0.98)){
			phi1.SetXiint(xiint*1.001);
			myfile<<"Xiint noveles!"<<endl;
			a=0;
			myfile<<"a="<<a<<", aktualis ido="<<zerotime<<" Megyek EFixerhez!"<<endl;
			return;
		}
		else if(actualE>(initE*1.02)){
			phi1.SetXiint(xiint*0.999);
			myfile<<"Xiint csokkentes!"<<endl;
			a=0;
			myfile<<"a="<<a<<", aktualis ido="<<zerotime<<" Megyek EFixerhez!"<<endl;
			return;
		}
		else{
			myfile<<"Energia rendben!"<<endl;
		}

		a=1;
		myfile<<"a="<<a<<", aktualis ido="<<zerotime<<" Megyek PhiTesterhez!"<<endl;
}

/*MFixerConstE segéd: tömegértékek összehasonlítása, módosítja m2-t*/
void MFeedBack(int &a,
	std::ofstream &myfile, std::ofstream &myfile1,std::ofstream &myfile2,std::ofstream &myfile3,std::ofstream &myfile4,
	double &zerotime,const double &maxtime,const double &dt,Model& phi1, 
	const double & mwanted, const double & mcalcul){
		double m2=phi1.GetM2();
		double difference=abs(mcalcul-mwanted)/mwanted;

		myfile<<"Hello: MFeedBack"<<endl;

		if(mcalcul>mwanted){//hûteni kell
			difference+=1;
			myfile<<"Hutottem! m2 növelve ";
		}
		else{
			difference=1-difference;
			myfile<<"Futottem! m2 csokkentve ";
		}
		phi1.SetM2(m2*difference);
		a=0;
		myfile<<"a="<<a<<", aktualis ido="<<zerotime<<" Megyek EFixerhez!"<<endl;
}

/*Tömeg fixáló fõprogram*/
void MFixer(string opendir,string dir, double randphi, double randpi, bool noise){
	
	//Modell létrehozása
	Model phi1(0.1,50,50,50,1,randphi,randpi,opendir.c_str());
	const char* infofile=(dir+"info.dat").c_str();
	
	
	double dt,maxtime,starttime,zerotime,actual;
	double errorlimit,Mwanted,Mcalc;
	double xi,gamma;
	char swal;
	char isConstEc;
	bool isConstE=false;
	int meastime,counter=0,a=0; //a: swithcer
	phi1.ModelInit();

	bool sign;
	if(phi1.MomentMeas()[0]>0){
		sign=true;
	} //phimean elõjele FeedB2 elõtt
	else{
		sign=false;
	}
	dt=phi1.Getdt();
	maxtime=phi1.GetMaxtime();
	meastime=phi1.GetMeasTime();

	cout<<"Ha korabbi merest folytatunk, mi a kezdeti idopont? (Random esetben 0)"<<endl;
	cin>>zerotime;

	actual=zerotime;

	cout<<"Mi az elvart tomegertek?"<<endl;
	cin>>Mwanted;

	cout<<"A beallitas milyen pontos legyen? (%)"<<endl;
	cin>>errorlimit;
	
	cout<<"Mennyi legyen xi,gamma?"<<endl;
	cin>>xi>>swal>>gamma;
	phi1.SetGamma(gamma);
	phi1.SetXiint(xi);

	cout<<"Konstans energián mérjünk?<y/n>"<<endl;
	cin>>isConstEc;
	if(isConstEc=='y'){
		isConstE=true;
	}
	const double Einit=phi1.CountSumEnergy();
	

	//info fájl nyitása
	std::ofstream myfile;
	myfile.open((dir+"info.dat").c_str(),std::ios::out);
	myfile<<"Racsallando: "<<phi1.GetGrid()<<endl;
	myfile<<"Tomegnegyzet: "<<phi1.GetM2()<<endl;
	myfile<<"Lambda: "<<phi1.GetL()<<endl;
	myfile<<"Idolepes: "<<phi1.Getdt()<<endl;
	myfile<<"Meresgyakorisag: "<<phi1.GetMeasTime();
	myfile<<"Meres vege: "<<phi1.GetMaxtime()<<endl;
	myfile<<"Randphi: "<<phi1.Getrandphi()<<endl;
	myfile<<"Randpi: "<<phi1.Getrandpi()<<endl;
	myfile<<"Elvart tomeg: "<<Mwanted<<endl;
	myfile<<"Tomegpontossag elvart hatara: "<<errorlimit<<"%"<<endl;
	myfile<<"Kezdeti xiint: "<<phi1.Getxiint()<<endl;
	myfile<<"(Kezdeti) gamma: "<<phi1.GetGamma()<<endl;
	myfile<<"Kezdeti energia: "<<Einit<<endl;

	errorlimit/=100;

	cout<<"A meresinicializalas kesz, kezdem a munkat."<<endl;

	//Fájlba írjuk miket mérünk
	myfile<<"Merem az osszenergia idofuggeset."<<endl;
	myfile<<"Merem phi varhato erteket es szorasat az ido fuggvenyeben."<<endl;
	myfile<<"Merem a homerseklet idofuggeset."<<endl;
	
	//kezdeti terek kiírása
	cout<<"A kezdeti teret nem irom ismet fajlba."<<endl;
	phi1.FFTModel();
	cout<<"Kezdeti FFT ter fajlba irasa..."<<endl;
	phi1.ToFileFFT((dir+"mivanbennemFFT_init.dat").c_str());
	cout<<"...kesz."<<endl;

	//kimeneti fájlok nyitása

	//idõ-energia
	std::ofstream myfile1;
	myfile1.precision(6);
	myfile1.open((dir+"Time_SumEnergy.dat").c_str());

	//idõ-hõmérséklet
	std::ofstream myfile2;
	myfile2.precision(6);
	myfile2.open((dir+"Ido_Homerseklet.dat").c_str());

	//idõ-momentumok
	std::ofstream myfile3;
	myfile3.precision(6);
	myfile3.open((dir+"time_phimean.dat").c_str());

	std::ofstream myfile4;
	myfile4.precision(6);
	myfile4.open((dir+"time_phi2nd.dat").c_str());
	
	//intervall hosszok és absz.eltérések
	std::ofstream myfile5;
	myfile5.precision(6);
	myfile5.open((dir+"intervalldata.dat").c_str());
	

	int count=0;

	vector<double> mom;
	//mérések és idõfejlõdés

	myfile<<endl;
	myfile<<"Hello! Indul az MFixer!"<<endl;

	do{
	switch(a){
	
	case 0:
		EFixer(a,myfile,myfile1,myfile2,myfile3,myfile4,actual,maxtime,dt,phi1);
		if((a==11)&&(!(isConstE))){ //ha nem konst energián akarunk mérni, nem kell meghívni EChanger-t
			a=1;
		}
		break;
	case 1:
		PhiTester(a,myfile,myfile1,myfile2,myfile3,myfile4,actual,maxtime,dt,phi1);
		break;

	case 2:
		Mcalc=MCalc(a,myfile,myfile1,myfile2,myfile3,myfile4,myfile5,actual,maxtime,dt,phi1,Mwanted,errorlimit);
		break;
	
	case 3:
		if(isConstE){
			MFeedBack(a,myfile,myfile1,myfile2,myfile3,myfile4,actual,maxtime,dt,phi1,Mwanted,Mcalc);
		}
		else{
			FeedBack(a,myfile,myfile1,myfile2,myfile3,myfile4,actual,maxtime,dt,phi1,Mwanted,Mcalc);
		}
		break;
	case 10:
			if(phi1.MomentMeas()[0]>0){
			sign=true;
		} //phimean elõjele FeedB2 elõtt
		else{
			sign=false;
		}
		FeedB2(a,myfile,myfile1,myfile2,myfile3,myfile4,actual,maxtime,dt,phi1,sign);
		break;
	case 11:
		EChanger(a,myfile,myfile1,myfile2,myfile3,myfile4,actual,maxtime,dt,phi1,Einit);
		break;
	
	}//switch vége
	}while(!(a>3 && a<10));

	myfile<<"****************"<<endl;
	myfile<<"*Viszlát MFixer*"<<endl;
	myfile<<"****************"<<endl;

	myfile<<"Végsõ gamma: "<<phi1.GetGamma()<<endl;
	myfile<<"Végsõ xiint: "<<phi1.Getxiint()<<endl;

	phi1.ToFile((dir+"valoster_evolved.dat").c_str());
	phi1.ToFileFFT((dir+"reciprokter_evolved.dat").c_str());

	myfile.close();
	myfile1.close();
	myfile2.close();
	myfile3.close();
	myfile4.close();
	myfile5.close();
}




/*Tömegfixáló indítása*/
void MFixerStartUp(){

	string opendir;
	string dir;
	double randphi;
	double randpi;
	bool noise=true;
	

	cout<<"Udvozol a TomegFixer Meres Program!\n"<<endl;
	cout<<"Ha random kezdeti feltetelt szeretnel, ird: RANDOM"<<endl;
	cout<<"Kulonben ird be a beolvasando fajlnevet eleresi uttal egyutt."<<endl;
	cout<<"Ide varom a valaszt: ";
	cin>>opendir;
	cout<<"Kerem a celkonyvtarat (win: vegen perjellel): "<<endl;
	cin>>dir;

	cout<<"Zajos meres!!!"<<endl;
	//cin>>noi;

	cout<<"Mi legyen randphi?"<<endl;
	cin>>randphi;
	cout<<"Mi legyen randpi?"<<endl;
	cin>>randpi;

	//ha random kezdeti feltételt szeretnénk
	if(!(opendir.compare("RANDOM"))){
		
		Model phirandom(0.1,50,50,50,1,randphi,randpi);
		phirandom.RandomFill();
		cout<<"Fajlba irom a generalt teret..."<<endl;
		phirandom.ToFile((dir+"randomgeninit.dat").c_str());
		cout<<"...fajlba iras kesz."<<endl;
		opendir=(dir+"randomgeninit.dat");
	}

	//mostantól jön a mérés.
	MFixer(opendir,dir,randphi,randpi,noise);


}

/*T12 és GTK mérése okosabban - EvolveModel nélkül*/
void DirectCombo(string opendir, string dir, double randphi, double randpi,bool noise){
	//Modell létrehozása
	Model phi1(0.1,50,50,50,1,randphi,randpi,opendir.c_str());
	const char* infofile=(dir+"info.dat").c_str();
	
	unsigned xmax=phi1.GetXmax();
	unsigned ymax=phi1.GetYmax();
	unsigned zmax=phi1.GetZmax();
	double dt,maxtime,zerotime;
	double xi,gamma;
	int meastime,measmin,timepoint,tinc,counter=0;
	int bins;
	char swal;
	cout<<"Maxtime legyen oszthato dt-vel!"<<endl;
	phi1.ModelInit();

	dt=phi1.Getdt();
	maxtime=phi1.GetMaxtime();
	meastime=phi1.GetMeasTime();

	cout<<"Ha korabbi merest folytatunk, mi a kezdeti idopont? (Random esetben 0)"<<endl;
	cin>>zerotime;

	cout<<"Mekkora legyen a minimális ablakmeret? (dt egysegben) "<<endl;
	cin>>measmin;

	cout<<"Mennyivel noveljuk az ablakmeretet (dt egysegben)? "<<endl;
	cin>>tinc;

	cout<<"Hany kulonbozo ablakot vegyunk fel? "<<endl;
	cin>>timepoint;

	cout<<"Bin-ek száma? "<<endl;
	cin>>bins;

	cout<<"Mennyi legyen xi,gamma?"<<endl;
	cin>>xi>>swal>>gamma;
	phi1.SetGamma(gamma);
	phi1.SetXiint(xi);

	//info fájl nyitása
	std::ofstream myfile;
	myfile.open((dir+"info.dat").c_str(),std::ios::out);
	myfile<<"Racsallando: "<<phi1.GetGrid()<<endl;
	myfile<<"Tomegnegyzet: "<<phi1.GetM2()<<endl;
	myfile<<"Lambda: "<<phi1.GetL()<<endl;
	myfile<<"Idolepes: "<<phi1.Getdt()<<endl;
	myfile<<"Meresgyakorisag: "<<phi1.GetMeasTime()<<endl;
	myfile<<"Meres vege: "<<phi1.GetMaxtime()<<endl;
	myfile<<"Randphi: "<<phi1.Getrandphi()<<endl;
	myfile<<"Randpi: "<<phi1.Getrandpi()<<endl;
	myfile<<"Min. ablakmeret: "<<measmin<<endl;
	myfile<<"Ablaknoveles: "<<tinc<<endl;
	myfile<<"Ablakszam: "<<timepoint<<endl;
	myfile<<"Xiint: "<<phi1.Getxiint()<<endl;
	myfile<<"Gamma: "<<phi1.GetGamma()<<endl;

	myfile.close();
	cout<<"A meresinicializalas kesz, kezdem a munkat."<<endl;

	//kezdeti terek kiírása
	cout<<"A kezdeti teret nem irom ismet fajlba."<<endl;
	phi1.FFTModel();
	cout<<"Kezdeti FFT ter fajlba irasa..."<<endl;
	phi1.ToFileFFT((dir+"mivanbennemFFT_init.dat").c_str());
	cout<<"...kesz."<<endl;

	
	//a mérésekhez szükséges:

	//nyers adatoknak
	int xdim=phi1.GetXmax()*phi1.GetYmax()*(phi1.GetZmax()/2+1);
	array3<double> gtk(timepoint,xdim,6,sizeof(Complex));
	//kmeas-t nullára inicializáljuk
	for(int l=0;l<timepoint;l++){	
	for(int i=0;i<xdim;i++){
			for(int j=0;j<6;j++){
				gtk(l,i,j)=0;
			}//yfor
		}//xfor
	}//lfor

	//binelt adatoknak
		array3<double> gtk2(timepoint,xdim,7,sizeof(Complex));
	//kmeas-t nullára inicializáljuk
		for(int l=0;l<timepoint;l++){
		for(int i=0;i<xdim;i++){
			for(int j=0;j<7;j++){
				gtk2(l,i,j)=0;
			}//yfor
		}//xfor
		}//lfor
	

		std::ofstream myfile1;
		myfile1.precision(6);
		myfile1.open((dir+"GTKPropagator_Bin.dat").c_str());

		int *kx=new int[xmax];
		int *ky=new int[ymax];
		int *meastime2; //mérésgyakoriságok eltolása - vektor
		int *rownum2; //ennyiedik sort kell majd számolni
		vector<Model> Modelset1; //Modellek lesznek benne, melyeket az idõfejlõdés során kimentünk, ez az eltolás nélküli oszlop
		int maxDist; //legnagyobb ablak mérete
		int maxDistCounter=0; //Ennyiszer telt el maxDist-nyi idõ
				

		unsigned p=0;
		for(p;p<(xmax/2)+1;p++){
			kx[p]=p;
		}
		for(p;p<xmax;p++){
			kx[p]=(p-xmax);
		}

		p=0;
		for(p;p<(ymax/2)+1;p++){
			ky[p]=p;
		}
		for(p;p<ymax;p++){
			ky[p]=(p-ymax);
		}
		
		p=0;
		
		
		//pointerfoglalások ahány t-re +1 (t'-nek)
		meastime2= new int[timepoint+1];
		rownum2=new int[timepoint+1]; //a nulladik helyen az alapidõ
		rownum2[0]=0;
		//mérésgyakoriságok és belõlük rögtön counter
		meastime2[0]=0; //counter, a többi counteri
		for(int i=0;i<timepoint;i++){
			meastime2[i+1]=0-(measmin+(i*tinc));
			rownum2[i+1]=0; //kezdetben még egy tagot sem számoltunk
			}//for mérésidõk feltöltése

		maxDist=meastime2[timepoint];

		//idõfejlõdés
		for(double timer=zerotime;timer<maxtime;timer+=dt){

			//csak mérésgyakoriságonként töltünk fel
			for(int i=0;i<(timepoint+1);i++){
				//ellenõrizni kell, hogy a counter már elérte-e a nullát
				if(meastime2[i]>=0){
					if((meastime2[i]%meastime)==0){
						if(i==0){ //alapidõ mérés eltároljuk
							Modelset1.push_back(phi1);
							int size=Modelset1.size();
							
							//cout<<"letaroltunk egy modelt, ido: "<<timer<<endl;
							Modelset1[size-1].FFTModel(); //rögtön elõállítjuk a Fourier teret, piFFT-ben az van ami kell


							//A letarolt modellben a pifield tárolja a T12-t
							//T12 elõállítása
								for(unsigned I=0;I<xmax;I++){
									for (unsigned J=0;J<ymax;J++){
										for (unsigned K=0;K<zmax;K++){
											double element=Modelset1[size-1].T12(I,J,K);
											Modelset1[size-1].SetPiFieldElement(I,J,K,element);

										}//Kfor
									}//Jfor
								}//Ifor

							//	cout<<"Letaroltuk T12-ket"<<endl;
								//T12-k fourier trafója
								rcfft3d ForwardPhi(Modelset1[size-1].GetZmax(),Modelset1[size-1].GetPIField(),Modelset1[size-1].GetdataFFTfield());
								ForwardPhi.fft(Modelset1[size-1].GetPIField(),Modelset1[size-1].GetdataFFTfield());
								cout<<"FFT data field-et is letaroltuk!"<<timer<<endl;
						}
						else{ //ha valamelyik ablak értéke jön, akkor ki is értékelünk
							
							Complex temp=0;
							Complex temp2=0;
							//aktuális fourier tér
							phi1.FFTModel();
							int actrow=rownum2[i];
							array3<double> T12array0(xmax,ymax,zmax,Modelset1[0].GeTAlign3());
							//T12 elõállítása
								for(unsigned I=0;I<xmax;I++){
									for (unsigned J=0;J<ymax;J++){
										for (unsigned K=0;K<zmax;K++){

											T12array0[I][J][K]=phi1.T12(I,J,K);
			
										}//Kfor
									}//Jfor
								}//Ifor


								rcfft3d ForwardPhi(Modelset1[0].GetZmax(),T12array0,phi1.GetdataFFTfield());
								ForwardPhi.fft(T12array0,phi1.GetdataFFTfield());

							//szamolas
							
							//bejárjuk a k teret és számoljuk sum:T12(t+t',k)*T12(*t',k)-t
							//valamint phik2k-t
							int idx=0;
							for(unsigned I=0;I<xmax;I++){
								for(unsigned J=0;J<ymax;J++){
									for(unsigned kz=0;kz<((zmax/2)+1);kz++){
				
										// ha létezõ index
										if(idx<(int)gtk.Ny()){					

											temp=phi1.GetDataelement(I,J,kz)*conj(Modelset1[actrow].GetDataelement(I,J,kz))*(meastime*dt);
											gtk(i-1,idx,1)+=real(temp);//ReT12
											gtk(i-1,idx,2)+=imag(temp);//ImT12

											temp2=phi1.GetphiFFTelement(I,J,kz)*conj(Modelset1[actrow].GetphiFFTelement(I,J,kz))*(meastime*dt);
											gtk(i-1,idx,4)+=real(temp2);
											gtk(i-1,idx,5)+=imag(temp2);
											gtk(i-1,idx,0)+=(phi1.CountPhik2(I,J,kz)); //Phikabs2
											gtk(i-1,idx,3)+=(4*sin(M_PI*kx[I]/xmax)*sin(M_PI*kx[I]/xmax))+(4*sin(M_PI*ky[J]/ymax)*sin(M_PI*ky[J]/ymax))+(4*sin(M_PI*kz/zmax)*sin(M_PI*kz/zmax));
											idx++;
										}
										else{
											cout<<"HIBAS INDEX"<<endl;
										}

									}//kzfor
								}//kyfor
							}//kxfor


							//noveljuk a sorszámot
							rownum2[i]++;
						//	cout<<timer<<" "<<(i*tinc+measmin)<<endl;
						}//szamolas vege
		
					}//ha counteri mérés van
				}//ha counteri>=0
			}//for vajon counteri >=0?
//			cout<<timer<<endl;
			//counterek növelése
			for(int i=0;i<(timepoint+1);i++){
				meastime2[i]+=1;
			}//counterek növelése
			
			//fejlesztünk
		if(noise){
			phi1.RefreshNoiseField();
		}
		else{
			phi1.RefreshField();	
		}
			
		}//idõfejlõdés vége
cout<<"atlagolas jon"<<endl;
	//Átlagolás
		
	for(int l=0;l<timepoint;l++){	
	for(int i=0;i<xdim;i++){
				gtk(l,i,0)/=rownum2[l+1];
				gtk(l,i,1)/=(rownum2[l+1]*dt*meastime);
				gtk(l,i,2)/=(rownum2[l+1]*dt*meastime);
				gtk(l,i,4)/=(rownum2[l+1]*dt*meastime);
				gtk(l,i,5)/=(rownum2[l+1]*dt*meastime);
				gtk(l,i,3)/=rownum2[l+1];
		}//xfor
//cout<<"l "<<l<<endl;
	}//lfor
		



		int idx;
		double sitemax=4*sin(M_PI/2)*sin(M_PI/2)+4*sin(M_PI/2)*sin(M_PI/2)+4*sin(M_PI/2)*sin(M_PI/2);
		
		//kitöltjük az elsõ sort a site (határ)értékekkel
		for(int l=0;l<timepoint;l++){
		for(int n=0;n<xdim;n++){
			gtk2(l,n,0)=sitemax/bins*(n+1);
		}
		}//lfor

cout<<"binel"<<endl;		
		for(int l=0;l<timepoint;l++){
		for(int n=0;n<xdim;n++){

			idx=bins*gtk(l,n,3)/sitemax;
			gtk2(l,idx,1)+=gtk(l,n,0);//phi
			gtk2(l,idx,2)+=gtk(l,n,1);//re T12
			gtk2(l,idx,3)+=gtk(l,n,2);//im T12
			gtk2(l,idx,4)++; //counter
			gtk2(l,idx,5)+=gtk(l,n,4);//re PhikPhik
			gtk2(l,idx,6)+=gtk(l,n,5);//im PhikPhik

		}
//cout<<"l "<<l<<endl;
		}//lfor
		

	cout<<"...bin-elve rendeztuk. Atlagolas es fajlba iras kovetkezik..."<<endl;
	
	for(int l=0;l<timepoint;l++){
cout<<"lciklus "<<l<<"tp "<<timepoint<<endl;
	for(int n=0;n<xdim;n++){
		if(gtk2(l,n,4)!=0){
		gtk2(l,n,1)/=gtk2(l,n,4);
		gtk2(l,n,2)/=gtk2(l,n,4);
		gtk2(l,n,3)/=gtk2(l,n,4);
		gtk2(l,n,5)/=gtk2(l,n,4);
		gtk2(l,n,6)/=gtk2(l,n,4);
		myfile1<<measmin+l*tinc<<"\t"<<gtk2(l,n,0)<<"\t"<<gtk2(l,n,1)<<"\t"<<gtk2(l,n,2)<<"\t"<<gtk2(l,n,3)<<"\t"<<gtk2(l,n,5)<<"\t"<<gtk2(l,n,6)<<endl;
		}//if
		
	}//for
	}//lfor
	
	
	

	//kmeas-t nullára inicializáljuk
	for(int l=0;l<timepoint;l++){	
	for(int i=0;i<xdim;i++){
			for(int j=0;j<(4+2);j++){
				gtk(l,i,j)=0;
			}//yfor
		}//xfor
	}//lfor

	for(int l=0;l<timepoint;l++){
		for(int i=0;i<xdim;i++){
			for(int j=0;j<7;j++){
				gtk2(l,i,j)=0;
			}//yfor
		}//xfor
		}//lfor


	myfile1.close();

	
	cout<<"Minden ablakra kesz!\nFajlba irjuk a vegso valos es reciprokteret!"<<endl;
	phi1.ToFile((dir+"valoster_evolved.dat").c_str());
	phi1.ToFileFFT((dir+"reciprokter_evolved.dat").c_str());

	cout<<"...fajlbairas kesz.!"<<endl;
	delete [] kx;
	delete [] ky;
	delete [] meastime2;
	delete [] rownum2;
}

void DirectComboStartUp(){
	string opendir;
	string dir;
	double randphi;
	double randpi;
	bool noise=false;
	char noi;


	cout<<"Udvozol a Combo Meres Program!\n"<<endl;
	cout<<"Ha random kezdeti feltetelt szeretnel, ird: RANDOM"<<endl;
	cout<<"Kulonben ird be a beolvasando fajlnevet eleresi uttal egyutt."<<endl;
	cout<<"Ide varom a valaszt: ";
	cin>>opendir;
	cout<<"Kerem a celkonyvtarat (win: vegen visszaperjellel): "<<endl;
	cin>>dir;

	cout<<"Zajos meres? <y,n>"<<endl;
	cin>>noi;

	if(noi=='y'){
		noise=true;
	}

	cout<<"Mi legyen randphi?"<<endl;
	cin>>randphi;
	cout<<"Mi legyen randpi?"<<endl;
	cin>>randpi;
	

	//ha random kezdeti feltételt szeretnénk
	if(!(opendir.compare("RANDOM"))){
		
		Model phirandom(0.1,50,50,50,1,randphi,randpi);
		phirandom.RandomFill();
		cout<<"Fajlba irom a generalt teret..."<<endl;
		phirandom.ToFile((dir+"randomgeninit.dat").c_str());
		cout<<"...fajlba iras kesz."<<endl;
		opendir=(dir+"randomgeninit.dat");
	}

	//mostantól jön a mérés.
		DirectCombo(opendir,dir,randphi,randpi,noise);
}
//Fõprogram

int main(){
	int nr;
	try {
		
		cout<<"*******************************************"<<endl;
		cout<<"** Udvozol a Phi4 modell meresi program! **"<<endl;
		cout<<"*******************************************"<<endl;
		do{
		cout<<"**********"<<endl;
		cout<<"** Menu **"<<endl;
		cout<<"**********"<<endl;
		cout<<"\n\n"<<endl;
		cout<<"Tesztek:"<<endl;
		cout<<"1: FFT tesztelese"<<endl;
		cout<<"2: Random kezdoallapot tesztelese"<<endl;
		cout<<"3: Fajlbol feltoltes tesztelese"<<endl;
		cout<<"Meresek: "<<endl;
		cout<<"4: Cellaterfogat novelesenek vizsgalata"<<endl;
		cout<<"5: Meresek egyetlen idopontban"<<endl;
		cout<<"6: Meresek idofejlodessel"<<endl;
		cout<<"7: Randomint-hom fugges"<<endl;
		cout<<"8: G(t,k_) meres"<<endl;
		cout<<"9: G(t,k_) PHI4 meres"<<endl;
		cout<<"10: Egyszeri viszkozitas meres"<<endl;
		cout<<"13: 10+8 combo"<<endl;
		cout<<"14: 10+8 ComboExtra"<<endl;
		cout<<"15: MFixer"<<endl;
		cout<<"16: DirectCombo"<<endl;
		cout<<"11: Kilepes"<<endl;
		cout<<"\n\n"<<endl;
		cout<<"Valasztott tevekenyseg szama: "<<endl;

		cin >> nr;


		 switch (nr) {
		 case 1:
		  test_1();
		  cout<<"Vegeztunk, vissza a menube. Varok egy entert.\n\n\n\n\n\n\n"<<endl;
		 cin.ignore(); 
		cin.get();

		  break;

		 case 2:
		  test_2();
		  cout<<"Vegeztunk, vissza a menube. Varok egy entert.\n\n\n\n\n\n\n"<<endl;
		 cin.ignore(); 
		cin.get();

		  break;	

		 case 3:
		  test_3();
		  cout<<"Vegeztunk, vissza a menube. Varok egy entert.\n\n\n\n\n\n\n"<<endl;
		 cin.ignore(); 
		cin.get();

		  break;	

		 case 4:
		  test_4();
		  cout<<"Vegeztunk, vissza a menube. Varok egy entert.\n\n\n\n\n\n\n"<<endl;
		 cin.ignore(); 
		cin.get();

		  break;	

		 case 5:
		  test_5();
		  cout<<"Vegeztunk, vissza a menube. Varok egy entert.\n\n\n\n\n\n\n"<<endl;
		 cin.ignore(); 
		cin.get();

		  break;	
	
		  case 6:
			  test_6();
			  cout<<"Vegeztunk, vissza a menube. Varok egy entert.\n\n\n\n\n\n\n"<<endl;
		 cin.ignore(); 
		cin.get();

			 break;	
			case 7:
			  test_7();
			  cout<<"Vegeztunk, vissza a menube. Varok egy entert.\n\n\n\n\n\n\n"<<endl;
		cin.ignore(); 
		cin.get();

			 break;	
		case 8:
			  GTK();
			  cout<<"Vegeztunk, vissza a menube. Varok egy entert.\n\n\n\n\n\n\n"<<endl;
		cin.ignore(); 
		cin.get();

			 break;	
		case 9:
			  Gphi4k();
			  cout<<"Vegeztunk, vissza a menube. Varok egy entert.\n\n\n\n\n\n\n"<<endl;
		 cin.ignore(); 
		cin.get();

			 break;	
		case 10:
			 T12T12();
			  cout<<"Vegeztunk, vissza a menube. Varok egy entert.\n\n\n\n\n\n\n"<<endl;
		cin.ignore(); 
		cin.get();

			 break;	
		
		case 11:
			  break;

	case 12:
			  modeltest();
			  cout<<"Vegeztunk, vissza a menube. Varok egy entert.\n\n\n\n\n\n\n"<<endl;
		 cin.ignore(); 
		cin.get();

			 break;	
	case 14:
			 ComboE();
			  cout<<"Vegeztunk, vissza a menube. Varok egy entert.\n\n\n\n\n\n\n"<<endl;
		 cin.ignore(); 
		cin.get();

			 break;

	case 15:
			MFixerStartUp();
			  cout<<"Vegeztunk, vissza a menube. Varok egy entert.\n\n\n\n\n\n\n"<<endl;
		 cin.ignore(); 
		cin.get();

			 break;

	case 16:
			DirectComboStartUp();
			  cout<<"Vegeztunk, vissza a menube. Varok egy entert.\n\n\n\n\n\n\n"<<endl;
		 cin.ignore(); 
		cin.get();

			 break;
	default: cout<<"Nincs ilyen tevekenyseg!"<<endl;
   }//switch

		 
		}while(nr!=11);

		cout<<"*******************************************"<<endl;
		cout<<"**  Koszonjuk, hogy mert a programmal!   **"<<endl;
		cout<<"**               Viszlat!                **"<<endl;
		cout<<"*******************************************"<<endl;
		cin.ignore(); 
		cin.get();


	}catch(const char* a)
	 {
		cerr<<"Nagy a baj\n"<<a<<endl;
		cerr<<"Biztos jol adtal meg minden parametert?"<<endl;
		cerr<<"A program kilep! Varok egy entert!"<<endl;
		cin.ignore(); 
		cin.get();

	 }
	catch(...){
		cerr<<"Nagy a baj\n"<<endl;
		cerr<<"A program kilep! Varok egy entert!"<<endl;
		cin.ignore(); 
		cin.get();

	}
	#ifdef _CRTDBG_MAP_ALLOC
  _CrtDumpMemoryLeaks();  // ellenõrzi, hogy volt-e memóriaszivárgás
#endif
	  return 0;
}
