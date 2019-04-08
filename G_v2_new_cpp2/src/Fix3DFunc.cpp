#include "Fix3DFunc.h"
#include<math.h>
#include <stdlib.h>
#include<cmath>
#include"vectoralgebra.h"
using namespace std;
const double pi = 3.1415926535897;




double f_total_no_k(double *k,size_t dim, void* p){
	params_int3d_total_loc &params=*reinterpret_cast<params_int3d_total_loc *>(p);
	double kx,ky,kz,omegak,omegaPik,omegaPfk,integrand;
	double k4[4],k4cmi[4],k4cmf[4],kcmi2,kcmf2,cutoff;
	double firstLine,secondLine;
	double delta1,delta2,delta2a,delta2b;
	//complex<double> deltai,deltaf;
	double fact,jacobian,kp,k0,k0new,mpi;
	mpi=params.mpi;
	kp=k[0];
	k0new=kp/(1-kp);
	kx=k0new*sin(k[1])*cos(k[2]);
	ky=k0new*sin(k[1])*sin(k[2]);
	kz=k0new*cos(k[1]);
	omegak=pow(k0new*k0new+mpi*mpi,0.5);
	k4[0]=k4cmi[0]=omegak;
	k4[1]=k4cmi[1]=kx;
	k4[2]=k4cmi[2]=ky;
	k4[3]=k4cmi[3]=kz;
	//matrix_prod_vec( 4,params.Boost_to_cmi,k4, k4cmi);
	matrix_prod_vec( 4,params.Boost_to_cmf,k4, k4cmf);
	kcmi2=pow(k4cmi[1],2)+pow(k4cmi[2],2)+pow(k4cmi[3],2);
	kcmf2=pow(k4cmf[1],2)+pow(k4cmf[2],2)+pow(k4cmf[3],2);
	cutoff=exp(-params.alphabar*(kcmi2-params.qicm2)*(kcmf2-params.qfcm2));
	//cout<<params.alphabar<<" "<<params.qicm2<<endl;
	omegaPik=pow((params.Pi[0]-kx)*(params.Pi[0]-kx)+(params.Pi[1]-ky)*(params.Pi[1]-ky)+(params.Pi[2]-kz)*(params.Pi[2]-kz)+mpi*mpi,0.5);
	omegaPfk=pow((params.Pf[0]-kx)*(params.Pf[0]-kx)+(params.Pf[1]-ky)*(params.Pf[1]-ky)+(params.Pf[2]-kz)*(params.Pf[2]-kz)+mpi*mpi,0.5);
	delta1=1/(2*omegak*(pow(params.Ef-omegak,2)-pow(omegaPfk,2))*(pow(params.Ei-omegak,2)-pow(omegaPik,2)));
	delta2a=1/(2*omegaPfk*(pow(params.Ef+omegaPfk,2)-pow(omegak,2))*(pow(params.Ei-params.Ef-omegaPfk,2)-pow(omegaPik,2)));
	delta2b=1/(2*omegaPik*(pow(params.Ei+omegaPik,2)-pow(omegak,2))*(pow(params.Ef-params.Ei-omegaPik,2)-pow(omegaPfk,2)));
	delta2=delta2a+delta2b;
	/*calculating cutoff here..*/

	firstLine=(delta1+delta2)*(cutoff-1);//
    secondLine=delta2*cutoff;
	jacobian=sin(k[1])*pow(kp,2)/(pow(kp-1,4));
	//fact=1/pow(2*pi,3);
	integrand=firstLine-secondLine;
	//cout<<cutoff<<" "<<delta1<<" "<<delta2a<<" "<<delta2b<<endl;
	//cout<<firstLine<<" "<<secondLine<<endl;
	//cout<<"here.."<<integrand<<endl;
	return integrand*jacobian;
}


double ReYJ0(int J,int M,double theta){
	double final;
	if(J==0)
		return final= 1./pow(4*pi,0.5);
	else if(J==1)
		return final=1./2.*pow(3/pi,0.5)*cos(theta);
	else if(J==2)
		return final=1./4.*pow(5/pi,0.5)*(3*cos(theta)*cos(theta)-1);
	else if(J==3)
		return final=1./4.*pow(7/pi,0.5)*(5*pow(cos(theta),3)-3*cos(theta));
	else if (J==4)
		return final=3./16.*pow(1/pi,0.5)*(35*pow(cos(theta),4)-30*pow(cos(theta),2)+3);
	return 0;
}



void Deltas_Nsigmas(double Lambda, double kx,double ky, double kz, void* p,double* delta,
		double* Nsigma, double*Kr){
	double deltari,deltarf;
	double Nsigmaoi;
	double Nsigmaof;
	params_int3d_total_loc &params=*reinterpret_cast<params_int3d_total_loc *>(p);
	double k2=pow(kx,2)+pow(ky,2)+pow(kz,2);
	double omegak,omegaPik,omegaPfk;
	//cout<<kx<<" "<<ky<<" "<< kz <<";k2="<<k2<<" ;Lambda="<<Lambda<<endl;
	omegak=pow(k2+pow(Lambda,2),0.5);
	//cout<<kx<<" "<<ky<<" "<< kz <<";k2="<<k2<<" ;Lambda="<<Lambda<<"; omegak"<<omegak<<endl;
	omegaPik=pow(pow(params.Pi[0]-kx,2)+
			pow(params.Pi[1]-ky,2)+pow(params.Pi[2]-kz,2)+pow(Lambda,2),0.5);
	omegaPfk=pow(pow(params.Pf[0]-kx,2)+
			pow(params.Pf[1]-ky,2)+pow(params.Pf[2]-kz,2)+pow(Lambda,2),0.5);
	*delta=1/(2*omegak*(pow(params.Ef-omegak,2)-pow(omegaPfk,2))*(pow(params.Ei-omegak,2)-pow(omegaPik,2)));
	//cout<<"first"<<*delta<<endl;
	deltarf=1/(2*omegaPfk*(pow(params.Ef+omegaPfk,2)-pow(omegak,2))*(pow(params.Ei-params.Ef-omegaPfk,2)-pow(omegaPik,2)));
	deltari=1/(2*omegaPik*(pow(params.Ei+omegaPik,2)-pow(omegak,2))*(pow(params.Ef-params.Ei-omegaPik,2)-pow(omegaPfk,2)));
	double k4_vec_cmi[4],k4_vec_cmf[4],k2_cmi,k2_cmf,k4[4],k4i[4],k4f[4];
	double k4i_vec_cmi[4],k4i_vec_cmf[4],k4f_vec_cmi[4],k4f_vec_cmf[4];
	//cout<<Lambda<<"omegak"<<omegak<<";delta"<<*delta<<endl;
	//cout<<params.Ei<<" "<<params.Ef<<endl;
	//cout<<"delta0"<<*delta<<endl;

	//cout<<"Lambda"<<Lambda<<endl;

	k4[0]=omegak;
	k4i[0]=params.Ei+omegaPik;
	k4f[0]=params.Ef+omegaPfk;
	k4[1]=k4i[1]=k4f[1]=kx;
	k4[2]=k4i[2]=k4f[2]=ky;
	k4[3]=k4i[3]=k4f[3]=kz;
	//cout<<"--------"<<endl;
	//print_matrixh(params.Boost_to_cmf);
	//cout<<"--------"<<endl;


	matrix_prod_vec( 4,params.Boost_to_cmi,k4, k4_vec_cmi);
	matrix_prod_vec( 4,params.Boost_to_cmf,k4, k4_vec_cmf);

	matrix_prod_vec( 4,params.Boost_to_cmi,k4i, k4i_vec_cmi);
	matrix_prod_vec( 4,params.Boost_to_cmf,k4i, k4i_vec_cmf);

	matrix_prod_vec( 4,params.Boost_to_cmi,k4f, k4f_vec_cmi);
	matrix_prod_vec( 4,params.Boost_to_cmf,k4f, k4f_vec_cmf);






	double k2i_cmi,k2f_cmi,k2i_cmf,k2f_cmf;

	k2_cmi=pow(k4_vec_cmi[1],2)+pow(k4_vec_cmi[2],2)+pow(k4_vec_cmi[3],2);
	k2_cmf=pow(k4_vec_cmf[1],2)+pow(k4_vec_cmf[2],2)+pow(k4_vec_cmf[3],2);
	k2i_cmi=pow(k4i_vec_cmi[1],2)+pow(k4i_vec_cmi[2],2)+pow(k4i_vec_cmi[3],2);
	k2i_cmf=pow(k4i_vec_cmf[1],2)+pow(k4i_vec_cmf[2],2)+pow(k4i_vec_cmf[3],2);
	k2f_cmi=pow(k4f_vec_cmi[1],2)+pow(k4f_vec_cmi[2],2)+pow(k4f_vec_cmi[3],2);
	k2f_cmf=pow(k4f_vec_cmf[1],2)+pow(k4f_vec_cmf[2],2)+pow(k4f_vec_cmf[3],2);
	//cout<<"k2_cmi="<<k2_cmi<<endl;
	//cout<<"k2i_cmi="<<k2i_cmi<<"; k2i_cmf"<<k2i_cmf<<"; k2f_cmi="<<k2f_cmi<<"; k2f_cmf="<<k2f_cmf<<endl;

	double thetai=acos(k4_vec_cmi[3]/pow(k2_cmi,0.5));
	double thetaf=acos(k4_vec_cmf[3]/pow(k2_cmf,0.5));
	//cout<<"argument thetaf"<<k4_vec_cmf[3]/pow(k2_cmf,0.5)<<endl;
	//cout<<"h="<<k2_cmi<<" "<<k2_cmf<<" "<<thetai<<" "<<thetaf<<endl;

	double thetaoii=acos(k4i_vec_cmi[3]/pow(k2i_cmi,0.5));
	double thetaoif=acos(k4i_vec_cmf[3]/pow(k2i_cmf,0.5));

	//cout<<kz<<" "<<k4i[3]<<" "<<k4_vec_cmi[3]<<endl;
	//cout<<"1="<<thetaoii<<" "<<thetaoif<<endl;

	double thetaofi=acos(k4f_vec_cmi[3]/pow(k2f_cmi,0.5));
	double thetaoff=acos(k4f_vec_cmf[3]/pow(k2f_cmf,0.5));
	//cout<<acos(k4f_vec_cmf[3]/pow(k2f_cmf,0.5))<<endl;
	//cout<<thetaoff<endl;

     //cout<<"2="<<thetaofi<<" "<<thetaoff<<endl;

	*Nsigma=4.*pi*pow(k2_cmf,0.5)*ReYJ0(1,0,thetaf)*k4[0]*ReYJ0(1,0,thetai)*pow(k2_cmi,0.5);
	//cout<<"Lambda="<<Lambda<<"; omegak="<<omegak<<"; k4[0]="<<k4[0]<<endl;
	//cout<<"delta="<<*delta<<endl;
	//cout<<"k2_cmi="<<k2_cmi<<"; k2_cmf="<<k2_cmf<<"; thetai="<<thetai<<"; thetaf="<<thetaf<<endl;
	//cout<<"Nsigma \omegak="<<*Nsigma<<"; kx,ky,kz="<<kx<<" "<<ky<<" "<<kz<<endl;
	//cout<<"thetai="<<thetai<<"; thetaf="<<thetaf<<endl;
	//cout<<pow(k2_cmf,0.5)<<" "<<ReYJ0(1,0,thetaf)<<" "<<k4[0]<<" "<<ReYJ0(1,0,thetai)<<" "<<pow(k2_cmi,0.5)<<endl;
	Nsigmaoi=4*pi*pow(k2i_cmi,0.5)*ReYJ0(1,0,thetaoif)*k4i[0]*ReYJ0(1,0,thetaoii)*pow(k2i_cmf,0.5);
	Nsigmaof=4*pi*pow(k2f_cmi,0.5)*ReYJ0(1,0,thetaoff)*k4f[0]*ReYJ0(1,0,thetaofi)*pow(k2f_cmf,0.5);
	//cout<<Nsigmaoi<<" "<<Nsigmaof<<endl
	//exit(EXIT_FAILURE);

	*Kr=Nsigmaoi*deltari+Nsigmaof*deltarf;
	//cout<<Lambda<<"delta0="<<*delta<<" "<<*Nsigma<<endl;



}

void print_matrixh(double (*m)[4]){
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			cout<<m[i][j]<<" ";
		}
		cout<<endl;
	}
}




double f_total_one_k(double *k,size_t dim, void* p){
	params_int3d_total_loc &params=*reinterpret_cast<params_int3d_total_loc *>(p);
	double kx,ky,kz,omegak,integrand;
	double k4[4],k4cmi[4],k4cmf[4],kcmi2,kcmf2,cutoff;
	//double delta1,delta2,delta2a,delta2b;
	//complex<double> deltai,deltaf;
	double fact,jacobian,kp,k0new,mpi=params.mpi;
	double Lambda0=mpi,Lambda1=3.*mpi, Lambda2=6.*mpi;
	double qicm20,qfcm20,qicm21,qfcm21,qicm22,qfcm22;
	double c0=1.,c1=-35./27.,c2=8./27.;
	kp=k[0];
	k0new=kp/(1-kp);
	kx=k0new*sin(k[1])*cos(k[2]);
	ky=k0new*sin(k[1])*sin(k[2]);
	kz=k0new*cos(k[1]);
	omegak=pow(pow(kx,2)+pow(ky,2)+pow(kz,2)+pow(mpi,2),0.5);

	k4[0]=omegak;
	k4[1]=kx;
	k4[2]=ky;
	k4[3]=kz;
	matrix_prod_vec( 4,params.Boost_to_cmi,k4, k4cmi);
	matrix_prod_vec( 4,params.Boost_to_cmf,k4, k4cmf);
	qicm20=k4cmi[0]*k4cmi[0]/4.-pow(Lambda0,2);
    qfcm20=k4cmf[0]*k4cmf[0]/4.-pow(Lambda0,2);

    qicm21=k4cmi[0]*k4cmi[0]/4.-pow(Lambda1,2);
    qfcm21=k4cmf[0]*k4cmf[0]/4.-pow(Lambda1,2);

    qicm22=k4cmi[0]*k4cmi[0]/4.-pow(Lambda2,2);
    qfcm22=k4cmf[0]*k4cmf[0]/4.-pow(Lambda2,2);






	kcmi2=pow(k4cmi[1],2)+pow(k4cmi[2],2)+pow(k4cmi[3],2);
	kcmf2=pow(k4cmf[1],2)+pow(k4cmf[2],2)+pow(k4cmf[3],2);
	cutoff=exp(-params.alphabar*(kcmi2-params.qicm2)*(kcmf2-params.qfcm2));
	/*double cutoff1,cutoff2;
	cutoff1=exp(-params.alphabar*(kcmi2-qicm21)*(kcmf2-qfcm21));
	cutoff2=exp(-params.alphabar*(kcmi2-qicm22)*(kcmf2-qfcm22));*/

	//cout<<"kcmi2="<<kcmi2<<endl;

	double delta0,delta1,delta2,Kr0,Kr1,Kr2;
	double Nsigma0ok,Nsigma1ok,Nsigma2ok;
	Deltas_Nsigmas(Lambda0,kx,ky,kz, p,&delta0,&Nsigma0ok,&Kr0);
	Deltas_Nsigmas(Lambda1,kx,ky,kz, p,&delta1,&Nsigma1ok,&Kr1);
	Deltas_Nsigmas(Lambda2,kx,ky,kz, p,&delta2,&Nsigma2ok,&Kr2);
	//cout<<Lambda0<<"delta0="<<delta0<<" "<<Nsigma0ok<<endl;
	//cout<<Lambda1<<"delta1="<<delta1<<"; Nomegak= "<<Nsigma1ok<<endl;
	//cout<<Lambda2<<"delta2="<<delta2<<"; Nomegak= "<<Nsigma2ok<<endl;



	double firstLine0,firstLine1,firstLine2,firstLine;
	double secondLine1,secondLine2,secondLine;
	double thirdLine;
	firstLine0=c0*(delta0*Nsigma0ok+Kr0);
	firstLine1=c1*(delta1*Nsigma1ok+Kr1);
	firstLine2=c2*(delta2*Nsigma2ok+Kr2);
	firstLine=(cutoff-1.)*(firstLine0+firstLine1+firstLine2);

	secondLine1=c1*delta1*Nsigma1ok;
	secondLine2=c2*delta2*Nsigma2ok;
	secondLine=-cutoff*(secondLine1+secondLine2);

	//secondLine=-cutoff1*secondLine1+cutoff2*secondLine2;
	//cout<<secondLine1<<" "<<secondLine2<<endl;
	//exit(EXIT_FAILURE);

	thirdLine=-cutoff*(c0*Kr0+c1*Kr1+c2*Kr2);
	//cout<<"k[1]"<<k[1]<<endl;
	integrand=firstLine+secondLine+thirdLine;//firstLine+secondLine+thirdLine;
	//+secondLine+thirdLine;//firstLine;//+secondLine+thirdLine;//+secondLine+thirdLine;//secondLine+thirdLine;//firstLine+secondLine+thirdLine;//secondLine;//+secondLine+thirdLine;

	jacobian=sin(k[1])*pow(kp,2)/(pow(kp-1.,4));
	//exit(EXIT_FAILURE);
	return integrand*jacobian;
}

double f_one_k_firstLine(double *k,size_t dim, void* p){
	params_int3d_total_loc &params=*reinterpret_cast<params_int3d_total_loc *>(p);
	double kx,ky,kz,omegak,integrand;
	double k4[4],k4cmi[4],k4cmf[4],kcmi2,kcmf2,cutoff;
	//double delta1,delta2,delta2a,delta2b;
	//complex<double> deltai,deltaf;
	double fact,jacobian,kp,k0new,mpi=params.mpi;
	double Lambda0=mpi,Lambda1=3.*mpi, Lambda2=6.*mpi;
	double c0=1.,c1=-35./27.,c2=8./27.;
	kp=k[0];
	k0new=kp/(1-kp);
	kx=k0new*sin(k[1])*cos(k[2]);
	ky=k0new*sin(k[1])*sin(k[2]);
	kz=k0new*cos(k[1]);
	omegak=pow(pow(kx,2)+pow(ky,2)+pow(kz,2)+pow(mpi,2),0.5);

	k4[0]=omegak;
	k4[1]=kx;
	k4[2]=ky;
	k4[3]=kz;
	matrix_prod_vec( 4,params.Boost_to_cmi,k4, k4cmi);
	matrix_prod_vec( 4,params.Boost_to_cmf,k4, k4cmf);


	kcmi2=pow(k4cmi[1],2)+pow(k4cmi[2],2)+pow(k4cmi[3],2);
	kcmf2=pow(k4cmf[1],2)+pow(k4cmf[2],2)+pow(k4cmf[3],2);
	cutoff=exp(-params.alphabar*(kcmi2-params.qicm2)*(kcmf2-params.qfcm2));

	double delta0,delta1,delta2,Kr0,Kr1,Kr2;
	double Nsigma0ok,Nsigma1ok,Nsigma2ok;
	Deltas_Nsigmas(Lambda0,kx,ky,kz, p,&delta0,&Nsigma0ok,&Kr0);
	Deltas_Nsigmas(Lambda1,kx,ky,kz, p,&delta1,&Nsigma1ok,&Kr1);
	Deltas_Nsigmas(Lambda2,kx,ky,kz, p,&delta2,&Nsigma2ok,&Kr2);
	double firstLine0,firstLine1,firstLine2,firstLine;
	double secondLine1,secondLine2,secondLine;
	double thirdLine;
	firstLine0=c0*(delta0*Nsigma0ok+Kr0);
	firstLine1=c1*(delta1*Nsigma1ok+Kr1);
	firstLine2=c2*(delta2*Nsigma2ok+Kr2);
	firstLine=(cutoff-1.)*(firstLine0+firstLine1+firstLine2);
	//cout<<"k[1]"<<k[1]<<endl;
	integrand=firstLine;//+secondLine+thirdLine;//firstLine;//+secondLine+thirdLine;//+secondLine+thirdLine;//secondLine+thirdLine;//firstLine+secondLine+thirdLine;//secondLine;//+secondLine+thirdLine;

	jacobian=sin(k[1])*pow(kp,2)/(pow(kp-1,4));
	return integrand*jacobian;
}

void Deltas_Nsigmas_2D(double Lambda, double kx,double ky, double kz, double k_par, void* p,double* delta,
		double* Nsigma, double*Kr){
	double deltari,deltarf;
	double Nsigmaoi;
	double Nsigmaof;
	params_int2d_total_loc &params=*reinterpret_cast<params_int2d_total_loc *>(p);
	double k2=pow(kz,2)+pow(k_par,2);
	double omegak,omegaPik,omegaPfk;
	//cout<<kx<<" "<<ky<<" "<< kz <<";k2="<<k2<<" ;Lambda="<<Lambda<<endl;
	omegak=pow(k2+pow(Lambda,2),0.5);
	//cout<<kx<<" "<<ky<<" "<< kz <<";k2="<<k2<<" ;Lambda="<<Lambda<<"; omegak"<<omegak<<endl;
	omegaPik=pow(pow(k_par,2)+pow(params.Pi[2]-kz,2)+pow(Lambda,2),0.5);
	omegaPfk=pow(pow(k_par,2)+pow(params.Pf[2]-kz,2)+pow(Lambda,2),0.5);
	*delta=1/(2*omegak*(pow(params.Ef-omegak,2)-pow(omegaPfk,2))*(pow(params.Ei-omegak,2)-pow(omegaPik,2)));
	//cout<<"first"<<*delta<<endl;
	deltarf=1/(2*omegaPfk*(pow(params.Ef+omegaPfk,2)-pow(omegak,2))*(pow(params.Ei-params.Ef-omegaPfk,2)-pow(omegaPik,2)));
	deltari=1/(2*omegaPik*(pow(params.Ei+omegaPik,2)-pow(omegak,2))*(pow(params.Ef-params.Ei-omegaPik,2)-pow(omegaPfk,2)));
	double k4_vec_cmi[4],k4_vec_cmf[4],k2_cmi,k2_cmf,k4[4],k4i[4],k4f[4];
	double k4i_vec_cmi[4],k4i_vec_cmf[4],k4f_vec_cmi[4],k4f_vec_cmf[4];
	//cout<<Lambda<<"omegak"<<omegak<<";delta"<<*delta<<endl;

	//cout<<"Lambda"<<Lambda<<endl;

	k4[0]=omegak;
	k4i[0]=params.Ei+omegaPik;
	k4f[0]=params.Ef+omegaPfk;
	k4[1]=k4i[1]=k4f[1]=kx;
	k4[2]=k4i[2]=k4f[2]=ky;
	k4[3]=k4i[3]=k4f[3]=kz;
	//cout<<"--------"<<endl;
	//print_matrixh(params.Boost_to_cmf);
	//cout<<"--------"<<endl;


	matrix_prod_vec( 4,params.Boost_to_cmi,k4, k4_vec_cmi);
	matrix_prod_vec( 4,params.Boost_to_cmf,k4, k4_vec_cmf);

	matrix_prod_vec( 4,params.Boost_to_cmi,k4i, k4i_vec_cmi);
	matrix_prod_vec( 4,params.Boost_to_cmf,k4i, k4i_vec_cmf);

	matrix_prod_vec( 4,params.Boost_to_cmi,k4f, k4f_vec_cmi);
	matrix_prod_vec( 4,params.Boost_to_cmf,k4f, k4f_vec_cmf);






	double k2i_cmi,k2f_cmi,k2i_cmf,k2f_cmf;

	k2_cmi=pow(k_par,2)+pow(k4_vec_cmi[3],2);
	k2_cmf=pow(k_par,2)+pow(k4_vec_cmf[3],2);
	k2i_cmi=pow(k_par,2)+pow(k4i_vec_cmi[3],2);
	k2i_cmf=pow(k_par,2)+pow(k4i_vec_cmf[3],2);
	k2f_cmi=pow(k_par,2)+pow(k4f_vec_cmi[3],2);
	k2f_cmf=pow(k_par,2)+pow(k4f_vec_cmf[3],2);
	//cout<<"k2_cmi="<<k2_cmi<<endl;
	//cout<<"k2i_cmi="<<k2i_cmi<<"; k2i_cmf"<<k2i_cmf<<"; k2f_cmi="<<k2f_cmi<<"; k2f_cmf="<<k2f_cmf<<endl;

	double thetai=acos(k4_vec_cmi[3]/pow(k2_cmi,0.5));
	double thetaf=acos(k4_vec_cmf[3]/pow(k2_cmf,0.5));
	//cout<<"h="<<k2_cmi<<" "<<k2_cmf<<" "<<thetai<<" "<<thetaf<<endl;

	double thetaoii=acos(k4i_vec_cmi[3]/pow(k2i_cmi,0.5));
	double thetaoif=acos(k4i_vec_cmf[3]/pow(k2i_cmf,0.5));

	//cout<<kz<<" "<<k4i[3]<<" "<<k4_vec_cmi[3]<<endl;
	//cout<<"1="<<thetaoii<<" "<<thetaoif<<endl;

	double thetaofi=acos(k4f_vec_cmi[3]/pow(k2f_cmi,0.5));
	double thetaoff=acos(k4f_vec_cmf[3]/pow(k2f_cmf,0.5));
	//cout<<acos(k4f_vec_cmf[3]/pow(k2f_cmf,0.5))<<endl;
	//cout<<thetaoff<endl;

     //cout<<"2="<<thetaofi<<" "<<thetaoff<<endl;

	*Nsigma=4*pi*pow(k2_cmf,0.5)*ReYJ0(1,0,thetaf)*k4[0]*ReYJ0(1,0,thetai)*pow(k2_cmi,0.5);
	//cout<<*Nsigma<<endl;
	//cout<<pow(k2_cmf,0.5)<<" "<<ReYJ0(1,0,thetaf)<<" "<<k4[0]<<" "<<ReYJ0(1,0,thetai)<<" "<<pow(k2_cmi,0.5)<<endl;
	Nsigmaoi=4*pi*pow(k2i_cmi,0.5)*ReYJ0(1,0,thetaoif)*k4i[0]*ReYJ0(1,0,thetaoii)*pow(k2i_cmf,0.5);
	Nsigmaof=4*pi*pow(k2f_cmi,0.5)*ReYJ0(1,0,thetaoff)*k4f[0]*ReYJ0(1,0,thetaofi)*pow(k2f_cmf,0.5);
	//cout<<Nsigmaoi<<" "<<Nsigmaof<<endl;
	*Kr=Nsigmaoi*deltari+Nsigmaof*deltarf;
	//cout<<Lambda<<"delta0="<<*delta<<" "<<*Nsigma<<endl;



}

double f_total_one_k_2D(double *k,size_t dim, void* p){
	params_int2d_total_loc &params=*reinterpret_cast<params_int2d_total_loc *>(p);
	double kx,ky,kz,k_par,omegak,integrand;
		double k4[4],k4cmi[4],k4cmf[4],kcmi2,kcmf2,cutoff;
		//double delta1,delta2,delta2a,delta2b;
		//complex<double> deltai,deltaf;
		double fact,jacobian,kp,k0new,mpi=params.mpi;
		double Lambda0=mpi,Lambda1=3.*mpi, Lambda2=6.*mpi;
		double c0=1.,c1=-35./27.,c2=8./27.;
		kp=k[0];
		k0new=kp/(1-kp);
		kx=k0new*sin(k[1])*cos(pi/3.);
		ky=k0new*sin(k[1])*sin(pi/3.);
		k_par=k0new*sin(k[1]);
		kz=k0new*cos(k[1]);
		omegak=pow(pow(kz,2)+pow(k_par,2)+pow(mpi,2),0.5);

	k4[0]=omegak;
	k4[1]=kx;
	k4[2]=ky;
	k4[3]=kz;
	matrix_prod_vec( 4,params.Boost_to_cmi,k4, k4cmi);
	matrix_prod_vec( 4,params.Boost_to_cmf,k4, k4cmf);


	kcmi2=pow(k4cmi[1],2)+pow(k4cmi[2],2)+pow(k4cmi[3],2);
	kcmf2=pow(k4cmf[1],2)+pow(k4cmf[2],2)+pow(k4cmf[3],2);
	cutoff=exp(-params.alphabar*(kcmi2-params.qicm2)*(kcmf2-params.qfcm2));

	double delta0,delta1,delta2,Kr0,Kr1,Kr2;
	double Nsigma0ok,Nsigma1ok,Nsigma2ok;
	Deltas_Nsigmas_2D(Lambda0,kx,ky,kz,k_par,p,&delta0,&Nsigma0ok,&Kr0);
	Deltas_Nsigmas_2D(Lambda1,kx,ky,kz,k_par,p,&delta1,&Nsigma1ok,&Kr1);
	Deltas_Nsigmas_2D(Lambda2,kx,ky,kz,k_par,p,&delta2,&Nsigma2ok,&Kr2);
	//cout<<Lambda0<<"delta0="<<delta0<<" "<<Nsigma0ok<<endl;
	/*cout<<Lambda1<<"delta1="<<delta1<<" "<<Nsigma1ok<<endl;
	cout<<Lambda2<<"delta2="<<delta2<<" "<<Nsigma2ok<<endl;*/


	double firstLine0,firstLine1,firstLine2,firstLine;
	double secondLine1,secondLine2,secondLine;
	double thirdLine;
	firstLine0=c0*(delta0*Nsigma0ok+Kr0);
	firstLine1=c1*(delta1*Nsigma1ok+Kr1);
	firstLine2=c2*(delta2*Nsigma2ok+Kr2);
	firstLine=(cutoff-1.)*(firstLine0+firstLine1+firstLine2);

	secondLine1=c1*delta1*Nsigma1ok;
	secondLine2=c2*delta2*Nsigma2ok;
	secondLine=-cutoff*(secondLine1+secondLine2);
	//cout<<"here"<<endl;

	thirdLine=-cutoff*(c0*Kr0+c1*Kr1+c2*Kr2);
	//cout<<"k[1]"<<k[1]<<endl;
	integrand=secondLine;//+secondLine+thirdLine;//firstLine;//+secondLine+thirdLine;//+secondLine+thirdLine;//secondLine+thirdLine;//firstLine+secondLine+thirdLine;//secondLine;//+secondLine+thirdLine;

	jacobian=sin(k[1])*pow(kp,2)/(pow(kp-1,4));
	return integrand*jacobian;
}


double f_one_k_2ndLine(double *k,size_t dim, void* p){
	params_int3d_total_loc &params=*reinterpret_cast<params_int3d_total_loc *>(p);
	double kx,ky,kz,k_par,omegak,integrand;
	double k4[4],k4cmi[4],k4cmf[4],kcmi2,kcmf2,cutoff;
	//double delta1,delta2,delta2a,delta2b;
	//complex<double> deltai,deltaf;
	double fact,jacobian,kp,k0new,mpi=params.mpi;
	double Lambda0=mpi,Lambda1=3.*mpi, Lambda2=6.*mpi;
	double c0=1.,c1=-35./27.,c2=8./27.;
	kp=k[0];
	k0new=kp/(1-kp);
	kx=k0new*sin(k[1])*cos(pi/3.);
	ky=k0new*sin(k[1])*sin(pi/3.);
	k_par=k0new*sin(k[1]);
	kz=k0new*cos(k[1]);
	omegak=pow(pow(kz,2)+pow(k_par,2)+pow(mpi,2),0.5);

	k4[0]=omegak;
	k4[1]=kx;
	k4[2]=ky;
	k4[3]=kz;
	matrix_prod_vec( 4,params.Boost_to_cmi,k4, k4cmi);
	matrix_prod_vec( 4,params.Boost_to_cmf,k4, k4cmf);


	kcmi2=pow(k4cmi[3],2)+pow(k_par,2);
	kcmf2=pow(k4cmf[3],2)+pow(k_par,2);
	cutoff=exp(-params.alphabar*(kcmi2-params.qicm2)*(kcmf2-params.qfcm2));

	double delta1,delta2,Kr1,Kr2;
	double Nsigma1ok,Nsigma2ok;
	Deltas_Nsigmas(Lambda1,kx,ky,kz, p,&delta1,&Nsigma1ok,&Kr1);
	Deltas_Nsigmas(Lambda2,kx,ky,kz, p,&delta2,&Nsigma2ok,&Kr2);

	double secondLine;
	secondLine=-cutoff*(c1*delta1*Nsigma1ok+c2*delta2*Nsigma2ok);

	integrand=secondLine;//+secondLine+thirdLine;//firstLine;//+secondLine+thirdLine;//+secondLine+thirdLine;//secondLine+thirdLine;//firstLine+secondLine+thirdLine;//secondLine;//+secondLine+thirdLine;

	jacobian=sin(k[1])*pow(kp,2)/(pow(kp-1,4));
	return integrand*jacobian;
}



