#include "Geval.h"
#include<iostream>
#include<stdlib.h>
#include"vectoralgebra.h"
#include<string>
using namespace std;
#include "Geval.h"
#include "Gvariables.h"
#include "FixFunc.h"
#include "Fix3DFunc.h"
#include <math.h>
#include<cmath>
#include<complex.h>
#include<numeric>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
const double pi = 3.1415926535897;
double eps=1e-8;
const   complex<double> I(0.0,1.0);

void g_munu(double *(metric_low)[4]){
for (int  i=0;i<4;i++){
	for (int j=0;j<4;j++){
		if(i==j){
			if(i==0)
				metric_low[i][j]=+1.;
			else
				metric_low[i][j]=-1.;
		}
		else
			metric_low[i][j]=0;

	}
}

}

void G_eval::Set_Boost_to_cm(){
	Boost_to_cmi[0][0]=gammai;
	Boost_to_cmi[0][1]=-betai_abs*ni[0]*gammai;
	Boost_to_cmi[0][2]=-betai_abs*ni[1]*gammai;
	Boost_to_cmi[0][3]=-betai_abs*ni[2]*gammai;
	Boost_to_cmi[1][0]=-gammai*betai_abs*ni[0];
	Boost_to_cmi[1][1]=1+(gammai-1)*pow(ni[0],2);
	Boost_to_cmi[1][2]=(gammai-1)*ni[0]*ni[1];
	Boost_to_cmi[1][3]=(gammai-1)*ni[0]*ni[2];
	Boost_to_cmi[2][0]=-gammai*betai_abs*ni[1];
	Boost_to_cmi[2][1]=-gammai*betai_abs*ni[0]*ni[1];
	Boost_to_cmi[2][2]=1+(gammai-1)*ni[1]*ni[1];
	Boost_to_cmi[2][3]=(gammai-1)*ni[0]*ni[1];
	Boost_to_cmi[3][0]=-gammai*betai_abs*ni[2];
	Boost_to_cmi[3][1]=(gammai-1)*ni[0]*ni[2];
	Boost_to_cmi[3][2]=(gammai-1)*ni[1]*ni[2];
	Boost_to_cmi[3][3]=1+(gammai-1)*ni[2]*ni[2];

	Boost_to_cmf[0][0]=gammaf;
	Boost_to_cmf[0][1]=-betaf_abs*nf[0]*gammaf;
	Boost_to_cmf[0][2]=-betaf_abs*nf[1]*gammaf;
    Boost_to_cmf[0][3]=-betaf_abs*nf[2]*gammaf;
    Boost_to_cmf[1][0]=-gammaf*betaf_abs*nf[0];
    Boost_to_cmf[1][1]=1+(gammaf-1)*pow(nf[0],2);
    Boost_to_cmf[1][2]=(gammaf-1)*nf[0]*nf[1];
    Boost_to_cmf[1][3]=(gammaf-1)*nf[0]*nf[2];
    Boost_to_cmf[2][0]=-gammaf*betaf_abs*nf[1];
    Boost_to_cmf[2][1]=-gammaf*betaf_abs*nf[0]*nf[1];
    Boost_to_cmf[2][2]=1+(gammaf-1)*nf[1]*nf[1];
    Boost_to_cmf[2][3]=(gammaf-1)*nf[0]*nf[1];
    Boost_to_cmf[3][0]=-gammaf*betaf_abs*nf[2];
    Boost_to_cmf[3][1]=(gammaf-1)*nf[0]*nf[2];
    Boost_to_cmf[3][2]=(gammaf-1)*nf[1]*nf[2];
    Boost_to_cmf[3][3]=1+(gammaf-1)*nf[2]*nf[2];
}

void G_eval:: Set_values(double Eii, double *Pii,double Efi, double *Pfi,double Li,double mpii){
	//cout<<"In set values"<<endl;
	Ei=Eii;
	Pi[0]=Pii[0];
	Pi[1]=Pii[1];
	Pi[2]=Pii[2];
	Ef=Efi;
	Pf[0]=Pfi[0];
	Pf[1]=Pfi[1];
	Pf[2]=Pfi[2];
	P4i[0]=Ei;
	P4i[1]=Pi[0];
	P4i[2]=Pi[1];
	P4i[3]=Pi[2];
	P4f[0]=Ef;
	P4f[1]=Pf[0];
	P4f[2]=Pf[1];
	P4f[3]=Pf[2];
	//define {\bf\beta}_i=+{\bf P}_i/E_i and similarly for {\bf \beta}_f.. the convention is opposite to the paper
	scal_prod_vec(3,1/Ei,Pi,betai);
	scal_prod_vec(3,1/Ef,Pf,betaf);
	betai_abs=pow(dot_product(betai,betai,3),0.5);
	betaf_abs=pow(dot_product(betaf,betaf,3),0.5);
	//cout<<Ef<<Pf[2]<<endl;
	//cout<<"betaf"<<betaf[0]<<betaf[1]<<betaf[2]<<endl;
	//cout<<"betai_abs"<<betai_abs<<endl;
	if(betai_abs!=0){
	scal_prod_vec(3,1/betai_abs,betai,ni);}
	else{
		ni[0]=0;
		ni[1]=0.;
		ni[2]=0.;
	}
	if(betaf_abs!=0){
	scal_prod_vec(3,1/betaf_abs,betaf,nf);}
	else{
		nf[0]=0;
		nf[1]=0.;
		nf[2]=0.;
	}
	//scal_prod_vec(3,1/betaf_abs,betaf,nf);
	//cout<<"ni versor="<<ni[0]<<ni[1]<<ni[2]<<endl;

	//cout<<"nf versor"<<nf[0]<<nf[1]<<nf[2]<<endl;
		//cout<<beta[0]<<beta[1]<<beta[2]<<endl;
	int length_betai=sizeof(betai)/sizeof(*betai);

     gammai=1/(pow(1-dot_product(betai,betai,length_betai),0.5));
     gammaf=1/(pow(1-dot_product(betaf,betaf,length_betai),0.5));
     //cout<<"if="<<gammai<<" "<<gammaf<<endl;
     cout<<betai[2]<<endl;

	Set_Boost_to_cm();
	matrix_prod_vec(4,Boost_to_cmi,P4i,P4icm);
	cout<<"P4icm  here"<<endl;
	cout<<P4icm[0]<<" "<<P4icm[1]<<" "<<P4icm[2]<<" "<<P4icm[3]<<endl;
	matrix_prod_vec(4,Boost_to_cmf,P4f,P4fcm);
	cout<<"P4fcm  here"<<endl;
	cout<<P4fcm[0]<<" "<<P4fcm[1]<<" "<<P4fcm[2]<<" "<<P4fcm[3]<<endl;

	L=Li;
	mpi=mpii;

	Eicm=P4icm[0];
	qicm2=Eicm*Eicm/4-mpi*mpi;


	Efcm=P4fcm[0];
	qfcm2=Efcm*Efcm/4-mpi*mpi;

	//cout<<"Eicm"<<Eicm<<" "<<mpi<<endl;
	//cout<<"Efcm"<<Efcm<<" "<<mpi<<endl;
	//cout<<qicm2<<" "<<gammai<<" "<<gammaf<<endl;

}

    void G_eval::  intGfunct(double* G_real, double* G_imag){
	double smooth_int,semian_int,result,sum;
	int nmax_old,nmax_new,nmax;
	double alpha=1.;
	double sum_old, sum_new,final_sum,semian_int_real,semian_int_imag;
	//cout<<vartest.test<<endl;
	//cout<<"test struct Gvariables"<<Gvariables.Ei<<endl;
	/*Initialize global variable L=lattice volume*/
	//L=LL;
	//mpi=mpii;
	/*cout<<"num iterations="<<num_iterations<<endl;*/
	/*cout<<"Calculating sum ..."<<endl;*/
	nmax_old=1;
	sum_old=1/(pow(L,3))*sum_over_k(alpha,nmax_old);
	nmax_new=nmax_old+3;
	sum_new=1/(pow(L,3))*sum_over_k(alpha,nmax_new);
	while(abs(sum_old-sum_new)>=1e-8){
		nmax_old=nmax_new;
		sum_old=sum_new;
		nmax_new=nmax_old+3;
		sum_new=sum_over_k(alpha,nmax_new);
	}
	final_sum=sum_new;
	/*cout<<"final sum="<<final_sum<<endl;
	cout<<"testing sum nmax="<<nmax_new<<endl;
	cout<<"Calculating semi analytical part .."<<endl;*/
	semi_an_int(&semian_int_real,&semian_int_imag);
	/*cout<<semian_int<<endl;
	cout<<"Calculating smooth 3D integral .."<<endl;*/
	/*Factor 1/(2*pi)^3 already included in the integral*/
	smooth_int=0;//integral_3d_total(alpha,"no_k");
	/*cout<<"smooth integral"<<endl;*/
	//!!!!!!! the line below is the correct one!!!!!!
	//*G_real=final_sum-smooth_int-semian_int_real;
	//!!!!!Line below to check the 3D integral!!!!!
	*G_real=smooth_int;
	*G_imag=-semian_int_imag;
	//return result;
}

void G_eval::test_fig8(){
	complex<double>int_1d;
	double Lambda0=1.,Lambda1=3.,Lambda2=6.;
	double c0=1.,c1=-35./27.,c2=8./27.;
	complex<double> IA0,IA1,IA2;
	IA_nu1nu2nu3(Lambda0,&IA0);
	IA_nu1nu2nu3(Lambda1,&IA1);
	IA_nu1nu2nu3(Lambda2,&IA2);
	cout<<"IA0="<<IA0<<endl;
	cout<<"IA1="<<IA1<<endl;
	cout<<"IA2="<<IA2<<endl;
	int_1d=c0*IA0+c1*IA1+c2*IA2;
	double alphabar=1./81.,final3d;
	string control="one_k";
	final3d=integral_3d_total(alphabar,control);
	double sumk,sumk5;
	int nmax=70,nmax5=5;
	sumk=sum_over_k_onek(alphabar,nmax);
	sumk5=sum_over_k_onek(alphabar,nmax5);
	cout<<"int 1d="<<int_1d<<endl;
	cout<<"int 3d="<<final3d<<endl;
	cout<<"sumk="<<sumk<<";nmax="<<nmax<<endl;
	cout<<"*****Sum-In="<<sumk-(final3d)<<endl;
	cout<<"*****Sum-In-sum5="<<sumk-(final3d)-sumk5<<endl;
	cout<<"*****final G="<<sumk-(real(int_1d)+final3d)<<endl;
	cout<<"*****final G with q="<<(sumk-(real(int_1d)+final3d))/(pow(qicm2,0.5)*pow(qfcm2,0.5))<<endl;
	cout<<"*****final G-sum5="<<sumk-(real(int_1d)+final3d)-sumk5<<endl;
	cout<<"qicm2="<<qicm2<<"qfcm2="<<qfcm2<<endl;





}

void G_eval::test_fig9(){
	complex<double>int_1d;
	double Lambda0=1.,Lambda1=3.,Lambda2=6.;
	double c0=1.,c1=-35./27.,c2=8./27.;
	complex<double> IA0,IA1,IA2;
	IA_nu1nu2nu3(Lambda0,&IA0);
	IA_nu1nu2nu3(Lambda1,&IA1);
	IA_nu1nu2nu3(Lambda2,&IA2);
	cout<<"IA0="<<IA0<<endl;
	cout<<"IA1="<<IA1<<endl;
	cout<<"IA2="<<IA2<<endl;
	int_1d=c0*IA0+c1*IA1+c2*IA2;
	//double alphabar=1./81.,final3d;
	//string control="one_k";
	//final3d=integral_3d_total(alphabar,control);
	double nmax=30,sumk;
	double alphabar=1./81.,final3d;
	string control="one_k";
	final3d=integral_3d_total(alphabar,control);
	sumk=sum_over_k_onek(alphabar,nmax);
	cout<<"int 1d="<<int_1d<<endl;
	cout<<"int 3d="<<final3d<<endl;
	cout<<"sum_over_k="<<sumk<<"; nmax="<<nmax<<endl;
	/*cout<<"*****final G="<<sumk-(real(int_1d)+final3d)<<endl;
	cout<<"qicm2="<<qicm2<<"qfcm2="<<qfcm2<<endl;*/
	complex<double> qicm,qfcm;
	qicm=pow(qicm2,0.5);
	qfcm=pow(qfcm2,0.5);
	cout<<"qicm2="<<qicm<<"qfcm2="<<pow(qfcm2,0.5)<<endl;





}

 void G_eval::test_int1d(double* int1d_real,double* int1d_imag,string control){
	 if(control=="no_k")
	 int1d_scal(int1d_real,int1d_imag);
	 else if (control=="one_k")
		 int1d_scal(int1d_real,int1d_imag);
 }

 void G_eval::test_int1d_fig8(complex<double>* int_1d){
	 double Lambda0=1.,Lambda1=3.,Lambda2=6.;
	 double c0=1.,c1=-35./27.,c2=8./27.;
	 complex<double> IA0,IA1,IA2;
	 IA_nu1nu2nu3(Lambda0,&IA0);
	 IA_nu1nu2nu3(Lambda1,&IA1);
	 IA_nu1nu2nu3(Lambda2,&IA2);
	 cout<<"IA0="<<IA0<<endl;
	 cout<<"IA1="<<IA1<<endl;
	 cout<<"IA2="<<IA2<<endl;
	 *int_1d=c0*IA0+c1*IA1+c2*IA2;


  }




 void G_eval::test_int3d_total(double* total, double alphabar,string control){
	 *total=integral_3d_total(alphabar,control);
	 //*total=integral_2d_total(alphabar,control);

 }

 void G_eval::test_sumk(double* sumk_real,double* sumk_imag){
 	 double alpha=1./81.,nmax=10.;
 	 *sumk_real=sum_over_k_onek(alpha,nmax);
 	 *sumk_imag=0.;
 }

double G_eval:: sum_over_k(double alpha,int nmax){
	double kx,ky,kz,k2,omegak,omegaPik,omegaPfk;
	double sumtot=0,den;
	int i,j,k;
	double k4[4],k4cmi[4],k4cmf[4],kcmi2,kcmf2,omegakcm;
	double cutoff,alphabar;
	alphabar=L/pow(2*pi,4)*alpha;
	for (i =-nmax;i<=nmax;i++){
			for(j=-nmax;j<=nmax;j++){
				for (k=-nmax;k<=nmax;k++){
					kx=i*2*pi/L;
					ky=j*2*pi/L;
					kz=k*2*pi/L;
					k2=pow(kx,2)+pow(ky,2)+pow(kz,2);
					omegak=pow(k2+mpi*mpi,0.5);
					k4[0]=omegak;
					k4[1]=kx;
					k4[2]=ky;
					k4[3]=kz;
					matrix_prod_vec(4,Boost_to_cmi,k4, k4cmi);
					matrix_prod_vec(4,Boost_to_cmf,k4, k4cmf);
					kcmi2=pow(k4cmi[1],2)+pow(k4cmi[2],2)+pow(k4cmi[3],2);
					kcmf2=pow(k4cmf[1],2)+pow(k4cmf[2],2)+pow(k4cmf[3],2);
					omegakcm=k4cmf[0];
					omegaPik=pow((Pi[0]-kx)*(Pi[0]-kx)+(Pi[1]-ky)*(Pi[1]-ky)+(Pi[2]-kz)*(Pi[2]-kz)+mpi*mpi,0.5);
					omegaPfk=pow((Pf[0]-kx)*(Pf[0]-kx)+(Pf[1]-ky)*(Pf[1]-ky)+(Pf[2]-kz)*(Pf[2]-kz)+mpi*mpi,0.5);
					den=1/(2*omegak*(pow(Ef-omegak,2)-pow(omegaPfk,2))*(pow(Ei-omegak,2)-pow(omegaPik,2)));
					cutoff=exp(-alphabar*(kcmi2-qicm2)*(kcmf2-qfcm2));
					sumtot=sumtot+den*cutoff;
					//sumtot=sumtot+den*expDamping(Ei,Ef,omegak,omegaPik,omegaPfk,alpha);
				}
			}
		}
	return sumtot;
}

double G_eval:: sum_over_k_onek(double alphabar,int nmax){
	double kx,ky,kz,k2,omegak,omegaPik,omegaPfk;
	double sumtot=0.,den;
	int i,j,k;
	double k4[4],k4cmi[4],k4cmf[4],kcmi2,kcmf2;
	double cutoff,num;
	double thetai,thetaf;
	for (i =-nmax;i<=nmax;i++){
			for(j=-nmax;j<=nmax;j++){
				for (k=-nmax;k<=nmax;k++){
					kx=i*2*pi/L;
					ky=j*2*pi/L;
					kz=k*2*pi/L;
					k2=kx*kx+ky*ky+kz*kz;
					omegak=pow(k2+mpi*mpi,0.5);
					k4[0]=omegak;
					k4[1]=kx;
					k4[2]=ky;
					k4[3]=kz;
					matrix_prod_vec(4,Boost_to_cmi,k4, k4cmi);
					matrix_prod_vec(4,Boost_to_cmf,k4, k4cmf);
					kcmi2=pow(k4cmi[1],2)+pow(k4cmi[2],2)+pow(k4cmi[3],2);
					kcmf2=pow(k4cmf[1],2)+pow(k4cmf[2],2)+pow(k4cmf[3],2);
					//omegakcm=k4cmf[0];
					/*cout<<"vec k3="<<k4[1]<<" "<<k4[2]<<" "<<k4[3]<<endl;
					cout<<"vec k43cmi="<<k4cmi[1]<<" "<<k4cmi[2]<<" "<<k4cmi[3]<<endl;
					cout<<"vec k43cmf="<<k4cmf[1]<<" "<<k4cmf[2]<<" "<<k4cmf[3]<<endl;
					cout<<"------------------"<<endl;*/
					omegaPik=pow(pow(Pi[0]-kx,2)+pow(Pi[1]-ky,2)+pow(Pi[2]-kz,2)+pow(mpi,2),0.5);
					omegaPfk=pow(pow(Pf[0]-kx,2)+pow(Pf[1]-ky,2)+pow(Pf[2]-kz,2)+pow(mpi,2),0.5);
					den=1/(2*omegak*(pow(Ef-omegak,2)-pow(omegaPfk,2))*(pow(Ei-omegak,2)-pow(omegaPik,2)));
					cutoff=exp(-alphabar*(kcmi2-qicm2)*(kcmf2-qfcm2));
					//if(kcmi2!=0){
					thetai=acos(k4cmi[3]/pow(kcmi2,0.5));
					thetaf=acos(k4cmf[3]/pow(kcmf2,0.5));
					num=4.*pi*pow(kcmi2,0.5)*pow(kcmf2,0.5)*ReYJ0(1,0,thetaf)*ReYJ0(1,0,thetai)*omegak;
					//cout<<"thetai"<<thetai<<"; kcmi2="<<pow(kcmi2,0.5)<<"; k4cmi[3]"<<k4cmi[3]<<endl;
					//cout<<"thetai"<<ReYJ0(1,0,thetai)<<"thetaf"<<ReYJ0(1,0,thetaf)<<"den="<<den<<endl;
					//cout<<"den="<<den<<endl;
					//cout<<"ratio="<<k4cmi[3]/pow(kcmi2,0.5)<<endl;
					/*cout<<"Y10"<<ReYJ0(1,0,thetai)<<";thetai="<<thetai<<";kcmi2="<<kcmi2<<endl;
					cout<<k4[1]<<k4[2]<<k4[3]<<endl;
					cout<<k4cmi[1]<<k4cmi[2]<<k4cmi[3]<<endl;
					cout<<"-------------"<<endl;*/
					sumtot=sumtot+num*den*cutoff;
				//}
					//cout<<"sumtot="<<sumtot<<"; num="<<num<<" "<<pow(kcmi2,0.5)<<" "<<pow(kcmf2,0.5)<<endl;
					//cout<<"Y10="<<ReYJ0(1,0,thetaf)<<" "<<ReYJ0(1,0,thetaf)<<" "<<omegak<<endl;
					//cout<<"sumtot="<<sumtot<<";kcmi2="<<kcmi2<<"; kcmf2="<<kcmf2<<"; H="<<cutoff<<endl;
					//sumtot=sumtot+den*expDamping(Ei,Ef,omegak,omegaPik,omegaPfk,alpha);
				}
			}
		}
	return sumtot/(pow(L,3));//*pow(qicm2*qfcm2,0.5));//*
}




/*Calculating 1dim integrals, real and imaginary part, making use
 * of the functions F_i(x) in the source file "FixFunc.cpp", the corresponding header
 * is "FixFunc.h". Integration performed with standard quadrature routines
 * of gsl.
 */
void G_eval::semi_an_int(double* semian_int_real,double* semian_int_imag){
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
	double error, result_real,result_imag;
	gsl_function Greal,Gimag;
	params1d_int params;
	params.Pi2=pow(Ei,2)-(pow(Pi[0],2)+pow(Pi[1],2)+pow(Pi[2],2));;
	params.Pf2=pow(Ef,2)-(pow(Pf[0],2)+pow(Pf[1],2)+pow(Pf[2],2));
	params.PiPf=Ei*Ef-Pi[0]*Pf[0]-Pi[1]*Pf[1]-Pi[2]*Pf[2];
	params.mpi=mpi;
	params.eps=1e-5;
	Greal.params=reinterpret_cast<void *>(&params);
	Gimag.params=reinterpret_cast<void *>(&params);
	Greal.function =&F1_real;
	Gimag.function=&F1_imag;
	gsl_integration_qags (&Greal,0.0,1.0,0,1e-8, 10000, w, &result_real, &error);
	gsl_integration_qags (&Gimag,0.0,1.0,0,1e-8, 10000, w, &result_imag, &error);
	*semian_int_real=result_real;
	*semian_int_imag=result_imag;

}



void G_eval:: int1d_scal(double * int1d_scal_real, double *int1d_scal_imag){
	/* no momentum insertions*/
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
	double error, result_real,result_imag;
	gsl_function Greal,Gimag;
	params1d_int params;
	params.Pi2=pow(Ei,2)-(pow(Pi[0],2)+pow(Pi[1],2)+pow(Pi[2],2));
	params.Pf2=pow(Ef,2)-(pow(Pf[0],2)+pow(Pf[1],2)+pow(Pf[2],2));
	params.PiPf=Ei*Ef-Pi[0]*Pf[0]-Pi[1]*Pf[1]-Pi[2]*Pf[2];
	params.mpi=mpi;
	params.eps=1./5000.;//1e-5;
	Greal.params=reinterpret_cast<void *>(&params);
	Gimag.params=reinterpret_cast<void *>(&params);
	Greal.function =&F1_real;
	Gimag.function=&F1_imag;
	gsl_integration_qags (&Greal,0.0,1.0,0,1e-8, 10000, w, &result_real, &error);
	gsl_integration_qags (&Gimag,0.0,1.0,0,1e-8, 10000, w, &result_imag, &error);
	*int1d_scal_real=result_real;
	*int1d_scal_imag=result_imag;

}

void G_eval::IA_nu1nu2nu3(double Lambda,complex<double>* IA){
	double I31_real, I31_imag, I32_real, I32_imag,I33_real,I33_imag;
	double I34_real, I34_imag,	I35_real, I35_imag, I36_real, I36_imag;
	/*one momentum insertion, to check with Mathematica, following conventions of
	 * covariant.pdf version 5, Eq.C37,C38,C39*/
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
	double error, result_11_real,result_11_imag,result_12_real,result_12_imag;
	gsl_function int1d_real,int1d_imag;
	params1d_int params;
	params.Pi2=pow(Ei,2)-(pow(Pi[0],2)+pow(Pi[1],2)+pow(Pi[2],2));;
	params.Pf2=pow(Ef,2)-(pow(Pf[0],2)+pow(Pf[1],2)+pow(Pf[2],2));
	params.PiPf=Ei*Ef-Pi[0]*Pf[0]-Pi[1]*Pf[1]-Pi[2]*Pf[2];
	params.mpi=Lambda;
	params.eps=1e-5;//1./500.;
	int1d_real.params=reinterpret_cast<void *>(&params);
	int1d_imag.params=reinterpret_cast<void *>(&params);

	int1d_real.function=&x3F1_real;
	int1d_imag.function=&x3F1_imag;
	gsl_integration_qags (&int1d_real,0.0,1.0,0,1e-8, 10000, w,
				              &I31_real, &error);
	gsl_integration_qags (&int1d_imag,0.0,1.0,0,1e-8, 10000, w,
					              &I31_imag, &error);

	//cout<<"I31"<<endl;
	int1d_real.function=&F4_real;
	int1d_imag.function=&F4_imag;
	gsl_integration_qags (&int1d_real,0.0,1.0,0,1e-6, 10000, w,
					              &I32_real, &error);
	gsl_integration_qags (&int1d_imag,0.0,1.0,0,1e-6, 10000, w,
						              &I32_imag, &error);
	//cout<<"I32"<<endl;

	int1d_real.function=&x2F2_real;
	int1d_imag.function=&x2F2_imag;
	gsl_integration_qags (&int1d_real,0.0,1.0,0,1e-6, 10000, w,
						              &I33_real, &error);
	gsl_integration_qags (&int1d_imag,0.0,1.0,0,1e-6, 10000, w,
							              &I33_imag, &error);
	//cout<<"I33"<<endl;

	int1d_real.function=&xF3_real;
	int1d_imag.function=&xF3_imag;
	gsl_integration_qags (&int1d_real,0.0,1.0,0,1e-7, 10000, w,
							              &I34_real, &error);
	gsl_integration_qags (&int1d_imag,0.0,1.0,0,1e-6, 10000, w,
								              &I34_imag, &error);
	//cout<<"I34"<<endl;

	int1d_real.function=&xF5_real;
	int1d_imag.function=&xF5_imag;
	gsl_integration_qags (&int1d_real,0.0,1.0,0,1e-8, 10000, w,
								              &I35_real, &error);
	gsl_integration_qags (&int1d_imag,0.0,1.0,0,1e-8, 10000, w,
									             &I35_imag, &error);

	int1d_real.function=&F6_real;
	int1d_imag.function=&F6_imag;
	gsl_integration_qags (&int1d_real,0.0,1.0,0,1e-8, 10000, w,
									              &I36_real, &error);
	gsl_integration_qags (&int1d_imag,0.0,1.0,0,1e-8, 10000, w,
										              &I36_imag, &error);

	/*out<<"I^{3,X}"<<endl;
	cout<<"I31 x3F1: real="<<I31_real<<" ;imag="<<I31_imag<<endl;
	cout<<"I32 F4: real="<<I32_real<<" ;imag="<<I32_imag<<endl;
	cout<<"I33 x2F2: real="<<I33_real<<" ;imag="<<I33_imag<<endl;
	cout<<"I34 xF3: real="<<I34_real<<" ;imag="<<I34_imag<<endl;
	cout<<"I35 xF5: real="<<I35_real<<" ;imag"<<I35_imag<<endl;
	cout<<"I36 F6: real="<<I36_real<<" ;imag"<<I36_imag<<endl;
	cout<<"------------"<<endl;*/
	complex<double> I000,I003,I303;
	I000=pow(Ef,3)*(I31_real+I*I31_imag)+pow(Ei,3)*(I32_real+I*I32_imag)
			+3.*pow(Ef,2)*Ei*(I33_real+I*I33_imag)+3.*pow(Ei,2)*Ef*(I34_real+I*I34_imag)+
			- 3./4.*Ef*(I35_real+I*I35_imag)- 3./4.*Ei*(I36_real+I*I36_imag);
	/*cout<<"test I000"<<endl;
	cout<<"I35 xF5: real="<<I35_real<<" ;imag"<<I35_imag<<" "<<3./4.*Ef<<endl;
	cout<<"I36 F6: real="<<I36_real<<" ;imag"<<I36_imag<<" "<<3./4.*Ei<<endl;*/
/*
	cout<<pow(Ef,3)*(I31_real+I*I31_imag)<<endl;
	cout<<pow(Ei,3)*(I32_real+I*I32_imag)<<endl;
	cout<<3.*pow(Ef,2)*Ei*(I33_real+I*I33_imag)<<endl;
	cout<<3.*pow(Ei,2)*Ef*(I34_real+I*I34_imag)<<endl;
	cout<<- 3./4.*Ef*(I35_real+I*I35_imag)<<endl;
	cout<<- 3./4.*Ei*(I36_real+I*I36_imag)<<endl;
	cout<<"------"<<endl;*/
	I003=-pow(Ef,2)*Pf[2]*(I31_real+I*I31_imag)-pow(Ei,2)*Pi[2]*(I32_real+I*I32_imag)+
			-(2.*Pf[2]*Ef*Ei+pow(Ef,2)*Pi[2])*(I33_real+I*I33_imag)+
			-(2.*Pi[2]*Ei*Ef+pow(Ei,2)*Pf[2])*(I34_real+I*I34_imag)+0.25*Pf[2]*(I35_real+I*I35_imag)+
			+0.25*Pi[2]*(I36_real+I*I36_imag);

	I303=Pf[2]*Pf[2]*Ef*(I31_real+I*I31_imag)+Pi[2]*Pi[2]*Ei*(I32_real+I*I32_imag)+
			+(2.*Pf[2]*Pi[2]*Ef+Pf[2]*Pf[2]*Ei)*(I33_real+I*I33_imag)+
			+(2.*Pf[2]*Pi[2]*Ei+Pi[2]*Pi[2]*Ef)*(I34_real+I*I34_imag)+
			+0.25*Ef*(I35_real+I*I35_imag)+0.25*Ei*(I36_real+I*I36_imag);
	complex<double> final,test;
	double betai_low=-betai[2],betaf_low=-betaf[2];
	//double betai_low=betai[2],betaf_low=betaf[2];
	//the one below reproduces max 1d integral..
	*IA=3.*(gammai*betai_low*gammaf*betaf_low*I000+gammai*gammaf*I303-gammai*gammaf*(betai_low+betaf_low)*I003);

	//I believe this one below is the correct one
	//*IA=3.*(gammai*betai_low*gammaf*betaf_low*I000+gammai*gammaf*I303+gammai*gammaf*(betai_low+betaf_low)*I003);



	/*formula below has been corrected*/
	//*IA=3.*(gammai*betai_low*gammaf*betaf_low*I000+gammai*gammaf*I303+gammai*gammaf*(betai_low+betaf_low)*I003);
	//test=I31_real+I*I31_imag;
	cout<<test<<endl;
	cout<<"betaiz"<<betai[2]<<"; betafz="<<betaf[2]<<endl;
	cout<<"gammai="<<gammai<<";gammaf="<<gammaf<<endl;
	cout<<"I000="<<I000<<"; I003="<<I003<<"; I303="<<I303<<endl;
	cout<<"final"<<final<<endl;
	/*cout<<"final"<<final<<endl;*/


	gsl_integration_workspace_free (w);


}



void
display_results (char *title, double result, double error)
{
  printf ("%s ==================\n", title);
  printf ("result = % .6f\n", result);
  printf ("sigma  = % .6f\n", error);
}

/*total integral without momentum insertions*/


void print_matrix(double (*m)[4]){
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			cout<<m[i][j]<<" ";
		}
		cout<<endl;
	}
}

void G_eval:: set_3d_integral_param(double alphabar, void*p ){
	params_int3d_total_loc &params=*reinterpret_cast<params_int3d_total_loc *>(p);
	params.eps=0.;//1e-7;
	params.Ei=Ei;
	params.Pi[0]=Pi[0];
	params.Pi[1]=Pi[1];
	params.Pi[2]=Pi[2];
	params.Ef=Ef;
	params.Pf[0]=Pf[0];
	params.Pf[1]=Pf[1];
	params.Pf[2]=Pf[2];
	params.alphabar=alphabar;
	params.mpi=mpi;
	params.qicm2=qicm2;
	params.qfcm2=qfcm2;
	for(int i=0;i<4;i++){
			  for(int j=0;j<4;j++){
			  params.Boost_to_cmi[i][j]=Boost_to_cmi[i][j];
			  params.Boost_to_cmf[i][j]=Boost_to_cmf[i][j];
			  }
		  }
}

double G_eval::integral_3d_total(double alphabar,string control ){
	  double res, err;
	  params_int3d_total_loc params;
	  set_3d_integral_param( alphabar, &params );
	  //cout<<"params pass"<<params.Pf[2]<<endl;
	  cout<<"print boost to cm"<<endl;
	  print_matrix(params.Boost_to_cmi);
	  print_matrix(params.Boost_to_cmf);


	  cout<<"alphabar here="<<alphabar<<endl;
	  double xl[3] = { 0, 0, 0 };
	  double xu[3] = { 0.9999, pi, 2*pi-.00001 };
	  /*I seem to reproduce MAx figure only for 0.99,0.999, 0.9999 upper  extremum..*/

	  const gsl_rng_type *T;
	  gsl_rng *r;
	  double (*integrand_3d)(double*,size_t , void*);
	  cout<<"k's in the numerator?";
	  if (control=="no_k"){
		  cout<<control<<endl;
		  integrand_3d=&f_total_no_k;}
	  else if (control=="one_k"){
		  cout<<control<<endl;
		 integrand_3d=&f_total_one_k ;
		  //integrand_3d=&f_one_k_firstLine;
		  //integrand_3d=&f_one_k_2ndLine;
		 }
	  else{
		  cout<<control<<endl;
		  integrand_3d=&f_total_no_k;}

	  gsl_monte_function F = {integrand_3d, 3, 0} ;
	  //gsl_monte_function F = {&f_total_one_k, 3, 0} ;
	  F.params=&params;

	  //cout<<"F address"<<&F<<"; integrand_3d address"<<&f_total_one_k<<endl;



	  size_t calls = 500000;

	  gsl_rng_env_setup ();

	  T = gsl_rng_default;
	  r = gsl_rng_alloc (T);
	  {
	      gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);

	      gsl_monte_vegas_integrate (&F, xl, xu, 3, 10000, r, s,
	                                 &res, &err);

	      //display_results ("vegas warm-up", res, err);

	      printf ("converging...\n");

	      do
	        {
	          gsl_monte_vegas_integrate (&F, xl, xu, 3, calls/5, r, s,
	                                     &res, &err);
	          printf ("result = % .6f sigma = % .6f "
	                  "chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
	        }
	      while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

	      display_results ("vegas final", res, err);

	      gsl_monte_vegas_free (s);
	    }

	    gsl_rng_free (r);
	    //cout<<params.qicm2<<" "<<params.qfcm2<<endl;
	    return res/pow(2*pi,3);///(pow(qicm2,0.5)*pow(qfcm2,0.5));
}

double G_eval::integral_2d_total(double alphabar,string control ){
	  double res, err;
	  params_int2d_total_loc params;
	  set_3d_integral_param( alphabar, &params );
	  //cout<<"params pass"<<params.Pf[2]<<endl;
	  cout<<"print boost to cm"<<endl;
	  print_matrix(params.Boost_to_cmi);

	  cout<<"alphabar here="<<alphabar<<endl;
	  double xl[3] = { 0, 0 };
	  double xu[3] = { 1., pi };
	  /*I seem to reproduce MAx figure only for 0.99,0.999, 0.9999 upper  extremum..*/

	  const gsl_rng_type *T;
	  gsl_rng *r;

	  gsl_monte_function F = {&f_total_one_k_2D, 2, 0} ;
	  F.params=&params;

	  size_t calls = 1000000;

	  gsl_rng_env_setup ();

	  T = gsl_rng_default;
	  r = gsl_rng_alloc (T);
	  {
	      gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);

	      gsl_monte_vegas_integrate (&F, xl, xu, 2, 10000, r, s,
	                                 &res, &err);

	      //display_results ("vegas warm-up", res, err);

	      printf ("converging...\n");

	      do
	        {
	          gsl_monte_vegas_integrate (&F, xl, xu, 2, calls/5, r, s,
	                                     &res, &err);
	          printf ("result = % .6f sigma = % .6f "
	                  "chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
	        }
	      while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

	      display_results ("vegas final", res, err);

	      gsl_monte_vegas_free (s);
	    }

	    gsl_rng_free (r);
	    //cout<<params.qicm2<<" "<<params.qfcm2<<endl;
	    return res/pow(2*pi,2);///(pow(qicm2,0.5)*pow(qfcm2,0.5));
}

