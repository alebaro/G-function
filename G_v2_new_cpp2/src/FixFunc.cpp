#include <iostream>
#include<complex.h>
#include <math.h>
#include<cmath>
using namespace std;
#include"FixFunc.h"
const double pi = 3.1415926535897;
const   complex<double> I(0.0,1.0);

void Fi_common_op(complex<double> *finalyp,
		complex<double> *finalym,double x, void* p){
	params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	double A, A2, B;
	complex<double> deltap;
	A=1-2*x*params.PiPf/params.Pi2;
	A2=pow(1-2*x*params.PiPf/params.Pi2,2);
	B=4*(-params.mpi*params.mpi/params.Pi2+x*(1-x)*params.Pf2/params.Pi2);
	deltap=decltype(deltap)(A2+B,params.eps);
	*finalyp=0.5*(A+sqrt(deltap));
	*finalym=0.5*(A-sqrt(deltap));
}

void L_peps_common_op(complex<double> *L_eps,complex<double> f,double x,void* p){
	params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex <double> first,second,third,Re_f,Im_f;
	Re_f=real(f);
	Im_f=imag(f);
	first=log(abs((1-x-f)/f));
	second=atan((1-x-Re_f)/(Im_f+params.eps));
	third=atan((Re_f)/(Im_f+params.eps));
	*L_eps=first+I*second+I*third;
}

void L_meps_common_op(complex<double> *L_eps,complex<double> f,double x,void* p){
	params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex <double> first,second,third,Re_f,Im_f;
	Re_f=real(f);
	Im_f=imag(f);
	first=log(abs((1-x-f)/f));
	second=atan((1-x-Re_f)/(Im_f-params.eps));
	third=atan((Re_f)/(Im_f-params.eps));
	*L_eps=first+I*second+I*third;
}

void F1_common_op(complex<double>* final,double x, void*p){
	params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> finalyp,finalym,fact1,fact2,Lp_eps,Lm_eps;
	Fi_common_op(&finalyp,&finalym,x,&params);
	L_peps_common_op(&Lp_eps,finalyp,x,&params);
	L_meps_common_op(&Lm_eps,finalym,x,&params);
	fact1=1.0/(pow(4*pi,2)*params.Pi2*(finalyp-finalym));
	fact2=Lp_eps-Lm_eps;
	*final=fact1*fact2;
}

double F1_real(double x, void* p){

	params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> final;
	F1_common_op(&final,x,&params);
	return real(final);
}

double F1_imag(double x, void* p){

	params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> final;
	F1_common_op(&final,x,&params);
	return imag(final);
}


double xF1_real(double x, void* p){

	params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> final;
	F1_common_op(&final,x,&params);
	return x*real(final);

}

double xF1_imag(double x, void* p){
	params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> final;
	F1_common_op(&final,x,&params);
	return x*imag(final);
}


double x3F1_real(double x, void* p){

	params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> final;
	F1_common_op(&final,x,&params);
	return pow(x,3)*real(final);
}

double x3F1_imag(double x, void* p){

	params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> final;
	F1_common_op(&final,x,&params);
	return pow(x,3)*imag(final);
}


void F2_common_op(complex<double>* final,double x, void*p){
	params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> finalyp,finalym,fact1,fact2,Lp_eps,Lm_eps;
	Fi_common_op(&finalyp,&finalym,x,&params);
	L_peps_common_op(&Lp_eps,finalyp,x,&params);
	L_meps_common_op(&Lm_eps,finalym,x,&params);
	fact1=1.0/(pow(4*pi,2)*params.Pi2*(finalyp-finalym));
	fact2=finalyp*Lp_eps-finalym*Lm_eps;
	*final=fact1*fact2;
}

double F2_real(double x, void* p){

	params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> final;
	F2_common_op(&final,x,&params);
	return real(final);
}

double F2_imag(double x, void* p){
	params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> final;
	F2_common_op(&final,x,&params);
	return imag(final);
}


double x2F2_real(double x, void* p){

	params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> final;
	F2_common_op(&final,x,&params);
	return real(final)*pow(x,2);

}

double x2F2_imag(double x, void* p){

	params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> final;
	F2_common_op(&final,x,&params);
	return imag(final)*pow(x,2);
}


void F3_common_op(complex<double>* final,double x, void*p){
	params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> finalyp,finalym,fact1,fact2,fact3,Lp_eps,Lm_eps;
	Fi_common_op(&finalyp,&finalym,x,&params);
	L_peps_common_op(&Lp_eps,finalyp,x,&params);
	L_meps_common_op(&Lm_eps,finalym,x,&params);
	fact1=1.0/(pow(4*pi,2)*params.Pi2);
	fact2=1-x;
	fact3=(pow(finalyp,2)*Lp_eps-pow(finalym,2)*Lm_eps)/(finalyp-finalym);
	//*final=fact1*(fact2+fact3);
	*final=fact1*fact3;

}

double F3_real(double x, void* p){

	params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> final;
	F3_common_op(&final,x,&params);
	return real(final);
}


double F3_imag(double x, void* p){

	params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> final;
	F3_common_op(&final,x,&params);
	return imag(final);

}


double xF3_real(double x, void* p){

	params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> final;
	F3_common_op(&final,x,&params);
	return x*real(final);
}


double xF3_imag(double x, void* p){

	params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> final;
	F3_common_op(&final,x,&params);
	return x*imag(final);

}

void F4_common_op(complex<double>* final,double x, void*p){
	params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> finalyp,finalym,fact1,fact2,fact3,Lp_eps,Lm_eps;
	//cout<<"F4_common_op"<<params.eps<<endl;
	Fi_common_op(&finalyp,&finalym,x,&params);
	L_peps_common_op(&Lp_eps,finalyp,x,&params);
	L_meps_common_op(&Lm_eps,finalym,x,&params);
	fact1=1.0/(pow(4*pi,2)*params.Pi2);
	//fact2=(1.-x)*(1.-x+2.*finalyp+2.*finalym)/2.;
	fact2=(1.-x)*(finalyp+finalym);
	fact3=(pow(finalyp,3)*Lp_eps-pow(finalym,3)*Lm_eps)/(finalyp-finalym);
	*final=fact1*(fact2+fact3);
}

double F4_real(double x, void* p){
	    params1d_int &params=*reinterpret_cast<params1d_int *>(p);
		complex<double> final;
		F4_common_op(&final,x,&params);
		//void L_peps_common_op(complex<double> *L_eps,complex<double> f,double x,void* p);
		return real(final);

}

double F4_imag(double x, void* p){
    params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> final;
	F4_common_op(&final,x,&params);
	//void L_peps_common_op(complex<double> *L_eps,complex<double> f,double x,void* p);
	return imag(final);
}

void F5_common_op(complex<double>* final,double x, void*p){
	params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> finalyp,finalym,fact1,fact2,fact3,fact4,Lp_eps,Lm_eps;
	Fi_common_op(&finalyp,&finalym,x,&params);
	L_peps_common_op(&Lp_eps,finalyp,x,&params);
	L_meps_common_op(&Lm_eps,finalym,x,&params);
	fact1=-1.0/(8*pow(pi,2));
	fact2=(1-x-finalym)*log(1-x-finalym);
	fact3=(1-x-finalyp)*log(1-x-finalyp);
	fact4=finalym*log(-finalym)+finalyp*log(-finalyp);
	*final=fact1*(fact2+fact3+fact4);
}

double F5_real(double x, void* p){
	    params1d_int &params=*reinterpret_cast<params1d_int *>(p);
		complex<double> final;
		F5_common_op(&final,x,&params);
		return real(final);

}

double F5_imag(double x, void* p){
    params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> final;
	F5_common_op(&final,x,&params);
	return imag(final);

}

double xF5_real(double x, void* p){
    params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> final;
	F5_common_op(&final,x,&params);
	return x*real(final);

}

double xF5_imag(double x, void* p){
    params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> final;
	F5_common_op(&final,x,&params);
	return x*imag(final);

}

void F6_common_op(complex<double>* final,double x, void*p){
	params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> finalyp,finalym,fact1,fact2,fact3,fact4,fact5,Lp_eps,Lm_eps;
	Fi_common_op(&finalyp,&finalym,x,&params);
	L_peps_common_op(&Lp_eps,finalyp,x,&params);
	L_meps_common_op(&Lm_eps,finalym,x,&params);
	fact1=-1.0/(pow(4*pi,2));
	//fact2=-(1-x)*(1-x+finalym+finalyp);
	fact2=-(1-x)*(finalym+finalyp);
	fact3=(pow(1-x,2)-pow(finalym,2))*log(1-x-finalym);
	fact4=(pow(1-x,2)-pow(finalyp,2))*log(1-x-finalyp);
	fact5=pow(finalym,2)*log(-finalym)+pow(finalyp,2)*log(-finalyp);
	*final=fact1*(fact2+fact3+fact4+fact5);
}

double F6_real(double x, void* p){
	    params1d_int &params=*reinterpret_cast<params1d_int *>(p);
		complex<double> final;
		F6_common_op(&final,x,&params);
		return real(final);
}

double F6_imag(double x, void* p){
    params1d_int &params=*reinterpret_cast<params1d_int *>(p);
	complex<double> final;
	F6_common_op(&final,x,&params);
	return imag(final);
}



