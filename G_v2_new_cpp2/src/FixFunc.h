#ifndef FixFunc_H
#define FixFunc_H
using namespace std;
#include<complex.h>
#include <math.h>
#include<cmath>

struct params1d_int{
	double Pi2,Pf2,PiPf,mpi,eps;
};

void Fi_common_op(complex<double> *finalyp,
		complex<double> *finalym,double x, void* p);

void L_peps_common_op(complex<double> *L_eps,complex<double> f,double x,void* p);
void L_meps_common_op(complex<double> *L_eps,complex<double> f,double x,void* p);


void F1_common_op(complex<double>* final,double x, void*p);
double F1_real(double x, void* p);/*to adjust*/
double F1_imag(double x, void* p);/*to adjust*/
double xF1_real(double x, void* p);/*to adjust*/
double xF1_imag(double x, void* p);/*to adjust*/
double x3F1_real(double x, void* p);/*to adjust*/
double x3F1_imag(double x, void* p);/*to adjust*/

void F2_common_op(complex<double>* final,double x, void*p);
double F2_real(double x, void* p);/*to adjust*/
double F2_imag(double x, void* p);/*to adjust*/
double x2F2_real(double x, void* p);/*to adjust*/
double x2F2_imag(double x, void* p);/*to adjust*/

void F3_common_op(complex<double>* final,double x, void*p);
double F3_real(double x, void* p);/*to adjust*/
double F3_imag(double x, void* p);/*to adjust*/
double xF3_real(double x, void* p);/*to adjust*/
double xF3_imag(double x, void* p);/*to adjust*/

void F4_common_op(complex<double>* final,double x, void*p);
double F4_real(double x, void* p);/*to adjust*/
double F4_imag(double x, void* p);/*to adjust*/

void F5_common_op(complex<double>* final,double x, void*p);
double F5_real(double x, void* p);/*to check*/
double F5_imag(double x, void* p);/*to check*/
double xF5_real(double x, void* p);/*to check*/
double xF5_imag(double x, void* p);/*to check*/

void F6_common_op(complex<double>* final,double x, void*p);
double F6_real(double x, void* p);/*to adjust*/
double F6_imag(double x, void* p);/*to adjust*/


#endif
