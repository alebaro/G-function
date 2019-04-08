#ifndef FIX3DFUNC_H_
#define FIX3DFUNC_H_
#include<iostream>
#include<math.h>
#include<cmath>
using namespace std;

struct params_int3d_total_loc{
	double Ei,Ef,mpi,eps,alphabar;
	double Boost_to_cmi[4][4],Boost_to_cmf[4][4];
	double Pi[3],Pf[3],qicm2,qfcm2;
};

struct params_int2d_total_loc{
	double Ei,Ef,mpi,eps,alphabar;
	double Boost_to_cmi[4][4],Boost_to_cmf[4][4];
	double Pi[3],Pf[3],qicm2,qfcm2;
};


double f_total_no_k(double *k,size_t dim, void* p);
double ReYJ0(int J,int M,double theta);
void Deltas_Nsigmas(double Lambda, double kx,double ky, double kz, void* p,double* delta,double* deltari,double* deltarf,
		double* Nsigma, double*Kr);
void print_matrixh(double (*m)[4]);
double f_total_one_k(double *k,size_t dim, void* p);
double f_one_k_firstLine(double *k,size_t dim, void* p);
double f_one_k_2ndLine(double *k,size_t dim, void* p);
void Deltas_Nsigmas_2D(double Lambda,  double kz,double k_par, void* p,double* delta,double* deltari,double* deltarf,
		double* Nsigma, double*Kr);
double f_total_one_k_2D(double *k,size_t dim, void* p);
//double f_total_one_k_2D(double *k,size_t dim, void* p);



#endif /* FIX3DFUNC_H_ */
