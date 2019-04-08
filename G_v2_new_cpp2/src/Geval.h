#ifndef GEVAL_H_
#define GEVAL_H_
#include<string>
#include <math.h>
#include<cmath>
#include<complex.h>
#include<stdlib.h>
using namespace std;
class G_eval {
	double Ei,Ef,L,mpi,beta_abs,P4i[4],P4f[4],ni[3],nf[3];
	double Eicm,Efcm,P4icm[4],P4fcm[4],qicm2,qfcm2;
	double betai[3],betaf[3],betai_abs,betaf_abs,gammai,gammaf;
	double Pi[3],Pf[3],Boost_to_cmi[4][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
	double Boost_to_cmf[4][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
	public:
	    void Set_Boost_to_cm();
	    void Set_values(double Eii, double *Pii,double Efi, double *Pfi,double Li,double mpii);
		void intGfunct(double* G_real,double* G_imag );
		void test_fig8();
		void test_fig9();
		void test_int1d(double* int1d_real,double* int1d_imag,string control);
		void test_int1d_fig8(complex<double>* int_1d);
		void test_int3d_total(double* total,double alphabar,string control);
		void test_sumk(double* sumk_real,double* sumk_imag);
		double sum_over_k(double alpha,int nmax);
		double sum_over_k_onek(double alphabar,int nmax);
		void semi_an_int(double* semian_int_real,double* semian_int_imag);
		void integral_test(int num_lorentz_indexes);
		void int1d_scal(double * int1d_scal_real, double *int1d_scal_imag);
		void IA_nu1nu2nu3(double Lambda,complex<double>* IA);
		//double integral( double alpha);
		void set_3d_integral_param(double alpha, void* p);
		double integral_3d_total(double alpha,string control );
		double integral_2d_total(double alphabar,string control );
		//double integral_total_one_k(double alpha );

};

#endif /* GEVAL_H_ */
