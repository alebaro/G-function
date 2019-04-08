//============================================================================
// Name        : full_G.cpp
// Authors     : AB RB MH FO DW
// Version     : v2
// Description : Code for evaluating G function. Output is the G function for lf=li=0.
//               Functions for 1dim integral in the file FixFunc.cpp (header
//               FixFunc.h), at the moment we have: F1,xF1,F2,F3 (real and imaginary).
// Execution   : Download IDE Eclipse Oxygen at http://www.eclipse.org/downloads/packages/release/Oxygen/M2,
//               download gsl library for mac at http://macappstore.org/gsl/,
//               add gsl at the project, build project, run project.
//============================================================================

#include <iostream>
#include <iostream>
#include<fstream>
#include<math.h>
#include<complex.h>
#include<string>
using namespace std;
//#include"GsClass.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
const double pi = 3.1415926535897;
#include"Geval.h"
#include"Gvariables.h"

void test_complex_square_roots(){
	double epsilon=1e-7;
	double a=1.2;
	complex <double> z1,z2,temp(-a,epsilon);
	z1=sqrt(temp);
	cout<<z1<<endl;
	//cout<<abs(z1)<<endl;
	cout<<real(z1)<<endl;
	cout<<imag(z1)<<endl;
	cout<<sqrt(z1)<<endl;
	cout<<log(10)<<endl;
	z2=decltype(z2)(1.,3.);
	cout<<z2<<endl;
	double b=2.01;
	complex<double> temph(1.,3.);
	cout<<"test complex product"<<pow(temph,2)<<endl;
}


void plot_Fig5_left_panel(){
	double G_real,G_imag,int3d_real,int3d_imag,Eicm,Efcm;
	Eicm=3.4;
	Efcm=2.1;
		Gvariables Gvariables;
		Gvariables.mpi=1;
		Gvariables.L=6.;
		double unit_mom=2*pi/Gvariables.L;
		//Gvariables.alpha=1.0;
		Gvariables.Pi[0]=0.;
		Gvariables.Pi[1]=0.;
		Gvariables.Pi[2]=0.;
		Gvariables.Pf[0]=0.;
		Gvariables.Pf[1]=0.;
		Gvariables.Pf[2]=unit_mom;
		Gvariables.Ei=pow(Eicm*Eicm+Gvariables.Pi[0]*Gvariables.Pi[0]+
						Gvariables.Pi[1]*Gvariables.Pi[1]+Gvariables.Pi[2]*Gvariables.Pi[2],0.5);
	    Gvariables.Ef=pow(Efcm*Efcm+Gvariables.Pf[0]*Gvariables.Pf[0]+
				Gvariables.Pf[1]*Gvariables.Pf[1]+Gvariables.Pf[2]*Gvariables.Pf[2],0.5);
		//test_complex_square_roots();
		G_eval Gclass;//(variables);
		Gclass.Set_values(Gvariables.Ei,Gvariables.Pi,
				Gvariables.Ef,Gvariables.Pf,Gvariables.L,Gvariables.mpi);
		string control="no_k";
		Gclass.test_int3d_total(&int3d_real,1./pow(3.,4.),control);
		//cout<<"Eicm="<<Eicm<<"Efcm="<<Efcm<<"; In of Fig.5 left panel now="<<int3d_real<<endl;
		cout<<"Eicm="<<Eicm<<"Efcm="<<Efcm<<"; In of Fig.5 left panel now="<<int3d_real<<endl;
}


void plot_Fig5_right_panel(){
	double G_real,G_imag,int3d_real,int3d_imag,Eicm,Efcm;
	Eicm=3.4;
	Efcm=2.5;
		Gvariables Gvariables;
		Gvariables.mpi=1;
		Gvariables.L=6.;
		double unit_mom=2*pi/Gvariables.L;
		//Gvariables.alpha=1.0;
		Gvariables.Pi[0]=0.;
		Gvariables.Pi[1]=0.;
		Gvariables.Pi[2]=unit_mom;
		Gvariables.Pf[0]=0.;
		Gvariables.Pf[1]=0.;
		Gvariables.Pf[2]=unit_mom;
		Gvariables.Ei=pow(Eicm*Eicm+Gvariables.Pi[0]*Gvariables.Pi[0]+
						Gvariables.Pi[1]*Gvariables.Pi[1]+Gvariables.Pi[2]*Gvariables.Pi[2],0.5);
	    Gvariables.Ef=pow(Efcm*Efcm+Gvariables.Pf[0]*Gvariables.Pf[0]+
				Gvariables.Pf[1]*Gvariables.Pf[1]+Gvariables.Pf[2]*Gvariables.Pf[2],0.5);
		//test_complex_square_roots();
		G_eval Gclass;//(variables);
		Gclass.Set_values(Gvariables.Ei,Gvariables.Pi,
				Gvariables.Ef,Gvariables.Pf,Gvariables.L,Gvariables.mpi);
		cout<<"Ei,Ef="<<Gvariables.Ei<<" "<<Gvariables.Ef<<endl;
		string control="one_k";
		Gclass.test_int3d_total(&int3d_real,1./pow(3.,4),control);
		//cout<<"Eicm="<<Eicm<<"Efcm="<<Efcm<<"; In of Fig.5 left panel now="<<int3d_real<<endl;
		cout<<"Eicm="<<Eicm<<"Efcm="<<Efcm<<"; In of Fig.5 right panel now 3d int="<<int3d_real<<endl;





}

void plot_Fig6_up_panel(){
	double G_real,G_imag,int3d_real,int3d_imag,Eicm,Efcm;
	Eicm=2.05;
	Efcm=2.08;
		Gvariables Gvariables;
		Gvariables.mpi=1;
		Gvariables.L=6.;
		double unit_mom=2*pi/Gvariables.L;
		//Gvariables.alpha=1.0;
		Gvariables.Pi[0]=0.;
		Gvariables.Pi[1]=0.;
		Gvariables.Pi[2]=0.;
		Gvariables.Pf[0]=0.;
		Gvariables.Pf[1]=0.;
		Gvariables.Pf[2]=unit_mom;
		Gvariables.Ei=pow(Eicm*Eicm+Gvariables.Pi[0]*Gvariables.Pi[0]+
						Gvariables.Pi[1]*Gvariables.Pi[1]+Gvariables.Pi[2]*Gvariables.Pi[2],0.5);
	    Gvariables.Ef=pow(Efcm*Efcm+Gvariables.Pf[0]*Gvariables.Pf[0]+
				Gvariables.Pf[1]*Gvariables.Pf[1]+Gvariables.Pf[2]*Gvariables.Pf[2],0.5);
		//test_complex_square_roots();
		G_eval Gclass;//(variables);
		Gclass.Set_values(Gvariables.Ei,Gvariables.Pi,
				Gvariables.Ef,Gvariables.Pf,Gvariables.L,Gvariables.mpi);
		string control="no_k";
		double int1d_real,int1d_imag;
		Gclass.test_int1d(& int1d_real,& int1d_imag, control);
		cout<<"Eicm="<<Eicm<<"Efcm="<<Efcm<<"; In of Fig.6 upper panel now="<<int1d_real<<endl;
		cout<<"Eicm="<<Eicm<<"Efcm="<<Efcm<<"; In of Fig.6 lower panel now="<<int1d_imag<<endl;



}


void plot_Fig8_right_panel(){
	double G_real,G_imag,int3d_real,int3d_imag,Eicm,Efcm;
	Eicm=2.25;
	double s;
	s=0.6;
	Efcm=s*(Eicm-2.)+2.;
	Gvariables Gvariables;
	Gvariables.mpi=1.;
	Gvariables.L=6.;
	double unit_mom=2*pi/Gvariables.L;
		//Gvariables.alpha=1.0;
	Gvariables.Pi[0]=0.;
	Gvariables.Pi[1]=0.;
	Gvariables.Pi[2]=unit_mom;
	Gvariables.Pf[0]=0.;
	Gvariables.Pf[1]=0.;
	Gvariables.Pf[2]=unit_mom;
	Gvariables.Ei=pow(Eicm*Eicm+Gvariables.Pi[0]*Gvariables.Pi[0]+
						Gvariables.Pi[1]*Gvariables.Pi[1]+Gvariables.Pi[2]*Gvariables.Pi[2],0.5);
	Gvariables.Ef=pow(Efcm*Efcm+Gvariables.Pf[0]*Gvariables.Pf[0]+
				Gvariables.Pf[1]*Gvariables.Pf[1]+Gvariables.Pf[2]*Gvariables.Pf[2],0.5);
		//test_complex_square_roots();
	G_eval Gclass;//(variables);
	Gclass.Set_values(Gvariables.Ei,Gvariables.Pi,
				Gvariables.Ef,Gvariables.Pf,Gvariables.L,Gvariables.mpi);
	string control="no_k";
	complex<double> int_1d;
	cout<<"Eicm="<<Eicm<<"Efcm="<<Efcm<<endl;
	Gclass.test_int1d_fig8( &int_1d);
	cout<<"Eicm="<<Eicm<<"Efcm="<<Efcm<<"; In of Fig.8 1d real Inu1nu2nu3="<<real(int_1d)<<endl;
	cout<<"Eicm="<<Eicm<<"Efcm="<<Efcm<<"; In of Fig.8 1d imag Inu1nu2nu3="<<imag(int_1d)<<endl;
	Gclass.test_int3d_total(&int3d_real,1./pow(3.,4),control);
	//cout<<"Eicm="<<Eicm<<"Efcm="<<Efcm<<"; In of Fig.5 left panel now="<<int3d_real<<endl;
	cout<<"Eicm="<<Eicm<<"Efcm="<<Efcm<<"; In of Fig.8 right panel now 3d int="<<int3d_real<<endl;
	cout<<"final Integral is="<<int_1d+int3d_real<<endl;
	double sumk_real,sumk_imag;
	Gclass.test_sumk(& sumk_real,& sumk_imag);
	cout<<"sumk"<<sumk_real<<endl;
	cout<<"*******full G= "<<real(int_1d)+int3d_real+sumk_real<<endl;


}

void plot_Fig8_right_panel_final(){
	double G_real,G_imag,int3d_real,int3d_imag,Eicm,Efcm;
	Eicm=2.3;
	double s=0.6;
	Efcm=s*(Eicm-2.)+2.;
	//Eicm=3.2;//2.16364;//2.34545; //2.25455;
	//Efcm=2.6;//3.70909;//2.61818;//3.52727;
	Gvariables Gvariables;
	Gvariables.mpi=1.;
	Gvariables.L=6.;
	double unit_mom=2*pi/Gvariables.L;
		//Gvariables.alpha=1.0;
	Gvariables.Pi[0]=0.;
	Gvariables.Pi[1]=0.;
	Gvariables.Pi[2]=unit_mom;
	Gvariables.Pf[0]=0.;
	Gvariables.Pf[1]=0.;
	Gvariables.Pf[2]=unit_mom;
	Gvariables.Ei=pow(Eicm*Eicm+Gvariables.Pi[0]*Gvariables.Pi[0]+
						Gvariables.Pi[1]*Gvariables.Pi[1]+Gvariables.Pi[2]*Gvariables.Pi[2],0.5);
	Gvariables.Ef=pow(Efcm*Efcm+Gvariables.Pf[0]*Gvariables.Pf[0]+
				Gvariables.Pf[1]*Gvariables.Pf[1]+Gvariables.Pf[2]*Gvariables.Pf[2],0.5);
		//test_complex_square_roots();
	G_eval Gclass;//(variables);
	Gclass.Set_values(Gvariables.Ei,Gvariables.Pi,
				Gvariables.Ef,Gvariables.Pf,Gvariables.L,Gvariables.mpi);
	cout<<endl;
	complex<double>int_1d;
		double Lambda0=1.,Lambda1=3.,Lambda2=6.;
		double c0=1.,c1=-35./27.,c2=8./27.;
		complex<double> IA0,IA1,IA2;
		Gclass.IA_nu1nu2nu3(Lambda0,&IA0);
		Gclass.IA_nu1nu2nu3(Lambda1,&IA1);
		Gclass.IA_nu1nu2nu3(Lambda2,&IA2);
		cout<<"IA0="<<IA0<<endl;
		cout<<"IA1="<<IA1<<endl;
		cout<<"IA2="<<IA2<<endl;
		int_1d=c0*IA0+c1*IA1+c2*IA2;
		cout<<"1d integral"<<int_1d<<";Eicm="<<Eicm<<"; Efcm="<<Efcm<<endl;
	Gclass.test_fig8();
	cout<<"2\pi/L="<<unit_mom<<endl;
	cout<<" vec{Pi}="<<Gvariables.Pi[0]<<" "<<Gvariables.Pi[1]<<" "<<Gvariables.Pi[2]<<endl;
	cout<<" vec{Pf}="<<Gvariables.Pf[0]<<" "<<Gvariables.Pf[1]<<" "<<Gvariables.Pf[2]<<endl;

	cout<<"Eicm="<<Eicm<<"; Efcm="<<Efcm<<"; s="<<s<<endl;




}

void plot_Fig9_right_panel_final(){
	double G_real,G_imag,int3d_real,int3d_imag,Eicm,Efcm;
	Eicm=2.1;
	double s;
	s=0.6;
	Efcm=s*(Eicm-2.)+2.;
	Gvariables Gvariables;
	Eicm=3.4;
	Efcm=2.05;
	Gvariables.mpi=1.;
	Gvariables.L=6.;
	double unit_mom=2*pi/Gvariables.L;
		//Gvariables.alpha=1.0;
	Gvariables.Pi[0]=0.;
	Gvariables.Pi[1]=0.;
	Gvariables.Pi[2]=0.;
	Gvariables.Pf[0]=0.;
	Gvariables.Pf[1]=0.;
	Gvariables.Pf[2]=unit_mom;
	Gvariables.Ei=pow(Eicm*Eicm+Gvariables.Pi[0]*Gvariables.Pi[0]+
						Gvariables.Pi[1]*Gvariables.Pi[1]+Gvariables.Pi[2]*Gvariables.Pi[2],0.5);
	Gvariables.Ef=pow(Efcm*Efcm+Gvariables.Pf[0]*Gvariables.Pf[0]+
				Gvariables.Pf[1]*Gvariables.Pf[1]+Gvariables.Pf[2]*Gvariables.Pf[2],0.5);
		//test_complex_square_roots();
	G_eval Gclass;//(variables);
	Gclass.Set_values(Gvariables.Ei,Gvariables.Pi,
				Gvariables.Ef,Gvariables.Pf,Gvariables.L,Gvariables.mpi);
	cout<<endl;
	Gclass.test_fig9();
	cout<<"2\pi/L="<<unit_mom<<endl;
		cout<<" vec{Pi}="<<Gvariables.Pi[0]<<" "<<Gvariables.Pi[1]<<" "<<Gvariables.Pi[2]<<endl;
		cout<<" vec{Pf}="<<Gvariables.Pf[0]<<" "<<Gvariables.Pf[1]<<" "<<Gvariables.Pf[2]<<endl;
	cout<<"Eicm="<<Eicm<<"; Efcm="<<Efcm<<"; s="<<s<<endl;




}

int main() {
	//test_complex_square_roots();
	/*plot_Fig5_left_panel();
	plot_Fig5_right_panel();*/
	//plot_Fig8_right_panel_final();
	//plot_Fig9_right_panel_final();

	plot_Fig8_right_panel_final();
	cout<<"----"<<endl;
	//plot_Fig6_up_panel();


	/*double G_real,G_imag;
	Gvariables Gvariables;
	Gvariables.Ei=1.8;
	Gvariables.Ef=3.999;
	Gvariables.mpi=1;
	Gvariables.L=4.;
	double unit_mom=2*pi/Gvariables.L;*/
	//Gvariables.alpha=1.0;
	/*Gvariables.Pi[0]=0.;
	Gvariables.Pi[1]=0.;
	Gvariables.Pi[2]=0.;
	Gvariables.Pf[0]=0.;
	Gvariables.Pf[1]=0.;
	Gvariables.Pf[2]=unit_mom;
	//test_complex_square_roots();
	G_eval Gclass;//(variables);
	Gclass.Set_values(Gvariables.Ei,Gvariables.Pi,
			Gvariables.Ef,Gvariables.Pf,Gvariables.L,Gvariables.mpi);
	Gclass.intGfunct(&G_real,&G_imag);
	cout<<"here"<<endl;
	cout<<"Ei="<<Gvariables.Ei<<"; Ef="<< Gvariables.Ef<<"; real G="<<G_real<<"; imag G="<<G_imag<<endl;
	double int1d_real,int1d_imag,int1d_11_real,int1d_11_imag;
	double int1d_12_real,int1d_12_imag;
	double Ei,Ef, Efhcm;
	double Efcm,Eicm,int3d_real,int3d_imag,sumk_real,sumk_imag;
	Gclass.test_int1d(&int1d_real,&int1d_imag);
	Gclass.test_int3d(&int3d_real,&int3d_imag);
	Gclass.test_sumk(&sumk_real,&sumk_imag);
	cout<<"-------------------------------------"<<endl;
	cout<<"Ei="<<Gvariables.Ei<<"; Ef="<<Gvariables.Ef<<endl;
	cout<<"int1d real="<<int1d_real<<"; int1d_imag="<< int1d_imag<<endl;
	cout<<"int3d real="<<int3d_real<<"; int3d_imag="<< int3d_imag<<endl;
	cout<<"sumk_real="<<sumk_real<<"; sumk_imag="<< sumk_imag<<endl;
	cout<<"---------------------------------------"<<endl;
	Eicm=1.9;
	Ei=pow(Eicm*Eicm+Gvariables.Pi[0]*Gvariables.Pi[0]+Gvariables.Pi[1]*Gvariables.Pi[1]+Gvariables.Pi[2]*Gvariables.Pi[2],0.5);
	cout<<"Ei-Eicm="<<Eicm-Ei<<endl;
	Efcm=1.8;
	cout<<"test"<<endl;
	ofstream outputs_3dreal_Efcm,outputs_1dreal_Efcm,outputs_1dimag_Efcm;
	outputs_3dreal_Efcm.open("int3dreal_Ecm.txt");
	outputs_1dreal_Efcm.open("int1dreal_Ecm.txt");
	outputs_1dimag_Efcm.open("int1dimag_Ecm.txt");*/
	//outputs_1dreal_Efcm.open("int1d_F5_real_Ecm.txt");
	//outputs_1dimag_Efcm.open("int1_F5_imag_Ecm.txt");
	//double F5_real,F5_imag,total_int;

		/*	for (int j=1;j<=20;j++){
				Efhcm=Efcm+j/10.;
				Ef=pow(Efhcm*Efhcm+Gvariables.Pf[0]*Gvariables.Pf[0]+Gvariables.Pf[1]*Gvariables.Pf[1]+Gvariables.Pf[2]*Gvariables.Pf[2],0.5);
				Gclass.Set_values(Ei,Gvariables.Pi,Ef,Gvariables.Pf,Gvariables.L,Gvariables.mpi);
			    cout<<"-----"<<endl;
				//Gclass.intGfunct(&G_real,&G_imag);
			    Gclass.test_int3d(&int3d_real,&int3d_imag);
			    Gclass.test_int1d(&int1d_real,&int1d_imag);
				//Gclass.test_int1d_v2(& int1d_11_real, & int1d_11_imag,
				    //     & int1d_12_real, & int1d_12_imag);
				//Gclass.test_int3d_total(&total_int);
				//Gclass.int1d_F5_test(& F5_real,& F5_imag);
				//Ef=pow(Efhcm*Efhcm+4*pi*pi/(Gvariables.L*Gvariables.L),0.5);
			    //Efcm=pow(Efh*Efh-4*pi*pi/(Gvariables.L*Gvariables.L),0.5);
				outputs_1dreal_Efcm<<Efhcm<<" "<<int1d_real<<endl;
				//outputs_1dimag_Efcm<<Efhcm<<" "<<int1d_imag<<endl;
				outputs_3dreal_Efcm<<Efhcm<<" "<<int3d_real<<endl;
				cout<<"In"<<int3d_real<<endl;

			    cout<<"Eicm="<<Eicm<<" Efcm="<<Efhcm<<"; num3d="<<int3d_real<<endl;
			}
			outputs_3dreal_Efcm.close();
			outputs_1dreal_Efcm.close();
			outputs_1dimag_Efcm.close();

	//int num_indexes=1;
	//Gclass.integral_test(num_indexes);

	//ofstream outputs_real,outputs_imag,outputs_real_Ecm,outputs_imag_Ecm,test;
	/*outputs_real.open("G_real3D_v3.txt");
	outputs_imag.open("G_imag3D_v3.txt");
	outputs_real_Ecm.open("G_real3D_Ecm_v3.txt");
	outputs_imag_Ecm.open("G_imag3D_Ecm_v3.txt");
	test.open("test.txt");
	double Ei,Ef,Eih, Efh;
	double Efcm;
	Ei=2.015;
	Ef=2.013;
	outputs_real.precision(5);
	outputs_imag.precision(5);
	//Holds only for the above specific case!!!!! with Pi=[0,0,0] and Pf=[0,0,2\pi/L]
	// outputs in two separate files for the real and imaginary part
	// To note that the class definition should be changed, in particular set values should take only
	//Pi,Pf,mpi,L while Ei anf Ef should be arguments of the method intGfunction (which contains both the sum and
	//the integral).
	for (int i=1;i<=10;i++){
		for (int j=1;j<=10;j++){
			Eih=Ei+i/50.;
			Efh=Ef+j/50.;
			Gclass.Set_values(Eih,Gvariables.Pi,Efh,Gvariables.Pf,Gvariables.L,Gvariables.mpi);
		    Gclass.intGfunct(&G_real,&G_imag);
		    Efcm=pow(Efh*Efh-4*pi*pi/(Gvariables.L*Gvariables.L),0.5);
		    test<<"Ei="<<Eih<<"Efcm="<<Efcm<<endl;
		    outputs_real<<" "<<Eih<<" "<<Efh<<" "<<G_real<<endl;
		    outputs_imag<<" "<<Eih<<" "<<Efh<<" "<<G_imag<<endl;
		    outputs_real_Ecm<<" "<<Eih<<" "<<Efcm<<" "<<G_real<<endl;
		    outputs_imag_Ecm<<" "<<Eih<<" "<<Efcm<<" "<<G_imag<<endl;
		}
	}
	//cout<<"Ei="<<Eih<<"; Ef="<< Efh<<"; real G="<<G_real<<"; imag G="<<G_imag<<endl;
	 // outputs<<" "<<Eih<<" "<< Efh<<" "<<G_real<<" "<<G_imag<<endl;
			    //outputs<<Eih<<" "<<Efh<<" "<<G_real<<" "<<G_imag<<endl;
	test.close();
	outputs_real.close();
	outputs_imag.close();
	outputs_real_Ecm.close();
	outputs_imag_Ecm.close();
	//outputs<<Gvariables.Ei<<"  "<< Gvariables.Ef<<" "<<G_real<<" "<<G_imag;*/

	return 0;
}
