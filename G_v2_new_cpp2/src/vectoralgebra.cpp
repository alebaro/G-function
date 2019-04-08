/*
 * vectoralgebra.cpp
 *
 *  Created on: Dec 4, 2018
 *      Author: alessandrobaroni
 */

#include "vectoralgebra.h"
#include<numeric>
#include<math.h>
#include<cmath>
#include<iostream>
#include<stdlib.h>
using namespace std;

double dot_3d_product(double*a, double*b){
	double final=a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	return final;
}

double dot_product(double *a,double* b,int length){
//	cout<<begin(a);
	return inner_product(a,a+length,b,0.0);
}


void matrix_prod_vec( int length,double (*m)[4],double* p, double* final){
 double sum;
 for(int i=0;i<4;i++){
	 sum=0;
	 for(int j=0;j<4;j++){
		 sum+=m[i][j]*p[j];
	 }
	 final[i]=sum;
 }
}

void scal_prod_vec(int size,double a,double*init,double* final){
	for (int i=0;i<size;i++){
		final[i]=init[i]*a;
	}
}


