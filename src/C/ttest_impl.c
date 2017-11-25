/*********************************************************************
 Copyright 2003 Bristol-Myers Squibb. All rights reserved.

	 Author: Isaac M. Neuhaus
	 Phone: (609) 818 3196
	 e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

	 Several subrutines were adapted from Numerical Recipes in C, W H
	 Press et al., 1988

 $Id: ttest.c,v 1.9 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

/* in the precalc versions you can leave datas as NULL but check twice */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#define EXTERN
#include "stats.h"
#include "ttest_impl.h"
#ifndef SQR
#define	SQR(A) ((A)*(A))
#endif
//#define	TTEST_DEBUG
#undef	TTEST_DEBUG

/*
Unequal sample sizes, equal variance
This test is used only when it can be assumed that the two distributions have the same variance.
*/	
void
ttest(double data1[], int n1, double data2[], int n2, double *t, double *df,
		double *prob, double *m1, double *v1, double *m2, double *v2)
{

	double var1,var2,svar,ave1,ave2;
	
	avevar(data1,n1,&ave1,&var1);
	avevar(data2,n2,&ave2,&var2);
#ifdef TTEST_DEBUG
		printf ("Data1 (%i) => ", n1);
		print_array_d(data1, n1);
		printf ("\n");
		printf ("Data2 (%i) => ", n2);
		print_array_d(data2, n2);
		printf ("\n");
#endif
	if( m1 != NULL ) *m1 = ave1;
	if( v1 != NULL ) *v1 = var1;
	if( m2 != NULL ) *m2 = ave2;
	if( v2 != NULL ) *v2 = var2;
	*df=n1+n2-2;
	svar=((n1-1)*var1+(n2-1)*var2)/(*df);
	*t=(ave1-ave2)/sqrt(svar*(1.0/n1+1.0/n2));
	*prob=betai(0.5*(*df),0.5,(*df)/((*df)+(*t)*(*t)));
}
/* any of m1, v1, m2, v2 may contain a pre-calculated statistic which we can reuse
   if have_ flag is TRUE (1)
*/
void
ttest_precalc(double data1[], int n1, double data2[], int n2, double *t, double *df,
		double *prob,
	double *m1,
	double *v1, char have_precalc_1,
	double *m2,
	double *v2, char have_precalc_2
){
	double var1,var2,svar,ave1,ave2;
	if( have_precalc_1 == 1 ){
		ave1 = *m1; var1 = *v1;
	} else {
		avevar(data1,n1,&ave1,&var1);
		if( m1 != NULL ) *m1 = ave1;
		if( v1 != NULL ) *v1 = var1;
	}
	if( have_precalc_2 == 1 ){
		ave2 = *m2; var2 = *v2;
	} else {
		avevar(data2,n2,&ave2,&var2);
		if( m2 != NULL ) *m2 = ave2;
		if( v2 != NULL ) *v2 = var2;
	}
#ifdef TTEST_DEBUG
		printf ("Data1 (%i) => ", n1);
		print_array_d(data1, n1);
		printf ("\n");
		printf ("Data2 (%i) => ", n2);
		print_array_d(data2, n2);
		printf ("\n");
#endif
	*df=n1+n2-2;
	svar=((n1-1)*var1+(n2-1)*var2)/(*df);
	*t=(ave1-ave2)/sqrt(svar*(1.0/n1+1.0/n2));
	*prob=betai(0.5*(*df),0.5,(*df)/((*df)+(*t)*(*t)));
}


/*
Unequal sample sizes, unequal variance
This test, also known as Welch's t-test, is used only when the two population variances are assumed to be different
(the two sample sizes may or may not be equal) and hence must be estimated separately.
*/
void
tutest(double data1[], int n1, double data2[], int n2, double *t, double *df,
		 double *prob, double *m1, double *v1, double *m2, double *v2)
{

	double var1,var2,ave1,ave2;
	
	avevar(data1,n1,&ave1,&var1);
	avevar(data2,n2,&ave2,&var2);
#ifdef TTEST_DEBUG
		printf ("Data1 (%i) => ", n1);
		print_array_d(data1, n1);
		printf ("\n");
		printf ("Data2 (%i) => ", n2);
		print_array_d(data2, n2);
		printf ("\n");
#endif
	if( m1 != NULL ) *m1 = ave1;
	if( v1 != NULL ) *v1 = var1;
	if( m2 != NULL ) *m2 = ave2;
	if( v2 != NULL ) *v2 = var2;
	*t=(ave1-ave2)/sqrt(var1/n1+var2/n2);
	*df=SQR(var1/n1+var2/n2)/(SQR(var1/n1)/(n1-1)+SQR(var2/n2)/(n2-1));
	*prob=betai(0.5*(*df),0.5,(*df)/((*df)+SQR(*t)));

}
void
tutest_precalc(double data1[], int n1, double data2[], int n2, double *t, double *df,
		 double *prob,
	double *m1,
	double *v1, char have_precalc_1,
	double *m2,
	double *v2, char have_precalc_2
){
	double var1,var2,svar,ave1,ave2;
	if( have_precalc_1 == 1 ){
		ave1 = *m1; var1 = *v1;
	} else {
		avevar(data1,n1,&ave1,&var1);
		if( m1 != NULL ) *m1 = ave1;
		if( v1 != NULL ) *v1 = var1;
	}
	if( have_precalc_2 == 1 ){
		ave2 = *m2; var2 = *v2;
	} else {
		avevar(data2,n2,&ave2,&var2);
		if( m2 != NULL ) *m2 = ave2;
		if( v2 != NULL ) *v2 = var2;
	}
#ifdef TTEST_DEBUG
		printf ("Data1 (%i) => ", n1);
		print_array_d(data1, n1);
		printf ("\n");
		printf ("Data2 (%i) => ", n2);
		print_array_d(data2, n2);
		printf ("\n");
#endif
	*t=(ave1-ave2)/sqrt(var1/n1+var2/n2);
	*df=SQR(var1/n1+var2/n2)/(SQR(var1/n1)/(n1-1)+SQR(var2/n2)/(n2-1));
	*prob=betai(0.5*(*df),0.5,(*df)/((*df)+SQR(*t)));
}

/* Paired samples t-tests typically consist of a sample of matched pairs of similar units, or one group of units that has been
	 tested twice (a "repeated measures" t-test). A typical example of the repeated measures t-test would be where subjects are
	 tested prior to a treatment, say for high blood pressure, and the same subjects are tested again after treatment with a
	 blood-pressure lowering medication.
*/
void
tptest(double data1[], double data2[], int n, double *t, double *df,
		 double *prob, double *m1, double *v1, double *m2, double *v2,
		double *covaria
){

	int j;
	double var1,var2,ave1,ave2,sd,cov=0.0;
	
	avevar(data1,n,&ave1,&var1);
	avevar(data2,n,&ave2,&var2);
#ifdef TTEST_DEBUG
		printf ("Data1 (%i) => ", n);
		print_array_d(data1, n);
		printf ("\n");
		printf ("Data2 (%i) => ", n);
		print_array_d(data2, n);
		printf ("\n");
#endif
	if( m1 != NULL ) *m1 = ave1;
	if( v1 != NULL ) *v1 = var1;
	if( m2 != NULL ) *m2 = ave2;
	if( v2 != NULL ) *v2 = var2;
	for (j=0;j<n;j++)
	cov += (data1[j]-ave1)*(data2[j]-ave2);
	cov /= *df=n-1;
	if( covaria != NULL ) *covaria = cov;
	sd=sqrt((var1+var2-2.0*cov)/n);
	*t=(ave1-ave2)/sd;
	*prob=betai(0.5*(*df),0.5,(*df)/((*df)+(*t)*(*t)));
}
void
tptest_precalc(double data1[], double data2[], int n, double *t, double *df,
		 double *prob,
	double *m1,
	double *v1, char have_precalc_1,
	double *m2,
	double *v2, char have_precalc_2,
	double *covaria, char have_precalc_cov
){
	double var1,var2,svar,ave1,ave2;
	if( have_precalc_1 == 1 ){
		ave1 = *m1; var1 = *v1;
	} else {
		avevar(data1,n,&ave1,&var1);
		if( m1 != NULL ) *m1 = ave1;
		if( v1 != NULL ) *v1 = var1;
	}
	if( have_precalc_2 == 1 ){
		ave2 = *m2; var2 = *v2;
	} else {
		avevar(data2,n,&ave2,&var2);
		if( m2 != NULL ) *m2 = ave2;
		if( v2 != NULL ) *v2 = var2;
	}

	int	j;
	double sd,cov=0.0;
#ifdef TTEST_DEBUG
		printf ("Data1 (%i) => ", n);
		print_array_d(data1, n);
		printf ("\n");
		printf ("Data2 (%i) => ", n);
		print_array_d(data2, n);
		printf ("\n");
#endif
	if( have_precalc_cov ){
		cov = *covaria;
	} else {
		for (j=0;j<n;j++)
		cov += (data1[j]-ave1)*(data2[j]-ave2);
		cov /= *df=n-1;
		if( covaria != NULL ) *covaria = cov;
	}
	sd=sqrt((var1+var2-2.0*cov)/n);
	*t=(ave1-ave2)/sd;
	*prob=betai(0.5*(*df),0.5,(*df)/((*df)+(*t)*(*t)));
}



