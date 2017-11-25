/*********************************************************************
 Copyright 2003 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net 

 $Id: var_metrics.c,v 1.9 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

#include <math.h>
#include <stdlib.h>
#define EXTERN
#include "var_metrics.h" 

double
var_euclidean(VAR_VECTOR *v1,
	      VECTOR *v2)
{
  double distance = 0.0;
  int i;
  double *v1m = v1->mean;
  double *v1s = v1->stdev;
  double *v2p = v2->data;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  for (i = 0; i < v1->n; i++) {
    double stdv = *v1s++;
    double d = stdv > 0.0 ? ((*v1m++) - (*v2p++)) / stdv : 0.0;
    distance += d * d;
  }
  return sqrt(distance / v1->n);

}

double
var_ieuclidean(VAR_VECTOR *v1,
	       VECTOR *v2)
{
  double distance = 0.0;
  int i;
  double *v1m = v1->mean;
  double *v1s = v1->stdev;
  double *v2p = v2->data;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  for (i = 0; i < v1->n; i++) {
    double stdv = *v1s++;
    double d = stdv > 0.0 ? (((*v1m++) * -1) - (*v2p++)) / stdv : 0.0;
    distance += d * d;
  }
  return sqrt(distance / v1->n);

}

double
var_minkowski(VAR_VECTOR *v1,
	      VECTOR *v2,
	      double power)
{
  double distance = 0.0;
  int i;
  double *v1m = v1->mean;
  double *v1s = v1->stdev;
  double *v2p = v2->data;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  for (i = 0; i < v1->n; i++) {
    double stdv = *v1s++;
    double d = stdv > 0.0 ? ((*v1m++) - (*v2p++)) / stdv : 0.0;
    distance += pow(fabs(d), power);
  }
  return (pow(distance, 1/power)) / v1->n;

}
 
double
var_iminkowski(VAR_VECTOR *v1,
	       VECTOR *v2,
	       double power)
{
  double distance = 0.0;
  int i;
  double *v1m = v1->mean;
  double *v1s = v1->stdev;
  double *v2p = v2->data;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  for (i = 0; i < v1->n; i++) {
    double stdv = *v1s++;
    double d = stdv > 0.0 ? (((*v1m++) * -1) - (*v2p++)) / stdv : 0.0;
    distance += pow(fabs(d), power);
  }
  return (pow(distance, 1/power)) / v1->n;

}

double
var_manhattan(VAR_VECTOR *v1,
	      VECTOR *v2)
{
  double distance = 0.0;
  int i;
  double *v1m = v1->mean;
  double *v1s = v1->stdev;
  double *v2p = v2->data;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  for (i = 0; i < v1->n; i++) {
    double stdv = *v1s++;
    double d = stdv > 0.0 ? fabs((*v1m++) - (*v2p++)) / stdv : 0.0;
    distance += d;
  }
  return  distance / v1->n;

}

double
var_imanhattan(VAR_VECTOR *v1,
	       VECTOR *v2)
{
  double distance = 0.0;
  int i;
  double *v1m = v1->mean;
  double *v1s = v1->stdev;
  double *v2p = v2->data;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  for (i = 0; i < v1->n; i++) {
    double stdv = *v1s++;
    double d = stdv > 0.0 ? fabs(((*v1m++) * -1) - (*v2p++)) / stdv : 0.0;
    distance += d;
  }
  return  distance / v1->n;

}

double
var_maximum(VAR_VECTOR *v1,
	    VECTOR *v2)
{
  double distance = 0.0;
  double value = 0.0;
  int i;
  double *v1m = v1->mean;
  double *v1s = v1->stdev;
  double *v2p = v2->data;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  for (i = 0; i < v1->n; i++) {
    double stdv = *v1s++;
    value = stdv > 0.0 ? fabs((*v1m++) - (*v2p++)) / stdv : 0.0;
    if (value > distance) {
      distance = value;  
    }
  }
  return  distance;

}

double
var_imaximum(VAR_VECTOR *v1,
	     VECTOR *v2)
{
  double distance = 0.0;
  double value   = 0.0;
  int i;
  double *v1m = v1->mean;
  double *v1s = v1->stdev;
  double *v2p = v2->data;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  for (i = 0; i < v1->n; i++) {
    double stdv = *v1s++;
    value = stdv > 0.0 ? fabs(((*v1m++) * -1) - (*v2p++)) / stdv : 0.0;
    if (value > distance) {
      distance = value;  
    }
  }
  return  distance;

}
 
double
var_pearson(VAR_VECTOR *v1,
	    VECTOR *v2)
{

  double sumX    = 0.0;
  double sumY    = 0.0;
  double sumX2   = 0.0;
  double sumY2   = 0.0;
  double sumXY   = 0.0;
  double SXX     = 0.0;
  double SYY     = 0.0;
  double SXY     = 0.0;
  int i;
  double *v1m = v1->mean;
  /* double *v1s = v1->stdev; */
  double *v2p = v2->data;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  for (i = 0; i < v1->n; i++) {
    sumX    += *v1m;
    sumY    += *v2p;
    sumX2   += (*v1m) * (*v1m);
    sumY2   += (*v2p) * (*v2p);
    sumXY   += (*v1m++) * (*v2p++);
  }
  SXX  = sumX2 - ((sumX * sumX) / v1->n);
  SYY  = sumY2 - ((sumY * sumY) / v2->n);
  SXY  = sumXY - (sumX * sumY) / v1->n;
  
  if (SXX * SYY == 0) {
    return 0.0;
  } else {
    return (SXY / (sqrt (SXX * SYY)));
  }

}
 
double
var_ipearson(VAR_VECTOR *v1,
	     VECTOR *v2)
{

  double sumX    = 0.0;
  double sumY    = 0.0;
  double sumX2   = 0.0;
  double sumY2   = 0.0;
  double sumXY   = 0.0;
  double SXX     = 0.0;
  double SYY     = 0.0;
  double SXY     = 0.0;
  int i;
  double *v1m = v1->mean;
  /* double *v1s = v1->stdev; */
  double *v2p = v2->data;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  for (i = 0; i < v1->n; i++) {
    double invp = (*v1m++ * -1);
    sumX           += invp;
    sumY           += (*v2p);
    sumX2          += invp * invp;
    sumY2          += (*v2p) * (*v2p);
    sumXY          += invp * (*v2p++);
  }
  SXX  = sumX2 - ((sumX * sumX) / v1->n);
  SYY  = sumY2 - ((sumY * sumY) / v2->n);
  SXY  = sumXY - (sumX * sumY) / v1->n;
  
  if (SXX * SYY == 0) {
    return 0.0;
  } else {
    return (SXY / (sqrt (SXX * SYY)));
  }

}

double
var_spearman(VAR_VECTOR *v1,
	     VECTOR *v2)
{

  int i = v1->n;
  int l = 0;
  double *v1m = v1->mean;
  /* double *v1s = v1->stdev; */
  double *v2p = v2->data;
  double d, zd, probd, rs, probrs;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  
  spearman_corr(v1m, v2p, i, l, &d, &zd, &probd, &rs, &probrs);
  
  return rs;
  
}

double
var_ispearman(VAR_VECTOR *v1,
	      VECTOR *v2)
{

  int i = v1->n;
  int l = 0;
  double *v1m = v1->mean;
  /* double *v1s = v1->stdev; */
  double *v2p = v2->data;
  double *v3p;
  double d, zd, probd, rs, probrs;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  
  v3p = ARRAY_ALLOC(v1->n, double);

  for (i = 0; i < v1->n; i++) {
    v3p[i] = *v1m++ * -1;
  }
  
  spearman_corr(v3p, v2p, i, l, &d, &zd, &probd, &rs, &probrs);
  
  return rs;

}

double
var_lspearman(VAR_VECTOR *v1,
	      VECTOR *v2)
{

  int i = v1->n;
  int l = 1;
  double *v1m = v1->mean;
  /* double *v1s = v1->stdev; */
  double *v2p = v2->data;
  double d, zd, probd, rs, probrs;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  
  spearman_corr(v1m, v2p, i, l, &d, &zd, &probd, &rs, &probrs);
  
  return rs;
  
}

double
var_ilspearman(VAR_VECTOR *v1,
	       VECTOR *v2)
{

  int i = v1->n;
  int l = 1;
  double *v1m = v1->mean;
  /* double *v1s = v1->stdev; */
  double *v2p = v2->data;
  double *v3p;
  double d, zd, probd, rs, probrs;
  
  if (v1->n != v2->n) {
    fprintf(stderr, "Vector size mismatch %d != %d\n", v1->n, v2->n);
    exit(1);
  }
  
  v3p = ARRAY_ALLOC(v1->n, double);

  for (i = 0; i < v1->n; i++) {
    v3p[i] = *v1m++ * -1;
  }
  
  spearman_corr(v3p, v2p, i, l, &d, &zd, &probd, &rs, &probrs);
  
  return rs;

}
