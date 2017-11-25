/*********************************************************************
 Copyright 2006 Bristol-Myers Squibb. All rights reserved.

   Author: Isaac M. Neuhaus
   Phone: (609) 818 3196
   e-mail: isaac.neuhaus@bms.com or neuhausi@comcast.net

   Several subrutines were adapted from Numerical Recipes in C, W H
   Press et al., 1988

 $Id: ttest.h,v 1.3 2007/02/28 21:58:29 neuhausi Exp $
**********************************************************************/

#ifndef ttest_impl_h
#define ttest_impl_h 

/* in the precalc versions you can leave datas as NULL but check twice */

/*
Unequal sample sizes, equal variance
This test is used only when it can be assumed that the two distributions have the same variance.
*/
void
ttest(double data1[], int n1, double data2[], int n2,
      double *t, double *df, double *prob,
      double *m1, double *v1, double *m2, double *v2);

/*
Unequal sample sizes, equal variance
This test is used only when it can be assumed that the two distributions have the same variance.
BUT
 any of m1, v1, m2, v2 may contain a pre-calculated statistic which we can reuse
   if have_ flag is TRUE (1)
*/
void
ttest_precalc(double data1[], int n1, double data2[], int n2, double *t, double *df,
                double *prob,
        double *m1,
        double *v1, char have_precalc_1,
        double *m2,
        double *v2, char have_precalc_2
);
/*
Unequal sample sizes, unequal variance
This test, also known as Welch's t-test, is used only when the two population variances are assumed to be different
(the two sample sizes may or may not be equal) and hence must be estimated separately.
*/
void
tutest(double data1[], int n1, double data2[], int n2,
       double *t, double *df, double *prob,
       double *m1, double *v1, double *m2, double *v2);
void
tutest_precalc(double data1[], int n1, double data2[], int n2,
       double *t, double *df, double *prob,
        double *m1,
        double *v1, char have_precalc_1,
        double *m2,
        double *v2, char have_precalc_2);
/* Paired samples t-tests typically consist of a sample of matched pairs of similar units, or one group of units that has been
   tested twice (a "repeated measures" t-test). A typical example of the repeated measures t-test would be where subjects are
   tested prior to a treatment, say for high blood pressure, and the same subjects are tested again after treatment with a
   blood-pressure lowering medication.
*/
void
tptest(double data1[], double data2[], int n,
       double *t, double *df, double *prob,
       double *m1, double *v1, double *m2, double *v2, double *covaria);
void
tptest_precalc(double data1[], double data2[], int n,
       double *t, double *df, double *prob,
        double *m1,
        double *v1, char have_precalc_1,
        double *m2,
        double *v2, char have_precalc_2,
        double *covaria, char have_precalc_covaria);
#endif
