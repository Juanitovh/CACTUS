/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   f_g_3CtrlPoints.h
 * Author: jppl
 *
 * Created on February 11, 2021, 4:00 PM
 */

#ifndef F_G_3CTRLPOINTS_H
#define F_G_3CTRLPOINTS_H
#include "strand_class.h"
using namespace std;

double F_CostFunc_Curve_3CntrlPoints(double * x , double * y , double *z); 
double F_OneStrand_3CtrlPoints(Strand strand,  double (*cost_function)(double* , double* , double *) );
double F_AllStrands_3CtrlPoints(vector<Strand> strands , double (*cost_function)(double* , double* , double *) );

void G_CostFunc_Curve_3CntrlPoints(double * x, double * y  , double *z, double * gx , double * gy , double *gz, int flag_ignore );
void G_OneStrand_3CtrlPoints(Strand strand  ,void (*grad_cost_function) (double * , double *   , double *, double *  , double *  , double *, int  ));
void G_AllStrands_3CtrlPoints(vector<Strand> strands ,void (*grad_cost_function) (double * , double *   , double *, double *  , double *  , double *, int  ));

#endif /* F_G_3CTRLPOINTS_H */
