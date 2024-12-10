/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   f_g_direction_strand.h
 * Author: jppl
 *
 * Created on April 24, 2021, 5:59 PM
 */

#ifndef F_G_DIRECTION_STRAND_H
#define F_G_DIRECTION_STRAND_H

#include "vector_funcions.h"
#include "strand_class.h"


double F_cost_angle_strand(double * x , double * y , double * direction_strand);
double F_OneStrand_2CtrlPoints(Strand strand,  double (*cost_function)(double* , double* , double *) );
double F_AllStrands_2CtrlPoints( vector<Strand> strands  ,double (*cost_function)(double *, double* , double*) );

void G_CostFunc_angle_strand(double * x, double * y  ,  double * gx , double * gy , double *direction_strand,  int flag_ignore );
void G_OneStrand_3CtrlPoints(Strand strand  ,void (*grad_cost_function) (double * , double *   ,  double *  , double *  , double *,  int  ));
void G_AllStrands_3CtrlPoints(vector<Strand> strands ,void (*grad_cost_function) (double * , double *   ,  double *  , double *  , double *,  int ));


#endif /* F_G_DIRECTION_STRAND_H */
