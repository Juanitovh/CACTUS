/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   f_g_2CtrlPoints.h
 * Author: jppl
 *
 * Created on February 11, 2021, 2:40 PM
 */

#ifndef F_G_2CTRLPOINTS_H
#define F_G_2CTRLPOINTS_H
#include "strand_class.h"

double  F_CostFunc_Lenght_2CntrlPoints(double * control0 , double * control1 , double L );
double F_CostFunc_Lenght_2CntrlPoints_hooke(double * control0 , double * control1  , double L);

double F_OneStrand_2CtrlPoints(Strand strand ,   double (*cost_function)(double* , double* , double)  );
double F_AllStrands_2CtrlPoints( vector<Strand> strands  ,double (*cost_function)(double *, double* , double) );

void G_CostFunc_Lenght_2CntrlPoints_hooke( double* xi , double * yj , double * gxi , double * gyj , int flag_ignore_right , double L_len );
void G_CostFunc_Lenght_2CntrlPoints( double* xi , double * yj , double * gxi , double * gyj , int flag_ignore_right , double L );
void G_OneStrand_2CtrlPoints( Strand strand ,  void (*gradient_cost_function ) (double*  , double *  , double *  , double *  , int  , double) );
void G_AllStrands_2CtrlPoints(vector<Strand>  strands, void (*gradient_cost_function ) (double*  , double *  , double *  , double *  , int , double ) );


#endif /* F_G_2CTRLPOINTS_H */
