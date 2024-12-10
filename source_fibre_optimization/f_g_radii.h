/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   f_g_radii.h
 * Author: jppl
 *
 * Created on February 11, 2021, 5:46 PM
 */

#ifndef F_G_RADII_H
#define F_G_RADII_H
#include"strand_class.h"
double f_cost_rads(vector<Strand> strands);
double f_cost_rads_half(vector<Strand> strands);


void grad_cost_rads(vector<Strand> strands);
void grad_cost_rads_half(vector<Strand> strands);


double  F_CostFunc_Rads_2CntrlPoints(double * control0 , double * control1 );


void G_CostFunc_Rads_2CntrlPoints( double* xi , double * yj , double * gxi , double * gyj , int flag_ignore_right );


#endif /* F_G_RADII_H */
