/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   f_g_compression.h
 * Author: jppl
 *
 * Created on February 12, 2021, 11:36 AM
 */

#ifndef F_G_COMPRESSION_H
#define F_G_COMPRESSION_H
#include "strand_class.h"



double f_compress_point(double * point , double *mass_gravity , double **bounding_box);
double f_compress_box(double * point , double *mass_gravity , double **bounding_box);
double F_CostFunc_compression(vector<Strand> strands , double *mass_gravity, double **bounding_box ,double (*cost_function)(double* , double* , double **));

void  g_compress_point(double * point , double *g_point  ,  double *mass_gravity , double **bounding_box);
void  g_compress_box(double * point , double *g_point  ,  double *mass_gravity , double **bounding_box);
void g_CostFunc_compression(vector<Strand> strands , double *mass_gravity , double **bounding_box, void (*cost_function)(double * , double* , double* , double **));
#endif /* F_G_COMPRESSION_H */
