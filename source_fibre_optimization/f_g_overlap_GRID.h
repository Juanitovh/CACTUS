/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   f_g_overlap_GRID.h
 * Author: jppl
 *
 * Created on February 5, 2021, 7:49 PM
 */

#ifndef F_G_OVERLAP_GRID_H
#define F_G_OVERLAP_GRID_H


#include "Grid_space.h"

#include <unordered_map>
#define pdd pair<double, double>
#define pii pair<int,int>
// #define map_ii_vd map<pii , vector<double>>
#define map_ii_vd map<pii , array<double, 8>>
// #define map_ii_vd std::unordered_map<pii, std::array<double, 8>, std::hash<pii>>




// #define map_ii_vd unordered_map<pii, vector<double>>

// #define map_ii_vd std::unordered_map<pii, std::vector<double>, std::hash<pii>>


extern double scale_overlap_gs;

// double clap(double a);

// constexpr inline double clap(double a);
void  min_distance_lines(double * __restrict__ p0 , double * __restrict__ p1 , double * __restrict__ q0  , double * __restrict__ q1 , double res[2]);
double f_overlap_two_points_cuadratic(double * control0 , double * control1 );
double f_overlap_two_capsules_cuadratic(double * p0 , double * p1  , double * q0 , double *q1);
double f_overlap_all_strands_grid(v4p grid , vector<Strand> strands  );


void grad_ovelap_two_capsules_cuadratic(double * __restrict__ p0 , double * __restrict__ p1 , double * __restrict__ q0 , double * __restrict__ q1 , double * __restrict__ gp0 , double * __restrict__ gp1 , double * __restrict__  gq0 , double * __restrict__ gq1 , bool ignore_p0 , bool ignore_q0 ,int size_cube, bool compress);
void grad_ovelap_two_points_cuadratic(double * xi , double * yj , double * gxi , double * gyj ,  int size_cube );

void g_overlap_all_strands_grid(v4p grid , vector<Strand> strands , double * grad_useless  , int n_total_params , bool compress);

#endif /* F_G_OVERLAP_GRID_H */
