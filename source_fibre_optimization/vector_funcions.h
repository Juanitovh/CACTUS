#ifndef VECTOR_FUNCIONS_H
#define VECTOR_FUNCIONS_H



#include<iostream>
#include <vector>
#include<math.h>
#include <ctime>
#include <string>
#include <fstream>



#include <stdio.h>
#include <stdlib.h>

using namespace std;

#define vf vector<double>
#define vvf vector<vector<double>>

void parallel_fill_n(vf& v, size_t count, int value);
void parallel_fill_n(double * v, size_t count, int value);

float float_rand( float min, float max );


inline double norm(double *x , int n)
{
    return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}



inline double norm(double * x , double * x2 , int n)
{
    return sqrt((x[0]-x2[0])*(x[0]-x2[0]) + (x[1]-x2[1])*(x[1]-x2[1]) + (x[2]-x2[2])*(x[2]-x2[2]));
}



void subarray(double * v, int a , int b , double * out);


void subarray_point(double * v , int a, int b  , double ** out	);
void copy_pointer(double * v , int a, double ** out);

double maximum(double a, double b);



void print_v(double * v , int n);


double* linspace(double start , double end , int n);


double ** linspace_nd(double * start , double *end , int n , int dimension);


double * flatten(double **  v, int axis , int m , int n);
double * flatten(vector<vector<double>>  v, int axis , int m , int n);

double ** create_matrix( int num_rows , int num_cols);


void clear_intelligent(double ** v);


void clear_arr(double * v);

void clear_arr(double ** v , int n);


void print_params(double * v , int n);


void print_vnd(double **v , int n , int m);



double dot_prod(double *v1 , double *v2 , int n);
// inline double dot_prod3d(double *v1, double *v2);

inline double dot_prod3d(double *v1, double *v2) {
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}



void suma_vec(double * v1 ,double *v2 , double alfa, double *res , int n);


void resta_alfav(double * v1 ,vector<double> v2, double alfa, int n);
void resta_alfav(double * v1 ,double *v2 , double alfa, int n);

vf operator*( double alfa , const vf& v);
vf operator*(const vf& v, const vf& v2);
vf sqrt_v(vf x);
vf operator+(const vf& v, const vf& v2);
vf operator-(const vf& v, const vf& v2);
vf operator/(const vf& v, const vf& v2);
double norm(vf x , vf x2);
double norm(vf x);
#endif
