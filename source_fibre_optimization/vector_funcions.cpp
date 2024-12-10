#include<iostream>
#include <vector>
#include<math.h>
#include <ctime>
#include <string>
#include <fstream>

using namespace std;

#include <omp.h>
#include "omp_config.h"


#include <stdio.h>
#include <stdlib.h>

#include "vector_funcions.h"


#define vf vector<double>
#define vvf vector<vector<double>>

int nsize = 10000;

void parallel_fill_n(vf &v, size_t count, int value) {
    #pragma omp parallel for num_threads(GLOBAL_NUM_THREADS_vector) schedule(static, nsize)
    for (size_t i = 0; i < count; ++i) {
        v[i] = value;
    }
}


void parallel_fill_n(double * v, size_t count, int value) {
    #pragma omp parallel for num_threads(GLOBAL_NUM_THREADS_vector) schedule(static,nsize)
    for (size_t i = 0; i < count; ++i) { 
        v[i] = value;
    }
}

double double_rand( double min, double max )
{
    double scale = rand() / (double) RAND_MAX; /* [0, 1.0] */
    return min + scale * ( max - min );      /* [min, max] */
}


// double norm(double *x , int n)
// {
//     double res = 0;
//     for(int i = 0; i < n ; i++)
//         res +=(x[i]*x[i]);
//     return sqrt(res);
// }





// double norm(double * x , double * x2 , int n)
// {
//     double res = 0;
//     //cout<<"entering hell"<<endl;
//     for(int i = 0; i < n ; i++)
//     {
//         //	cout<<x[i] <<" - "<<x2[i]<<endl;
//         res +=(x[i]-x2[i])*(x[i]-x2[i]);
//     }
//     return sqrt(res);
// }



void subarray(double * v, int a , int b , double * out)
{
    for(int i = a; i < b; i++)
        out[i-a] = v[i];
}

void subarray_point(double * v , int a, int b  , double ** out	)
{
    *out = v+a;
}

void copy_pointer(double * v , int a, double ** out	)
{
    //integer "a" is the displacement between the pointers
    *out = v+a;
}

////////////////////////////////

double maximum(double a, double b)
{
    if (a>b)
        return a;
    return b;
}




void print_v(double * v , int n)
{
    for(int i=0 ; i < n ; i++)
        cout<< v[i] <<" ";
    cout<<endl;
}







double* linspace(double start , double end , int n)
{
    double* pasos = new double[n];
    double step = (end- start )/(n-1);

    for(int i = 0 ; i < n ; i++)
        pasos[i] = start + step*i;
    return pasos;
}




double ** linspace_nd(double * start , double *end , int n , int dimension)
{

    double ** tot = new double*[dimension] ;
    for(int i = 0 ; i < dimension ; i ++)
    {
        tot[i] = linspace(start[i] ,end[i] , n);
    }
    return  tot;
}


double * flatten(double **  v, int axis , int m , int n)
    //with 0 == concatenate
{
    double * new_n  = new double [n*m];
    if (axis == 0)
    {
        for (int i = 0; i < m ; i++)
            for(int j =0 ; j <n ; j++)
                new_n[i*n+j]= v[i][j];

    }
    else if (axis == 1)
    {
        for(int j =0 ; j <n ; j++)
            for (int i = 0; i < m ; i++)
                new_n[j*  + i]= v[i][j];

    }
    return new_n;
}


double * flatten(vector<vector<double>>  v, int axis , int m , int n)
    //with 0 == concatenate
{
    double * new_n  = new double [n*m];
    if (axis == 0)
    {
        for (int i = 0; i < m ; i++)
            for(int j =0 ; j <n ; j++)
                new_n[i*n+j]= v[i][j];

    }
    else if (axis = 1)
    {
        for(int j =0 ; j <n ; j++)
            for (int i = 0; i < m ; i++)
                new_n[j*m + i]= v[i][j];

    }
    return new_n;
}


///////////////////

double ** create_matrix( int num_rows , int num_cols)
{
    double** matrix;

    matrix = new double*[num_rows];
    matrix[0] = new double[num_rows * num_cols];

    for(int i = 1; i < num_rows; i++) {
        matrix[i] = matrix[i-1] + num_cols;
    }
    return matrix;
}

void clear_intelligent(double ** v)
{
    delete [] v[0];
    delete [] v; 
}


void clear_arr(double * v)
{
    delete [] v;
}

void clear_arr(double ** v , int n)
{
    for(int i = 0 ; i < n ; i++)
        delete [] v[i];
    delete [] v;
}

void print_params(double * v , int n)
{
    for(int i = 0 ; i <n ; i++)
    {
        if(i%4==0)
            cout<<" [ ";
        cout <<v[i] <<" ";
        if(i%4==3)
        {
            cout<<" ] ";
        }
    }
}

void print_vnd(double **v , int n , int m)
{
    for(int i =0 ; i< n ; i++)
    {
        for(int j =0 ; j < m ; j++)
            cout<< v[i][j] <<" ";
        cout<<endl;
    }
}


// double dot_prod3d(double *v1, double *v2)
// {
//     return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]; 
// }

// inline double dot_prod3d(const double* v1, const double* v2) {
//     return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]; 
// }


double dot_prod(double *v1 , double *v2 , int n)
{
    double res = 0; 
    for(int i = 0 ; i < n ; i++)
        res +=v1[i]*v2[i];
    return res;
}

void suma_vec(double * v1 ,double *v2 , double alfa, double *res , int n)
{

#pragma omp parallel for num_threads(GLOBAL_NUM_THREADS_vector) schedule(static,nsize)
    for(int i = 0 ; i < n ; i++)
    {
        //cout<<v1[i]<<" -*- "<<alfa*v2[i] <<" , ";
        res[i] =v1[i] +  alfa*v2[i];
    }
    //cout<<endl<<endl;

}


void resta_alfav(double * v1 ,vector<double> v2, double alfa, int n)

{
#pragma omp parallel for num_threads(GLOBAL_NUM_THREADS_vector) schedule(static,nsize)
    for(int i = 0 ; i < n ; i++)
    {
        //cout<<v1[i]<<" -*- "<<alfa*v2[i] <<" , ";
        v1[i] -= alfa*v2[i];
    }
    //cout<<endl<<endl;

}

void resta_alfav(double * v1 ,double *v2 , double alfa, int n)

{
#pragma omp parallel for num_threads(GLOBAL_NUM_THREADS_vector) schedule(static,nsize)
    for(int i = 0 ; i < n ; i++)
    {
        //cout<<v1[i]<<" -*- "<<alfa*v2[i] <<" , ";
        v1[i] -= alfa*v2[i];
    }
    //cout<<endl<<endl;

}


////////////////







vf operator*( double alfa , const vf& v)
{
    vector<double> v2 (v.size());

#pragma omp parallel for num_threads(GLOBAL_NUM_THREADS_vector) schedule(static,nsize)
    for(int i = 0 ; i < v.size() ; i++)
        v2[i] = alfa*v[i];

    return v2;
}

vf operator*(const vf& v, const vf& v2)
{
    vf v3 (v.size());

    // if( v.size() != v2.size())[[unlikely]]
    // {
    //     cout<< " Vector dimensions don't match"<<endl;
    // }

#pragma omp parallel for num_threads(GLOBAL_NUM_THREADS_vector) schedule(static,nsize)
    for(int i = 0 ; i < v.size() ; i++)
        v3[i] = v2[i]*v[i];

    return v3;
}

vf sqrt_v(vf x)
{
    vf res ((int) x.size());

#pragma omp parallel for num_threads(GLOBAL_NUM_THREADS_vector) schedule(static,nsize)
    for(int i = 0 ; i < x.size() ; i++)
    {
        res[i] = sqrt(x[i]);
    }
    return res;
}

vf operator+(const vf& v, const vf& v2)
{
    vf v3 (int(v.size()));

    // if( v.size() != v2.size())
    // {
    //     cout<< " Vector dimensions don't match"<<endl;
    // }

#pragma omp parallel for num_threads(GLOBAL_NUM_THREADS_vector) schedule(static,nsize)
    for(int i = 0 ; i < v.size() ; i++)
        v3[i] = v2[i]+v[i];

    return v3;
}

vf operator-(const vf& v, const vf& v2)
{
    vf v3 (int(v.size()));

    // if( v.size() != v2.size())
    // {
    //     cout<< " Vector dimensions don't match"<<endl;
    // }

#pragma omp parallel for num_threads(GLOBAL_NUM_THREADS_vector) schedule(static,nsize)
    for(int i = 0 ; i < v.size() ; i++)
        v3[i] = v[i]- v2[i];

    return v3;
}

vf operator/(const vf& v, const vf& v2)
{
    vf v3 (int(v.size()));

    // if( v.size() != v2.size())
    // {
    //     cout<< " Vector dimensions don't match"<<endl;
    // }

#pragma omp parallel for num_threads(GLOBAL_NUM_THREADS_vector) schedule(static,nsize)
    for(int i = 0 ; i < v.size() ; i++)
        v3[i] = v[i]/v2[i];

    return v3;
}

double norm(vf x)
{
    double res = 0;
    for(int i = 0; i < x.size() ; i++)
        res +=(x[i]*x[i]);
    return sqrt(res);
}

double norm(vf x , vf x2)
{
    double res = 0;

#pragma omp parallel for num_threads(GLOBAL_NUM_THREADS_vector) schedule(static,nsize)
    for(int i = 0; i < x.size() ; i++)
        res +=(x[i]-x2[i])*(x[i]-x2[i]);
    return sqrt(res);
}
