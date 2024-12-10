/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   f_g_3CtrlPoints.cpp
 * Author: jppl
 * 
 * Created on February 11, 2021, 4:00 PM
 */

#include "f_g_3CtrlPoints.h"

#include"vector_funcions.h"

double F_CostFunc_Curve_3CntrlPoints(double * x , double * y , double *z)
{
    double eps = 1e-16;
    double r1 = x[3];
    double r2 = y[3];
    double r3 = z[3];

    //double zeros[3] = {0,0,0};
    double nv0 = norm(y,x,3);
    double nv1 = norm(y,z,3);

    double dot_v0v1 = 0;
    for(int i =0 ; i < 3 ; i++)
        dot_v0v1 += (z[i] - y[i])*(y[i] - x[i]);

    double base = (1- dot_v0v1/(nv0*nv1 + eps));

    double pre = 3.0/ (r1+r2+r3);

    return base*base*pre;
}


double F_OneStrand_3CtrlPoints(Strand strand,  double (*cost_function)(double* , double* , double *) )
{
    int n_control = strand.n_control;
    double  *control0 , *control1 ,*control2;
    int pos_start;


    //cout<<"dist are cool"<<endl;

    double costo = 0;
    for(int i = 1; i < n_control-1  ; i++)
    {
        pos_start = 4*i;
        copy_pointer(strand.params , pos_start-4   ,&control0);
        copy_pointer(strand.params , pos_start        ,&control1);
        copy_pointer(strand.params , pos_start+4  ,&control2);

        costo += cost_function(control0 , control1, control2);
    }
    //copy_pointer(strand.params , 0  ,&control0);
    //copy_pointer(strand.params , 4 ,&control1);
    //costo += cost_function(strand.start , control0 , control1);


    //copy_pointer(strand.params , 4*n_control-8 ,&control0);
    //copy_pointer(strand.params , 4*n_control-4 ,&control1);
    //costo += cost_function(control0 , control1 , strand.end);
    //cout <<" costo 1 strand" << costo<<endl;
    return costo;
}

double F_AllStrands_3CtrlPoints(vector<Strand> strands , double (*cost_function)(double* , double* , double *) )
{

    double res = 0;
    for(int i = 0; i < strands.size() ; i++)
    {

        if (strands[i].n_control != 1)
            res += F_OneStrand_3CtrlPoints(strands[i]   ,cost_function);

    }
    return res;
}

////////////////

void G_CostFunc_Curve_3CntrlPoints(double * x, double * y  , double *z, double * gx , double * gy , double *gz, int flag_ignore )
{
    double eps =1e-16;
    double r1 = x[3];
    double r2 = y[3];
    double r3 = z[3];

    //double zeros[3] = {0,0,0};
    double nv0 = norm(y,x,3);
    double nv1 = norm(y,z,3);

    double dot_v0v1 = 0;
    for(int i =0 ; i < 3 ; i++)
        dot_v0v1 += (z[i] - y[i])*(y[i] - x[i]);

    double base = 2*( dot_v0v1/(nv0*nv1 + eps) -1 );

    double pre =  3.0/(r1+r2+r3);



    //double grad_r = .0*   2*pre/3.0*base*base;
    //cout<<"second part"<<endl;

    //gy[3] += grad_r;
    //cost functions triple_checked_ apparently in version with overleaf
    for(int i = 0 ; i < 3 ; i++)
        gy[i]  += pre*base* (  nv1*nv0*(z[i]+x[i] -2*y[i])   -   dot_v0v1 *( nv1*(y[i]-x[i])/(nv0+ eps) + nv0*(y[i]-z[i])/(nv1 + eps)  )      )/ (nv0*nv0*nv1*nv1 + eps);
    //gy[i]  += pre*base* (  +nv1*nv0*(z[i]+x[i] -2*y[i])   -   dot_v0v1 *    ( nv1*(y[i]-x[i])/(nv0+ eps) + nv0*(y[i]-z[i])/(nv1 + eps)  )      )/ (nv0*nv0*nv1*nv1 + eps);

    if(flag_ignore != -1) [[likely]]
    {
        //do x
        //	gx[3] += grad_r;
        for(int i = 0 ; i < 3 ; i++)
            gx[i]  += pre*base* (  nv1*nv0*(-z[i] + y[i]) - dot_v0v1*(x[i] - y[i])*nv1/(nv0+eps)   )/ (nv0*nv0*nv1*nv1 + eps);

    }	

    if (flag_ignore !=1) [[likely]]
    {
        // do z
        //gz[3] += grad_r;
        for(int i = 0 ; i < 3 ; i++)
            gz[i]  += pre*base* (  nv1*nv0*(y[i]-x[i]) -  dot_v0v1 *(z[i] - y[i])*nv0/(nv1+eps)   )/ (nv0*nv0*nv1*nv1 + eps);
    }

}

void G_OneStrand_3CtrlPoints(Strand strand  ,void (*grad_cost_function) (double * , double *   , double *, double *  , double *  , double *, int  ))
{
    int n_control = strand.n_control;
    int pos_starti   ;
    double * x , *y , *z;
    double  * gx , *gy , *gz;
    int flag_ignore;
    for(int i = 2 ; i < n_control-2 ; i++)
    {
        if (i ==2) [[unlikely]]
            flag_ignore = -1;
        else if (i == (n_control-3)) [[unlikely]]
            flag_ignore =1;
        else [[likely]]
            flag_ignore = 0;
        pos_starti = 4*i;
        copy_pointer(strand.params , pos_starti -4,  &x);
        copy_pointer(strand.g_params , pos_starti-4 , &gx);

        copy_pointer(strand.params , pos_starti , &y);
        copy_pointer(strand.g_params , pos_starti, &gy);

        copy_pointer(strand.params , pos_starti+4, &z);
        copy_pointer(strand.g_params , pos_starti+4, &gz);

        grad_cost_function(x,y,z,gx,gy,gz , flag_ignore);
        //	cout<<"subcomp3 - "<<j<<endl;
    }
    //copy_pointer(strand.params , 0 , &x);
    //copy_pointer(strand.g_params , 0 , &gx);

    //copy_pointer(strand.params , 4 , &y);
    //copy_pointer(strand.g_params , 4 , &gy);


    //grad_cost_function(strand.start ,x,y ,gx , gx,gy , 0);

    //copy_pointer(strand.params  , n_control*4-8 , &y);
    //copy_pointer(strand.g_params , n_control*4-8 , &gy);

    //copy_pointer(strand.params , n_control*4-4 , &z);
    //copy_pointer(strand.g_params , n_control*4-4 , &gz);
    //grad_cost_function(y , z, strand.end, gy ,gz ,gz , 2);


}


void G_AllStrands_3CtrlPoints(vector<Strand> strands ,void (*grad_cost_function) (double * , double *   , double *, double *  , double *  , double *, int  ))
{
    for(int i = 0 ; i < strands.size() ; i ++)
    {	
        if (strands[i].n_control != 1) [[likely]]
            G_OneStrand_3CtrlPoints(strands[i]  , grad_cost_function );
    }
}
