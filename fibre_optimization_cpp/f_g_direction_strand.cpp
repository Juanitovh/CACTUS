/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   f_g_direction_strand.cpp
 * Author: jppl
 * 
 * Created on April 24, 2021, 5:59 PM
 */

#include "f_g_direction_strand.h"
double F_cost_angle_strand(double * x , double * y , double * direction_strand)
{
    double eps = 1e-16;
    double r1 = x[3];
    double r2 = y[3];


    double nv0 = norm(y,x,3);
    double d0 = norm(direction_strand,3);


    double dot_v0v1 = 0;
    for(int i =0 ; i < 3 ; i++)
    {
        dot_v0v1 += (direction_strand[i] )*(y[i] - x[i]);
    }
    dot_v0v1 = abs(dot_v0v1);

    double base = (1- dot_v0v1/(nv0*d0 + eps));

    double pre = 2.0/ (r1+r2);

    return pre*base*base;
}



double F_OneStrand_2CtrlPoints(Strand strand,  double (*cost_function)(double* , double* , double *) )
{
    int n_control = strand.n_control;
    double  *control0 , *control1 ;
    int pos_start;


    //cout<<"dist are cool"<<endl;

    double costo = 0;
    for(int i = 1; i < n_control-2  ; i++)
    {
        pos_start = 4*i;
        copy_pointer(strand.params , pos_start        ,&control0);
        copy_pointer(strand.params , pos_start+4  ,&control1);

        costo += cost_function(control0 , control1, strand.direction);
    }

    return costo;
}


double F_AllStrands_2CtrlPoints( vector<Strand> strands  ,double (*cost_function)(double *, double* , double*) )
{
    double res = 0;
    #pragma omp parallel for schedule(dynamic,20) reduction(+:res)
    for(int i = 0; i < strands.size() ; i++)
    {		
        if (strands[i].n_control != 1) [[likely]]
            res += F_OneStrand_2CtrlPoints(strands[i]  ,  cost_function);
    }
    return res;
}


/////////////


void G_CostFunc_angle_strand(double * x, double * y  ,  double * gx , double * gy , double *direction_strand,  int flag_ignore )
{
    //flag_ignore  : -1, ignore left, +1 ignore right.
    double eps =1e-8;
    double r1 = x[3];
    double r2 = y[3];

    double zeros[3] = {0,0,0};
    double nxy = norm(y,x,3);

    double d0 = norm(direction_strand , 3);

    double dot_v0v1 = 0;
    for(int i =0 ; i < 3 ; i++)
        dot_v0v1 += (direction_strand[i] )*(y[i] - x[i]);

    int signo =0 ;
    if (dot_v0v1 < 0 )
    {
        signo = -1 ;
    }
    else 
    {
        signo = 1; 
    }
    dot_v0v1 = abs(dot_v0v1);
    double base = 2*(1- dot_v0v1/(nxy*d0 + eps));

    double pre = 2.0/(r1+r2);



    //double grad_r = .0*   2*pre/3.0*base*base;
    //cout<<"second part"<<endl;

    //gy[3] += grad_r;

    if(flag_ignore != -1) [[likely]]
        for(int i = 0 ; i < 3 ; i++)
            gx[i]  += pre* base*(  -nxy*d0*(direction_strand[i])*signo   +   dot_v0v1 *    ( d0*(x[i]-y[i])/(nxy+ eps)  )      )/ (nxy*nxy*d0*d0 + eps);

    if(flag_ignore != 1) [[likely]]
    {
        //do x
        //	gx[3] += grad_r;
        for(int i = 0 ; i < 3 ; i++)
            gy[i]  += pre* base*(  -nxy*d0*(direction_strand[i])*signo  +   dot_v0v1 *    ( d0*(y[i]-x[i])/(nxy+ eps)  )      )/ (nxy*nxy*d0*d0 + eps);

    }	


}

void G_OneStrand_3CtrlPoints(Strand strand  ,void (*grad_cost_function) (double * , double *   ,  double *  , double *  , double *,  int  ))
{
    int n_control = strand.n_control;
    int pos_starti   ;
    double * x , *y ;
    double  * gx , *gy ;
    int flag_ignore;
    for(int i = 1 ; i < n_control-2 ; i++)
    {
        if (i ==1) [[unlikely]]
            flag_ignore = -1;
        else if (i == (n_control-3)) [[unlikely]]
            flag_ignore =1;
        else [[likely]]
            flag_ignore = 0;
        pos_starti = 4*i;
        copy_pointer(strand.params , pos_starti,  &x);
        copy_pointer(strand.g_params , pos_starti, &gx);

        copy_pointer(strand.params , pos_starti  +4     , &y);
        copy_pointer(strand.g_params , pos_starti +4, &gy);


        grad_cost_function(x,y,gx,gy ,strand.direction, flag_ignore);
        //	cout<<"subcomp3 - "<<j<<endl;
    }


}


void G_AllStrands_3CtrlPoints(vector<Strand> strands ,void (*grad_cost_function) (double * , double *   ,  double *  , double *  , double *,  int ))
{
    for(int i = 0 ; i < strands.size() ; i ++)
    {	
        if (strands[i].n_control != 1) [[likely]]
            G_OneStrand_3CtrlPoints(strands[i]  , grad_cost_function );
    }
}
