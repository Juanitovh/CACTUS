/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   f_g_2CtrlPoints.cpp
 * Author: jppl
 * 
 * Created on February 11, 2021, 2:40 PM
 */

#include "f_g_2CtrlPoints.h"
#include <iostream>
#include"vector_funcions.h"
#include "strand_class.h"
#include <functional>
using namespace std;


double  F_CostFunc_Lenght_2CntrlPoints(double * control0 , double * control1 )
{
    double eps = 1e-16;
    double r1 = control0[3];
    double r2 = control1[3];

    double dista = norm(control0 , control1 , 3);

    double res = maximum( dista/(eps + r1+ r2)  - 2, 0);
    //	cout <<" ---* " <<res<<endl;
    // cout<<"radas --- "<< r1 <<"  " <<r2 <<"   " << dista  <<" ---  "<< res <<endl;
    //cout <<" ---- " <<res <<endl;
    return res*res;
}



double  F_CostFunc_Lenght_2CntrlPoints_hooke(double * control0 , double * control1 , double L_len )
{
    double eps = 1e-16;
    double r1 = control0[3];
    double r2 = control1[3];

    double k = 2/(1/r1+1/r2) * 10;

    double dista = norm(control0 , control1 , 3);

    double res  = (dista - L_len);
    if (res<0) res = 0;
    //	cout <<" ---* " <<res<<endl;
    // cout<<"radas --- "<< r1 <<"  " <<r2 <<"   " << dista  <<" ---  "<< res <<endl;
    //cout <<" ---- " <<res <<endl;
    return res*res*k/2.0;
}


double F_OneStrand_2CtrlPoints(Strand strand ,   double (*cost_function)(double* , double* , double)  )
{
    int n_control = strand.n_control;

    double  *control0 , *control1;
    int pos_start;

    double costo = 0;
    for(int i = 0; i < n_control-1  ; i++)
    {
        pos_start = 4*i;
        copy_pointer(strand.params , pos_start        ,&control0);
        copy_pointer(strand.params , pos_start+4  ,&control1);

        //this is for the generic function of two points
        costo += cost_function(control0 , control1 , strand.L_segment);
    }
    //copy_pointer(strand.params , 0   ,&control0);
    //copy_pointer(strand.params , 4*n_control-4  ,  &control1);

    //costo += cost_function(control0 , strand.start,strand.L_segment);
    //costo += cost_function(control1 , strand.end ,strand.L_segment);
    //cout <<" costo un strando "<<costo<<endl; 
    return costo;
}


double F_AllStrands_2CtrlPoints( vector<Strand> strands  ,double (*cost_function)(double *, double* , double) )
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


//////////////////////////////////////////////////////////////




void G_CostFunc_Lenght_2CntrlPoints( double* xi , double * yj , double * gxi , double * gyj , int flag_ignore_right , double L)
{
    double r1,r2;
    double eps = 1e-16;
    r1 = xi[3];
    r2= yj[3];
    //esta resta si es importante
    double xy_norm = norm(xi , yj , 3);
    //xy_norm = maximum(xy_norm , 0.0001);

    double maxi = maximum(xy_norm/(r1+r2+ eps) -2  , 0);

    double grad_r =-xy_norm/(eps + (r1+r2)*(r1+r2)) *.01 ;
    if(maxi > 0)
    {

        //cout<<"second part"<<endl;
        //gxi[3] += grad_r*2*maxi;

        //if(flag_ignore_right != 1)
        //	gyj[3] += grad_r*2*maxi;

        for(int i = 0 ; i < 3 ; i++)
            gxi[i]  +=  2*maxi*(xi[i] - yj[i])/(eps + (r1 +r2)* xy_norm);


        if(flag_ignore_right != 1)
        {
            for(int i = 0 ; i < 3 ; i++)
                gyj[i]  += -2*maxi*(xi[i] - yj[i])/(eps + (r1 +r2)* xy_norm);

        }		

    }
    //else 
    //	maxi = maxi/1000.0;




}

void G_CostFunc_Lenght_2CntrlPoints_hooke( double* xi , double * yj , double * gxi , double * gyj , int flag_ignore , double L_len )
{
    // flag_ignore -1 ignore left , +1 ignore right , 0 ignore none
    double r1,r2;
    double eps = 1e-16;
    r1 = xi[3];
    r2= yj[3];
    double k = 2/(1/r1+1/r2) *10;
    //esta resta si es importante
    double xy_norm = norm(xi , yj , 3);

    double res = (xy_norm - L_len);
    //xy_norm = maximum(xy_norm , 0.0001);

    //cout<<"second part"<<endl;
    //gxi[3] += grad_r*2*maxi;

    //if(flag_ignore_right != 1)
    //	gyj[3] += grad_r*2*maxi;

    if (res <0)
        return;

    if (flag_ignore != -1)[[likely]]
    {
        for(int i = 0 ; i < 3 ; i++)
            gxi[i]  +=  k*res*(xi[i] - yj[i])/(eps +  xy_norm);
    }


    if(flag_ignore != 1) [[likely]]
    {
        for(int i = 0 ; i < 3 ; i++)
            gyj[i]  += k*res*(-xi[i] + yj[i])/(eps + xy_norm);

    }	



}




void G_OneStrand_2CtrlPoints( Strand strand ,  void (*gradient_cost_function ) (double*  , double *  , double *  , double *  , int , double) )
{
    int pos_starti   ;
    int n_control = strand.n_control;
    double * x_actual , *y_actual;
    double  * gx_actual , *gy_actual;
    //cout<<" -- strand ii " <<endl;
    int flag_ignore=0;
    for(int i = 1 ; i < n_control-2 ; i++)
    {

        if (i == 1 ) [[unlikely]]
            flag_ignore =-1;
        else if (i == (n_control-3) ) [[unlikely]]
            flag_ignore =1;
        else  [[likely]]
            flag_ignore =0 ;
        pos_starti = 4*i;
        copy_pointer(strand.params , pos_starti , &x_actual);
        copy_pointer(strand.g_params , pos_starti , &gx_actual);

        copy_pointer(strand.params , pos_starti+4 , &y_actual);
        copy_pointer(strand.g_params , pos_starti+4, &gy_actual);


        gradient_cost_function(x_actual ,y_actual, gx_actual , gy_actual, flag_ignore , strand.L_segment);
    }
    //copy_pointer(strand.params , 0 , &x_actual);
    //copy_pointer(strand.g_params , 0 , &gx_actual);
    //gradient_cost_function(x_actual ,strand.start, gx_actual , gy_actual, 1,strand.L_segment);

    //copy_pointer(strand.params , n_control*4-4 , &x_actual);
    //copy_pointer(strand.g_params , n_control*4-4, &gx_actual);
    //gradient_cost_function(x_actual ,strand.end, gx_actual , gy_actual, 1 ,strand.L_segment);

}


void G_AllStrands_2CtrlPoints(vector<Strand>  strands, void (*gradient_cost_function ) (double*  , double *  , double *  , double *  , int , double  ) )
{
    // #pragma omp parallel for schedule(dynamic,20) reduction(+:res)
    for(int i = 0 ; i < strands.size() ; i ++)
    {		

        if (strands[i].n_control != 1) [[likely]]
            G_OneStrand_2CtrlPoints(strands[i] , gradient_cost_function );		
    }
}
