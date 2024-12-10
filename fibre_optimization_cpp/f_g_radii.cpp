/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   f_g_radii.cpp
 * Author: jppl
 * Created on February 11, 2021, 5:46 PM
 */

#include "f_g_radii.h"

#include"vector_funcions.h"
#include "strand_class.h"
using namespace std;

double f_cost_rads(vector<Strand> strands)
{
    double res = 0;
    double aux;
    for(int i = 0 ; i < strands.size() ; i++)
    {
        for(int j = 0 ; j < strands[i].n_control ; j++ )
        {
            aux = strands[i].params[4*j+3] - strands[i].original_rad[j];
            res += aux*aux;
        }
    }
    return res;

}

double f_cost_rads_half(vector<Strand> strands)
{
    double res = 0 ;
    double aux;
    for(int i = 0 ; i < strands.size() ; i++)
    {
        for(int j = 0 ; j < strands[i].n_control ; j++ )
        {
            aux = maximum( strands[i].params[4*j+3] / strands[i].original_rad[j]  -6.0/5.0,0 );
            res += aux*aux;
        }
    }
    return res;

}

void grad_cost_rads(vector<Strand> strands)
{
    for(int i = 0 ; i < strands.size() ; i++)
    {
        for(int j = 0 ; j < strands[i].n_control ; j++ )
        {
            if (strands[i].n_control ==1)
                continue;
            if (strands[i].params[4*j+3]  < 0)
                strands[i].params[4*j+3]  = strands[i].original_rad[j]*.5;
            strands[i].g_params[4*j+3 ] += 2*( strands[i].params[4*j+3] - strands[i].original_rad[j]);
        }
    }
}

void grad_cost_rads_half(vector<Strand> strands)
{
    double aux ;
    for(int i = 0 ; i < strands.size() ; i++)
    {
        for(int j = 0 ; j < strands[i].n_control ; j++ )
        {
            aux = maximum( strands[i].params[4*j+3] / strands[i].original_rad[j]  -6.0/5.0,0 );
            strands[i].g_params[4*j+3 ] += -2*aux/strands[i].original_rad[j];
        }
    }
}


double  F_CostFunc_Rads_2CntrlPoints(double * control0 , double * control1 )
{
    //double eps = 1e-16;
    double r1 = control0[3];
    double r2 = control1[3];
    double res = (r1-r2);

    //double dista = norm(control0 , control1 , 3);

    //double res = maximum( dista/(eps + r1+ r2)  - .5, 0);
    //	cout <<" ---* " <<res<<endl;
    //     cout<<"radas --- "<< r1 <<"  " <<r2 <<"   " << dista  <<" ---  "<< res <<endl;
    //cout <<" ---- " <<res <<endl;
    return res*res;
}


void G_CostFunc_Rads_2CntrlPoints( double* xi , double * yj , double * gxi , double * gyj , int flag_ignore_right )
{
    double r1,r2 ,base ;
    //double eps = 1e-16;
    r1 = xi[3];
    r2= yj[3];
    //esta resta si es importante
    base = 2*(r1-r2);
    gxi[3] += base;
    if(flag_ignore_right != 1)[[unlikely]]
    {
        gyj[3] += -base;
    }
}
