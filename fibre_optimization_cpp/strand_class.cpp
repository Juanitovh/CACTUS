/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   strand_class.cpp
 * Author: jppl
 * 
 * Created on February 4, 2021, 12:49 PM
 */
#include <iostream>
#include <algorithm>
#include "strand_class.h"

#include"vector_funcions.h"





Strand::Strand()
{
    //	start =NULL;
    //end = NULL;

}

void Strand::initiate(vector<vector<double>> v ,  double * global_params , double *global_g_params )
{
    n_control = v.size();
    // flattening vector<vector> from dimensions  nx4
    double * strand0 = flatten(v, 0, v.size(), 4);

    copy_pointer(global_params  , 0 , &params);
    copy_pointer(global_g_params  , 0 , &g_params);

    copy(strand0 , strand0 + n_control*4, params); // copy from std

    original_rad =  new double [n_control];
    for (int i = 0; i <  n_control; i++) {
        original_rad[i] = params[i * 4 + 3];
    }

    set_fixed_endpoints( );

    delete [] strand0;
}

void Strand::set_fixed_endpoints( )
{


    direction = new double[3];

    double *s0 , *s1;
    copy_pointer(params , 0, &s0);
    copy_pointer(params , n_control * 4 - 4, &s1);

    for (int i = 0 ; i <3 ; i++)
        direction[i] = s1[i] - s0[i];


    double dista = norm(s0, s1, 3);
    L_segment = dista/n_control*1.5;


    /*

       start = new double [4];
       end = new double [4];

     *
     copy_pointer(params , 0, &s0)
     copy_pointer(params , 4, &s1);
     for(int k = 0 ; k<4 ; k++)
     start[k] = s0[k] - s1[k] ;

     double norma = norm(start , 3);
     for(int k = 0 ; k<3 ; k++)
     start[k] = start[k]/norma*s0[3]*4 + s0[k];
     start[3] = s0[3]*.9;


     copy_pointer(params , n_control * 4 - 4, &s0);
     copy_pointer(params , n_control * 4 - 8, &s1);
     for(int k = 0 ; k<4 ; k++)
     end[k] = s0[k] - s1[k] ;

     norma = norm(end , 3);
     for(int k = 0 ; k<3 ; k++)
     end[k] = end[k]/norma*s0[3]*4 + s0[k];
     end[3] = s0[3]*.9;

     double dista = norm(start, end , 3);

*/
}


void Strand::clean_gradient()
{
    // fill_n(g_params , 4*n_control , 0);
     parallel_fill_n(g_params , 4*n_control , 0);
}


void Strand::delete_strand()
{

    //if (start != NULL)
    //{
    //	delete [] start;
    //	delete [] end;
    //}
    //delete [] params;
    //delete [] g_params;
    delete [] original_rad;

}
