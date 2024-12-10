/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   strand_class.h
 * Author: jppl
 *
 * Created on February 4, 2021, 12:49 PM
 */



#ifndef STRAND_CLASS_H
#define STRAND_CLASS_H


#include <iostream>
#include <math.h>
#include<vector>
using namespace std;

class Strand
{
    public:
        double *params;      // The parameters of the control points
        double *g_params; //gradient of  the control points
        int n_control;       //number of control points
        double L_segment;
        double *direction;
        
        
        double *original_rad;

        bool is_fixed;       // determine if the object will be fixed in end, start
        //double * start;     
        //double *end;
    
   

        Strand(); 
        void delete_strand();
        void initiate(vector<vector<double>> v ,  double * global_params , double *global_g_params );
        void set_fixed_endpoints();
        
        void clean_gradient();
 

};

#endif /* STRAND_CLASS_H */
