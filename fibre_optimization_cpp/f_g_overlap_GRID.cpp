/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   f_g_overlap_GRID.cpp
 * Author: jppl
 * 
 * Created on February 5, 2021, 7:49 PM
 */

#include "f_g_overlap_GRID.h"

#include "Grid_space.h"
#include"vector_funcions.h"
#include <algorithm> 

#include <cstring>
#include <omp.h>
#include "omp_config.h"
#include<set>
#include <map>
//double scale_overlap_gs;


// double clap(double a)
// {
//
//     return min(max(0.0,a) , 1.0);
// }

constexpr double clap(double a) {
    return (a < 0.0) ? 0.0 : (a > 1.0) ? 1.0 : a;
}
void min_distance_lines(double* __restrict__ p0, double* __restrict__ p1, double* __restrict__ q0, double* __restrict__ q1, double res[2]) {
    double lx = p0[0] - q0[0];
    double ly = p0[1] - q0[1];
    double lz = p0[2] - q0[2];

    double vx = p1[0] - p0[0];
    double vy = p1[1] - p0[1];
    double vz = p1[2] - p0[2];

    double bx = q1[0] - q0[0];
    double by = q1[1] - q0[1];
    double bz = q1[2] - q0[2];

    double bb = bx*bx + by*by + bz*bz;
    double vv = vx*vx + vy*vy + vz*vz;
    double bv = bx*vx + by*vy + bz*vz;
    double lb = lx*bx + ly*by + lz*bz;
    double lv = lx*vx + ly*vy + lz*vz;

    double det = bb * vv - bv * bv;

    if (det != 0) [[likely]]
    {

        double inv_det = 1.0 / det;
        res[0] = clap((-bb * lv + lb * bv) * inv_det);
        res[1] = clap((lb * vv - bv * lv) * inv_det);
    } else {
        if (bb == 0) {
            res[0] = vv ? clap(-lv / vv) : 0;
            res[1] = 0;
        } else {
            res[0] = 0;
            res[1] = vv ? 0 : clap(lb / bb);
        }
    }
}



void min_distance_lines3(double *p0, double *p1, double *q0, double *q1, double res[2]) {
    double l[3], v[3], b[3];

    for (int i = 0; i < 3; ++i) {
        l[i] = p0[i] - q0[i];
        v[i] = p1[i] - p0[i];
        b[i] = q1[i] - q0[i];
    }

    double bb = dot_prod3d(b, b);
    double vv = dot_prod3d(v, v);
    double bv = dot_prod3d(b, v);
    double lb = dot_prod3d(l, b);
    double lv = dot_prod3d(l, v);

    double det = bb * vv - bv * bv;
    double det_r = lb * vv - bv * lv;
    double det_t = -bb * lv + lb * bv;

    if (det != 0) {
        res[0] = det_t / det;
        res[1] = det_r / det;
    } else {
        if (bb == 0) {
            res[1] = 0;
            res[0] = (vv == 0) ? 0 : clap(-lv / vv);
        } else {
            res[0] = 0;
            res[1] = (vv == 0) ? clap(lb / bb) : 0;
        }
        return;
    }

    double t = res[0];
    double r = res[1];
    double clap_t = clap(t);
    double clap_r = clap(r);

    if (clap_r != r || clap_t != t) {
        res[0] = clap_t;
        res[1] = clap_r;
        return;
        double clap_r1 = clap((clap_t * bv + lb) / bb);
        double clap_t2 = clap((clap_r * bv - lv) / vv);

        double aux1 = clap_t * clap_t * vv + clap_r1 * clap_r1 * bb + 2 * (clap_t * lv - clap_r1 * lb - clap_r1 * clap_t * bv);
        double aux2 = clap_t2 * clap_t2 * vv + clap_r * clap_r * bb + 2 * (clap_t2 * lv - clap_r * lb - clap_r * clap_t2 * bv);

        if (aux1 > aux2) {
            res[0] = clap_t2;
            res[1] = clap_r;
        } else {
            res[0] = clap_t;
            res[1] = clap_r1;
        }
    }
}

void min_distance_lines2(double *p0 , double *p1 , double *q0  , double *q1 , double res[2])
{

    //	double *l , *v , *b;
    //	l = new double[3];
    //	v = new double[3];
    //b = new double[3];


    double l[3] , v[3] , b[3];

    for (int i = 0 ; i < 3;  i ++)
    {
        l[i] = p0[i] - q0[i];
        v[i] = p1[i] - p0[i];
        b[i] = q1[i] - q0[i];
    }
    double bb =dot_prod3d(b,b);
    double vv =dot_prod3d(v,v);
    double bv =dot_prod3d(b,v);
    double lb =dot_prod3d(l,b);
    double lv =dot_prod3d(l,v);

    double det = bb* vv- bv*bv;
    double det_r = lb*vv - bv*lv;
    double det_t = -bb*lv + lb*bv;

    double t=.5;
    double r=.5; 


    if (det != 0 )
    {
        // t = det_t / det;
        res[0] = det_t/det;
        // r = det_r/det;
        res[1] = det_r/det;

    }
    else
    {

        if (bb == 0 && vv!=0)
        {
            // r = 0;
            // t = (r*bv-lv)/vv;
            // t =  clap(t);
            res[1] = 0;
            res[0] = clap((-lv)/vv);

        }
        else if (bb != 0 && vv == 0)
        {
            // t = 0;
            // r = (t*bv+lb)/bb;
            // r = clap(r);
            res[0] = 0;
            res[1] = clap((lb)/bb);
        }
        else
        {
            // t = 0;
            // r= 0;
            res[0] = 0;
            res[1] = 0;
        }

        // res[0] = t; 
        // res[1] = r;
        return;

    }


    double clap_t1 = clap(t);
    double clap_r2 = clap(r);

    if ( clap_r2 != r || clap_t1 != t )
    {
        double clap_r1 = clap((clap_t1*bv+lb)/bb ) ;
        double clap_t2 = clap((clap_r2*bv-lv)/vv ) ;

        double aux1 = clap_t1*clap_t1*vv + clap_r1*clap_r1*bb + 2*(clap_t1*lv - clap_r1*lb - clap_r1*clap_t1*bv );
        double aux2 = clap_t2*clap_t2*vv + clap_r2*clap_r2*bb + 2*(clap_t2*lv - clap_r2*lb - clap_r2*clap_t2*bv );

        if (aux1 > aux2)
        {

            res[0] = clap_t2; 
            res[1] = clap_r2;
        } 
        else
        {

            res[0] = clap_t1; 
            res[1] = clap_r1;
        }


    }
    else{
        res[0] = t; 
        res[1] = r;
    }

}



double f_overlap_two_capsules_cuadratic(double * p0 , double * p1  , double * q0 , double *q1)
{
    double eps = 1e-16;
    double tr[2];
    min_distance_lines(p0,p1,q0,q1, tr);
    double t ,r;
    t = tr [0] ; 
    r = tr[1];
    //delete [] tr;
    double r1 = p0[3]*(1-t) + p1[3]*(t);
    double r2 = q0[3]*(1-r) + q1[3]*(r);

    q0[3] = abs(q0[3]);
    q1[3] = abs(q1[3]);
    p0[3] = abs(p0[3]);
    p1[3] = abs(p1[3]);


    double control0[4] , control1[4];

    for (int i = 0 ; i < 3 ; i++)
    {
        control0[i] = p0[i] + (p1[i] - p0[i])*t;
        control1[i] = q0[i] + (q1[i] - q0[i])*r;
    }

    double dista = norm(control0 , control1 , 3);


    double p0p1 = norm(p0 , p1 , 3);
    double q0q1 = norm(q0 , q1 , 3);

    double mulp =1;
    double mulq =1;
    if (p0p1 == 0)
    {
        mulq *= 2;
        p0p1 = 2*p0[3];
    }

    if (q0q1 == 0)
    {
        mulp  *= 2;
        q0q1 = 2*q0[3];
    }
    if(mulp == mulq & mulp ==2)
    {
        mulp=10;
        mulq=10;
    }
    p0p1 *=mulp;
    q0q1 *=mulq;



    double res = (1.05-   (dista) /(r1+r2  + eps));
    res = maximum(res , 0 );
    //	cout <<" ---* " <<res<<endl;

    //cout <<" ---- " <<res <<endl;
    double fin = res *res*p0p1*q0q1*r1*r2;
    if (fin < 0  )
    {
        cout <<"puta vida" << res*res << " " << p0p1 << " " << q0q1 << " " << r1 << " " <<" "<< r2 <<endl;
    }
    return fin  ;//*res ; //sqrt(res);
}




double f_overlap_all_strands_grid(v4p grid , vector<Strand> strands  )
{
    double res = 0;
    //#pragma omp parallel 
#pragma omp parallel for collapse(3) schedule(dynamic,5) reduction(+:res)
    for(int x = 0; x < grid.size() ; x++)
    {
        for(int y = 0 ; y < grid[0].size() ; y++)
        {
            for(int z = 0 ; z < grid[0][0].size() ; z++)
            {
                v1p  &current_stack = grid[x][y][z];
                //cout << "stack size "<< current_stack.size() << " - " << x <<" - " << y <<" - " << z <<endl;
                double *p0 , *p1, *q0 , *q1;

                for (int i = 0 ; i < current_stack.size()  ; i++)
                {
                    pii &ball1 = current_stack[i]; 
                    copy_pointer(strands[ball1.first].params ,ball1.second*4         , &p0 );
                    if (ball1.second != -1) [[likely]]
                        copy_pointer(strands[ball1.first].params ,ball1.second*4 + 4 , &p1 );
                    else [[unlikely]]
                        copy_pointer(strands[ball1.first].params ,ball1.second*4     , &p1 );
                    for(int j = 0 ; j < i ; j++)
                    {
                        pii &ball2 = current_stack[j]; 
                        //cout<< ball1.first << " - " <<  ball2.first <<endl;
                        if ( ball1.first != ball2.first)
                        {
                            //							cout << x << ""
                            //cout <<"comparando different balls"<<endl;
                            copy_pointer(strands[ball2.first].params ,ball2.second*4       , &q0);
                            if (ball2.second != -1) [[likely]]
                                copy_pointer(strands[ball2.first].params ,ball2.second*4+4 , &q1);
                            else  [[unlikely]]
                                copy_pointer(strands[ball2.first].params ,ball2.second*4    , &q1);
                            double aux = f_overlap_two_capsules_cuadratic(p0 , p1 , q0 , q1);
                            res += aux;
                        }
                    }
                }
            }
        }
    }
    //delete [] strand0;
    //delete [] strand1;
    return res;
}


///////////////////////////////////////////////////////////////////




void grad_ovelap_two_capsules_cuadratic(double * __restrict__ p0 , double * __restrict__ p1 , double * __restrict__ q0 , double * __restrict__ q1 , double * __restrict__ gp0 , double * __restrict__ gp1 , double * __restrict__  gq0 , double * __restrict__ gq1 , int ignore_p0 , int ignore_q0 , int size_cube, bool compress)
{
    double eps = 1e-16;
    double tr[2];
    min_distance_lines(p0,p1,q0,q1 , tr);
    double t ,r;
    t = tr [0] ; 
    r = tr[1];
    //delete [] tr;
    double r0 = p0[3]*(1-t) + p1[3]*(t);
    double r1 = q0[3]*(1-r) + q1[3]*(r);

    double p0p1 = norm(p0 , p1 , 3);
    double q0q1 = norm(q0 , q1 , 3);

    double mulp =1;
    double mulq =1;
    // if (p0p1 == 0)
    // {
    //     mulq *= 2;
    //     p0p1 = 2*p0[3];
    // }
    //
    // if (q0q1 == 0)
    // {
    //     mulp  *= 2;
    //     q0q1 = 2*q0[3];
    // }
    // if(mulp == mulq & mulp ==2)
    // {
    //     mulp=10;
    //     mulq=10;
    // }
    // p0p1 *=mulp;
    // q0q1 *=mulq;


    double control0[4] , control1[4];

    for (int i = 0 ; i < 3 ; i++)
    {
        control0[i] = p0[i] + (p1[i] - p0[i])*t;
        control1[i] = q0[i] + (q1[i] - q0[i])*r;
    }

    double dista = norm(control0 , control1 , 3);

    double res = (1.05-   (dista) /(r0+r1  + eps));


    //cout<<"second part"<<endl;

    //#pragma omp atomic
    //if (res > 0 && compress)
    // if (0)
    // {	
    //     double grad_r0 = ( 2*res * dista / (   (r1+r0)*(r1+r0)  + eps  ) *r0 + res*res )* p0p1 * q0q1*  1e-3;
    //     double grad_r1 = ( 2*res * dista / (   (r1+r0)*(r1+r0)  + eps  ) *r1 + res*res )* p0p1 * q0q1*  1e-3;
    //
    //     if (ignore_p0 != -1)
    //         gp0[3] += grad_r0 * (1-t)  ; 		
    //     if (ignore_p0 != 1)
    //         gp1[3]  += grad_r0 * (t)  ; 		
    //
    //     if (ignore_q0 == -1)
    //         gq0[3] += grad_r1 * (1-r) ; 		
    //     if (ignore_q0 == 1)
    //         gq1[3]  += grad_r1 * (r)  ; 	
    //
    //
    // }
    if (res <= 0    ) //&& compress 
    {
        return ;
        // t = r = .5; 
        // r0 = p0[3]*(1-t) + p1[3]*(t);
        // r1 = q0[3]*(1-r) + q1[3]*(r);
        // for (int i = 0 ; i < 3 ; i++)
        // {
        //     control0[i] = p0[i] + (p1[i] - p0[i])*t;
        //     control1[i] = q0[i] + (q1[i] - q0[i])*t;
        // }
        // dista = norm(control0 , control1 , 3);
        // res = (1.05-   (dista) /(r0+r1  + eps));		
        // res  =-res/(5000000 *size_cube*size_cube );
        // //				res = -res*res*res*res;
        // if (res >1 )
        //     cout <<" VVVVVV" << endl;
        //
        //
        // //		res = 0; 
        // if  (res < 1 ) 
        //     res = -res*res;
    }
    //else if (res < 0  && compress == false)
    //	res = 0;

    //res  =res/(500000 *size_cube*size_cube );


    //p0 ,p1 , 
    //q0 , q1 

    double middle_p[3] , middle_q[3];
    for (int i = 0 ; i < 3 ; i++)
    {
        middle_p[i] = p0[i]*.5 +  p1[i]*.5;
        middle_q[i] = q0[i]*.5 +  q1[i]*.5;
    }

    double norm_middle_p = norm( middle_p , 3);
    double norm_middle_q = norm( middle_q , 3);

    for(int i = 0 ; i < 3 ; i++)
    {
        middle_p[i] /= norm_middle_p;
        middle_q[i] /= norm_middle_q;
    }



    double current_gp0[3] , current_gp1[3] , current_gq0[3] , current_gq1[3];


    float condition_number = .0;
    double pre_derivative =  -2*res/((r0+r1)*dista +eps);
    double aux;


    if (ignore_p0 != -1) [[likely]]
    {
        for(int i = 0 ; i < 3 ; i++)
        {
            aux = p0p1*q0q1*pre_derivative * ( p0[i]*(1+t*t -2*t) + p1[i]*(t-t*t) + q0[i]*(-1+t+r-t*r) + q1[i]*(-r+t*r) + condition_number)  ;
            aux = aux + q0q1*res*res*(p0[i] - p1[i] + condition_number)/(p0p1 + eps);
            gp0[i] += aux*r1*r0;
        }
        // current_gp0[i] = aux*r1*r0;
    }
    if (ignore_p0 != 1) [[likely]]
    {

        for(int i = 0 ; i < 3 ; i++)
        {
            aux =	pre_derivative*( p0[i]*(t-t*t) + p1[i]*(t*t) + q0[i]*(-t +t*r) + q1[i] *( -t*r) + condition_number);
            aux = aux + q0q1*res*res*(p1[i] - p0[i] + condition_number )/(p0p1 + eps);
            gp1[i]  += aux*r1*r0;
            // current_gp1[i] = aux*r1*r0;
        }
    }


    if (ignore_q0 != -1) [[likely]]
    {

        for(int i = 0 ; i < 3 ; i++)
        {
            aux = pre_derivative * ( q0[i]*(1+r*r -2*r) + q1[i]*(r-r*r) + p0[i]*(-1+t+r-t*r) + p1[i]*(-t+t*r) + condition_number)  ;
            aux = aux + p0p1*res*res*(q0[i] - q1[i] + condition_number)/(q0q1 + eps);
            gq0[i] += aux*r1*r0;
            // current_gq0[i] = aux*r1*r0;
        }
    }
    if (ignore_q0 != 1 ) [[likely]]
    {

        for(int i = 0 ; i < 3 ; i++)
        {
            aux =pre_derivative*( q0[i]*(r-r*r ) + q1[i]*(r*r) + p0[i]*(-r +t*r) + p1[i] *( -t*r) + condition_number);
            aux = aux + p0p1*res*res*(q1[i] - q0[i] + condition_number )/(q0q1 + eps);
            gq1[i]  += aux*r1*r0;
            // current_gq1[i] = aux*r1*r0;
        }
    }



}

// // return ;
// double projection;
// double expand_force  = 0 ; //.15;
// if (ignore_p0 != -1) [[likely]]
// {
//     projection = dot_prod(current_gp0 , middle_p , 3);
//     if (projection > 0)
//         projection = 0;
//
//     for (int i = 0 ; i < 3 ; i++)
//         gp0[i] += current_gp0[i]  + projection*middle_p[i]*expand_force;
// }
//
// if (ignore_p0 != 1) [[likely]]
// {
//     projection = dot_prod(current_gp1 , middle_p , 3);
//     if (projection > 0)
//         projection = 0;
//
//     for (int i = 0 ; i < 3 ; i++)
//         gp1[i] += current_gp1[i]  + projection*middle_p[i]*expand_force;
// }
//
// if (ignore_q0 != -1) [[likely]]
// {
//     projection = dot_prod(current_gq0 , middle_q , 3);
//     if (projection > 0)
//         projection = 0;
//
//     for (int i = 0 ; i < 3 ; i++)
//         gq0[i] += current_gq0[i]  + projection*middle_q[i]*expand_force;
// }
//
// if (ignore_q0 != 1) [[likely]]
// {
//     projection = dot_prod(current_gq1 , middle_q , 3);
//     if (projection > 0)
//         projection = 0;
//
//     for (int i = 0 ; i < 3 ; i++)
//         gq1[i] += current_gq1[i]  + projection*middle_q[i]*expand_force;
// }
//


// }



//	if (seto.size() != 0)
//		cout << "seto " << seto.size();


//Certainly! Here's the code with std::unordered_map using the same #define structure:

// cpp
// Copy code
#include <unordered_map>
#include <vector>
#include <functional>
#include <cstring>
#include <array>

#define pii pair<int,int>
// #define map_ii_vd map<pii , vector<double>>

#define map_ii_vd map<pii , array<double, 8>>
// #define map_ii_vd std::unordered_map<pii, std::vector<double>, std::hash<pii>>
//
// namespace std {
//     template <> struct hash<pii> {
//         size_t operator()(const pii& p) const {
//             return p.first * 1000033 + p.second;
//         }
//     };
// }
//
//
// map_ii_vd strands_in_box(const std::vector<pii>& current_stack)
// {
//     map_ii_vd updates;
//
//     std::vector<double> empty(8, 0.0);
//     for (const pii& ball1 : current_stack)
//     {
//         updates[ball1] = empty;
//     }
//
//     return updates;
// }

// #define map_ii_vd std::unordered_map<pii, std::array<double, 8>, std::hash<pii>>

// namespace std {
//     template <> struct hash<pii> {
//         size_t operator()(const pii& p) const {
//             return p.first * 1000033 + p.second;
//         }
//     };
// }

map_ii_vd strands_in_box(const std::vector<pii>& current_stack) {
    map_ii_vd updates;

    std::array<double, 8> empty = {};
    // vector<double> empty(8, 0.0);

    for (const pii& ball1 : current_stack) {
        updates[ball1] = empty;
    }

    return updates;
}

double ** strands_in_box_arr(int size) {

    // allocate the memory in an array of size  ( 8 size)
    cout <<" &*&*&*%^&*%^&*%^& " << size << endl  ;
    double **updates = new double*[size];
    updates[0] = new double[size * 8];
    memset(updates[0], 0, size * 8 * sizeof(double));


    for(int i = 1; i < size; ++i)
        updates[i] = updates[i-1] + 8;

    return updates;
}

double ** declare_buffer_upates(int size) {

    // allocate the memory in an array of size  ( 8 size)
    // double **updates = new double*[size];
    // updates[0] = new double[size * 8];
    // // memset(updates[0], 0, size * 8 * sizeof(double));
    //
    // for(int i = 1; i < size; ++i)
    //     updates[i] = updates[i-1] + 8;

    double ** updates = new double*[size];

    for(int i = 0; i < size; ++i)
        updates[i] = new double[8];


    return updates;
}

void empty_buffer_upates(double ** updates , int size)
{
    // memset(updates[0], 0, size * 8 * sizeof(double));
    for(int i = 0; i < size; ++i)
        memset(updates[i], 0, 8 * sizeof(double));
}

double ** resize_buffer_upates(double ** updates2 , int size , int size_old)
{
    // delete[] updates2[0];
    // delete[] updates2;
    //
    // double ** updates = new double*[size];
    //
    // updates = new double*[size];
    // updates[0] = new double[size * 8];
    //
    // for(int i = 1; i < size; ++i)
    //     updates[i] = updates[i-1] + 8;
    //
    // return updates;

    for(int i = 0; i < size_old; ++i)
        delete[] updates2[i];

    delete[] updates2;

    double ** updates = new double*[size];

    for(int i = 0; i < size; ++i)
        updates[i] = new double[8];

    return updates;
}


void free_buffer_updates(double ** updates )
{
    delete[] updates[0];
    delete[] updates;
}

// map_ii_vd strands_in_box(v1p current_stack)
// {
//     map_ii_vd updates ; 
//
//     vector<double> empty (8 , 0.0);
//
//      vector<double> empty(8, 0.0);
//     for (const pii& ball1 : current_stack)
//     {
//         updates[ball1] = empty;
//     }
//
//
//     return updates;
//
// }



void g_overlap_all_strands_grid(v4p grid , vector<Strand> strands , double * grad_useless , int n_total_params  , bool compress )
{
    //double res = 0;
    //#pragma omp parallel 

    // cout<<"starting gradient"<<endl;
#pragma omp parallel for collapse(2)  schedule(dynamic,5)  num_threads(GLOBAL_NUM_THREADS)  //reduction(+:grad_useless[:n_total_params])
    for(int x = 0; x < grid.size() ; x++)
    {
        for(int y = 0 ; y < grid[0].size() ; y++)
        {

            // cout<<"thread paralelzation"<<endl;
            int size_local_buffer = 50;
            int old_size = size_local_buffer;
            double ** updates =  declare_buffer_upates( size_local_buffer );
            // cout<<"declare buffer"<<endl;
            for(int z = 0 ; z < grid[0][0].size() ; z++)
            {

                v1p  &current_stack = grid[x][y][z];
                int my_size = current_stack.size();

                if (my_size > size_local_buffer)
                {
                    // cout << "resizing buffer " << endl;
                    while (my_size > size_local_buffer)
                        size_local_buffer *= 2;
                    // cout << "new size " << size_local_buffer << endl;
                    updates = resize_buffer_upates(updates , size_local_buffer, old_size);
                    old_size = size_local_buffer;

                }
                // cout<<"empty buffer"<<endl;
                empty_buffer_upates(updates , my_size);
                // cout<<"emptied buffer"<<endl;

                // map_ii_vd updates =  strands_in_box(current_stack );
                //cout << "stack size "<< current_stack.size() << " - " << x <<" - " << y <<" - " << z <<endl;
                double *p0 , *p1 , *gp0 , *gp1 , *q0 , *q1 , * gq0 , * gq1;

                for (int i = 0 ; i < current_stack.size()  ; i++)
                {
                    pii &ball1 = current_stack[i]; 

                    copy_pointer(strands[ball1.first].params ,ball1.second*4 , &p0 );
                    // copy_pointer(&updates[ball1][0]          , 0             , &gp0 );
                    copy_pointer(&updates[i][0]          , 0             , &gp0 );

                    if (ball1.second != -1) [[likely]]
                    {

                        copy_pointer(strands[ball1.first].params  ,ball1.second*4 +4 , &p1 );
                        // copy_pointer(&updates[ball1][4]           , 0                , &gp1 );
                        copy_pointer(&updates[i][4]           , 0                , &gp1 );
                    }
                    else [[unlikely]]
                    {

                        copy_pointer(strands[ball1.first].params  ,ball1.second*4  , &p1 );
                        // copy_pointer(&updates[ball1][0]           , 0              , &gp1 );
                        copy_pointer(&updates[i][0]           , 0              , &gp1 );
                    }
                    for(int j = 0 ; j < i ; j++)
                    {
                        pii &ball2 = current_stack[j]; 
                        //cout<< ball1.first << " - " <<  ball2.first <<endl;

                        if ( ball1.first != ball2.first)
                        {
                            int ignore_p=0, ignore_q =0; 
                            if (ball1.second == 2) [[unlikely]]
                                ignore_p = -1;
                            else if (ball1.second == strands[ball1.first].n_control-4  || ball1.second ==-1) [[unlikely]]
                                ignore_p = 1;
                            else [[likely]]
                                ignore_p = 0;

                            if (ball2.second == 2) [[unlikely]]
                                ignore_q = -1;
                            else if (ball2.second == strands[ball2.first].n_control-4  || ball2.second ==-1) [[unlikely]]
                                ignore_q = 1;
                            else  [[likely]]
                                ignore_q = 0;

                            //							cout << x << ""
                            //cout <<"comparando different balls"<<endl;
                            copy_pointer(strands[ball2.first].params  ,ball2.second*4 , &q0);
                            // copy_pointer(&updates[ball2][0]           ,0              , &gq0);
                            copy_pointer(&updates[j][0]           ,0              , &gq0);

                            if (ball2.second != -1) [[likely]]
                            {
                                copy_pointer(strands[ball2.first].params ,ball2.second*4 +4 , &q1);
                                // copy_pointer(&updates[ball2][4]          ,0                 , &gq1);
                                copy_pointer(&updates[j][4]          ,0                 , &gq1);
                            }
                            else [[unlikely]]
                            {

                                copy_pointer(strands[ball2.first].params ,ball2.second*4    , &q1);
                                // copy_pointer(&updates[ball2][0]          ,0                 , &gq1);
                                copy_pointer(&updates[j][0]          ,0                 , &gq1);

                            }
                            grad_ovelap_two_capsules_cuadratic(p0 , p1 , q0 ,q1,  gp0 , gp1 , gq0 , gq1 , ignore_p,ignore_q , current_stack.size() , compress);
                        }
                    }
                }
                // cout<<"computing gradient done "<<endl;
                int size_stack = current_stack.size();
                // map_ii_vd::iterator it;
                pii balli (0,0);
                // vector<double> gradi;
                // array<double, 8> gradi = {};
                double *gradi;

                pii it;
                // for( it=updates.begin() ; it != updates.end() ; it++)
                // for( it=current_stack.begin() ; it != current_stack.end() ; it++)
                // cout<< " emptying stack!! almost there!"<<endl;
                for( int i = 0 ; i < size_stack ; i++)
                {
                    // it = current_stack[i];
                    // balli = it->first;
                    // gradi = it->second;

                    balli = current_stack[i];
                    copy_pointer(&updates[i][0] , 0 , &gradi);

                    if (balli.second == -1) [[unlikely]]
                        balli.second =0;

                    for (int a  = 0 ; a <8 ; a++)
                    {
#pragma omp atomic update
                        strands[balli.first].g_params[balli.second*4 + a] += gradi[a];
                        // strands[balli.first].g_params[balli.second*4 + a] += 0 ;//gradi[a];
                    }
                }


            }

            // cout<<"freening memory"<<endl;
            free_buffer_updates(updates );
            // cout<<"freed memory"<<endl;
        }
    }
    //delete [] strand0;
    //delete [] strand1;
}
