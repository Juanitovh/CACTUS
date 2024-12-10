/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Grid_space.cpp
 * Author: jppl
 * 
 * Created on February 4, 2021, 6:00 PM
 */

#include <math.h> 
#include<algorithm>
#include <set>
#include <omp.h>
#include "omp_config.h"
#include "Grid_space.h"
#include"vector_funcions.h"



void get_bounding_box(double * params, int n_params ,double box[2][3] , double offset )
{


    for(int i =0 ; i < 3 ; i ++)
    {
        box[0][i] = params[i];
        box[1][i]  = params[i];
    }

    for(int i= 0 ; i < n_params/4 ; i++)
    {
        for(int k = 0 ; k < 3 ; k ++)
        {
            box[0][k] = min (box[0][k] , params[4*i+k] - params[4*i+3] - offset);
            box[1][k]  = max(box[1][k] , params[4*i+k] + params[4*i+3] + offset);

        }
    }

}


void get_index_inGrid(double * point , double  bounding_box [2][3] , double step , int position [3] )
{
    for(int i  = 0 ; i < 3  ; i++)
    {
        position[i] = floor( ( point[i] - bounding_box[0][i])/step   ); 
        //cout << position[i] <<" ";
    }
    //cout<<endl;
}

void get_index_inGrid_BB(double * point , double  bounding_box [2][3] , double step , int position [2][3] )
{
    for(int i  = 0 ; i < 3  ; i++)
    {
        position[0][i] = floor( ( point[i] - bounding_box[0][i] - point[3])/step   ); 
        position[1][i] = floor( ( point[i] - bounding_box[0][i] + point[3])/step   ); 
        //cout << position[i] <<" ";
    }
    //cout<<endl;
}




v4p create_grid(int box_dim[3])
{
    v4p grid (box_dim[0] , v3p ( box_dim[1], v2p (box_dim[2] , v1p (0)  ) ) );

    return grid;
}

void clear_grid(v4p &grid  , int box_dim[3])
{
    for(int x = 0 ;  x< box_dim[0] ; x++)
    {
        for(int y =0 ; y < box_dim[1] ; y ++)
        {
            for(int z = 0 ; z < box_dim[2] ; z++)
            {
                grid[x][y][z].clear();
            }
        }
    }
}


bool  get_27_neighbour_position (int   index[3] , int  new_index[3] , int dim_grid[3] , int n_neghbour)
{
    int xyz[3];
    xyz[0] = ( n_neghbour%3) -1;
    xyz[1] = ( (n_neghbour/3)  %3) -1;
    xyz[2] = ( (n_neghbour/9)%3) -1;


    for(int i =0 ; i < 3 ; i++)
    {
        new_index[i] = index[i] + xyz[i];
        if ( new_index[i] <0 ||  new_index[i] >= dim_grid[i])
        {
            //cout<<"super error index grid!!"<< new_index[0]<<" " << new_index[1] << " "<<new_index[2] <<endl;
            new_index[0] = new_index[1] = new_index[2] = -100;
            // cout<<"error index outside grid!"<<endl;
            return false ;
        }
    }
    return true;
}

bool point_inside_cube(double point[3] , double origin_box[3] , double step )
{
    for(int i = 0 ; i < 3 ; i++)
    {
        if ( point[i] < origin_box[i] || point[i] > (origin_box[i] + step))
        {
            return false;
        }
    }
    return true;
}

bool sphere_inside_cube(double sphere[4] , double origin_box[3] , double step )
{
    double aux_point[3];
    for(int i = 0 ; i < 3 ; i++)
        aux_point[i] = sphere[i] - sphere[3];

    bool aux = point_inside_cube(aux_point , origin_box , step);

    for(int i = 0 ; i < 3 ; i++)
        aux_point[i] = sphere[i] + sphere[3];
    bool aux2 = point_inside_cube(aux_point , origin_box , step);

    if (aux && aux2)
        return true;

    return false;

}


bool sphere_cube_intersection(double origin_box[3] , double step , double sphere[4])
{
    double mid_point[3];
    double v_dir[3];
    for(int i = 0 ; i < 3 ; i++)
    {
        mid_point[i] = origin_box[i] + step/2;
        v_dir[i] = mid_point[i] - sphere[i];
    }
    //normalize v_dir
    double v_dir_norm = 0;
    for(int i = 0 ; i < 3 ; i++)
        v_dir_norm += v_dir[i]*v_dir[i];
    v_dir_norm = sqrt(v_dir_norm);
    for(int i = 0 ; i < 3 ; i++)
        v_dir[i] /= v_dir_norm;

    double closest_point[3];
    for(int i = 0 ; i < 3 ; i++)
    {
        closest_point[i] = sphere[i] + v_dir[i]*sphere[3];
    }
    if ( point_inside_cube(closest_point , origin_box , step) )
    {
        return true;
    }
    else
    {
        return false;
    }

}




void get_8_vertices_cube(double p0[3], double step , double next_vertex[3] , int n_neighbour )
{
    next_vertex[0] = p0[0] +  ( n_neighbour%2)     *step;
    next_vertex[1] = p0[1] +  ( (n_neighbour/2)%2) *step;
    next_vertex[2] = p0[2] +  ( (n_neighbour/4)%2) *step;

}


bool cube_intersection(double p0[3], double step_p0 , double p1[3], double step_p1)
{

    double next_vertex[3];

    for(int nn =0 ; nn < 8 ; nn++)
    {
        int flag_dim=0;
        get_8_vertices_cube(p1, step_p1 , next_vertex , nn);
        //cout <<"n_neighbour "<< next_vertex[0] << " " <<next_vertex[1] <<" " <<next_vertex[2]<<endl;
        for(int k = 0 ; k < 3 ; k++)
        {
            if( p0[k] <= next_vertex[k] && next_vertex[k] <= p0[k] + step_p0  )
                flag_dim++;

        }
        if (flag_dim ==3)
            return true;
    }

    return false;
}

void point_interpolate(double *p0 , double *p1 , int step , int n , double pi[4])
{
    for (int i = 0; i < 4 ; i++)
    {
        pi[i] = p0[i]+ (p1[i] - p0[i])*step/((float)n);
    }
}


void fill_grid(v4p &grid ,const vector<Strand>& strands , double bounding_box[2][3] , int dim_grid[3] , double step_grid)
{

    double start, end;
    // int n_threads = 36 ; //omp_get_num_threads();


    vector <v4p> all_grids  ( GLOBAL_NUM_THREADS ,  v4p  (dim_grid[0] , v3p ( dim_grid[1], v2p (dim_grid[2] , v1p (0)  ) ) ) );


    start = omp_get_wtime();


#pragma omp parallel for  schedule(dynamic,10) num_threads(GLOBAL_NUM_THREADS)
    for(int i = 0 ; i < strands.size() ; i++)
    {
        int id_thread = omp_get_thread_num();

        double *p0;
        double *p1;
        //double *point = new double[4];
        double point[4];
        int  pos_grid[3];
        int  pos_grid_new[3];
        int  pos_grid_bb[2][3];

        double box_sphere[3];
        double box_voxel[3];
        vector<int> pos_v (3);

        // Strand &s = strands[i];

        const Strand &s = strands[i];
        //for (int k= 1 ; k < s.n_control -2 ; k++)
        // cout<<"error here" <<endl;

        if (s.n_control == 1)[[unlikely]]
        {
            continue;
        }

        for (int k= 2 ; k < s.n_control -3  ; k++)
        {
            copy_pointer(s.params , 4*k  , &p0);
            copy_pointer(s.params , 4*k +4  , &p1);

            // for interpolate points
            pii sphere (i,k);


            vector<vector<int>>future_set ;
            // int n_step =  2+ norm(p0 , p1 , 3)/( .5*p0[3] + .5*p1[3] );
            int n_step =  2+ norm(p0 , p1 , 3)/(2 * sqrt(step_grid));
            // int n_step =  2+ norm(p0 , p1 , 3)/( .2*p0[3] + .2*p1[3] );
            for (int step = 0 ; step <= n_step; step ++)
            {
                point_interpolate(p0 , p1 , step , n_step , point);
                // get_index_inGrid( point ,   bounding_box , step_grid , pos_grid );
                get_index_inGrid_BB( point ,   bounding_box , step_grid , pos_grid_bb );

                for (int bb_a = pos_grid_bb[0][0] ; bb_a <= pos_grid_bb[1][0] ; bb_a ++)
                {
                    for(int bb_b = pos_grid_bb[0][1] ; bb_b <= pos_grid_bb[1][1] ; bb_b ++)
                    {
                        for (int bb_c = pos_grid_bb[0][2] ; bb_c <= pos_grid_bb[1][2] ; bb_c ++)
                        {
                            pos_v[0] = bb_a;
                            pos_v[1] = bb_b;
                            pos_v[2] = bb_c;
                            future_set.push_back(pos_v);
                        }
                    }
                }

            }


            set<vector<int>> seto (future_set.begin(), future_set.end());
            for (set<vector<int>>::iterator it=seto.begin(); it!=seto.end(); ++it)
            {
                pos_v = *it; 
                all_grids[id_thread][ pos_v[0] ][ pos_v[1] ][ pos_v[2] ].push_back(sphere);
            }

            // cout <<"--after all_grids" << endl;
            future_set.clear();
            seto.clear();
        }
        // cout <<"--seto futuro clean" << endl;
    }


    // cout<<"--end for"<<endl;
    // cout <<"merging grids" << endl;


// #pragma omp parallel for collapse(3) num_threads(GLOBAL_NUM_THREADS)
    for (int current_thread =  0 ; current_thread < GLOBAL_NUM_THREADS ; current_thread++ )
    {

        for(int x = 0; x <  dim_grid[0]; x++)
        {
            for(int y = 0 ; y < grid[0].size() ; y++)
            {
                for(int z = 0 ; z < grid[0][0].size() ; z++)
                {
                    v1p  &current_stack = all_grids[current_thread][x][y][z];

                    for (int ii = 0 ; ii < current_stack.size()  ; ii++)
                    {

                        pii &ball1 = current_stack[ii]; 

                        grid[x][ y][z ].push_back(ball1);

                    }
                }
            }
        }

    }

    for (int current_thread =  0 ; current_thread < GLOBAL_NUM_THREADS ; current_thread++ )
        clear_grid(all_grids[current_thread]  ,  dim_grid);


    // cout <<"--end merging grids" << endl;

    end = omp_get_wtime();
    //cout <<" -- total time filling " << end - start << endl;
}




/// case for single ball
// else
// {
//     //cout <<"zero case" << endl;
//     copy_pointer(s.params , 0  , &p0);
//
//     point_interpolate(p0 , p0 , 0 , 2 , point);
//     // for interpolate points
//     pii sphere (i,-1);
//
//
//     vector<vector<int>>future_set ;
//     get_index_inGrid( point ,   bounding_box , step_grid , pos_grid );
//
//     box_sphere[0] = point[0] -  point[3];
//     box_sphere[1] = point[1]   - point[3];
//     box_sphere[2] = point[2]  - point[3];
//
//     for (int n = 0 ; n < 27 ; n++)
//     {
//         bool flag_27 = get_27_neighbour_position(pos_grid , pos_grid_new , dim_grid, n);
//         if (flag_27 == false)
//             continue;
//         box_voxel[0] = pos_grid_new[0]*step_grid + bounding_box[0][0];
//         box_voxel[1] = pos_grid_new[1] *step_grid+ bounding_box[0][1];
//         box_voxel[2] = pos_grid_new[2]*step_grid+ bounding_box[0][2];
//
//         if (    cube_intersection(box_sphere , 2*point[3] , box_voxel , step_grid ) ||  cube_intersection( box_voxel , step_grid , box_sphere , 2*point[3] ) )
//         {
//             pos_v[0] = pos_grid_new[0] ;
//             pos_v[1] = pos_grid_new[1] ;
//             pos_v[2] = pos_grid_new[2] ;
//             future_set.push_back(pos_v);
//
//         }
//     }
// }

