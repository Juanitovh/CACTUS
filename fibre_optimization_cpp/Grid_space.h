/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Grid_space.h
 * Author: jppl
 *
 * Created on February 4, 2021, 6:00 PM
 */

#ifndef GRID_SPACE_H
#define GRID_SPACE_H

#include <vector>
#include <iostream>
#include <array>
#include "strand_class.h"
using namespace std;

#define pii pair<int,int>
#define v1p vector<pii>
#define v2p vector<v1p>
#define v3p vector<v2p>
#define v4p vector<v3p>


void get_bounding_box(double * params , int n_params,double box[2][3] , double max_rad);
void get_index_inGrid(double * point , double  box[2][3]  , double step , int position [3] );
v4p create_grid(int box_dim[3]);
void clear_grid(v4p &grid  , int box_dim[3]);

// void fill_grid(v4p &grid , vector<Strand>  strands, double bounding_box[2][3] , int dim_grid[3] , double step_grid);
void fill_grid(v4p &grid ,const vector<Strand>& strands , double bounding_box[2][3] , int dim_grid[3] , double step_grid);
bool  are_equal(int * p1 , int *p2);

bool  get_27_neighbour_position (int   index[3] , int  new_index[3] , int dim_grid[3] , int n_neghbour);

bool cube_intersection(double p0[3], double step_p0 , double p1[3], double step_p1);
void get_8_vertices_cube(double p0[3], double step , double next_vertex[3] , int n_neighbour );

#endif /* GRID_SPACE_H */
