/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   read_strands_txt.h
 * Author: jppl
 *
 * Created on January 12, 2021, 1:54 PM
 */

#ifndef READ_STRANDS_TXT_H
#define READ_STRANDS_TXT_H

#include"strand_class.h"
using namespace std;

vector<string> split(string v , char delimiter);

string get_rad_from_name (string name);

vector <string> read_listdir(string name );

bool is_nfg_strand(string name);

vector <vector<double>> read_nfg_cylinder_file(string path , string name);

vector<vector<vector<double> >> read_nfg_cylinders_in(string path);

vector <vector<vector<double>>>  read_generic_cylinder_list(string path , float *box_size);


bool is_3DR_strand(string name);

vector <vector<double>> read_3DR_cylinder_file(string path , string name);

vector<vector<vector<double> >> read_3DR_cylinders_in(string path);

void write_generic_cylinder_list(vector<Strand> strands , string name , float box_size);
#endif /* READ_STRANDS_TXT_H */
