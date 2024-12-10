/*
   i To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   read_strands_txt.cpp
 * Author: jppl
 * 
 * Created on January 12, 2021, 1:54 PM
 */



#include <iostream>
#include <fstream>
#include <string>

#include <iostream>
#include <algorithm>
#include <vector>
#include <string>



#include <dirent.h>
#include <sys/types.h>

using namespace std;

#include <vector>

#include "read_strands_txt.h"


//reading files nfg style

vector<string> split(string v , char delimiter)
{
    int start=0 , end = 0;
    int n = v.size();
    vector<string> res;
    while(end <= n)
    {
        if (v[end] == delimiter || end == n)
        {
            string aux = v.substr(start , end - start);
            res.push_back(aux); 
            start = end+1;
            end = start;
            continue;
        }
        end++;

    }
    return res;
}

string get_rad_from_name (string name)
{
    auto res = split(name , 'r');
    string s2 = res[2];
    s2 = s2.substr( 0 , s2.size() - 4);
    return s2;
}



vector <string> read_listdir(string name )
{
    DIR *dr;
    struct dirent *en;
    vector<string> files ;
    dr = opendir(name.c_str()); //open all directory
    if (dr) 
    {
        while ((en = readdir(dr)) != NULL) {
            files.push_back(en->d_name);
            //cout<<" \n"<<en->d_name;
            //print all directory name
        }
        closedir(dr); //close all directory
    }
    return files;
}

bool is_nfg_strand(string name)
{
    string base = "strand";

    for(int i = 0 ; i < min(base.size() , name.size()) ; i++)
    {
        if (base[i] != name[i])
        {
            return false;
        }
    }
    return true;
}


vector <vector<double>> read_nfg_cylinder_file(string path , string name)
{
    ifstream myfile (path + '/' + name);
    vector<vector<double>> strands;

    //  	cout <<" -- "<<get_rad_from_name(name)<<endl;
    double radis = stof(get_rad_from_name(name));

    vector<double> si (4);
    if (myfile.is_open())
    {
        double x;
        int i = 0;

        si[3] = radis;
        while ( myfile >> x )
        {
            si[i] = x; 

            if (i ==2)
            {

                strands.push_back(si);
            }
            i = (i+1)%3; 

        }
        myfile.close();

    }
    return strands;
}


vector<vector<vector<double> >> read_nfg_cylinders_in(string path)
{
    vector<vector<vector<double>>> strands;
    vector <string> files = read_listdir(path);
    for(int i = 0 ; i < files.size() ; i++)
    {
        if (is_nfg_strand(files[i]))
        {
            vector<vector<double>> si = read_nfg_cylinder_file(path , files[i]);
            //cout<<  files[i] << " - " << si.size() <<endl;
            strands.push_back(si);
        }
    }
    return strands;

}


/// reading files my_style


bool is_3DR_strand(string name)
{
    string base = "fibre";

    for(int i = 0 ; i < min(base.size() , name.size()) ; i++)
    {
        if (base[i] != name[i])
        {
            return false;
        }
    }
    return true;
}


vector <vector<double>> read_3DR_cylinder_file(string path , string name)
{
    ifstream myfile (path + '/' + name);
    vector<vector<double>> strands;

    //  	cout <<" -- "<<get_rad_from_name(name)<<endl;

    vector<double> si (4);
    if (myfile.is_open())
    {
        double x;
        int i = 0;


        while ( myfile >> x )
        {
            si[i] = x;
            if (i ==3)
                strands.push_back(si);

            i = (i+1)%4; 
        }
        myfile.close();

    }
    return strands;
}


vector<vector<vector<double> >> read_3DR_cylinders_in(string path)
{
    vector<vector<vector<double>>> strands;
    vector <string> files = read_listdir(path);
    for(int i = 0 ; i < files.size() ; i++)
    {
        if (is_3DR_strand(files[i]))
        {
            vector<vector<double>> si = read_3DR_cylinder_file(path , files[i]);
            //cout<<  files[i] << " - " << si.size() <<endl;
            strands.push_back(si);
        }
    }
    return strands;

}


vector <vector<vector<double>>>  read_generic_cylinder_list(string path , float *box_size)
{

    ifstream infile; 
    infile.open(path);
    float box_lenght ;
    infile>>box_lenght;
    *box_size = box_lenght;
    // gola


    int n_cyilinders;
    infile >> n_cyilinders; 
    vector <vector<vector<double>>> cylinder_list;
    for ( int i = 0 ; i < n_cyilinders ; i++)
    {
        int n_control;
        infile >> n_control;
        vector<vector<double>> cylinder;
        for(int k = 0 ; k < n_control; k++)
        {
            vector<double> punto(4);
            for(int m =0 ; m < 4 ; m++)
                infile>>punto[m];
            cylinder.push_back(punto);
        }
        cylinder_list.push_back(cylinder); 
    }

    infile.close();
    return cylinder_list;
}


void write_generic_cylinder_list(vector<Strand> strands , string name , float box_lenght)
{
    ofstream myfile (name);
    if (myfile.is_open())
    {
        myfile<< box_lenght  <<endl;
        myfile<<strands.size()  <<endl;
        for(int i = 0  ; i < strands.size() ; i++ )
        {
            myfile<<strands[i].n_control<<endl;
            for (int k = 0 ; k < strands[i].n_control ; k++)
            {
                for(int m = 0 ; m < 4 ; m++)
                {
                    myfile<<strands[i].params[k*4 + m]<<" ";
                }		
                myfile<<endl;
            }

        }
        myfile.close();
    }
    else cout << "Unable to open file";

}
