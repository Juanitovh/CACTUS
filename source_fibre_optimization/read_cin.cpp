#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include "read_cin.h"
#include <cstdlib>
#include <dirent.h>
#include <stdio.h> 
#include <algorithm>

using namespace std;


vector<string> list_files(string path , string extension)
{
  vector<string> files_ex;
  DIR *d;
  struct dirent *dir;
  d = opendir(path.c_str());
  if (d) {
    while ((dir = readdir(d)) != NULL) {
        string file = dir->d_name;
        if (file.find(extension) != string::npos) {
          
            // cout << file << endl;
            files_ex.push_back(file);
        }
      // printf("%s\n", dir->d_name);
    }
    closedir(d);
  }
  sort(files_ex.begin(), files_ex.end());
  return files_ex;
}

vector<string> split_string(string s, string del )
{
    vector<string> tokens;
    int start = 0;

    int end = s.find(del);
    while (end != -1) {
        // cout << s.substr(start, end - start) << endl;
        tokens.push_back(s.substr(start, end - start));
        start = end + del.size();
        end = s.find(del, start);
    }
    // cout << s.substr(start, end - start);
    tokens.push_back(s.substr(start, end - start));
    return tokens;
}

void create_folder(string folder_name)
{
      string command = "mkdir " + folder_name;
    const int dir_err = system(command.c_str());
    if (-1 == dir_err)
    {
        printf("Error creating directory!n");
    }
}

int ** create_matrix_int(int sizeY , int sizeX)
{
  int **ary = new int*[sizeY];
  for(int i = 0; i < sizeY; ++i) {
    ary[i] = new int[sizeX];
  }
  return ary;
}

void clear_matrix(int ** matrix , int sizeY , int sizeX)
{
  for(int i = 0 ; i < sizeY ; i++)
    delete [] matrix[i];
  delete [] matrix;
}

int str_dist(string s, string t)
{
  ulong len_s = s.length();
  ulong len_t = t.length();

  /* base case: empty strings */
  if (len_s == 0) return int(len_t);
  if (len_t == 0) return int(len_s);

  if(len_s == 1 && len_t ==1)
    return s[0] != t[0];

  //Eigen::MatrixXd costos(len_s,len_t);
  int ** costos = create_matrix_int((int) len_s ,  (int) len_t);

  for(unsigned i = 0 ; i < s.size(); i++){
    for (unsigned j = 0 ; j < t.size(); j++){
      costos[i][j] = 0;
      costos[0][j] = j;
    }
    costos[i][0] = i;
  }

  int cost;

  for(unsigned i = 1 ; i < s.size(); i++){
    for (unsigned j = 1 ; j < t.size(); j++){
      /* test if last characters of the strings match */
      if (s[i] == t[j])
        cost = 0;
      else
        cost = 1;

      /* return minimum of delete char from s, delete char from t, and delete char from both */
      costos[i][j] =  min(min( costos[i-1][j] + 1,  costos[i][j-1] + 1) , costos[i-1][j-1] + cost);
    }
  }
  int res = costos[s.length()-1 ][ t.length() -1];
  clear_matrix(costos , len_s , len_t);
  return res ;
}

parameters::parameters()
{

}

void parameters::read_conf_file(string file)
{
  ifstream file_in;
  file_in.open(file);
  string aux;
  string line ;
  while (getline ( file_in, line ))
    {
      if (line[0] == '#' || (line[0] == line[1] == '/'))
      {
        cout <<"-----jumping line"<<endl;
        continue;
      }
      stringstream in (line);
      in >> aux;

      if (str_dist(aux , "path") < 1)
        in >> path;

      if (str_dist(aux , "path_out") < 1)
        in >> path_out;

      if (str_dist(aux , "initial_conditions") < 3)
        in >> initial_conditions;

      //j								if (str_dist(aux , "initial_conditions") < 1)
      //									in>> scale_rad_smoth;

    }

}

