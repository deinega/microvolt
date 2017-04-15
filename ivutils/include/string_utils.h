#ifndef _STRING_UTILS_H
#define _STRING_UTILS_H

///\en @file string_utils.h \brief Some usefull functions for work with strings

#include <string.h>
#include <string>
#include <vector>
#include "logexc.h"

using namespace std;

///\en calculate substring of string str similar to string::substr
/// used because at some architectures string::substr is unstable
string my_substr(const string &str, size_t pos, size_t n);

//# define SUBSTR(a,b,c) a.substr(b,c)
#define SUBSTR(a,b,c) my_substr(a,b,c)

///\en delete spaces and tabulation at the beginning and the end of the string
void trim(string &str);

///\en format a string 
const char *fmt(const char *format,...);

///\en returns string corresponding to value
const char *val2str(int val);
const char *val2str(double val);

void str2val(const char *s, int &a); /// extract integer number from string
void str2val(const char *s, double &a); /// extract double number from string
void str2val(const char *s, float &a); /// extract double number from string
void str2val(const char *s, string &a); /// just copy string

///\en separate list of strings s to individual strings which will be put into the vector
/// if single_char=0 then d is separator used in list of strings s
/// if single_char=1 then d is string that contains all possible separators,
/// and each separator can be only one character
/// f. e. separate_list("one,two",",",0) returns ("one","two")
/// f. e. separate_list("one,two","e,t",0) returns ("on","wo")
/// f. e. separate_list("one,two.tree",".,",0) returns ()
/// f. e. separate_list("one,two.tree",".,",1) returns ("one","two,"three")
vector<string> separate_list(string s, const string d=",", int single_char=0);

///\en separates list of strings s to individual strings, which are estimated using function ato
/// and saved in variables a1, a2, ...
template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9>
void extract_list(string s,T1 *a1,T2 *a2,T3 *a3,T4 *a4,T5 *a5,T6 *a6,T7 *a7,T8 *a8,T9 *a9){
  vector<string> list=separate_list(s);
  str2val(list[0].c_str(),*a1);
  if(list.size()>1 && a2){
    str2val(list[1].c_str(),*a2);
    if(list.size()>2 && a3){
      str2val(list[2].c_str(),*a3);
      if(list.size()>3 && a4){
        str2val(list[3].c_str(),*a4);
        if(list.size()>4 && a5){
          str2val(list[4].c_str(),*a5);
          if(list.size()>5 && a6){
            str2val(list[5].c_str(),*a6);
            if(list.size()>6 && a7){
              str2val(list[6].c_str(),*a7);
              if(list.size()>7 && a8){
                str2val(list[7].c_str(),*a8);
                if(list.size()>8 && a9){
                  str2val(list[8].c_str(),*a9);
                }
              }
            }
          }
        }
      }
    }
  }
}

template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8>
void extract_list(string s,T1 *a1,T2 *a2,T3 *a3,T4 *a4,T5 *a5,T6 *a6,T7 *a7,T8 *a8){
  extract_list(s,a1,a2,a3,a4,a5,a6,a7,a8,(T1 *)NULL);
}

template<class T1,class T2,class T3,class T4,class T5,class T6,class T7>
void extract_list(string s,T1 *a1,T2 *a2,T3 *a3,T4 *a4,T5 *a5,T6 *a6,T7 *a7){
  extract_list(s,a1,a2,a3,a4,a5,a6,a7,(T1 *)NULL);
}

template<class T1,class T2,class T3,class T4,class T5,class T6>
void extract_list(string s,T1 *a1,T2 *a2,T3 *a3,T4 *a4,T5 *a5,T6 *a6){
  extract_list(s,a1,a2,a3,a4,a5,a6,(T1 *)NULL);
}

template<class T1,class T2,class T3,class T4,class T5>
void extract_list(string s,T1 *a1,T2 *a2,T3 *a3,T4 *a4,T5 *a5){
  extract_list(s,a1,a2,a3,a4,a5,(T1 *)NULL);
}

template<class T1,class T2,class T3,class T4>
void extract_list(string s,T1 *a1,T2 *a2,T3 *a3,T4 *a4){
  extract_list(s,a1,a2,a3,a4,(T1 *)NULL);
}

template<class T1,class T2,class T3>
void extract_list(string s,T1 *a1,T2 *a2,T3 *a3){
  extract_list(s,a1,a2,a3,(T1 *)NULL);
}

template<class T1,class T2>
void extract_list(string s,T1 *a1,T2 *a2){
  extract_list(s,a1,a2,(T1 *)NULL);
}

template<class T1>
void extract_list(string s,T1 *a1){
  extract_list(s,a1,(T1 *)NULL);
}

///\en returns string with used amount of bytes (if bytes<1024), Kb (if bytes<1024*1024)
/// Mb (if bytes<1024^3) or Gb (otherwise)
string size_in_bytes(long long bytes);

///\en returns string in format hours:minutes:seconds
string stime(int sec);

///\en returns rotating stroke |, /, -, \ corresponding to value rot%4
char rot4(int rot);

///\en adds / (UNIX) or \ (Windows) to the directory path
void correct_directory_name(string &dir);

///\en checks if directory exists
bool check_directory(const string &dir);

///\en checks if path is absolute
bool abs_path(const string &path);

///\en Same as fopen() but prepends directory name (if not NULL) to the file name before opening
FILE *open_file(const char *fname, const char *mode, const string &directory="");

/**\en reads table from file to array *ptr
 records number of raws to N1 and number of columns to N2.
 if column_raw, numbers are recorded to ptr corresponding their order in file
 (fastest index corresponds to raws),
 otherwise fastest index corresponds to columns
 if skip_first then first ras is ignored (it can be used if this is text raw) */
/**\ru считывает таблицу из файла в массив *ptr
 в N1 записывается количество строк в файле, в N2 - количество столбцов столбцов.
 если column_raw = true, то данные в массив ptr записываются по порядку
 (быстрый индекс пробегает по строкам),
 в противном случае быстрый индекс пробегает по столбцам
 в случае skip_first первая строчка при считывании пропускается (нужно в том случае, если она текстовая) */
template<class T>
int read_table(const char *fname, T **ptr, int &N1, int &N2, int column_raw=1, int skip_first=0);

///\en records array ptr to the file
/// first is first text raw
template<class T>
int write_table(const char *fname, T *ptr, int N1, int N2, int column_raw=1, const char *first=NULL, int gnu_space=-1);

///\en get interpolated value from table
/// arg is argument, data corresponding to arrument is in j_arg column
/// j_arg is column with values
template<class T>
T get_table_value(T arg, T *ptr, int N1, int N2, int j_arg, int j_val, int column_raw=1, int ord=1);

/**\en ptr is table in format x - y - z - some-variables (less then 3 coordinates is possible, for example just x - y)
  arg corresponds to columns numbers for coordinates x,y,z
  function calculates grid parameters p1 (left lower vertex), p2 (right top point), size sz 
  and order ord for coordinates, from fastest to slowest*/
template<class T>
int table_grid(T *ptr, int N1, int N2, int arg_num, const int *arg, T *p1, T *p2, int *sz, int *ord, int column_raw=1);
/*
/// thin out given table ptr (with N1 rows and N2 columns) and record the result to thptr
/// arg_num - number of dimensions
/// arg - number of corresponding columns
/// thin - thinning coefficient (1 - no thinning out)
int thin_table(double *ptr, int N1, int N2, double **thptr, int &N, int arg_num, const int *arg, const int *thin, 
  bool column_raw=true, int *start=NULL);

/// arg - array of arg_num dimension with column numbers corresponding to coordinates
/// carg - from arg
/// arg value that defines project plane
/// N - new N1
int project_table(double *ptr, int N1, int N2, double **thptr, int &N, int arg_num, const int *arg, 
int carg, double val, bool column_raw=true);
*/
template <class inp_it>
size_t string_scan_t(const char *str,char *format,char *delim, inp_it beg, inp_it end){
  char *buf=new char[strlen(str)+1];
  strcpy(buf,str);
  typename iterator_traits<inp_it>::value_type val;
  char *c1=strtok(buf,delim);
  size_t c=0;
  while(c1 && beg!=end){
    int res=sscanf(c1,format,&val);
    if(res>0){
      c++;
      if(beg!=end)
        *beg++=val;
      else
        break;
    }
    c1=strtok(NULL,delim);
  }
  delete [] buf;
  return c;
}

///\en struct for use in average_files
struct av_comm_t{
  string text;
  int flag;
};

///\en avrerages the contents of space-separated text files containing numbers into the output
/// maximal number of columns: 100
/// maxinal line length: 10000 symbols
/// '#' symbols are treated as comments and transferred to output from the first valid file in sequence  
template < class filename_it>
int average_files(filename_it beg, filename_it end, const char *output){
  vector<av_comm_t> comments; 
  char buff[10000];
  double values[100];
  int ncol=0, nfiles=0; 
  vector< vector<double> > gtable;
 
  bool fix_comm=false;
  FILE *f=NULL;
  for(;beg!=end;++beg){
    vector< vector<double> > table(ncol);
    if(f)
      fclose(f);
    if(!fix_comm)
      comments.clear();
    f=fopen((*beg).c_str(),"rt");
    if(!f){
      LOGMSG(vblWARN,fmt("average_files: can't open the file '%s', skipped!",(*beg).c_str()),0);
      continue;
    }
    int line=0;
    size_t iline=0;
    bool skip=false;
    while(fgets(buff,4000,f)){
      ++line;
      char *comm=strstr(buff,"#");
      av_comm_t s;
      s.flag=-1;
      if(comm){
        s.text=comm;
        *comm=0;
      }
      if(!fix_comm)
        comments.push_back(s);
      
      int num=string_scan_t(buff,"%lf"," ",values,values+100);
      if(num<=0)// completely failed, next string
        continue;
      if(!ncol){
        ncol=num;
        table.resize(ncol);
      }
      else if(num<ncol){
        LOGMSG(vblWARN,fmt("average_files: unexpected number of entries at line %d in '%s', file skipped!",line,(*beg).c_str()),0);
        skip=true;
        break;
      }
      ++iline;
      for(int i=0;i<ncol;i++){
        if(table[i].size()<iline)
          table[i].resize(iline,values[i]);
        else 
          table[i][iline-1]+=values[i];
      }
      if(!fix_comm)
        comments.back().flag=1;
    }
    if(skip)
      continue;
    
    if(!gtable.size())
      gtable=table;
    else if(gtable[0].size()>table[0].size()){
      LOGMSG(vblWARN,fmt("average_files: number of lines in '%s' is too small, file skipped!",(*beg).c_str()),0);
      continue;
    }
    else{ // summing up
      for(size_t i=0;i<gtable.size();i++)
        for(size_t j=0;j<gtable[i].size();j++)
          gtable[i][j]+=table[i][j];
    }
    fix_comm=true; // fixing comments
    nfiles++;
  }
  if(f)
    fclose(f);
  if(!fix_comm)
    return LOGERR(0,fmt("average_files: no valid files with column data found!"),0);
  
  f=fopen(output,"wt");
  if(!f)
    return LOGERR(-1,fmt("average_files: can't open the file '%s' for writing",output),0);
  // retaining comment structure of the first valid file
  size_t iline=0;
  for(size_t i=0;i<comments.size();i++){
    if(comments[i].flag>0){
      for(size_t j=0;j<gtable.size();j++)
        fprintf(f,"%g ",gtable[j][iline]/nfiles);
      ++iline;
    }
    fprintf(f,"%s",comments[i].text.c_str());
    if(!comments[i].text.size())
      fprintf(f,"\n");
  }
  fclose(f);
  return 1;
}

///\en Writes a set of functions given by [\a fbeg, \a fend) into a multicolumn ASCII file
///    evaluated over a set of arguments [\a abeg, \a aend).
///    Arguments and functions must be compatible with the printf formats
///    supplied in \a arg_format and \a func_format respectively.
///    Arguments are printed in the first column, f(arg) values in the next columns.
///    Header is printed at the beginning,  tailer at end of file (with no extra \n after them).
///    Mode \a mode is supplied to fopen second argument (can be "wt" or "at").
///    \return >0 if OK.
template <class arg_it, class func_it>
int write_ascii(const char *fname, arg_it abeg, arg_it aend, 
               func_it fbeg, func_it fend, const char *header="",
               const char *arg_format="%g",const char *func_format="%g", const char *tab=" ",
               const char *mode="wt", const char *tailer=""){
  FILE *f=fopen(fname,mode);
  if(!f)
    return LOGERR(-1,fmt("write_ascii: can't open file '%s' intended for writing (mode: %s)!",fname,mode),0);
  fprintf(f,"%s",header);
  for(;abeg!=aend;++abeg){
    fprintf(f,arg_format,(*abeg));
    for(func_it fit=fbeg;fit!=fend;++fit){
      fprintf(f,"%s",tab);
      fprintf(f,func_format,(*fit)(*abeg));
    }
    fprintf(f,"\n");
  }
  fprintf(f,"%s",tailer);
  fclose(f);
  return 1;
}

#endif
