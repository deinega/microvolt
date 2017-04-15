#include "string_utils.h"
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef UNIX
#include <sys/types.h>
#include <dirent.h>
#else
#ifdef _MSC_VER
#include "windows.h"
#endif
#endif

string my_substr(const string &str, size_t pos, size_t n){
  if(n<=0) return "";
  if(pos+n>str.length())
    n=str.length()-pos;
  char *ch = new char[n+1];
  for(size_t i=0; i<n; i++)
    ch[i] = pos+i<str.length() ? str[pos+i] : 0;
  ch[n]=0;
  string my(ch);
  delete[]ch;
  return my;
}

void trim(string &str) {
  int i1=0;
  for (; i1<(int)str.length(); i1++) {
    if (str[i1]!=' ' && str[i1]!=char(9)) break;
  }
  if (i1==(int)str.length()) {
    str="";
  }
  else {
    int i2=(int)str.length()-1;
    for (; i2>=0; i2--) {
      if (str[i2]!=' ' && str[i2]!=char(9)) break;
    }
    str=SUBSTR(str, i1, i2-i1+1);
  }
}

const char *fmt(const char *format,...){
  va_list args;
  va_start(args,format);
  static char buff[2048];
  vsnprintf(buff,2048,format,args);
  return buff;
}

const char *val2str(int val){
  return fmt("%d",val);
}

const char *val2str(double val){
  return fmt("%g",val);
}

void str2val(const char *s, int &a){a=atoi(s);}
void str2val(const char *s, double &a){a=atof(s);}
void str2val(const char *s, float &a){a=float(atof(s));}
void str2val(const char *s, string &a){a=s;}

vector<string> separate_list(string s, const string d, int single_char){
  vector<string> list;
  int f;
  do{
    if(!single_char) // searching for the whole string
      f=(int)s.find(d);
    else{ //searching for different characters of d
      f=-1;
      for(int i=0;i<(int)d.length();i++){
        int fi=(int)s.find(d[i]);
        if(fi>=0 && (f<0 || fi<f))
          f=fi;
      }
    }
    if(f<0)
      break;
    list.push_back(SUBSTR(s,0,int(f)));
    s=SUBSTR(s,int(f+1),s.length()-int(f+1));
    //f=(int)s.find(d);
  }while(1);

  list.push_back(s);
  return list;
}

string size_in_bytes(long long bytes){
  char sz[100];
  long long Kb=1024;
  if(bytes<Kb)
    sprintf(sz,"%d b",(int)bytes);
  else if(bytes<Kb*Kb)
    sprintf(sz,"%.2f Kb",float(bytes/1024.));
  else if(bytes<Kb*Kb*Kb)
    sprintf(sz,"%.2f Mb",float((bytes/1024.)/1024.));
  else
    sprintf(sz,"%.2f Gb",float(((bytes/1024.)/1024.)/1024.));
  return string(sz);
}

string stime(int sec){
  int s=sec%60;
  int m=(sec/60)%60;
  int h=sec/(60*60);
  char res[100];
  sprintf(res,"%d:%02d:%02d",h,m,s);
  return res;
}

char rot4(int rot){
  int r=rot%4;
  if(r==0)
    return '|';
  if(r==1)
    return '/';
  if(r==2)
    return '-';
  if(r==3)
    return '\\';
  return ' ';
}

void correct_directory_name(string &dir) {
  size_t len=dir.length();
  if (len && dir[len-1]!='\\' && dir[len-1]!='/') {
#ifdef UNIX
    dir+='/';
#else
    dir+='\\';
#endif
  }
}

bool check_directory(const string &dir) {
# ifdef UNIX
  DIR *D=opendir(dir.c_str());
  if (D==NULL)
    return false;
  closedir(D);
  return true;
# else
#ifndef VS2005
  TCHAR Buffer[MAX_PATH];
  DWORD dwRet = GetCurrentDirectory(MAX_PATH, Buffer);
  if (dwRet==0)
    return false;
  if (dwRet>MAX_PATH)
    return false;
  if(!SetCurrentDirectory((LPCSTR)dir.c_str()))
    return false;
  if(!SetCurrentDirectory(Buffer))
    return false;
#endif
  return true;
# endif
}

bool abs_path(const string &path){
#ifdef UNIX
  if(path.length()>1)
    if(path[0]=='/')
      return true;
#else
  if(path.length()>3)
    if(path[1]==':')
      return true;
#endif
  return false;
}

FILE *open_file(const char *fname, const char *mode, const string &directory){
  string name=directory+string(fname);
  FILE *f=fopen(name.c_str(), mode);
  return f;
}

template<class T>
int read_table(const char *fname, T **ptr, int &N1, int &N2, int column_raw, int skip_first){
  FILE *f=NULL;
  if((f=fopen(fname, "r"))==NULL)
    return LOGERR(-1,fmt("File %s is not found\n",fname),0);

  const int num=200000;
  char str[num];
  char str_ex[num];
  char *res;
  // count strings number
  N1=0;
  while(!feof(f)){
    res=fgets(str,num,f);
    if(!res)
      break;
//    if(strlen(str)+1>num)
//      return LOGERR(-1,fmt("Too lagrge line sizes in file %s\n",fname),0);
    int count=0;
    for(size_t i=0;i<strlen(str);i++){
      if(str[i]!=' ' && str[i]!=char(9) && str[i]!=char(10)){
        count=1;
        break;
      }
    }
    if(count){
      N1++;
      if((!skip_first && N1==1) || (skip_first && N1==2))
        strcpy(str_ex,str);
    }
  }
  if(skip_first)
    N1--;
  if(N1<=0)
    return LOGERR(-1,fmt("File %s is empty\n",fname),0);
  // count columns number
  N2=0;
  bool del=true;
  for(size_t i=0;i<strlen(str_ex);i++){
    if(str_ex[i]==',' || str_ex[i]==' ' || str_ex[i]==char(9) || str_ex[i]==char(10)){
      if(del==false){
        del=true;
      }
    }
    else if(del==true){
      del=false;
      N2++;
    }
  }
  fclose(f);

  if(N2==0)
    return LOGERR(-1,fmt("Table could not be read from file %s\n",fname),0);

  *ptr = new T[N1*N2];
  f=fopen(fname, "r");
  if(skip_first)
    res=fgets(str,num,f);

  int ind=0;
  int strn=0;
  while(!feof(f)){
    res=fgets(str,num,f);
    if(!res)
      break;
    int count=0;
    for(size_t i=0;i<strlen(str);i++){
      if(str[i]!=' ' && str[i]!=char(9) && str[i]!=char(10)){
        count=1;
        break;
      }
    }
    if(!count)
      continue;

    char dbl[num];
    int dbli=0,j=0;
    bool del=true;
    for(size_t i=0;i<strlen(str)+1;i++){
      if(str[i]==',' || str[i]==' ' || str[i]==char(9) || str[i]==char(10) || str[i]=='\n' || i==strlen(str)){
        if(del==false){
          del=true;
          j++;
          dbl[dbli++]='\n';
          dbl[dbli++]=0;
          int index=ind;
          if(!column_raw){
            int c=ind%N2, r=ind/N2;
            index=N1*c+r; // changed order
//            index=N2*r+c; // current order
          }
          (*ptr)[index]=T(atof(dbl));
          char snan[5]={'N','a','N',10,0};
//          int nan=strcmp(dbl,"NaN");
          int nan=strcmp(dbl,snan);
          if(nan==0)
//            (*ptr)[index]=1e50;
            (*ptr)[index]=0;
          else
            (*ptr)[index]*=1;
          if(j>N2 || ((*ptr)[index]==0 && nan!=0 && (dbl[0]!='0' && ((dbl[0]!='.' && dbl[0]!='-') || dbl[1]!='0')))){
            delete[](*ptr);
            fclose(f);
            return LOGERR(-1,fmt("Table in file %s contains not only numbers\n",fname),0);
          }
          ind++;
        }
      }
      else if(del==true){
        del=false;
        dbli=0;
      }
      if(del==false){
        dbl[dbli++]=str[i];
      }
    }
    if(j!=N2){
      delete[](*ptr);
      fclose(f);
      return LOGERR(-1,fmt("Table could not be read from file %s\n",fname),0);
    }
    strn++;
  }

  fclose(f);
  return 1;
}

template<class T>
int write_table(const char *fname, T *ptr, int N1, int N2, int column_raw, const char *first, int gnu_space){
  FILE *f=NULL;
  if((f=fopen(fname, "w"))==NULL)
    return LOGERR(-1,fmt("File %s can not be created\n",fname),0);

  if(first)fprintf(f,"%s\n",first);

  T gval=0;
  if(gnu_space>=0){
    gval = column_raw ? ptr[gnu_space] : ptr[gnu_space*N1];
  }

  for(int i=0;i<N1;i++){
    if(gnu_space>=0 && i>0){
      T gval1 = column_raw ? ptr[i*N2+gnu_space] : ptr[gnu_space*N1+i];
      if(gval1<gval)
        fprintf(f,"\n");
      gval=gval1;
    }
    for(int j=0;j<N2;j++){
      double val = column_raw ? ptr[i*N2+j] : ptr[j*N1+i];
      fprintf(f,"%g",val);
      if(j<N2-1)
        fprintf(f,"%c",9);
    }
    if(i<N1-1)
      fprintf(f,"\n");
  }

  fclose(f);
  return 1;
}

template<class T>
T get_table_value(T arg, T *ptr, int N1, int N2, int j_arg, int j_val, int column_raw, int ord){
  T first = column_raw ? ptr[j_arg] : ptr[j_arg*N1];
  T last = column_raw ? ptr[j_arg+(N1-1)*N2] : ptr[j_arg*N1+N1-1];
  if(ord && arg <= first || !ord && arg >= first)
    return column_raw ? ptr[j_val] : ptr[j_val*N1];
  else if (ord && arg >= last || !ord && arg <= last)
    return column_raw ? ptr[j_val+(N1-1)*N2] : ptr[j_val*N1+N1-1];

  for(int i=0;i<N1-1;i++){
    int ind = ord ? i+1 : N1-1-(i+1);
    T next = column_raw ? ptr[j_arg+ind*N2] : ptr[j_arg*N1+ind];
    if(next>arg){
      T narg[2];
      T val[2];
      for(int j=0;j<2;j++){
        ind = ord ? i+j : N1-1-(i+j);
        narg[j] = column_raw ? ptr[j_arg+ind*N2] : ptr[j_arg*N1+ind];
        val[j] = column_raw ? ptr[j_val+ind*N2] : ptr[j_val*N1+ind];
      }
      return (val[0]*(narg[1]-arg)+val[1]*(arg-narg[0]))/(narg[1]-narg[0]);
    }
  }
  return 0;
}

template<class T>
int table_grid(T *ptr, int N1, int N2, int arg_num, const int *arg, T *p1, T *p2, int *sz, int *ord, int column_raw){

  for(int i=0;i<arg_num;i++){
    if(p1)
      p1[i] = column_raw ? ptr[arg[i]] : ptr[arg[i]*N1];
    if(p2)
      p2[i] = column_raw ? ptr[arg[i]+(N1-1)*N2] : ptr[arg[i]*N1+(N1-1)];
  }

  if(!sz && !ord)
    return 1;

  int *ord_val = ord ? ord : new int[arg_num];
  for(int i=0;i<arg_num;i++)
    ord_val[i]=-1;
  int ord_ind=0;
  for(int ni=0;ni<N1;ni++){ // find the order (which coordinate is fastest, which coordinate is slowest)
    for(int i=0;i<arg_num;i++){
      if(ord_val[i]<0){
        T val = column_raw ? ptr[arg[i]+ni*N2] : ptr[arg[i]*N1+ni];
        if(val!=p1[i])
          ord_val[i]=ord_ind++;
      }
    }
  }

  if(!sz)
    return 1;

  const int max=3;
  T val[max];

  for(int i=0;i<arg_num;i++){
    sz[i]=1;
    val[i]=-1e10;
  }
  sz[ord[0]]=0;

  ord_ind=0; // we will find sz for fastest coordinate first, and so on
  for(int i=0;i<N1;i++){ // find size sz
    while(true){
      T v = column_raw ? ptr[arg[ord[ord_ind]]+i*N2] : ptr[arg[ord[ord_ind]]*N1+i];
      if(v<val[ord_ind]){
        val[ord_ind]=v;
        ord_ind++;
        if(ord_ind>=arg_num-1){
          sz[ord[ord_ind]]=N1;
          for(int j=0;j<arg_num-1;j++)
            sz[ord[ord_ind]]/=sz[ord[j]];
          i=N1; // break from for
          break;
        }
        continue;
      }
      else if(v>val[ord_ind]){
        val[ord_ind]=v;
        sz[ord[ord_ind]]++;
      }
      break;
    }
  }
  if(!ord)delete []ord_val;
  return 1;
}
/*
int thin_table(double *ptr, int N1, int N2, double **thptr, int &N, int arg_num, const int *arg, const int *thin, 
bool column_raw, int *start){

  const int max=3;
  int sz[max]={0,1,1};
  double val[max];
  if(table_grid(ptr,N1,N2,arg_num,arg,NULL,NULL,sz,column_raw)<0)
    return -1;
  int fill=0;

  N=1;
  for(int j=0;j<arg_num;j++){
    N*=int(ceil(double(sz[j])/thin[j]));
//    val[j] = column_raw ? ptr[arg[fill]] : ptr[arg[fill]*N1];
    val[j] = -1e10;
  }
  *thptr = new double[N*N2];

  int st[max]={0,0,0};
  int ord[max];
  for(int i=0;i<max;i++){
    if(start)
      st[i]=start[i];
    ord[i]=st[i];
  }
  ord[0]-=1;
  int wr=0,wrmax=0;
  for(int i=0;i<max;i++){
    if(!ord[i])
      wr|=1<<i;
    wrmax|=1<<i;
  }

  for(int i=0,i1=0;i<N1;i++){
    fill=0;
    while(true){
      double v = column_raw ? ptr[arg[fill]+i*N2] : ptr[arg[fill]*N1+i];
      if(fill==2)
        int a=1;
      if(v<val[fill]){
        val[fill]=v;
        ord[fill]=st[fill];
        if(!ord[fill])
          wr|=1<<fill;
        fill++;
        continue;
      }
      if(v>val[fill]){
        val[fill]=v;
        ord[fill]++;
        if(ord[fill]%thin[fill])
          wr&=~(1<<fill);
        else
          wr|=1<<fill;
      }
      break;
    }
    if(wr==wrmax){
      for(int j=0;j<N2;j++){
        int index = column_raw ? j+i1*N2 : N*j+i1;
        (*thptr)[index] = column_raw ? ptr[j+i*N2] : ptr[j*N1+i];
      }
      i1++;
    }
  }
  return 0;
}

int project_table(double *ptr, int N1, int N2, double **thptr, int &N, int arg_num, const int *arg,
int carg, double val, bool column_raw){
  double v = -1e10;
  double coef=0;
  int ind=0;
  for(int i=0;i<N1;i++){
    double v1 = column_raw ? ptr[arg[carg]+i*N2] : ptr[arg[carg]*N1+i];
    if(v1>=val){
      coef = (v1-val)/(v1-v);
      break;
    }
    if(v1!=v){
      v=v1;
      ind++;
    }
  }
  int thin[3]={1,1,1};
  thin[carg]=int(1e6);
  int start[3]={0,0,0};
  start[carg]-=ind;
  return thin_table(ptr,N1,N2,thptr,N,arg_num,arg,thin,column_raw,start);
}
*/

template int read_table(const char *fname, double **ptr, int &N1, int &N2, int column_raw, int skip_first);
template int write_table(const char *fname, double *ptr, int N1, int N2, int column_raw, const char *first, int gnu_space);
template double get_table_value(double arg, double *ptr, int N1, int N2, int j_arg, int j_val, int column_raw, int ord);
template int table_grid(double *ptr, int N1, int N2, int arg_num, const int *arg, double *p1, double *p2, int *sz, int *ord, int column_raw);

template int read_table(const char *fname, float **ptr, int &N1, int &N2, int column_raw, int skip_first);
template int write_table(const char *fname, float *ptr, int N1, int N2, int column_raw, const char *first, int gnu_space);
template float get_table_value(float arg, float *ptr, int N1, int N2, int j_arg, int j_val, int column_raw, int ord);
template int table_grid(float *ptr, int N1, int N2, int arg_num, const int *arg, float *p1, float *p2, int *sz, int *ord, int column_raw);
