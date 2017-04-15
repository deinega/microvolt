#ifndef _TABLE_FUNCTION_H
#define _TABLE_FUNCTION_H

#include <string.h>
#include <stdio.h>

enum FTYPES {XYDATA,YDATA};
enum BTYPES {OUT_ZERO,OUT_LAST};

typedef double realtype;
#define REALFRM "%lf" 

FILE *Err_fopen(const char *file,const char *mode);
int fgetline(FILE *f,char *str,int len_tresh);

// returns the coefficient d
// fav = f(-)*d+f(+)*(1-d)

realtype DxScale(realtype dx,realtype x,int &k);

// compare function to sort realtype array indexed by integer
// order list by qsort
// if CmpAction==1 sorts in ascending order
// if CmpAction==-1 sorts in descending order
extern realtype *CmpArray;
extern int CmpAction;
int CompIndexed(const void *c1,const void *c2);

// Sets array in order according to 'order' list
void SetInOrder(unsigned *order,unsigned n, void *array, size_t size);

// reads one- or two- column file into the arrays
// allocates memory in yy and (if needed) xx.
// Returns the number of entries
int ReadTable(const char *file,realtype* &xx,realtype* &yy,int col1,int col2);

int ReadTableDirect(FILE *f,realtype* &xx,realtype* &yy,int col1,int col2,int ws);

class TableFunction{
public:
//private:
 void sort();
 static TableFunction* Fcur;
//public:
//protected:
 char alloc;
 realtype *xx;
 realtype *yy;
public:
 enum FTYPES ftype;
  // public:
 friend realtype TabFunc(realtype x);
 friend void SetFunc(TableFunction *f);

 int n;
 enum BTYPES btype;
 realtype xstart,xend;
 TableFunction(){ xx=yy=NULL; alloc=0;n=0; }
 TableFunction(realtype x1, realtype x2, int dim);
 TableFunction(int dim,realtype *arr);
 void reset(int dim,realtype *arrx,realtype *arry,int alloc_=0);
 TableFunction(int dim,realtype *arrx,realtype *arry,int alloc_=0);
 void reset(const char *file,int col1=-1,int col2=-1);
 TableFunction(const char *file,int col1=-1,int col2=-1);
 TableFunction(FILE *f, int col1, int col2);
// void write(char *file,realtype yscale=1.,TableMapFunc mfunc=NULL);
 virtual ~TableFunction(){
  if((alloc&0x01))delete [] yy;
  if((alloc&0x02))delete [] xx;
  xx=yy=NULL;
  n=0;
 }

 // trapecial functions
 virtual realtype integral();
 // returns the first x where integral reaches val
 // if not, returns xend
 // works ONLY for positively defined functions!!!
 virtual realtype integral_reaches(realtype val);

 // integrates current function, result is intergral(x)
 // scales the value with the given scale factor
 virtual TableFunction& integrate(realtype scale=1.);

 /* avoid using this function */
 virtual TableFunction& operator=(const TableFunction& F){
   if((alloc&0x01))delete[] yy;
   if((alloc&0x02))delete[] xx;
   //eprintf("ehe\n");
   if(F.alloc){
     printf("TableFunction: equating allocated item!\n");
   }
   memcpy((void*)this,(void*)&F,sizeof(TableFunction));
   alloc=0;
   return *this;
 }

 TableFunction(TableFunction &F){
   //eprintf("oo\n");
   if(F.alloc){
     printf("TableFunction: copying allocated item!\n");
   }
   memcpy((void*)this,(void*)&F,sizeof(TableFunction));
   alloc=0;
   F.alloc=0;
 }


  // copy creation
 virtual TableFunction& operator<<(const TableFunction &F);



 virtual TableFunction& operator*=(TableFunction &F);
 virtual TableFunction& operator*=(realtype f(realtype));
 virtual TableFunction& operator*=(realtype c);
 virtual TableFunction& operator/=(TableFunction &F);
 virtual TableFunction& operator+=(TableFunction &F);
 virtual TableFunction& operator+=(realtype f(realtype));
 virtual TableFunction& operator+=(realtype c);
 virtual TableFunction& operator-=(TableFunction &F);
 virtual TableFunction& operator-=(realtype f(realtype));
 virtual TableFunction& operator-=(realtype c);

 virtual TableFunction& operation(TableFunction &F,realtype fop(realtype,realtype));  

 virtual realtype *arg_table(){
   return xx;
 }
 virtual realtype *func_table(){
   return yy;
 }

 realtype x(int i){
  if(ftype==YDATA)return xstart+(xend-xstart)*i/(realtype)(n-1);
  return xx[i];
 }
 realtype &y(int i){ return yy[i];}
 realtype xscale(realtype x1,realtype x2);
 int   nsteps,stpk,stpi,insert_stat;
 realtype  stpx,yold;
 void  insert_begin(int steps);
 realtype insert_next(realtype y);
 virtual realtype operator()(realtype x)const;

 // for compatibility with emPhotons FIX THIS !!!!!!
 virtual realtype operator()(realtype x, realtype dx )const{
   return (*this)(x);
 }


 virtual realtype Dx(realtype x);
 virtual realtype Dy(realtype x); 
};

#endif
