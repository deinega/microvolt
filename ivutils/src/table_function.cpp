#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "table_function.h"
#include "logexc.h"

int press_key_=0;

void serror(const char *format,...){
 va_list vl;
 va_start(vl,format);
 vprintf(format,vl);
 if(press_key_){
  printf("Press ENTER to exit...\n");
  getchar();
 }
 exit(13);
}

void (*fatal_error)(const char *format, ...)= serror;

FILE* Err_fopen(const char *file,const char *mode){
 char str[30];
 FILE *f;
 f=fopen(file,mode);
 if(!f){
  sprintf(str,"Can't open file %s.\n",file);
  fatal_error(str);
 }
 return f;
}

int fgetline(FILE *f,char *str,int len_tresh=0){
 char ch;
 int i=0;
 ch=(char)fgetc(f);
 while(!feof(f)){
  if(ch=='\n')break;
  str[i++]=ch;
  if(i==len_tresh)break;
  ch=(char)fgetc(f);
 }
 str[i]=0;
 return i;
}

// returns the coefficient d
// fav = f(-)*d+f(+)*(1-d)

realtype DxScale(realtype dx,realtype x,int &k){
 k=(int)(x/dx);
 realtype r=x-dx*k;
 if(r<0)r=-r, k--;
 return r/dx;
}

realtype *CmpArray;
int CmpAction=1;

int CompIndexed(const void *c1,const void *c2){
 unsigned i1,i2;
 i1=*(int*)c1;
 i2=*(int*)c2;
 if(CmpArray[i1]<CmpArray[i2])return -CmpAction;
 if(CmpArray[i1]>CmpArray[i2])return  CmpAction;
 return 0;
}

void SetInOrder(unsigned *order,unsigned n, void *array, size_t size){
 unsigned i,k;
 void *tmp=malloc(size);
 if(!tmp)fatal_error("SetInOrder: Memory allocation error.\n");
 char *a1=(char *)array;

 for(i=0;i<n;i++){
  if(order[i]!=i){
   k=order[i];
   while(k<i)k=order[k];
   if(k!=i){
    memcpy(tmp,a1+size*i,size); // swapping
    memcpy(a1+size*i,a1+size*k,size);
    memcpy(a1+size*k,tmp,size);
   }
  }
 }
 free(tmp);
}

// reads the table from a file
// if whole_scan==1 reads only till the first failure
// col1 -- x column number
// col2 -- y column number
// if col1<=0, reads only y, x is assumed to be [0,1.] range scale
// if col2<=0 --> autodetect, first two file columns are used

int ReadTableDirect(FILE *InFile,realtype* &xx,realtype* &yy,int col1, int col2, int whole_scan=1){
 
 char str[1000];
 realtype x,y;
 int Count=0,InType=2,cur1,cur2,i,read_failed;
 char frm1[250]="", frm2[250]="";
 
 if(col2<=0){
   InType=0;
   col1=1;
   col2=2;
 }
 else if(col1<=0){
   col1=col2;
   col2=-1;
   InType=1;
 }
 
   
 for(i=0;i<col1-1;i++)strcat(frm1,"%*f ");
 strcat(frm1,"" REALFRM "");
 for(i=0;i<col2-1;i++)strcat(frm2,"%*f ");
 if(col2>=0)strcat(frm2,"" REALFRM "");

 //printf("frm: (%s, %s)\n",frm1,frm2);

 long initpos=ftell(InFile);
 i=0;
 do{
  fgetline(InFile,str,1000);

  cur1=sscanf(str,frm1,&x);
  if(cur1==EOF)cur1=0;
  cur2=sscanf(str,frm2,&y);
  if(cur2==EOF)cur2=0;
  //printf("C: %s <- %d %d\n",str,cur1,cur2);

  if(cur1 && InType<=(cur1+cur2)){
   if(!InType)InType=cur1+cur2;

   if(i){ //writing
    if(InType==1)yy[Count]=x;
    else{      
      xx[Count]=x;
      yy[Count]=y;
     
    }
   }
   Count++; //counting
   read_failed=0;
  }
  else read_failed=1;
  if(feof(InFile) || ((!whole_scan) && read_failed)){
   if(Count==0 || i!=0)break;

   if(yy) delete [] yy;
   yy=new realtype[Count];
   if(!yy)fatal_error("ReadTable: memory allocation error.\n");
   if(InType==2){
     if(xx) delete [] xx;
     xx=new realtype[Count];
     if(!xx)fatal_error("ReadTable: memory allocation error.\n");
   }
   i++;
   fseek(InFile,initpos,SEEK_SET);
   Count=0;
  }
 }while(i<2);

 return Count;
}

int ReadTable(const char *file,realtype* &xx,realtype* &yy, int col1=-1, int col2=-1){
 FILE *InFile;
 InFile=Err_fopen(file,"rt");
 int Count=ReadTableDirect(InFile,xx,yy,col1,col2,1);
 fclose(InFile);
 return Count;
}

TableFunction::TableFunction(int dim,realtype *arr=NULL){
 xx=NULL;
 yy=arr;
 n=dim;
 if(n<=0)fatal_error("TableFunction: Invalid array dimension: %d\n",n);
 alloc=0;
 if(yy==NULL){
   yy=new realtype[dim];
   if(!yy)fatal_error("TableFunction: memory allocation error.\n");
   alloc=0x1;
   for(int i=0;i<dim;i++)yy[i]=0.;
 }
 xstart=0.;
 xend=1.;
 ftype=YDATA;
 btype=OUT_ZERO;
}

TableFunction::TableFunction(realtype x1, realtype x2,int dim){
  TableFunction a(dim, NULL);
  *this=a;
  a.alloc=0;
  xstart=x1;
  xend=x2;
}


void TableFunction::reset(int dim,realtype *arrx,realtype *arry,int alloc_){
 n=dim;
 if(n<=0)fatal_error("TableFunction: Invalid array dimension: %d\n",n);
 if(alloc_==0){
   xx=arrx;
   yy=arry;
   alloc=0;
 }
 else{
   xx = new double[n];
   yy = new double[n];
   alloc=0x1 | 0x2;
   for(int i=0;i<n;i++){
     xx[i]=arrx[i];
     yy[i]=arry[i];
   }
 }
 sort();
 xstart=xx[0];
 xend=xx[n-1];
 ftype=XYDATA;
 btype=OUT_ZERO;
}

TableFunction::TableFunction(int dim,realtype *arrx,realtype *arry,int alloc_){
  reset(dim,arrx,arry,alloc_);
}

void TableFunction::reset(const char *file,int col1,int col2){
 if(alloc&0x1 && xx)
   delete [] xx;
 if(alloc&0x2 && yy)
   delete [] yy;

 xx=yy=NULL;
 n=ReadTable(file,xx,yy,col1,col2);
 if(n<=0)fatal_error("Can not initialize TableFunction from file %s\n",file);
 if(xx==NULL){
  ftype=YDATA;
  alloc=0x1;
 }
 else{
  ftype=XYDATA;
  alloc=0x3;
 }
 sort();
 xstart=xx[0];
 xend=xx[n-1];
 btype=OUT_ZERO;
}

TableFunction::TableFunction(const char *file,int col1,int col2):alloc(0){
  reset(file,col1,col2);
}

TableFunction::TableFunction(FILE *f, int col1=-1, int col2=-1){
 xx=yy=NULL;
 n=ReadTableDirect(f,xx,yy,col1,col2,0);
 if(n<=0)fatal_error("Can not initialize TableFunction from file\n");
 if(xx==NULL){
  ftype=YDATA;
  alloc=0x1;
 }
 else{
  ftype=XYDATA;
  alloc=0x3;
 }
 sort();
 xstart=xx[0];
 xend=xx[n-1];
 btype=OUT_LAST;
}


TableFunction& TableFunction::operator<<(const TableFunction &F){
  if((alloc&0x01))delete[] yy;
  if((alloc&0x02))delete[] xx;
  
  memcpy((void*)this,(void*)&F,sizeof(TableFunction));
  
  if(ftype==XYDATA)alloc=0x03;
  else alloc=0x01;
  
  int i;
  if(alloc&0x01){
    yy = new realtype[n];
    if(!yy)serror("TableFunction: MAE.\n");
    for(i=0;i<n;i++)yy[i]=F.yy[i];
  }
  if(alloc&0x02){
     xx = new realtype[n];
     if(!xx)serror("TableFunction: MAE.\n");
     for(i=0;i<n;i++)xx[i]=F.xx[i];
  }
  return *this;
}

/*
void TableFunction::write(char *file, realtype yscale, int (*mapfunc)(realtype *x, realtype *y)){
 int i;
 realtype x=0, y;
 FILE *f=Err_fopen(file,"wt");
 for(i=0;i<n;i++){
   switch(ftype){
   case XYDATA: 
     x=xx[i];
     break;
   case YDATA:
     x=xstart;
     if(n>1)x+=(-xstart+xend)*i/(n-1);
   }
   y=yy[i]*yscale;
   if(mapfunc){
     if(mapfunc(&x,&y)<=0)continue;
   }
   fprintf(f,"" REALFRM " " REALFRM "\n",x,y);
 }
 fclose(f);
}
*/   

realtype TableFunction::operator()(realtype x)const{
 
  //printf("" REALFRM " " REALFRM " " REALFRM "\n",x,xend,xstart);
 if(x<xstart){
  switch(btype){
   case OUT_ZERO: return 0.;
   case OUT_LAST: return yy[0];
  }
 }
 if(x>xend){
  switch(btype){
   case OUT_ZERO: return 0.;
   case OUT_LAST: return yy[n-1];
  }
 }
 

 realtype d=0;
 int i;

 if(xstart!=xend)d=(x-xstart)*(n-1)/(xend-xstart);
 i=(int)d;
 if(i>=n-1)return yy[n-1];

 //return yy[0]+x;

 if(ftype==YDATA){
  d=d-(realtype)i;
  d=*(yy+i)+(*(yy+i+1)-*(yy+i))*d;
  return d;
 }
 
 i++;
 while(xx[i]<x)i++;
 while(xx[i-1]>x)i--;

 d=xx[i]-xx[i-1];
 if(d!=0.)d=(x-xx[i-1])/d;
 return yy[i-1]+(yy[i]-yy[i-1])*d;
}


TableFunction& TableFunction::operator*=(TableFunction &F){
  int i;
  realtype dx=-(xstart-xend)/(n-1);
  realtype x;

  for(i=0;i<n;i++){
    if(ftype==YDATA)x=xstart+dx*i;
    else x=xx[i];

    //if(x>=F.xstart && x<= F.xend){
      yy[i]*=F(x);
    //}
  }
  return *this;
}


TableFunction& TableFunction::operator/=(TableFunction &F){
  int i;
  realtype dx=-(xstart-xend)/(n-1);
  realtype x;

  for(i=0;i<n;i++){
    if(ftype==YDATA)x=xstart+dx*i;
    else x=xx[i];

    //if(x>=F.xstart && x<= F.xend){
      realtype del=F(x);
      if(fabs(del)<1e-32){
        if(fabs(yy[i])>1e-32)yy[i]=1e32;
        else yy[i]=0;
      }
      else yy[i]/=del;
    //}
  }
  return *this;
}

TableFunction& TableFunction::operator*=(realtype F(realtype)){
  int i;
  realtype dx=-(xstart-xend)/(n-1);
  realtype x;

  for(i=0;i<n;i++){
    if(ftype==YDATA)x=xstart+dx*i;
    else x=xx[i];
    yy[i]*=F(x);
  }
  return *this;
}

TableFunction& TableFunction::operator+=(realtype F(realtype)){
  int i;
  realtype dx=-(xstart-xend)/(n-1);
  realtype x;

  for(i=0;i<n;i++){
    if(ftype==YDATA)x=xstart+dx*i;
    else x=xx[i];
    yy[i]+=F(x);
  }
  return *this;
}



TableFunction& TableFunction::operator+=(TableFunction &F){
  int i;
  realtype dx=-(xstart-xend)/(n-1);
  realtype x;

  for(i=0;i<n;i++){
    if(ftype==YDATA)x=xstart+dx*i;
    else x=xx[i];
    //if(x>=F.xstart && x<= F.xend){
      yy[i]+=F(x);
    //}
  }
  return *this;
}



TableFunction& TableFunction::operator+=(realtype c){
 
  int i;
  for(i=0;i<n;i++){
    yy[i]+=c;
  
  }
  return *this;
}

TableFunction& TableFunction::operator-=(realtype F(realtype)){
  int i;
  realtype dx=-(xstart-xend)/(n-1);
  realtype x;

  for(i=0;i<n;i++){
    if(ftype==YDATA)x=xstart+dx*i;
    else x=xx[i];
    yy[i]-=F(x);
  }
  return *this;
}



TableFunction& TableFunction::operator-=(TableFunction &F){
  int i;
  realtype dx=-(xstart-xend)/(n-1);
  realtype x;

  for(i=0;i<n;i++){
    if(ftype==YDATA)x=xstart+dx*i;
    else x=xx[i];
    //if(x>=F.xstart && x<= F.xend){
      yy[i]-=F(x);
    //}
  }
  return *this;
}



TableFunction& TableFunction::operation(TableFunction &F,
					realtype fop(realtype,realtype)){
  int i;
  realtype dx=-(xstart-xend)/(n-1);
  realtype x;

  for(i=0;i<n;i++){
    if(ftype==YDATA)x=xstart+dx*i;
    else x=xx[i];
    //if(x>=F.xstart && x<= F.xend){
      yy[i]=fop((*this)(x),F(x));
    //}
  }
  return *this;
}





TableFunction& TableFunction::operator-=(realtype c){

  int i;
  for(i=0;i<n;i++){
    yy[i]-=c;

  }
  return *this;
}



TableFunction& TableFunction::operator*=(realtype c){

  int i;
  for(i=0;i<n;i++){
    yy[i]*=c;

  }
  return *this;
}


realtype TableFunction::Dx(realtype x){

 if(ftype==YDATA){
   if(n<=1)return 0.;
   else return (xend-xstart)/(n-1);
 }

 if(x<xstart || x>=xend || n<2){
   return 0.;
 }

 realtype d=0;
 int i;
 if(xstart!=xend)d=(x-xstart)*(n-1)/(xend-xstart);
 i=(int)d;

 i++;
 while(xx[i]<x)i++;
 while(xx[i-1]>x)i--;
 d=xx[i]-xx[i-1];
 if(d!=0.)d=(x-xx[i-1])/d;
 i--;

 realtype dx1, dx2;

 if(i<=0)dx1=xx[1]-xx[0];
 else if(i>=n-1)dx1=xx[n-1]-xx[n-2];
 else dx1=(xx[i+1]-xx[i-1])/2;

 if(i>=n-2)dx2=xx[i+1]-xx[i];
 else dx2=(xx[i+2]-xx[i])/2;

 d=dx1+(dx2-dx1)*d;
 return d;
}

realtype TableFunction::Dy(realtype x){
 if(x<xstart || x>=xend || n<2){
   return 0.;
 }

 realtype d=0;
 int i;

 if(xstart!=xend)d=(x-xstart)*(n-1)/(xend-xstart);
 i=(int)d;

 if(ftype!=YDATA){
   i++;
   while(xx[i]<x)i++;
   while(xx[i-1]>x)i--;
   d=xx[i]-xx[i-1];
   if(d!=0.)d=(x-xx[i-1])/d;
   i--;
 }
 else{
   d=d-(realtype)i;
 }

 realtype dy1, dy2;

 if(i<=0)dy1=yy[1]-yy[0];
 else if(i>=n-1)dy1=yy[n-1]-yy[n-2];
 else dy1=(yy[i+1]-yy[i-1])/2;

 if(i>=n-2)dy2=yy[i+1]-yy[i];
 else dy2=(yy[i+2]-yy[i])/2;

 d=dy1+(dy2-dy1)*d;
 return d;
}




realtype TableFunction::xscale(realtype x1,realtype x2){
 realtype k/*,x0*/;
 if(x1>x2){k=x1;x1=x2;x2=k;}
 k=(x2-x1)/(xend-xstart);

 int i;
 switch(ftype){
  case XYDATA:
   for(i=0;i<n;i++)xx[i]=x1+k*(xx[i]-xstart);
  case YDATA:
   xstart=x1;
   xend=x2;
 }
 return k;
}

void TableFunction::insert_begin(int steps){
  nsteps=steps;
  stpk=0;
  stpx=0.;
  stpi=0;
  insert_stat=0;
  yold=0.;
}


realtype TableFunction::insert_next(realtype y){
   int dumm;
   realtype d, yr=yold;

   if(stpi==0){ // first step
     yy[stpi++]=y;
     stpx=((realtype)(nsteps-1))*stpi/(n-1); // where the next insert place is
   }
   else{ // scaling
     while((realtype)stpk>=stpx && stpi<n){
       d=DxScale(1.,stpx,dumm);
       yr=yold*(1-d)+y*d;

       yy[stpi]=yr;

       stpi++;
       stpx=((realtype)(nsteps-1))*stpi/(n-1);// where the next insert place is

     }
   }
   if(stpi>=n){
     insert_stat=0;
   }
   else{
     yold=y;
     stpk++; // where the next data comes to
     insert_stat=1;
   }
   return yr;
}


TableFunction * TableFunction::Fcur=NULL;
void SetFunc(TableFunction *f){ TableFunction::Fcur=f;}

realtype TabFunc(realtype x){
 if(TableFunction::Fcur)return (TableFunction::Fcur)->operator()(x);
 else return 0.;
}

realtype TableFunction::integral(){
  if(n<2)return 0.;
  int i;
  realtype sum=0;
  if(ftype==YDATA){
    realtype scale=(xend-xstart)/(n-1);
    for(i=0;i<n-1;i++)sum+=0.5f*(yy[i]+yy[i+1]);
    sum*=scale;
  }
  else{
    for(i=0;i<n-1;i++)sum+=0.5f*(yy[i]+yy[i+1])*(xx[i+1]-xx[i]);
  }
  return sum;
}

TableFunction& TableFunction::integrate(realtype scale){
  if(n<2)return *this;
  int i;
  realtype sum=0, left=0;
  if(ftype==YDATA){
    scale*=(xend-xstart)/(n-1);
    for(i=0;i<n;i++){
      sum+=0.75*yy[i]; // 3/4
      if(i>0)sum+=left/8.;
      else sum+=yy[i]/8.;
      if(i<n-1)sum+=yy[i+1]/8.;
      else sum+=yy[i]/8.;
      left=yy[i];
      yy[i]=sum*scale;
    }
  }
  else{
    double rx=1, lx=1;
    for(i=0;i<n;i++){
      if(i<n-1)rx=0.5*(xx[i+1]-xx[i]);
      if(i==0)lx=rx;
      sum+=0.75*yy[i]*(rx+lx); // 3/4
      if(i>0)sum+=left*lx/8.;
      else sum+=yy[i]*lx/8.;
      if(i<n-1)sum+=yy[i+1]*rx/8.;
      else sum+=yy[i]*rx/8.;
      left=yy[i];
      yy[i]=sum*scale;
      lx=rx;
    }
  }
  return *this;
}

realtype TableFunction::integral_reaches(realtype val){
  int i;
  realtype sum=0;
  if(ftype==YDATA){
    val/=(xend-xstart);
    for(i=0;i<n-1;i++){
      sum+=0.5f*(yy[i]+yy[i+1]);
      if(sum>=val)break;
    }
    if(i>=n-1)return xend;

    realtype sum0=sum-0.5f*(yy[i]+yy[i+1]);
    realtype dy=yy[i+1]-yy[i];

    realtype det=yy[i]*yy[i]+2*dy*(val-sum0);
    if(fabs(dy)<1e-32 || det <0)return xstart+(xend-xstart)*i/(n-1);
    realtype dx=(realtype)(sqrt(det)-yy[i])/dy;
    return xstart+(xend-xstart)*i/(n-1)+dx;
  }
  else{
    for(i=0;i<n-1;i++){
      sum+=0.5f*(yy[i]+yy[i+1])*(xx[i+1]-xx[i]);
      if(sum>=val)break;
    }
    if(i>=n-1)return xend;
    realtype sum0=sum-0.5f*(yy[i]+yy[i+1]);
    realtype dy=yy[i+1]-yy[i];

    realtype det=yy[i]*yy[i]+2*dy*(val-sum0);
    if(fabs(dy)<1e-32 || det <0)return xx[i];
    realtype dx=(realtype)(sqrt(det)-yy[i])/dy;
    return xx[i]+dx;
  }
}

void TableFunction::sort(){
 CmpArray=xx;
 unsigned *order=new unsigned[n];
 if(!order)fatal_error("TableFunction:sort: memory allocation error.\n");
 unsigned i;
 for(i=0;i<(unsigned)n;i++)order[i]=i;
 qsort(order,n,sizeof(int),CompIndexed);
 SetInOrder(order,n,xx,sizeof(realtype));
 SetInOrder(order,n,yy,sizeof(realtype));
 for(i=1;i<(unsigned)n;i++){
  if(fabs(xx[i-1]-xx[i])<1e-32){
   memmove(xx+i,xx+i+1,sizeof(int)*(n-(i+1)));
   n--;
  }
 }
 delete [] order;
}

