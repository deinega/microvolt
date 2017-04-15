/*****************************************************************************
 *
 *   Kinetic Technologies Ltd. 2006       
 *
 *   Authors	: Ilya Valuev, Alexei Deinega,  KINTECH, Moscow, Russia
 *
 *   Project	: FDTD-II
 *
 *   $Revision: 1.5 $
 *   $Date: 2013/12/04 22:24:35 $
 *   @(#) $Header: /home/dev/Photon/ivutils/src/ifmpi.cpp,v 1.5 2013/12/04 22:24:35 lesha Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/src/ifmpi.cpp,v $
$Revision: 1.5 $
$Author: lesha $
$Date: 2013/12/04 22:24:35 $
*/
/******************************************************************************
 * $Log: ifmpi.cpp,v $
 * Revision 1.5  2013/12/04 22:24:35  lesha
 * *** empty log message ***
 *
 * Revision 1.4  2013/12/03 22:50:47  lesha
 * *** empty log message ***
 *
 * Revision 1.3  2013/12/03 22:49:33  lesha
 * *** empty log message ***
 *
 * Revision 1.2  2013/12/03 22:47:10  lesha
 * time_MPI_transfers is included
 *
 * Revision 1.1  2012/12/06 02:24:05  lesha
 * *** empty log message ***
 *
 * Revision 1.4  2012/10/17 06:50:09  lesha
 * documentation
 *
 * Revision 1.3  2012/05/11 15:43:27  valuev
 * fixed box flux calculation, reverted to russian documentation
 *
 * Revision 1.2  2012/03/21 17:17:53  lesha
 * documentation
 *
 * Revision 1.1  2012/02/15 19:49:42  lesha
 * moved to ivutils
 *
*******************************************************************************/

#include <stdio.h>
#include "ifmpi.h"

#ifndef USE_MPI

/// for serial version apMPIHandler is just empty class
class apMPIHandler{};

apMPIComm::apMPIComm():myrank(0),np(1),init(0),is_exit(0){
  _impl=new apMPIHandler;
}

apMPIComm::~apMPIComm(){
  if(init)Finalize();
  delete _impl;
}

#endif


int apMPIComm::Init(int argc, char **argv){
  int res=1;
#ifdef USE_MPI
  res=_mpiinit(argc,argv);
#endif
  init=1;
  return res;
}

int apMPIComm::Finalize(){
#ifdef USE_MPI
  _mpifinalize();
#endif
  init=0;
  return 1;
}

int apMPIComm::get_myrank(){
  return myrank;
}

int apMPIComm::get_nproc(){
  return np;
}

void apMPIComm::set_exit(){
  is_exit=1;
}

int apMPIComm::check_exit(){
#ifdef USE_MPI
  _mpicheck_exit();
#endif
  return is_exit;
}

#ifdef USE_MPI

int apMPIComm::time_MPI_transfers(const char *fname, const int sz){

  double *time_vector = new double[np];
  time_vector[myrank]=0;
  double *time_table = myrank ? NULL : new double[np*np];
  double *buff = new double[sz];
  for(int ni=0;ni<np*np;ni++){
    int i=ni/np;
    int j=ni%np;
    if(i==j)
      continue;
    if(myrank==i){
      send_vec_type(buff,sz,j,ni);
    }
    else if(myrank==j){
      double t1=wtime();
      recv_vec_type(buff,sz,i,ni);
      double t2=wtime();
      time_vector[i]=t2-t1;
 //     time_vector[i]=myrank;
    }
    barrier();
  }
  gather_vectype(time_vector,np,time_table,np,0);
  delete []buff;
  delete []time_vector;
  FILE *f=NULL;
  if(myrank==0){
    f=fopen(fname,"wt");
    if(!f){
      delete []time_table;
      set_exit();
    }
  }
  if(check_exit()) return -1;
  if(myrank==0){
    for(int i=0;i<np;i++){
      for(int j=0;j<np;j++){
        fprintf(f,"%g", time_table[i*np+j]);
        fprintf(f,"\t");
      }
      fprintf(f,"\n");
    }
    fclose(f);
    delete []time_table;
  }
  return 1;
}

#endif



apMPIComm theComm;
