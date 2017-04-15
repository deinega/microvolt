/*****************************************************************************
 *
 *   Kinetic Technologies Ltd. 2006       
 *
 *   Authors	: Ilya Valuev, Alexei Deinega,  KINTECH, Moscow, Russia
 *
 *   Project	: FDTD-II
 *
 *   $Revision: 1.8 $
 *   $Date: 2013/12/05 06:33:48 $
 *   @(#) $Header: /home/dev/Photon/ivutils/src/implmpi.cpp,v 1.8 2013/12/05 06:33:48 lesha Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/src/implmpi.cpp,v $
$Revision: 1.8 $
$Author: lesha $
$Date: 2013/12/05 06:33:48 $
*/
/******************************************************************************
 * $Log: implmpi.cpp,v $
 * Revision 1.8  2013/12/05 06:33:48  lesha
 * *** empty log message ***
 *
 * Revision 1.7  2013/12/02 22:56:44  lesha
 * *** empty log message ***
 *
 * Revision 1.6  2013/12/02 22:30:37  lesha
 * mpi_gather is included
 *
 * Revision 1.5  2013/11/24 19:02:48  lesha
 * timers are made as double (not vec_type)
 *
 * Revision 1.4  2013/10/30 14:05:59  valuev
 * added in place reduce
 *
 * Revision 1.3  2013/08/05 16:19:08  lesha
 * allreduce_vec_type_maxloc is added
 *
 * Revision 1.2  2012/12/17 12:48:05  chorkov
 * work around: build on lomonosov
 *
 * Revision 1.1  2012/12/06 02:24:05  lesha
 * *** empty log message ***
 *
 * Revision 1.6  2012/10/17 06:50:09  lesha
 * documentation
 *
 * Revision 1.5  2012/10/01 15:23:28  lesha
 * *** empty log message ***
 *
 * Revision 1.4  2012/03/29 11:35:09  valuev
 * universal scatter_data() for emSourceWave
 *
 * Revision 1.3  2012/03/26 17:03:09  valuev
 * added emWave::scatter_data  (untested)
 *
 * Revision 1.2  2012/02/22 04:21:23  lesha
 * MPI_VEC_TYPE is used
 *
 * Revision 1.1  2012/02/15 19:49:42  lesha
 * moved to ivutils
 *
 * Revision 1.6  2012/02/15 19:43:57  lesha
 * emMPIComm: emdefs is excluded
 *
 * Revision 1.5  2010/10/23 10:55:08  lesha
 * *** empty log message ***
 *
 * Revision 1.4  2010/06/08 22:53:34  lesha
 * make UNIX compilable
 *
 * Revision 1.3  2010/04/19 14:26:47  lesha
 * comments by Yan
 *
 * Revision 1.2  2010/04/07 10:09:17  valuev
 * full MPI separation
 *
 * Revision 1.1  2010/04/07 09:48:12  valuev
 * Separated MPI interface from implementation
 *
*******************************************************************************/
#ifdef USE_MPI // mpi implementation: compiled only with parallel version

#include <mpi.h>
#include <stdio.h>
#include "implmpi.h"
#include "ifmpi.h"

apMPIComm::apMPIComm():myrank(0),np(1),init(0),is_exit(0){
  _impl=new apMPIHandler;
}

apMPIComm::~apMPIComm(){
  if(init)Finalize();
  delete _impl;
}

// функция-обёртка для запуска инициализаии MPI
int apMPIComm::_mpiinit(int argc, char **argv){
  int res=MPI_Init(&argc,&argv);
  impl()->set_all(MPI_COMM_WORLD);
  MPI_Comm_rank(impl()->get_all(),&myrank);
  MPI_Comm_size(impl()->get_all(),&np);
  if(myrank==0){
    printf("Started MPI run with %d processes.\n",np);
  }
  return res;
}

void apMPIComm::_mpifinalize(){
  MPI_Finalize();
}


void apMPIComm::_mpicheck_exit(){
  int res=is_exit;
  if(myrank)
    MPI_Send(&res,1,MPI_INT,0,myrank,impl()->get_all());
  else{
    if(res)
      fprintf(stderr,"fatal error on node with rank %d\n",0);
    for(int i=1;i<np;i++){
      MPI_Status status;
      MPI_Recv(&res,1,MPI_INT,i,i,impl()->get_all(),&status);
      if(res){
        is_exit=res;
        fprintf(stderr,"fatal error on node with rank %d\n",i);
      }
    }
  }
  MPI_Bcast(&is_exit,1,MPI_INT,0,impl()->get_all());
}

int apMPIComm::bcast_int(int *numrec, int sz, int rank){
  return MPI_Bcast(numrec,sz,MPI_INT,rank,impl()->get_all());
}
int apMPIComm::bcast_char(char *numrec, int sz, int rank){
  return MPI_Bcast(numrec,sz,MPI_CHAR,rank,impl()->get_all());
}
int apMPIComm::bcast_vec_type(vec_type *numrec, int sz, int rank){
  return MPI_Bcast(numrec,sz,MPI_VEC_TYPE,rank,impl()->get_all());
}
int apMPIComm::bcast_double(double *numrec, int sz, int rank){
  return MPI_Bcast(numrec,sz,MPI_DOUBLE,rank,impl()->get_all());
}

int apMPIComm::recv_vec_type(void *buff, int sz, int dest, int tag){
  MPI_Status status;
  return MPI_Recv(buff,sz,MPI_VEC_TYPE,dest,tag,impl()->get_all(),&status);
}

int apMPIComm::send_vec_type(void *buff, int sz, int src, int tag){
  return MPI_Send(buff,sz,MPI_VEC_TYPE,src,tag,impl()->get_all());
}
int apMPIComm::barrier(){
  return MPI_Barrier(impl()->get_all());
}

int apMPIComm::allreduce_vec_type_maxloc(void *b1, void *b2, int sz){
  return MPI_Allreduce(b1,b2,sz,MPI_VEC_TYPE_INT,MPI_MAXLOC,impl()->get_all());
}
int apMPIComm::reduce_vec_type_sum(void *b1, void *b2, int sz, int rank){
  return MPI_Reduce(b1,b2,sz,MPI_VEC_TYPE,MPI_SUM,rank,impl()->get_all());
}
/// MPI_IN_PLACE is not compiled at mvs100k
//int apMPIComm::reduce_vec_type_sum_in_place(void *b1,int sz, int rank){
//  bool on_root = (get_myrank() == rank);
//  return MPI_Reduce(on_root ? MPI_IN_PLACE : b1,on_root ? b1 : NULL,sz,MPI_VEC_TYPE,MPI_SUM,rank,impl()->get_all());
//}

double apMPIComm::wtime(){
  return MPI_Wtime();
}

int apMPIComm::dup_all(void **mycomm){
  MPI_Comm *dupcomm=new MPI_Comm;
  *mycomm=(void *)dupcomm;
  return MPI_Comm_dup(impl()->get_all(), dupcomm);
}

void apMPIComm::remove_comm(void *mycomm){
  MPI_Comm *mpicomm=(MPI_Comm *)mycomm;
  if(*mpicomm!=MPI_COMM_NULL)
    MPI_Comm_free(mpicomm);
  delete mpicomm;
}

int apMPIComm::comm_rank(void *mycomm, int *myrank){
  return MPI_Comm_rank(*(MPI_Comm *)mycomm,myrank);
}

int apMPIComm::create_type(size_t sz, void **typeptr){
  MPI_Datatype *dtype= new MPI_Datatype;
  *typeptr=(void *)dtype;
  MPI_Type_contiguous( sz, MPI_BYTE, dtype); 
  return MPI_Type_commit(dtype);
}

int apMPIComm::free_type(void *typeptr){
  int res= MPI_Type_free((MPI_Datatype *)typeptr); 
  delete (MPI_Datatype *)typeptr;
  return res;
}

/// Returns sizeof(MPI_Request) to span the requests
size_t apMPIComm::create_nblck_array(int size, void **reqptr, void **statptr){
  if(size){
    *reqptr = (void *)new MPI_Request[size];
    *statptr = (void *)new MPI_Status[size];
  }
  return sizeof(MPI_Request);
}

int apMPIComm::startall(int nsend, void *srequests){
  return MPI_Startall(nsend,(MPI_Request *)srequests);
}

int apMPIComm::waitall(int nsend, void *srequests, void *sstatuses){
  return MPI_Waitall(nsend,(MPI_Request *)srequests,(MPI_Status *)sstatuses);
}


void apMPIComm::free_nblck_array(int size, void **reqptr, void **statptr){
  MPI_Request *req=(MPI_Request *)*reqptr;
  MPI_Status *stat=(MPI_Status *)*statptr;
  for(int i=0;i<size;i++)
    MPI_Request_free(&req[i]);
  if(reqptr)
    delete [] req;
  if(statptr)
    delete [] stat;
  *reqptr=NULL;
  *statptr=NULL;
}

int apMPIComm::vec_type_send_init(void* outbuf, int outsz, int dest, int rank, void * /*MPI_Comm*/ mycomm, void * /*MPI_Request* */request){
  MPI_Request* mpireq = (MPI_Request *)request;
  return MPI_Send_init(outbuf,outsz,MPI_VEC_TYPE,dest,rank,*(MPI_Comm *)mycomm,mpireq);
}

int apMPIComm::vec_type_recv_init(void* inbuf, int insz, int dest, int rank, void * /*MPI_Comm*/ mycomm, void * /*MPI_Request* */request){
  MPI_Request* mpireq = (MPI_Request *)request;
  return MPI_Recv_init(inbuf,insz,MPI_VEC_TYPE,dest,rank,*(MPI_Comm *)mycomm,mpireq);
}

int apMPIComm::free_request(void *reqptr){
  return MPI_Request_free((MPI_Request *)reqptr); 
}

int apMPIComm::int_alltoall(void *buffs,int size,void *buffr,int size2,void *mycomm){
  return MPI_Alltoall(buffs,size,MPI_INT,buffr,size2, MPI_INT,*(MPI_Comm *)mycomm);
}

int apMPIComm::alltoallv(void* sbuf, int *svecsz, int *sdispl, void * /*MPI_Datatype*/tptr1, void* rbuf, int *rvecsz, int *rdispl, void * /*MPI_Datatype*/tptr2,void * /*MPI_Comm*/mycomm){
  return MPI_Alltoallv(sbuf,svecsz,sdispl,*(MPI_Datatype *)tptr1, rbuf, rvecsz, rdispl,*(MPI_Datatype *)tptr2,*(MPI_Comm *)mycomm);
}

int apMPIComm::gather_vectype(void *sbuf, int sz, void *rbuf, int rsz, int root){
  return MPI_Gather(sbuf,sz,MPI_VEC_TYPE,rbuf,rsz,MPI_VEC_TYPE,root,impl()->get_all());
}

int apMPIComm::scatterv_vectype(void *sendbuf, int *sendcnts, int *displs, void *recvbuf, int recvcnt, int root){
  return MPI_Scatterv(sendbuf,sendcnts,displs,MPI_VEC_TYPE,recvbuf,recvcnt,MPI_VEC_TYPE,root,impl()->get_all());
  //return MPI_Scatter(sendbuf,90000 /*sendcnts[0]*/,MPI_VEC_TYPE,recvbuf,90000 /*recvcnt*/,MPI_VEC_TYPE,root,impl()->get_all());
}

#endif
