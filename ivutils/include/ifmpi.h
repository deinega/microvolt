/*****************************************************************************
 *
 *   Kinetic Technologies Ltd. 2006       
 *
 *   Authors	: Ilya Valuev, Alexei Deinega,  KINTECH, Moscow, Russia
 *
 *   Project	: FDTD-II
 *
 *   $Revision: 1.8 $
 *   $Date: 2013/12/05 06:34:03 $
 *   @(#) $Header: /home/dev/Photon/ivutils/include/ifmpi.h,v 1.8 2013/12/05 06:34:03 lesha Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/include/ifmpi.h,v $
$Revision: 1.8 $
$Author: lesha $
$Date: 2013/12/05 06:34:03 $
*/
/******************************************************************************
 * $Log: ifmpi.h,v $
 * Revision 1.8  2013/12/05 06:34:03  lesha
 * *** empty log message ***
 *
 * Revision 1.7  2013/12/04 22:24:21  lesha
 * *** empty log message ***
 *
 * Revision 1.6  2013/12/03 22:47:10  lesha
 * time_MPI_transfers is included
 *
 * Revision 1.5  2013/12/02 22:30:37  lesha
 * mpi_gather is included
 *
 * Revision 1.4  2013/11/24 19:02:48  lesha
 * timers are made as double (not vec_type)
 *
 * Revision 1.3  2013/10/30 14:05:59  valuev
 * added in place reduce
 *
 * Revision 1.2  2013/08/05 16:19:08  lesha
 * allreduce_vec_type_maxloc is added
 *
 * Revision 1.1  2012/12/06 02:24:04  lesha
 * *** empty log message ***
 *
 * Revision 1.8  2012/10/17 06:50:09  lesha
 * documentation
 *
 * Revision 1.7  2012/09/28 21:47:30  lesha
 * *** empty log message ***
 *
 * Revision 1.6  2012/09/27 19:11:51  lesha
 * documentation
 *
 * Revision 1.5  2012/05/11 15:43:27  valuev
 * fixed box flux calculation, reverted to russian documentation
 *
 * Revision 1.4  2012/03/26 17:03:09  valuev
 * added emWave::scatter_data  (untested)
 *
 * Revision 1.3  2012/03/22 07:31:20  valuev
 * updated increments for iterators, prepared parallel/universal near-to-far (not working)
 *
 * Revision 1.2  2012/03/21 17:00:12  lesha
 * some comments
 *
 * Revision 1.1  2012/02/15 19:50:11  lesha
 * moved to ivutils
 *
 * Revision 1.28  2012/02/15 19:43:57  lesha
 * emMPIComm: emdefs is excluded
 *
 * Revision 1.27  2010/09/07 21:12:07  lesha
 * exception processing is corrected
 *
*******************************************************************************/
#ifndef _IFMPI_H
#define _IFMPI_H

/** @file ifmpi.h
    @brief Classes for mpi interface/ mpi stubs.
**/


#include "vector_3.h"

class apMPIHandler;

///\en interface to work with MPI. 
/// works in serial and MPI (if macro USE_MPI is defined) regimes.
/// In serial regime MPI function are not called.
///\ru класс интерфейса для работы с MPI
class apMPIComm {
  int init; // 1 if Init was called
  int myrank; // current procossor number
  int np; // processors number
  int is_exit; ///<\en 1 if some error happened 

  /// for parallel version apMPIHandler has a MPI communicator, for serial version this is just empty class
  apMPIHandler *_impl;
public:
  ///\en Accessor to MPI implementation.
  ///    The returned pointer should be used to acces MPI-specifiec functions
  ///    The apMPIHandler class definition is available only when including implmpi.h
  ///    that directly includes mpi.h
  apMPIHandler *impl(){
    return _impl;
  }

  apMPIComm();

  ///\en calls MPI_Init and initialize fields of the class
  ///\ru вызывает обертку для MPI инициализации
  int Init(int argc, char **argv);
  // calls MPI_Finalize
  int Finalize();

  // number of current processors
  int get_myrank();
  // processors number
  int get_nproc();

#ifdef USE_MPI
  ///\en interface for MPI functions, their realization can be found in implmpi.cpp
  ///\ru все эти функции находятся в файле implmpi.cpp
  int bcast_int(int *numrec, int sz, int rank);
  int bcast_char(char *numrec, int sz, int rank);
  int bcast_vec_type(vec_type *numrec, int sz, int rank);
  int bcast_double(double *numrec, int sz, int rank);
  int recv_vec_type(void *buff, int sz, int src, int tag);
  int send_vec_type(void *buff, int sz, int dest, int tag);
  int barrier();
  int allreduce_vec_type_maxloc(void *b1, void *b2, int sz);
  int reduce_vec_type_sum(void *b1, void *b2, int sz, int rank);
  ///\en reduces the sum at the same buffer
//  int reduce_vec_type_sum_in_place(void *b1, int sz, int rank);
  double wtime();
  int dup_all(void **mycomm);
  void remove_comm(void *mycomm);
  int comm_rank(void *mycomm, int *myrank);
  int create_type(size_t sz, void **typeptr);
  int free_type(void *typeptr);
  ///\en Returns sizeof(MPI_Request) to span the requests
  size_t create_nblck_array(int size, void **reqptr, void **statptr);
  int startall(int nsend, void *srequests);
  int waitall(int nsend, void *srequests, void *sstatuses);
  void free_nblck_array(int size, void **reqptr, void **statptr);
  int vec_type_send_init(void* outbuf, int outsz, int dest, int rank, void * /*MPI_Comm*/ mycomm, void * /*MPI_Request* */request);
  int vec_type_recv_init(void* inbuf, int insz, int dest, int rank, void * /*MPI_Comm*/ mycomm, void * /*MPI_Request* */request);
  int free_request(void *reqptr);
  int int_alltoall(void *buffs,int size,void *buffr,int size2,void *mycomm);
  int alltoallv(void* sbuf, int *svecsz, int *sdispl, void * /*MPI_Datatype*/tptr1, void* rbuf, int *rvecsz, int *rdispl, void * /*MPI_Datatype*/tptr2,void * /*MPI_Comm*/mycomm);
  int gather_vectype(void *sbuf, int sz, void *rbuf, int rsz, int root);
  int scatterv_vectype(void *sendbuf, int *sendcnts, int *displs, void *recvbuf, int recvcnt, int root);
protected:
  ///\en calls MPI_Init and initialize fields of the class
  ///\ru функция-обёртка для запуска инициализаии MPI
  int _mpiinit(int argc, char **argv);
  // calls MPI_Finalize
  void _mpifinalize();
  void _mpicheck_exit();
#endif
public:
  ///<\en set is_exit=1
  ///\ru установка флага is_exit единице
  void set_exit();
  ///<\en if is_exit=1 at some processor, set is_exit=1 at all processors. returns is_exit
  ///\ru проверка, зафиксирована ли на каком-либо процессоре ошибка. Если да - выход из программы
  int check_exit();
#ifdef USE_MPI
  ///\en record to the file fname table np x np (np is processors number)
  /// value in (i,j) element corresponds to the time required to transfer from i to j processor array of sz doubles 
  int time_MPI_transfers(const char *fname="time_mpi.d", const int sz=1000000);
#endif
  ~apMPIComm();
};

extern apMPIComm theComm; /// all MPI functions are called using this object

#endif
