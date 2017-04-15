/*****************************************************************************
 *
 *   Kinetic Technologies Ltd. 2010      
 *
 *   Authors	: Ilya Valuev, Alexei Deinega,  KINTECH, Moscow, Russia
 *
 *   Project	: FDTD-II
 *
 *   $Revision: 1.2 $
 *   $Date: 2013/08/05 16:19:08 $
 *   @(#) $Header: /home/dev/Photon/ivutils/include/implmpi.h,v 1.2 2013/08/05 16:19:08 lesha Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/include/implmpi.h,v $
$Revision: 1.2 $
$Author: lesha $
$Date: 2013/08/05 16:19:08 $
*/
/******************************************************************************
 * $Log: implmpi.h,v $
 * Revision 1.2  2013/08/05 16:19:08  lesha
 * allreduce_vec_type_maxloc is added
 *
 * Revision 1.1  2012/12/06 02:24:04  lesha
 * *** empty log message ***
 *
 * Revision 1.4  2012/10/17 06:50:09  lesha
 * documentation
 *
 * Revision 1.3  2012/03/21 17:00:12  lesha
 * some comments
 *
 * Revision 1.2  2012/02/22 04:21:23  lesha
 * MPI_VEC_TYPE is used
 *
 * Revision 1.1  2012/02/15 19:50:11  lesha
 * moved to ivutils
 *
 * Revision 1.3  2012/02/15 19:43:57  lesha
 * emMPIComm: emdefs is excluded
 *
 * Revision 1.2  2010/10/23 10:55:08  lesha
 * *** empty log message ***
 *
 * Revision 1.1  2010/04/07 09:48:12  valuev
 * Separated MPI interface from implementation
 *
*******************************************************************************/
# ifndef _IMPLMPI_H
# define _IMPLMPI_H

#include <mpi.h>
#ifndef SINGLE_PRECISION
#define MPI_VEC_TYPE MPI_DOUBLE
#define MPI_VEC_TYPE_INT MPI_DOUBLE_INT
#else
#define MPI_VEC_TYPE MPI_FLOAT
#define MPI_VEC_TYPE_INT MPI_FLOAT_INT
#endif

/// for parallel version apMPIHandler has a MPI communicator
class apMPIHandler{
  friend class apMPIComm;
  MPI_Comm all;
public:
  void set_all(MPI_Comm all_){
    all=all_;
  }
  MPI_Comm get_all() const{
    return all;
  }
};

# endif
