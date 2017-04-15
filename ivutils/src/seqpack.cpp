/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2006        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: ivutils
 *
 *   $Revision: 1.4 $
 *   $Date: 2013/09/10 23:08:24 $
 *   @(#) $Header: /home/dev/Photon/ivutils/src/seqpack.cpp,v 1.4 2013/09/10 23:08:24 lesha Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/src/seqpack.cpp,v $
$Revision: 1.4 $
$Author: lesha $
$Date: 2013/09/10 23:08:24 $
*/
/*s****************************************************************************
 * $Log: seqpack.cpp,v $
 * Revision 1.4  2013/09/10 23:08:24  lesha
 * compilation on icc 12.1.3
 *
 * Revision 1.3  2013/09/05 09:25:08  valuev
 * added default instantiation for index_pack<int>
 *
 * Revision 1.2  2013/09/04 15:22:56  valuev
 * restructure for 64-bits (added index_t)
 *
 * Revision 1.1  2012/12/06 02:24:05  lesha
 * *** empty log message ***
 *
 * Revision 1.15  2012/10/16 19:45:42  lesha
 * documentation
 *
 * Revision 1.14  2012/10/05 19:08:32  lesha
 * *** empty log message ***
 *
 * Revision 1.13  2012/06/16 01:06:05  valuev
 * sync with GridMD project
 *
 * Revision 1.7  2010/06/12 17:57:31  valuev
 * some workflow coding
 *
 * Revision 1.11  2009/02/27 00:40:58  lesha
 * 2 calles put_record_end bug is fixed
 *
 * Revision 1.10  2009/01/30 13:54:05  valuev
 * restructured as a library
 *
 * Revision 1.5  2008/04/09 17:31:25  valuev
 * Added bands and charges, modified wequil analyser
 *
 * Revision 1.4  2008/02/27 13:49:47  valuev
 * corrected variable delink from namespace
 *
 * Revision 1.3  2008/02/21 14:02:51  valuev
 * Added parametric methods
 *
 * Revision 1.2  2007/07/09 21:29:07  valuev
 * plasma with wave packets
 *
 * Revision 1.7  2007/02/20 10:26:12  valuev
 * added newlines at end of file
 *
 * Revision 1.6  2007/02/20 10:11:32  valuev
 * Debugging with g++ 2.96
 *
 * Revision 1.5  2007/02/19 23:50:09  valuev
 * New makefile style; got rid of some compiler warnings
 *
 * Revision 1.4  2006/12/03 22:46:22  lesha
 * Error in put_record_end is corrected
 *
 * Revision 1.3  2006/12/02 01:53:43  lesha
 * Error in put_record_end is corrected
 *
 * Revision 1.2  2006/11/24 20:17:31  valuev
 * Added CVS headers
 *
*******************************************************************************/

# include <string.h>
# include "seqpack.hpp"

//template class index_pack<ptrdiff_t>;
//template class index_pack<int>;

template int get_next_record(const vector<int> &ind, int &pos, int &incrcount, int single);
template int put_next_record(vector<int> &ind, int cur, int single);
template int put_record_end(vector<int> &ind,int single);

template int get_next_record(const vector<ptrdiff_t> &ind, ptrdiff_t &pos, ptrdiff_t &incrcount, ptrdiff_t single);
template int put_next_record(vector<ptrdiff_t> &ind, ptrdiff_t cur, ptrdiff_t single);
template int put_record_end(vector<ptrdiff_t> &ind,ptrdiff_t single);

