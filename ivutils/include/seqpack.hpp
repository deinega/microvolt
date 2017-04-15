/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2006        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: ivutils
 *
 *   $Revision: 1.1 $
 *   $Date: 2013/09/04 15:22:56 $
 *   @(#) $Header: /home/dev/Photon/ivutils/include/seqpack.hpp,v 1.1 2013/09/04 15:22:56 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/include/seqpack.hpp,v $
$Revision: 1.1 $
$Author: valuev $
$Date: 2013/09/04 15:22:56 $
*/
/*s****************************************************************************
 * $Log: seqpack.hpp,v $
 * Revision 1.1  2013/09/04 15:22:56  valuev
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
# include "seqpack.h"


/// extracts the next value at a given position 
/// pos must be less than or equal ind.size()-2
/// @returns the value or -1 if end of sequence or the postion is incorrect
/// modifies pos (the next reading position) and incrcount
template<class index_t>
int get_next_record(const vector<index_t> &ind, index_t &pos, index_t &incrcount, index_t single){
  do{
    index_t ic=pos;
    index_t cval0=ind[pos];
    //index_t *pind=(index_t *)(&ind[0]);
    if(cval0>=0){
      ic++;
      index_t cval=ind[ic];
      if(cval<0){
        if(incrcount<-cval){ 
          incrcount++;
        }
        else{  
          if(ind[ic+1]<0)pos=ic+2;
          else pos=ic+1;
          incrcount=0;
          continue;
        }
        if(ind[ic+1]<0){
          cval0+=(incrcount-1)*(-ind[ic+1]-1);
        }
        else{
          cval0+=(incrcount-1)*single;
        }
      }
      else{
        pos=ic;
        incrcount=0;
      }
      return cval0;
    }
    else{
      return -1; // end of record or incorrect position
    }
  }while(1);
  return -1;
}

/// specifies the next value to be recorded
/// only nonnegative values are possible
/// @retuns    0 if the value was not packed (just stored)
///            1 packed with selected increment
///            2 packed with other increment
///            3 packed as continued sequence
///           -1 tried to store after end-of sequence 
///           -2 tried to store negative value
template<class index_t>
int put_next_record(vector<index_t> &ind, index_t cur, index_t single){
  if(cur<0)return -2; //negative, not allowed
  index_t last=(index_t)ind.size()-1;
  if(last<=0){ // this is the first or the second element
    ind.push_back(cur);
  }
  else{
    if(ind[last]<0){ // incremented
      if(ind[last-1]<0){// the increment is not <single>
        if(ind[last-2]<0)return -1; // end of sequence ?
        index_t inc0=-ind[last]-1;
        index_t lval=ind[last-2]+inc0*(-ind[last-1]-1);
        if(cur-lval!=inc0)ind.push_back(cur); // not confroming to the increment
        else{
          --(ind[last-1]); // one more to the count
          return 3;
        }
      }
      else{ // the increment is <single>
        index_t lval=ind[last-1]+(-ind[last]-1)*single;
        if(cur-lval!=single)ind.push_back(cur); // not confroming to the increment
        else{
          --(ind[last]); // one more to the count
          return 3;
        }
      }
    }
    else{ // not incremented
      if(ind[last-1]<0){ // only one value in the buffer
        ind.push_back(cur); // wait till the next value comes
      }
      else{
        index_t inc0=ind[last]-ind[last-1];
        index_t inc=cur-ind[last];
        if(inc0==inc){ // increments coinside 
          if(inc==single){ // the increment is <single>, packing count only
            ind[last]=-3;
            return 1;
          }
          else if(inc>=0){ // only if they are non negative 
            // the increment is not <single>, packing increment and count
            ind[last]=-3;
            ind.push_back(-(inc+1));
            return 2;
          }
        }
        ind.push_back(cur); // store unpacked
      }
    }
  }
  return 1;
}

/// puts the end of record indicator 
template<class index_t>
int put_record_end(vector<index_t> &ind,index_t single){
  index_t last=(index_t)ind.size()-1;
  if (last>=0) {
    if (ind[last]<0){ // alrady in increment sequence
      if(last>0){
        if(ind[last-1]>=0){ // this was a single increment
          // modyfing it to double-increment
          ind.push_back(-(single+1));
        }
        else if(last>1 && ind[last-2]<0){
          return 1; // put_record_end called earlier
        }
      }
      else {
        if(last==0)
          return 1; // put_record_end called earlier
        ind.push_back(-1); // !!! this is impossible !
      }
      ind.push_back(-1);
    }
    else{
      ind.push_back(-1);
      ind.push_back(-1);
      ind.push_back(-1);
    }
  }
  else {
    ind.push_back(-1); // if the vector is empty
  }
  return 1;
}
