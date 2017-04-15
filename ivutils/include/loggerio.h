/*s***************************************************************************
 *
 *   Project	: ivutils
 *
 *   $Revision: 1.2 $
 *   $Date: 2013/03/05 14:14:44 $
 *   @(#) $Header: /home/dev/Photon/ivutils/include/loggerio.h,v 1.2 2013/03/05 14:14:44 belousov Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/include/loggerio.h,v $
$Revision: 1.2 $
$Author: belousov $
$Date: 2013/03/05 14:14:44 $
*/
/*s****************************************************************************
 * $Log: loggerio.h,v $
 * Revision 1.2  2013/03/05 14:14:44  belousov
 * *** empty log message ***
 *
 * Revision 1.1  2012/12/06 02:24:04  lesha
 * *** empty log message ***
 *
 * Revision 1.23  2012/09/25 20:11:48  lesha
 * *** empty log message ***
 *
 * Revision 1.22  2012/09/25 19:48:30  lesha
 * documentation
 *
 * Revision 1.21  2012/09/25 19:23:07  lesha
 * documentation
 *
 * Revision 1.20  2012/09/25 05:30:33  lesha
 * *** empty log message ***
 *
 * Revision 1.19  2011/09/27 12:36:59  biaks
 * some bags fixed for mingw config
 *
*******************************************************************************/
# ifndef LOGGERIO_H
# define LOGGERIO_H

# include <stdio.h>
# include <string.h>
# include <string>

using namespace std;

/**\en @file loggerio.h
       @brief Classes for storing data sequences in files. 
**/

/**\en class for storing slices of data_t values d:
 i            <-- sequence length (size_t)
 i            <-- number of slices (size_t)
 d d d d d d  <-- slice0
 d d d d d d  <-- slice1, etc.
 in a file
 the slices are of equal length
 header (i i) cab be skipped if flag if_header = false

 error codes for all functions:
              -1 out of range
              -2 invalid file state (not open/ already open)
              -3 invalid file (can't open)
              -4 wrong file format */

class SeqRecord {
  /// length of unit data in bytes
  size_t data_size;
  /// length of the data sequence
  size_t seqlen;
  /// current file
  FILE *file;
  /// the maximum number left for the last slice
  size_t rest;
  /// the current number of slices
  size_t numslices;
  /// presence of header at the begining of the file
  bool if_header;
public:
  /// constructor
  /// parameter dlen is the length of each slice
  SeqRecord(size_t ds=sizeof(double),size_t dlen=1,bool if_header_=false):
    data_size(ds),seqlen(dlen),file(NULL),if_header(if_header_){}

  /// returns -2 if file is open: can't change existing file
  int SetDataSize(size_t ds) {
    if(file)return -2;
    data_size=ds;
    return 1;
  }

  int SetSeqLen(size_t dlen) {
    if(file)return -2;
    seqlen=dlen;
    return 1;
  }

  /// opens a file
  /// if the mode is "w" then replaces the old file if it exists and starts new record
  /// if the mode is "a" then appends the old file if it exists (after checking format ) or starts new record if not
  /// if the mode is "r" analyses the given file (checks it format) and can perform GetData()
  /// if the mode is "ar" can do both reading and appending
  /// other combinations of r and w are not supported
  /// record will start from shift slice
  /// @return the current slice number, -2 if file is already open
  int OpenRecord(const string &name, const char *mode="w", const int shift=0) {
    if(file)return -2; // already open
    file=fopen(name.c_str(),"r");
    bool is_file=(file!=NULL);
    if (is_file)
      fclose(file);

    char real_mode[10];
    if (*mode=='r') {
      if (is_file)
        strcpy(real_mode, "rb"); // open existing file for reading
      else
        return -3;
    }
    else if (*mode=='w') {
      strcpy(real_mode, "w+b"); // open emty file for reading and writing (if file exists then it will be deleted)
      is_file=false;
    }
    else if (*mode=='a') {
      if (is_file)
        strcpy(real_mode, "r+b"); // open existing file for reading and writing
      else
        strcpy(real_mode, "w+b"); // open emty file for reading and writing (if file exists then it will be deleted)
    }
    else
      return -3;
    
    file=fopen(name.c_str(),real_mode);
    if(!file)
      return -3;

    rest=seqlen;

    if (is_file) {
      if (if_header) { // check seqlen and numslices
        size_t seqlen1;
        int res=(int)fread(&seqlen1,sizeof(size_t), 1, file);
        res+=(int)fread(&numslices, sizeof(size_t), 1, file);
        if(res!=2 || seqlen1!=seqlen){
          CloseRecord(); // wrong format
          return -4; 
        }
      }
      else {
        if(shift){
          numslices=shift;
          long long pos=(long long)(numslices)*(long long)(seqlen)*(long long)(data_size);
#ifndef _MSC_VER // gcc compiler
          fseeko64(file, pos, SEEK_SET);
# else
         _fseeki64(file, pos, SEEK_SET);
# endif
        }
        else{
          // checking the size
#ifndef _MSC_VER // gcc compiler
          fseeko64(file, 0, SEEK_END);
          long long fsize=ftello64(file);
# else
          _fseeki64(file, 0, SEEK_END);
          long long fsize=_ftelli64(file);
# endif
          numslices=(size_t)(fsize/(data_size*seqlen));
        }
      }
    }
    else {
      numslices=0;
      if (if_header) { // record header
        int res = fwrite(&seqlen,sizeof(size_t),1,file);
        res += fwrite(&numslices,sizeof(size_t),1,file);
//        fflush(file);
      }
    }
    return (int)numslices;
  }

  /// appends next arrlen records to a file from rec array
  /// if arrlen exceeds the maximum number left for current slice, appends the correct number of elements only
  /// DATA IS APPENDED BY COPYING FROM MEMORY (using fwrite) <- this behaviour may be changed in the future...
  /// returns the number of elements appended
  int AppendData(const void *rec, size_t arrlen=1) {
    if(!file)return -2;
    size_t count=(rest>arrlen) ? arrlen : rest;
/*    long long pos=((long long)(numslices+1)*(long long)(seqlen)-(long long)(rest))*(long long)(data_size);
    if(if_header)
      pos+=(long long)(2*sizeof(size_t));
# ifdef UNIX
    fseeko(file, pos, SEEK_SET);
# else
   _fseeki64(file, pos, SEEK_SET);
# endif*/
    count=fwrite(rec, data_size, count, file);
    rest-=count;
//    fflush(file);
    return (int)count;
  }

  /// switches to the next slice and fflushes the file
  /// @return the current slice number or <0 in case of error
  int NextSlice() {
    if(!file)return -2;
    if (!rest) {
      numslices++;
      if (if_header) {
        fseek(file, sizeof(size_t), SEEK_SET);
        int res = fwrite(&numslices, sizeof(size_t), 1, file);
        if(res!=1)
          return -3;
//        fflush(file);
      }
      rest=seqlen;
      return (int)numslices;
    }
    else
      return -1;
  }

  int Flush(){
    if(!file)return -2;
    return fflush(file);
  }

  /// closes the file
  /// @return the current slice number
  int CloseRecord() {
  /// doesn't remove unfinished slice, but this slice will be ignored next time because number of slices is not incremented
    if(file){
      fclose(file);
      file=NULL;
      return (int)numslices;
    }
    else return -2;
  }

  /// gets data from  slice  slicenum starting from beg till end
  /// all counters start from zero
  /// returns the number of elements read
  int GetData(void *rec, size_t slicenum, size_t beg=0, size_t end=0) {
    if(!file)return -2;
    if ((slicenum<numslices)&&(beg<seqlen)) {
      if (end>=seqlen)end=seqlen-1;
      long long pos=((long long)(slicenum)*(long long)(seqlen)+(long long)(beg))*(long long)(data_size);
      if(if_header)
        pos+=(long long)(2*data_size);
#ifndef _MSC_VER
      fseeko64(file, pos, SEEK_SET);
# else
     _fseeki64(file, pos, SEEK_SET);
# endif
      return (int)fread(rec, data_size, end-beg+1, file);
    }
    else
      return -1;
  }
};

# endif
