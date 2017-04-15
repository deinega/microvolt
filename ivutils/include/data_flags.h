#ifndef DATA_FLAGS_H
#define DATA_FLAGS_H

///\en @file data_flags.h \brief Some usefull bit flags

///\en Bit flags to specify which field argument values are needed in output
enum argTYPE{
  ///\en output time or frequency
  argt=0x1,
  ///\en output x coordinate
  argx=0x2,
  ///\en output y coordinate
  argy=0x4,
  ///\en output z coordinate
  argz=0x8,
  ///]en output all coordinates
  argr=argx|argy|argz,
  ///\en output all
  argAll=argt|argx|argy|argz,
  ///\en use wavelength instead of frequency
  ///\ru выводить длину волны, а не частоту
  argwv=0x10, 
};

///\en Bit modifiers to specify which part of the value is needed in case of complex number
enum outCOMPLEX_MODIFIERS{
  ///\en output absolute value
  outMOD=0x1,
  ///\en output real part
  outRE=0x2,
  ///\en output imaginary part
  outIM=0x4
};

///\en Specify change of what argument changes the sequence number of field output file for a set of space points (detector set)
enum outSWEEP{
  ///\en use one file for all data
  outMERGE=0,
  ///\en files numbered according to space points of the set
  outR=0x1,
  ///\en files numbered according to time (frequency) values 
  outT=0x2
};

///\en Argument change orders
enum emArgChangeOrders {
  ///\en time/frequency coordinate t(f) changes first
  ORD_TX=0,
  ///\en the same
  ORD_FX=0,
  ///\en space coordinate changes first: component ordering depend on the set, usually the order is (zyx) 
  ORD_XF=1,
  ORD_XT=1,
  ORD_RANDOM=2
};

///\en Argument (time or space) and direction of FFT 
enum FFTdir{
  TF_BACKWARD=0,
  XK_BACKWARD=0,
  TF_FORWARD=1,
  XK_FORWARD=2
};

///\en Box side descriptors
enum BOX_SIDES{
  BOX_BACK_X=0x1,
  BOX_BACK_Y=0x2,
  BOX_BACK_Z=0x4,
  BOX_FRONT_X=0x8,
  BOX_FRONT_Y=0x10,
  BOX_FRONT_Z=0x20,
  BOX_ALL=BOX_BACK_X|BOX_FRONT_X|BOX_BACK_Y|BOX_FRONT_Y|BOX_BACK_Z|BOX_FRONT_Z,
  BOX_OUTSIDE=0,
  BOX_INSIDE=0x40
};

///\en Format flags for tabular text files
enum TABLE_LOG{
  TABLE_TAB=0x1, ///<\en separate columns by tab (otherwise space)
  TABLE_HORIZONTAL=0x2, ///<\en print rows as a columns and columns as a rows
  TABLE_HEADER=0x4, ///<\en print header with column names at the first line
  TABLE_COMMENT=0x8, ///<\en comment header by '#'
  TABLE_APPEND=0x10 ///<\en append to the file
};

#endif
