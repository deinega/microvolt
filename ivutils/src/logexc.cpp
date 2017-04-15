//# ifdef USE_STDAFX
//# include "stdafx.h"
//# endif

#include "logexc.h"

stdfile_logger default_log("",0,stdout,stderr,0,vblALLBAD|vblMESS1,vblFATAL,1);

message_logger &message_logger::global(){
  if(!glogp){
    default_log.set_global(true);
  }
  return *glogp;
}

message_logger *message_logger::glogp=NULL;
