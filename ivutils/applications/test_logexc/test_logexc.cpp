#include "stdafx.h"
#include "logexc.h"
#include "string_utils.h"

/// to define some new exception, the following should be done:
/// 1. define some class (maybe derived from std::exception)
class my_exception: public logic_error{
public:
  my_exception(const char *str):logic_error(str){}
};

/// 2. optionally define log_exception_traits to specify vbLEVEL, name, and adding text
///    if not defined, the default traits (see logexc.h) will be used
template <>
struct log_exception_traits<my_exception>{
  /// exeption level according to the vbLEVELS: 
  static int level(const my_exception &signal){ return vblOERR; } 
  /// the string name of exception category
  static string name(const my_exception &signal){ return "my_exception";}
  /// adds words to what()
  static my_exception add_words(const my_exception &orig, const char *words){
    return my_exception(fmt("%s, %s",orig.what(),words));
  }
};

/// test function giving error message
int test(){
  return message(my_exception("Test Logic Error"),-2," at %s:%d",__FILE__,__LINE__);
}
 
/// test class giving error messages with local logger
struct funcclass{
  message_logger log;
  funcclass():log("funcclass"){
  }
  /// this function will never throw exceptions
  int func1(){
    log.set_throw(0);
    return log.message(vblOERR,-1,"error in func1");
  }
  /// this may throw
  int func2(int thr_ex=0){
    log.set_throw(thr_ex);
    return log.message(my_exception("error that may cause exception "),-1,"in func2");
  }
};


int main(){
  //1. using global logger
  // int error
  message(vblOERR,-1,"some error at %s:%d",__FILE__,__LINE__);	
  // int warning
  message(vblWARN,0,"some warning message");
  // int message (not reported by default)
  message(vblMESS4,0,"some  message");
  // changing log level
  message_logger::global().set_levels(vblALL);
  // int message (now it is reported)
  message(vblMESS4,0,"some  message");

  // now testing with exceptions
  do{
    // loggers have life time within a block but used globally
    // this will not throw exceptions, but return errorcode
    stdfile_logger log1("logger1",0,stdout,stderr,vblALLBAD,vblFATAL,1);
    int res=0;
    try{
      res=test();
      printf("got errorcode: %d\n",res);
    }
    catch(exception &e){
      printf("caught exception: %s\n",e.what());
    }
    // this will throw exception
    stdfile_logger log2("logger2",1,stdout,stderr,vblALLBAD,vblFATAL,1);
    try{
      res=test();
      printf("got errorcode: %d\n",res);
    }
    catch(exception &e){
      printf("caught exception: %s\n",e.what());
    }
  }while(0);

  /// now we automatically returned to the old global logger
  message(vblMESS4,0,"returned to old logger!\n");

  /// now we are testing local loggers
  funcclass mytest;
  try{
    mytest.func1();
    mytest.func2();
    mytest.func2(1);
  }
  catch(exception &e){
     printf("caught exception: %s\n",e.what());
  }
  /// now we are still with the old global logger
  // using macros
  LOGMSG(vblMESS4,"macro message without line report",0);
  LOGMSG(vblMESS4,"macro message with line report",1);
  LOGMSG(vblMESS4,"macro message with default line report",LINFO);
  LOGERR(-1,fmt("formatted error %s","message"),LINFO);
  message(vblMESS4,0,"finished!\n");

  return 0;
}

