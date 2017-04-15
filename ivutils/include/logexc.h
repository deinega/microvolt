#ifndef LOGEXC_H
#define LOGEXC_H

/// @file logexc.h \brief Interface to process exceptions and output messages

#include <stdarg.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <typeinfo>

#include "refobj.h"

/// this specifies whether to put file/line info in error messages
#ifndef LINFO
#define LINFO 1
#endif

using namespace std;

/// verbosity levels
enum vbLEVELS{
  vblNONE   = 0,  ///< completely silent
  vblFATAL  = 0x1,  ///< report fatal errors
  vblOERR   = 0x2, ///< report other errors
  vblERR    = vblFATAL|vblOERR, ///< report all errors
  vblWARN   = 0x4, ///< report warnings
  vblALLBAD = vblWARN|vblERR, ///< report all errors and warnings
  vblMESS1  = 0x8, ///< report messages of level 1 (important)
  vblMESS2  = 0x10, ///< report messages of level 2 (less important)
  vblMESS3  = 0x20, ///< report messages of level 3 (not important)
  vblMESS4  = 0x40, ///< report messages of level 4 (anything possible)
  vblALLMESS = vblMESS1|vblMESS2|vblMESS3|vblMESS4, ///< report all messages
  vblPROGR  = 0x80,           ///< report progress 
  vblALL    = 0xffff
};


/// traits structure to deduce exception level and name from exception class
/// by default all exceptions have vblFATAL level
template<class exc_t>
struct log_exception_traits{
  /// exeption level according to the vbLEVELS
  static int level(const exc_t &signal){ return vblFATAL; } 
  /// the string name of exception category
  static string name(const exc_t &signal){ return typeid(exc_t).name();}
  /// adds some more explanations to the description
  /// default behaviour: nothing done
  static exc_t add_words(const exc_t &orig, const char *words){
    return orig;
  }
};


/// integer exceptions have the level equal to their value
template<>
struct log_exception_traits<int>{
  /// exeption level according to the vbLEVELS
  static int level(const int &signal){ return signal; } 
  /// the string name of exception category
  static string name(const int &signal){ 
    if(signal&vblFATAL)
      return "fatal error";
    else if(signal&vblOERR)
      return "error";
    else if(signal&vblWARN)
      return "warning";
    else
      return "";
    /*
    else if(signal&vblALLMESS)
      return "message";
    else if(signal&vblPROGR)
      return "progress report";
    else
      return "integer exception";*/
  }
  /// default behaviour: nothing done
  static int add_words(const int &orig, const char *words){
    return orig;
  }
};

/// vbLEVELS exceptions act as integers
template<>
struct log_exception_traits<enum vbLEVELS>{
  static int level(const enum vbLEVELS &signal){ return log_exception_traits<int>::level(signal); } 
  static string name(const enum vbLEVELS &signal){ return log_exception_traits<int>::name(signal); } 
  static enum vbLEVELS add_words(const enum vbLEVELS &orig, const char *words){
    return orig;
  }
};


///\en
/// Logger class to control (computational) function behaviour when something requiring user attention has happened.
/// message(signal,errcode, text) is used to either throw an exception or return errorcode
/// At first, the the level of error is determined via log_exception_traits<>::level(signal)
/// For integer (enum) signals the level is the signal itself.
/// Then text is printed, if signal level is listed in output levels or (or in extra outlevels, if they are set) 
/// via log_text() function.
/// If level has vblERR bit, the behaviour is controlled by the flag specified in set_throw(flag):
/// flag=0:   nothing done;
/// flag=1:   calls add_words() for signal and throws signal;
/// flag=2:   throws pair<>(errcode, text);
/// flag=3:   throws errcode.
/// Then, if the level is listed in stop_levels (or in extra stop levels, if they are set), the program is aborted,
/// otherwise errcode is returned;
/// The function set_levels(out_levels,stop_levels) is used to specify bitflags for the levels which
/// require message output or/and program termination. Stop level has effect only when exceptions are not thrown.
/// The function extra_levels(eout_levels,estop_levels) is used to temporarily set the corresponding levels,
/// they are unset (the original levels are restored) by calling extra_levels(0,0).
///\ru
/// Логгер предназначен для обработки ошибок, возникающих во внутренних функциях библиотеки.
/// Первый механизм - это вызов исключений. Во внутренней функции кидается exeption и он ловится где-то снаружи.
/// Второй вариант - ты кидаешь код возврата. Тогда каждая функция выше должна анализировать код возврата.
class message_logger {
  // global message is a friend
  template<class exc_t>
  friend int message(const exc_t &signal, int errcode, const char *what, ...);
protected:
  string descriptor;
  int throw_ex;
  int outlevel, eoutlevel;
  int stoplevel, estoplevel;
  
  static message_logger *glogp;
  /// used to restore the previous global logger
  message_logger *prev, *next;
public:
  
  message_logger(const string &descriptor_="", int out_level=vblALLBAD|vblMESS1, 
                 int stop_level=vblFATAL, int throw_exceptions=0, int use_globally=0)
    :descriptor(descriptor_),prev(NULL),next(NULL){
    set_throw(throw_exceptions);
    set_levels(out_level,stop_level);
    extra_levels(0,0);
    if(use_globally){
      set_global(true);
    }
  }

  // clones this logger
  virtual message_logger *clone() const {
    message_logger *p= new message_logger(descriptor,outlevel,stoplevel,throw_ex,0);
    p->extra_levels(eoutlevel,estoplevel);
    return p;
  }
  
  /// returns a reference to global logger
  /// if not set, links with default message_logger
  static message_logger &global();
  
  /// sets/unsets this logger as the global logger
  int set_global(bool set){
    if(set){
      if(prev) // already set
        return -1;
      if(glogp)
        glogp->next=this;
      prev=glogp;
      glogp=this;
    }
    else{
      if(glogp!=this) // was not set as the global
        return -1; 
      glogp=prev;
      if(glogp)
        glogp->next=NULL;
      prev=NULL;
    }
    return 1;
  }
  
  virtual void set_throw(int throw_exceptions){
    throw_ex=throw_exceptions;
  }

  virtual int get_throw() const {
    return throw_ex;
  }

  virtual void set_levels(int out_level=vblALLBAD|vblMESS1, int stop_level=vblFATAL){
    outlevel=out_level;
    stoplevel=stop_level;
  }

  /// nonzero extra levels are applied instead of set ones
  virtual void extra_levels(int out_level=vblALLBAD|vblMESS1, int stop_level=vblFATAL){
    eoutlevel=out_level;
    estoplevel=stop_level;
  }
  
  template<class exc_t>
  int message(const exc_t &signal, int errcode, const char *what, ...){
    int level=log_exception_traits<exc_t>::level(signal);
    char buff[4048];
    if(level&(eoutlevel ? eoutlevel : outlevel)){ //needs to print a message
      va_list args;
      va_start(args,what);
      vsnprintf(buff,4048,what,args);
      log_text(level,log_exception_traits<exc_t>::name(signal).c_str(),buff);
    }
    int m_throw_ex= get_throw();
    if(level&vblERR){
      if(m_throw_ex==1) // throws exc_t exception 
        throw log_exception_traits<exc_t>::add_words(signal,buff);
      else if(m_throw_ex==2) // throws pair<>(int,const char*) exception 
        throw make_pair(errcode,what);
      else if(m_throw_ex==3) // throws int exception 
        throw errcode;
    } 
    if(level&(estoplevel ? estoplevel: stoplevel) ){ // needs to stop
      exit(errcode); 
    }
    return errcode;
  }

  virtual void log_text(int level, const char *messtype, const char *messtext){}


  /// checks that the deleted one is not in global logger chain 
  virtual ~message_logger(){
    if(prev){
      prev->next=next;
      if(next)
        next->prev=prev;
    }
    set_global(false);
  }
};

/// global message function
template<class exc_t>
int message(const exc_t &signal, int errcode, const char *what, ...){
  if(message_logger::glogp){
    va_list args;
    va_start(args,what);
    char buff[4048];
    vsnprintf(buff,4048,what,args);
    return message_logger::glogp->message(signal,errcode,buff);
  }
  else 
    return -1;
}

enum STDLOG_FORMAT {
  STDLOG_NEWLINE=0x1
};

/// message logger for which std and error streams may be specified
class stdfile_logger: public message_logger {
protected:
  FILE *fout, *ferr;
  int format;
  bool flush_flag;
public:
  stdfile_logger(const string &descriptor_="", int throw_exceptions=0, 
    FILE *out=stdout, FILE *err=stderr, int format=0,
    int out_level=vblALLBAD|vblMESS1,int stop_level=vblFATAL,int use_globally=0):
    message_logger(descriptor_,out_level,stop_level,throw_exceptions,use_globally),
    fout(NULL),ferr(NULL),flush_flag(true){
    set_out(out);
    set_err(err);
    set_format(format);
  }

  ~stdfile_logger(){
    if(fout && fout!=stdout)
      fclose(fout);
    if(ferr!=fout && ferr && ferr!=stderr)
      fclose(ferr);
  }

  // clones this logger
  virtual message_logger *clone() const {
    stdfile_logger *p= new stdfile_logger(descriptor,throw_ex,fout,ferr,format,outlevel,stoplevel,0);
    p->extra_levels(eoutlevel,estoplevel);
    p->flush_flag=flush_flag;
    return p;
  }

  virtual void set_out(FILE *out, int close_prev=0){
    if(close_prev && fout && fout!=stderr && fout !=stdout)
      fclose(fout);
    fout=out;
  }
  
  virtual void set_err(FILE *err, int close_prev=0){
    if(close_prev && ferr && ferr!=stderr && ferr !=stdout)
      fclose(ferr);
    ferr=err;
  }

  void set_format(int format_){
    format=format_;
  }

  virtual void log_text(int level, const char *messtype, const char *messtext){
    FILE *f= (level&vblALLBAD ? ferr : fout);
    if(!f)
      return;
    if(descriptor!="") // descriptor is used as header
      fprintf(f,"%s:\n",descriptor.c_str());
    if(string(messtype)!="")
      fprintf(f,"%s: ",messtype);
    fprintf(f,"%s",messtext);
    if(format&STDLOG_NEWLINE)
      fprintf(f,"\n");
    if(flush_flag)fflush(f);
  }
};

///\en This is a compound logger. It reveives messages and sends them to all loggers
/// in the list.
///\ru Это составной логгер. Он принимает сообшения и передает
/// их всем логгерам, которые находятся у него в списке.
class vector_logger: public refvector<message_logger>, public message_logger{
  virtual void log_text(int level, const char *messtype, const char *messtext){
    for(size_t i=0;i<size();i++){
      (*this)[i]->log_text(level,messtype,messtext);
    }
  }
public:
  vector_logger():refvector<message_logger>(1){}
  // if parent does not throw and
  // at least one of the child loggers throws exceptions then 1 is returned
  virtual int get_throw() const {
    if(throw_ex)
      return throw_ex;
    for(size_t i=0;i<size();i++){
      if((*this)[i]->get_throw())
        return 1;
    }
    return 0;
  }
  // clones this logger
  virtual message_logger *clone() const {
    vector_logger *p= new vector_logger();
    for(size_t i=0;i<size();i++){
      p->push_back((*this)[i]->clone());
    }
    return p;
  }
};

/// macros with common usage
#define LOGFATAL(code,text,lineinfo) ((lineinfo) ? ::message(vblFATAL,(code),"%s at %s:%d\n",(text),__FILE__,__LINE__) : \
                                   (::message(vblFATAL,(code),"%s",(text))) )


#define LOGERR(code,text,lineinfo) ((lineinfo) ? ::message(vblOERR,(code),"%s at %s:%d\n",(text),__FILE__,__LINE__) : \
                                   (::message(vblOERR,(code),"%s",(text))) )


#define LOGMSG(cat,text,lineinfo) ((lineinfo) ? ::message((cat),0,"%s at %s:%d\n",(text),__FILE__,__LINE__) : \
                                  (::message((cat),0,"%s",(text))) )






#endif
