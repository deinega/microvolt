/*****************************************************************************
 *
 *   Kinetic Technologies Ltd. 2006-2011       
 *
 *   Authors	: Ilya Valuev, Alexei Deinega,  KINTECH, Moscow, Russia
 *
 *   Project	: ivutils, FDTD-II
 *
 *   $Revision: 1.9 $
 *   $Date: 2013/11/24 19:02:48 $
 *   @(#) $Header: /home/dev/Photon/ivutils/include/component.h,v 1.9 2013/11/24 19:02:48 lesha Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/include/component.h,v $
$Revision: 1.9 $
$Author: lesha $
$Date: 2013/11/24 19:02:48 $
*/
/******************************************************************************
 * $Log: component.h,v $
 * Revision 1.9  2013/11/24 19:02:48  lesha
 * timers are made as double (not vec_type)
 *
 * Revision 1.8  2013/10/27 20:11:15  lesha
 * mult_nodes is moved as an argument to apInit
 *
 * Revision 1.7  2013/08/05 16:10:04  lesha
 * destructur is added to stdfile_logger
 *
 * Revision 1.6  2013/06/24 11:42:03  valuev
 * xfel reader
 *
 * Revision 1.5  2013/02/22 11:01:26  valuev
 * thin surface interface
 *
 * Revision 1.4  2013/01/24 16:58:46  valuev
 * removed virtual functions from base contour
 *
 * Revision 1.3  2012/12/17 22:57:30  lesha
 * *** empty log message ***
 *
 * Revision 1.2  2012/12/17 12:48:05  chorkov
 * work around: build on lomonosov
 *
 * Revision 1.1  2012/12/06 02:24:04  lesha
 * *** empty log message ***
 *
 * Revision 1.25  2012/11/01 14:20:25  lesha
 * 4 level warning is fixed
 *
 * Revision 1.24  2012/10/26 05:43:07  lesha
 * nothing important
 *
 * Revision 1.23  2012/10/04 17:38:51  lesha
 * make_vec is moved to utiltl
 *
 * Revision 1.22  2012/10/02 23:11:29  lesha
 * SpaceRegion is renamed to Region_3
 *
 * Revision 1.21  2012/10/01 15:15:27  lesha
 * open_file is modified
 *
 * Revision 1.20  2012/10/01 15:06:59  lesha
 * GetRAM is added
 *
 * Revision 1.19  2012/09/27 21:29:18  lesha
 * *** empty log message ***
 *
 * Revision 1.18  2012/09/27 19:08:40  lesha
 * documentation
 *
 * Revision 1.17  2012/09/25 05:20:39  lesha
 * gnudump is modified
 *
*******************************************************************************/
#ifndef _COMPONENT_H
#define _COMPONENT_H

/**\en @file component.h
    @brief Basic types and auxiliary templates for application construction
  
*/

#include <map>
#include <set>
#include <string>
#include <vector>
#include <time.h>

using namespace std;

#include "refobj.h"
#include "ifmpi.h"
#include "logexc.h"
#include "logger.h"
#include "string_utils.h"


# ifndef NO_GNUDUMP

#include "gnudump.h" // to do: separate the dumper from geometry

# else // define dumper stubs

#include "region_3.h" // needed for dumper

template<int N>
class RegDumper{
public:
  virtual bool Dump(vector<VecContour<N> > &sides, Box *lim=NULL) { return false;}
};

template<class reg_t>
RegDumper<reg_t::dimension> *GetRegionDumper(const reg_t *reg){
  return new RegDumper<reg_t::dimension>();
}


# endif

///\en calls message (defined in logexc.h) prepending its text with node number 
template<class exc_t>
int pmessage(const exc_t &signal, int errcode, const char *what, ...){
  va_list args;
  va_start(args,what);
  char buff[1024];
  sprintf(buff,"%d node: ",theComm.get_myrank());
  vsnprintf(buff+strlen(buff),1024,what,args);
  return ::message(signal,errcode,buff);
}

int apInitArgs(int argc, char **argv);

///\en when silent is true, only info is printed as vblMESS1, otherwise other starting parameters are specified
///\en when mult_nodes is true, all nodes produce log files proc_xxx.dat (xxx - node number),
/// otherwise only first node does it (proc.dat)
int apInitDir(string outdir, string rawdir, bool silent = false, bool mult_nodes = false);


///\en configuration settings
///\ru конфигурационные настройки
class apConfig{
  friend int apInitArgs(int argc, char **argv);
  friend int apInitDir(string outdir, string rawdir, bool silent, bool mult_nodes);

  typedef string (*_info)();

  static string base_info(){
    return "This is IVUTILS library, (c) 1996-2011 Ilya Valuev, Moscow\n";
  }

  string outdir; ///\en directory for application output files
  string rawdir; ///\en alternative directory for application output files, can be used for raw data (big files)

  ///\en accessible RAM (in bytes). can be used to optimize application, 
  /// f. e. to store extra data (which can not be placed in RAM) in binary files
  long long RAM;

  vector_logger ui_log; ///\en used logger to process exceptions and print output messages (see logexc.h)

  ///\en checks if format of directory idir is correct (or tries to fix it otherwise), and set it as a odir
  /// if directory idir doesn't exist, returns -1
  int set_dir(const string &idir, string &odir);

  ///\en set location of directory for application output files
  int SetOutputDir(const string &dir){
    return set_dir(dir,outdir);
  }

  ///\en set location of alternative directory for application output files
  int SetRawDataDir(const string &dir){
    return set_dir(dir,rawdir);
  }

public:

  vector<string> argv; ///\en command line arguments

  _info info; ///\en application name

  apConfig():RAM(200*1024*1024), info(&base_info){}

  string GetOutputDir() const{
    return outdir;
  }

  string GetRawDataDir() const{
    return rawdir;
  }

  void SetRAM(long long RAM_){
    RAM=RAM_;
  }

  long long GetRAM(){
    return RAM;
  }

  ///\en sets argv on all processors correspongingly to argc_ and argv_ on zero processor
  void SetSysArgs(int argc_=0, char **argv_=NULL, apMPIComm *comm=&theComm);

  ///\en Same as fopen() but prepends outdir directory name to the file name before opening
  FILE *open_file(const char *fname, const char *mode){
    return ::open_file(fname, mode, outdir);
  }
};


///\en Global configuration settings object. By default all components
///    use the pointer to this config object. Must be initialized 
///    at the start of application (in \ref apInit() function)
extern apConfig *theConfig, theConfigObj;

///\en Initialization of the applicaton
int apInit(int argc, char **argv, string outdir="", string rawdir="", bool silent = false, bool mult_nodes=false);

/**\en timer that can be started, stopped and measure its working time
 usually timer corresponds to some process, for example solving number of equations
 timer can be part of some timers hierarchy which corresponds to hierarchy of processes,
 for example solving of some equations is part of modeling of some physical process (upper level process),
 but it has as a part solving of one particular equation (lower level process).
 hierarchy of timers should correspond this hierarchy of processes.
 starting lower level timer causes starting upper level timers,
 so upper level timers should summate time measured by their lower level timers 
 */
class apTimer {
  vector<apTimer *> parentv; ///\en timers of upper level
  double time; ///\en measured working time
  double t0; ///\en start time
  int stlev; ///\en number of started timers of lower level plus this timer
  
public:

  apTimer(int started=0, apTimer *parent_=NULL):parentv(0),time(0),t0(0),stlev(0){
    if(parent_)bind(parent_);
    if(started)start();
  }

  ///\en return current time in seconds 
  /// if local==true, it will be time on current processor,
  /// otherwise time on zero processor
  static double gettime(bool local=true);

  ///\en attach timer p to the list of upper level timers parentv
  apTimer *bind(apTimer *p);
 
  ///\en if timer is not started, start this timer and all timers of upper level
  int start(double t0_=-1);

  ///\en if all lower level timers are stopped, stop timer 
  /// and inform upper level timers about this
  /// if force, stop timer and upper level timers
  int stop(int force=0, double tset_=-1); 

  ///\en updates time and t0. returns time
  double update(); 

  double *get_ptr(){
    return &time;
  }

  ~apTimer(){
    stop(1);
  }
};

/**\en collection of timers and data about memory used by different operations.
 time measured by timers and data about memory can be printed by Write and Usage */
/**\ru Хранилище таймеров, которые вызываются пользователем,
 и хранилище данных о памяти, задействованной разными объектами.
 Даные таймеров и памяти можно распечатать с помощью Write и Usage */
class apRunLog: public process_detail<double,double>{

  refvector<apTimer> timers; ///\en storage for timers
  map<string, size_t> names; ///\en timer's name - its number in vector timers
  int ut; ///\en if timers will be used

  ///\en information about what memory is reserved for some object
  struct memstat_t{
    string name; ///\en object name
    long long full; ///\en full unpacked memory in bytes 
    long long packed; ///\en actually reserved memory in bytes
    long long mpi[2]; ///\en memory reserved for MPI transfers (send and receive)
    memstat_t(const string &name_):name(name_),full(0),packed(0){
      mpi[0]=mpi[1]=0;
    }
  };
  vector<memstat_t> mems; ///\en storage for memstat_t
  map<string, size_t> memnames; ///\en memstat_t's name - its number in vector mems

  apConfig *config;

  ///\en just transfer values 
  virtual int process(vector<double> &vout, const double &r, int refnum=0){
    return (int)TransferValues(vout);
  }

  memstat_t &get_mem(const string &name);

public:

  apRunLog(apConfig *config_=theConfig):ut(1), config(config_), timers(1){}

  ~apRunLog(){
    for(size_t i=0;i<timers.size();i++)
      timers[i]->stop(1);
  }

  ///\en creates timer with given name, adds it to timers 
  /// and binds it to par_id timer (or doesn't bind to any timer if par_id=-1)
  int AddTimersField(const string &name, int par_id=-1);

  int start(int id){
    if(ut)return timers[id]->start();
    else return -1;
  }

  int stop(int id, int force=0){
    if(ut)return timers[id]->stop(force);
    else return -1;
  }

  double get_time(int id){
    if(ut)return *(timers[id]->get_ptr());
    else return -1;
  }

  void stop_all(){
    for(size_t i=0;i<timers.size();i++){
      timers[i]->stop(1);
    }
  }
  void update_all(){
    for(size_t i=0;i<timers.size();i++){
      timers[i]->update();
    }
  }

  ///\en writes a tabular file: timer name - measured working time
  int TimeUsage(const char *filename, int format);

  ///\en adds memory statistics for some operation described by memstat_t from mems
  /// (if it is not present in mems, creates new one)
  void add_memusage(const string &opname, size_t data_size, size_t packed_size=0);

  ///\en adds memory statistics (packed in packer) for some operation 
  template<class packer_t>
  void add_memusage(const string &opname,const packer_t *packer){
    add_memusage(opname,(size_t)packer->data_size(),(size_t)packer->packed_size());
  }

  ///\en adds mpi statistics for some operation 
  template<class packer_t>
  void add_mpiusage(const string &opname,int sendrcv,const packer_t *packer){
    memstat_t &mst=get_mem(opname);
    mst.mpi[sendrcv]+=packer->mpi_size(sendrcv);
  }

  ///\en logs recorded memory usage statistics
  void MemoryUsage(message_logger *mt=&message_logger::global()) const;
  
  ///\en logs recorded MPI statistics
  void MPIUsage(message_logger *mt=&message_logger::global()) const;

  void Usage(message_logger *mt=&message_logger::global()) const{
    MemoryUsage(mt);
    MPIUsage(mt);
  }
};


///\en Global logging variable. By default all components
///    use the pointer to this object. Must be initialized 
///    at the start of application (in \ref apInit() function)
extern apRunLog *theLog, theLogObj;

class apComponent;

/**\en Files manager.
       Opens file with given name, collects its descriptor in some storage,
       if file with this name is required, retrieve its descriptor from the storage
       and doesn't open file second time */
/**\ru Менеджер файлов, следящий за их открытием и закрытием.
       Открывает файл с заданным именем, хранит дескриптор на него в хранилище, 
       при повторном обращении к файлу извлекает дескриптор из хранилища
       и не допускает повторного открытия файла */
class file_manager{

  struct struc_t{
    FILE *f; ///\en file descriptor
    int count; ///\en how many times function open_file was called for this file
    int group_ind; ///\en index in some group
    struc_t():f(NULL),count(0),group_ind(0){}
  };

  typedef map<string,struc_t> map_t;
  map_t fmap; ///\en file map

public:

  ///\en gets file descriptor and open the file if it is not opened already
  /// records to count how many times function open_file was called for this file
  FILE *open_file(const string &fname, int *count=NULL, int group_ind=0, const char *mode="wt");

  ///\en closes files with given group_ind or all files if group_ind==-1
  void close_files(int group_ind=-1);
};


class apDump;
extern apDump *theDump;


/**\en Class for requesting dumps of objects and components (including geometric) into files.
 Manager files to manage file which are opened or closed
 */
class apDump{
  ///\en auxilliary structure to call non-member DumpComponent(comp_tt *, dumper_tt *, ...) function
  ///    from apDump::Dump for specific compononets based on typeid
  template <class dumper_tt,class comp_tt>
  struct comp_dumper_t{
    static void Dump(const apComponent *comp,apDump *dump=theDump, bool lim=false){
      DumpComponent((comp_tt *)comp,(dumper_tt *)dump,lim);
    }
  };

  template <class dumper_tt,class comp_tt>
  friend void AddDumpComponent();
protected:
  typedef void (*fdump)(const apComponent *, apDump *, bool);
  static map<string,fdump> fmap; // this class name.component name - function which dumps component
  
  apConfig *config;
  apMPIComm *comm;
  mngptr<Box> lim;
  file_manager fm; ///\en collection of files currently used for dumping

  string prefix; // for dump files
  int format;

  ///\en add prefix and component name as a prefix to given name
  string add_spec(const apComponent *comp,string name);

public:

  ///\en Constructor, all specializations are on by default
  apDump(apConfig *config_=theConfig, apMPIComm *comm_=&theComm): config(config_),comm(comm_),format(GNUPLOT|VTK){}
  
  virtual void SetLimitingBox(mngarg<Box> ptr){
    lim.reset(ptr);
  }

  void SetFormat(int format_){
    format=format_;
  }

  void SetPrefix(const string &prefix_) {
    prefix=prefix_;
  }
 
  ///\en find if this dump and component are presented in fmap and apply corresponding dump function
  virtual void Dump(const apComponent *comp, bool lim=false);
  
  int DumpRegion(const Region_3 *reg,const string &fname, bool if_lim=true);

  int DumpDirections(int n, const Vector_3 *origin, const Vector_3 *k, const string &fname, bool if_lim=true);

  template<class vset_t>
  int DumpVectorSet(const vset_t *vset,const VecTransform *trans,const apComponent *comp,const string &fname, bool if_lim=true){
    if(comm->get_myrank())return 0;
    int count;
    FILE *f=fm.open_file(config->GetOutputDir()+add_spec(comp,fname),&count);
    if(!f)
      return LOGERR(-1,fmt("emDump::DumpVectorSet: can't (re)open file '%s' for writing!\n",add_spec(comp,fname).c_str()),0); 
    ::DumpVectorSet(f,vset,trans,!count);
    fflush(f);
    return 1;
  }

  ~apDump(){
    fm.close_files();
  }

};

extern apDump *theDump, theDumpObj;

template <class dumper_tt,class comp_tt>
void AddDumpComponent(){
  string d=typeid(dumper_tt).name();
  string s=typeid(comp_tt).name(); // Microsoft bug: typeid causes invalid 
                             //memory leak detection (http://social.msdn.microsoft.com/forums/en-US/vclanguage/thread/90b17894-aad7-4c9d-b769-9816d0734450)
  apDump::fmap[d+"."+s]=apDump::comp_dumper_t<dumper_tt,comp_tt>::Dump;
}


/**\en class that controls general behaviour of all components of the library.
 it has pointers to global objects of classes apConfig (to use some global settings), 
 apDump (to dump information about itself), apRunLog (to measure its working time using timers), 
 apMPIComm (to use MPI commutator) */
/**\ru умеет рисовать себя, считать время своей работы с помощью таймеров, 
 открывать файлы, пользоваться заданным MPI коммутатором */
class apComponent{
protected:
  int dump; ///\en if information about this object will be dumped
  int ut; ///\en if timers will be used
  string name; ///\en name used in dump and runlog files
  typedef set<apComponent *> registry_t;
  static registry_t registry;

  mngptr<apConfig> config;
//  mngptr<apVirtDump> dumper;
  mngptr<apDump> dumper; // пока сделан emDump, поскольку в нем есть невиртуальные шаблонные функции (DumpContour)
  mngptr<apRunLog> log;
  mngptr<apMPIComm> comm;

  ///\en connect internal timers hierarchy from log to external timers p_id
  virtual void AddDefaultTimers(const vector<int> &p_id){}

  int AddTimersField(const string &tname, int par_id){
    return ut ? log->AddTimersField(name+"."+tname,par_id) : -1;
  }

  ///\en start timers
  void start(int id){
    if(ut)log->start(id);
  }

  ///\en stop timers
  void stop(int id, int force=0){
    if(ut)log->stop(id,force);
  }

public:

  apComponent():config(theConfig),dumper(theDump),log(theLog),comm(&theComm){
    // registering the component
    registry.insert(this);
    InitComponent();
  }
  
  void InitComponent(string name_="", int dump_=0, int ut_=0, const vector<int> &p_id=vector<int>()){
    name=name_;
    dump=dump_;
    ut=ut_;
    if(ut)
      AddDefaultTimers(p_id);
  }

  ///\en borrow settings from other component
  void InitComponent(const apComponent &other, const vector<int> &p_id=vector<int>()){
    config.reset(other.config.ptr());
    dumper.reset(other.dumper.ptr());
    log.reset(other.log.ptr());
    comm.reset(other.comm.ptr());
    InitComponent(other.name,other.dump,other.ut,p_id);
  }

  virtual void SetLimitingBox(mngarg<Box> ptr){
    dumper->SetLimitingBox(ptr);
  }

  ///\en dump some component using dumper
  virtual void DumpOther(const apComponent *comp, bool lim=false);

  ///\en dump some component using dumper and checking bit flags 
  virtual void DumpOther(int flags, const apComponent *comp, bool lim=false){
    if(dump&flags)
      DumpOther(comp,lim);
  }

  ///\en dump itself using dumper
  virtual void Dump(bool lim=false);

  void set_name(const string &name_){
    name=name_;
  }
  
  string get_name() const{
    return name;
  }

  apMPIComm *get_comm() const{
    return comm.ptr();
  }

  virtual ~apComponent(){
    // unregistering the component
    registry.erase(this);
  }
};

/**\en Base class for saving (reading) internal data to (from) the file
  using finction save_data (load_data).
  Derived classes inherit this functionality and override virtual function write_read */
/**\ru Базовый класс, осуществляющий запись внутренних данных в файл и ее выгрузку оттуда 
  с помощью функций save_data и load_data.
  Сам этот класс ничего не делает и служит лишь как предок для других классов,
  которые наследуют его функциональность и переопределяют виртуальную функцию write_read */
class restorer{
  static size_t fwrite_data(void *p, size_t size, size_t n, FILE *f){
    return fwrite(p,size,n,f);
  }
  static size_t fread_data(void *p, size_t size, size_t n, FILE *f){
    return fread(p,size,n,f);
  }
public:
  typedef size_t write_read_f(void *p, size_t size, size_t n, FILE *f);
  /*\en write/read data to/from the file 
   this function should be overrided in derived classes
   second argument is fwrite_data of fread_data correspindingly to
   if file is being written or read */
  /**\ru Функция для записи/считывания из файла, переопределяемая в классах-потомках.
   В write_read_f передается либо fwrite_data, либо fread_data, в зависимости от того, 
   идет ли запись в файл или считывание оттуда */
  virtual int write_read(FILE *f, write_read_f *fun){
    return 1;
  }
  ///\en save data to the file
  int save_data(FILE *f){
    return write_read(f,fwrite_data);
  }
  ///\en read data from the file
  int load_data(FILE *f){
    return write_read(f,fread_data);
  }
};

# endif
