/*****************************************************************************
 *
 *   Kinetic Technologies Ltd. 2006-2011       
 *
 *   Authors	: Ilya Valuev, Alexei Deinega,  KINTECH, Moscow, Russia
 *
 *   Project	: ivutils, FDTD-II
 *
 *   $Revision: 1.8 $
 *   $Date: 2013/11/24 19:07:54 $
 *   @(#) $Header: /home/dev/Photon/ivutils/src/component.cpp,v 1.8 2013/11/24 19:07:54 lesha Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/src/component.cpp,v $
$Revision: 1.8 $
$Author: lesha $
$Date: 2013/11/24 19:07:54 $
*/
/******************************************************************************
 * $Log: component.cpp,v $
 * Revision 1.8  2013/11/24 19:07:54  lesha
 * bug with timers is fixed
 *
 * Revision 1.7  2013/11/24 19:02:48  lesha
 * timers are made as double (not vec_type)
 *
 * Revision 1.6  2013/10/27 20:11:15  lesha
 * mult_nodes is moved as an argument to apInit
 *
 * Revision 1.5  2013/08/05 16:13:22  lesha
 * apInitDir is modified
 *
 * Revision 1.4  2013/06/24 11:42:03  valuev
 * xfel reader
 *
 * Revision 1.3  2013/01/26 01:44:06  lesha
 * nothing important
 *
 * Revision 1.2  2013/01/24 16:58:46  valuev
 * removed virtual functions from base contour
 *
 * Revision 1.1  2012/12/06 02:24:05  lesha
 * *** empty log message ***
 *
 * Revision 1.24  2012/11/14 03:55:57  lesha
 * *** empty log message ***
 *
 * Revision 1.23  2012/10/04 17:38:51  lesha
 * make_vec is moved to utiltl
 *
 * Revision 1.22  2012/10/03 00:51:21  lesha
 * VTK dumping is moved to gnudump
 *
 * Revision 1.21  2012/10/02 23:11:30  lesha
 * SpaceRegion is renamed to Region_3
 *
 * Revision 1.20  2012/09/27 19:08:44  lesha
 * documentation
 *
 * Revision 1.19  2012/09/25 01:12:46  lesha
 * documentation
 *
*******************************************************************************/

#include "component.h"
#include "data_flags.h"

FILE *file_manager::open_file(const string &fname, int *count, int group_ind, const char *mode){
  struc_t &s=fmap[fname];
  if(!s.f)
    s.f=::open_file(fname.c_str(),mode);
  if(count)
    *count=s.count;
  s.count++;
  s.group_ind=group_ind;
  return s.f; 
}

void file_manager::close_files(int group_ind){
  for(map_t::iterator it=fmap.begin(), e=fmap.end();it!=e;it++){
    if(it->second.f && (it->second.group_ind==group_ind || group_ind==-1)){
      fclose(it->second.f);
      it->second.f=NULL;
    }
    //it->second.count=0;
  }
}


int apConfig::set_dir(const string &idir, string &odir){
  odir=idir;
  if(odir.length()){
    correct_directory_name(odir);
    if(!check_directory(odir))
      return LOGERR(-1,fmt("Directory '%s' is inaccessible\n", odir.c_str()),0);
  }
  return 0;
}

void apConfig::SetSysArgs(int argc_, char **argv_, apMPIComm *comm){
  int argc=argc_;
#ifdef USE_MPI
  //MPI_Bcast(&argc,1,MPI_INT,0,comm->get_all());
  comm->bcast_int(&argc,1,0);
#endif
  for(int i=0;i<argc;i++){
    char argbuf[251];
    if(comm->get_myrank()==0){
      strncpy(argbuf,argv_[i],250);
    }
#ifdef USE_MPI
    //MPI_Bcast(argbuf,250,MPI_CHAR,0,comm->get_all());
    comm->bcast_char(argbuf,250,0);
#endif
    argv.push_back(argbuf);
  }
}




double apTimer::gettime(bool local){
#ifdef USE_MPI
  if(local)
    return theComm.wtime(); //MPI_Wtime();
  else{
    double sec;
    if(theComm.get_myrank()==0)
      sec=theComm.wtime(); //MPI_Wtime();
    //MPI_Bcast(&sec,1,MPI_VEC_TYPE,0,theComm.get_all());
    theComm.bcast_double(&sec,1,0);
    return sec;
  }
#else
  return (double)clock()/(double)CLOCKS_PER_SEC;
#endif
}

apTimer *apTimer::bind(apTimer *p){
  parentv.push_back(p);
  return this;
}

int apTimer::start(double t0_){
  if(!stlev){
    t0 = t0_<0 ? gettime() : t0_;
    for(size_t i=0;i<parentv.size();i++)
      parentv[i]->start(t0);
  }
  return stlev++;
}

int apTimer::stop(int force, double tset_){
  if(!stlev)return 0; // not started
  if(force)stlev=0;
  else stlev--;
  if(!stlev){
    double tset = tset_<0 ? gettime() : tset_;
    time+=tset-t0;
    for(size_t i=0;i<parentv.size();i++)
      parentv[i]->stop(force, tset);
  }
  return stlev;
}

double apTimer::update(){
  if(stlev){
    double t1=gettime();
    time+=t1-t0;
    t0=t1;
  }
  return time;
}

/*FILE *apRunLog::open_file(const char *fname, const char *mode){
  return config->open_file(fname,mode);
}*/

int apRunLog::AddTimersField(const string &name, int par_id){
  if(!ut)return -1;
  if(names.find(name)!=names.end()){
    return names[name];
  }
  apTimer *t=new apTimer();
  if(par_id>=0)t->bind(timers[par_id]);
  timers.push_back(t);
  AddField(t->get_ptr(),name);
  names[name]=(int)timers.size()-1;
  return (int)timers.size()-1;
}

int apRunLog::TimeUsage(const char *filename, int format){
  update_all();
  FILE *f;
  char tab=format&TABLE_TAB ? char(9) : ' ';
  if(format&TABLE_APPEND)
    format|=TABLE_HORIZONTAL;
  if(format&TABLE_APPEND){
    // checking if the file exists to write header only once
    f=config->open_file(filename,"rt");
    if(!f)format&=(~TABLE_APPEND);
    else fclose(f);
  }
//  f=open_file(filename,append? "at": "wt");
  f=config->open_file(filename, format&TABLE_APPEND? "at": "wt");
  if(f==NULL){
    message(vblWARN,-1,"could not open file %s\n",filename);
    theComm.set_exit();
  }
  if(theComm.check_exit())return -1;

  vector<double> vouta;
  // printing  descriptors

  if(format&TABLE_HORIZONTAL){
    if(!(format&TABLE_APPEND)){
      if(format&TABLE_COMMENT)
        fprintf(f,"#");
      vector<string> descra;
      int nfa=(int)GetFields(descra);
      for(int j=0;j<nfa;j++)fprintf(f,"%d %s%c",j+1,descra[j].c_str(),tab);
      fprintf(f,"\n");
    }
    int nfa=process(vouta,0);
    for(int j=0;j<nfa;j++)fprintf(f,"%g%c",vouta[j],tab);
    fprintf(f,"\n");
  }
  else {
    vector<string> descra;
    int nfa=(int)GetFields(descra);
    if (nfa!=process(vouta,0))
      return -1;
    char ch[100];
    for(int i=0;i<100;i++)
      ch[i]=' ';
    size_t max=0;
    for(size_t i=0;i<descra.size();i++){
      if(descra[i].length()>max)
        max=descra[i].length();
    }
    for(int j=0;j<nfa;j++) {
      if(vouta[j]==0)continue;
      string name=descra[j];
      ch[max-descra[j].length()]=0;
      name.append(ch);
      fprintf(f,"%d %s%c",j+1,name.c_str(),tab);
      fprintf(f,"%g",vouta[j]);
      fprintf(f,"\n");
      ch[max-descra[j].length()]=' ';
    }
  }
  fclose(f);
  return 1;
}

apRunLog::memstat_t &apRunLog::get_mem(const string &name){
  size_t num;
  if(memnames.find(name)!=memnames.end())
    num=memnames[name];
  else{
    num=mems.size();
    memnames[name]=num;
    mems.push_back(memstat_t(name));
  }
  return mems[num];
}

void apRunLog::add_memusage(const string &opname, size_t data_size, size_t packed_size){
  memstat_t &mst=get_mem(opname);
  mst.full+=data_size;
  mst.packed+=packed_size>0?packed_size:data_size;
}

/// logs recorded memory usage statistics
void apRunLog::MemoryUsage(message_logger *mt) const {
  mt->message(vblMESS1,0,"Memory usage statistics:\n");
  long long sum=0;
  for(vector<memstat_t>::const_iterator it=mems.begin(), e=mems.end();it!=e;++it){
    long long full=it->full;
    if(full){
      string sfull=size_in_bytes(full);
      long long packed=it->packed;
      if(packed<full){
        string spacked=size_in_bytes(packed);
        mt->message(vblMESS1,0,"%s: size %s, unpacked_size %s (packing %.3f)\n",it->name.c_str(),spacked.c_str(),sfull.c_str(),float(packed)/float(full));
      }
      else
        mt->message(vblMESS1,0,"%s: size %s\n",it->name.c_str(),sfull.c_str());
      sum+=full;
    }
  }
  string ssum=size_in_bytes(sum);
  mt->message(vblMESS1,0,"Total size %s\n",ssum.c_str());
}

/// logs recorded MPI statistics
void apRunLog::MPIUsage(message_logger *mt) const {
  if(theComm.get_nproc()==1)return;
  mt->message(vblMESS1,0,"MPI usage statistics:\n");
  for(vector<memstat_t>::const_iterator it=mems.begin(), e=mems.end();it!=e;++it){
    long long send=it->mpi[1];
    long long recv=it->mpi[0];
    if(send || recv){
      string ss=size_in_bytes(send);
      string sr=size_in_bytes(recv);
      mt->message(vblMESS1,0,"%s: send %s, receive %s\n",it->name.c_str(),ss.c_str(),sr.c_str());
    }
  }
}

apConfig theConfigObj, *theConfig=&theConfigObj;
apRunLog theLogObj, *theLog=&theLogObj;
apComponent::registry_t apComponent::registry;

int apInitArgs(int argc, char **argv){
  theComm.Init(argc,argv); // MPI intialization
  if(theComm.get_myrank()==0){
    for(int i=0;i<argc;i++){
      if(argv[i]==NULL){
        argc=i;
        break;
      }
    }
  }
  theConfig->SetSysArgs(argc,argv);

  return 0;
}

int apInitDir(string outdir, string rawdir, bool silent, bool mult_nodes){
  int res=0;
  // assign directory for output text files
  if(theConfig->SetOutputDir(outdir)<0)
    res=-1;
  // assign directory for output binary files
  if(theConfig->SetRawDataDir(rawdir)<0)
    res=-1;
  if(res<0)
    return -1;

  if(mult_nodes || !theComm.get_myrank()){
    char tmp[100];
    // file for processor(s) output messages
    if(mult_nodes)
      sprintf(tmp,"proc_%03d.dat",theComm.get_myrank()); // get processor rank
    else
      sprintf(tmp,"proc.dat");
    string myout=theConfig->outdir+tmp;
    FILE *pfile=fopen(myout.c_str(), "w");
    theConfig->ui_log.push_back(new stdfile_logger("",0,pfile,pfile));
  }
#if  !defined UNIX || !defined USE_MPI || defined LOG_STDOUT
  // substitute the previously set logger
  message_logger *prev=message_logger::global().clone();
  theConfig->ui_log.push_back(prev); /*new stdfile_logger("",0,stdout,stderr)*/
#endif
  // make this logger as a global (it will collect messages from LOGERR, LOGMSG)
  theConfig->ui_log.set_global(true);
  theConfig->ui_log.set_managed(1);
//  message_logger::global().set_throw(1);

  if(theComm.get_myrank()==0){
    string info=theConfig->info();
    message(vblMESS1,0,info.c_str());
  }

#ifdef USE_MPI
  if(!silent)
    message(vblMESS1,0,"Run is started with %d nodes. This is node %d\n",theComm.get_nproc(),theComm.get_myrank());
#endif

  int ac=(int)theConfig->argv.size();
  if(!silent){
    if(ac)
      message(vblMESS1,0,"There %s %d command line parameter%s\n", ac>1 ? "are" : "is", ac, ac>1? "s:" : ":");
    for(int i=0;i<ac;i++)
      message(vblMESS1,0,"%d %s\n",i+1,theConfig->argv[i].c_str());
    message(vblMESS1,0,"\n");
  }
  return 0;
}

int apInit(int argc, char **argv, string outdir, string rawdir, bool silent, bool mult_nodes){
  apInitArgs(argc,argv);
  return apInitDir(outdir, rawdir, silent, mult_nodes);
}
/*
///\en outdir and rawdir are numbers of command line parameters with correspondig directories
int apInit(int argc, char **argv, size_t outdir, size_t rawdir){
  apInitArgs(argc,argv);
  string sout,sraw;
  if(theConfig->argv.size()>outdir)
    sout=theConfig->argv[outdir];
  if(rawdir==size_t(-1))
    sraw=sout;
  return apInitDir(sout,sraw);
}
*/

void apComponent::DumpOther(const apComponent *comp, bool lim){
  dumper->Dump(comp,lim);
}

void apComponent::Dump(bool lim){
  dumper->Dump(this,lim);
}

string apDump::add_spec(const apComponent *comp,string name){
  if(prefix!="")
    name=prefix+"_"+name;
  string cn=comp->get_name();
  if(cn!="")
    name=cn+"_"+name;
  return name;
}

void apDump::Dump(const apComponent *comp, bool lim){
  string dname(typeid(*this).name()); // complete type name of the dumper class 
  string cname(typeid(*comp).name()); // complete type name of the component class 
  
  map<string,fdump>::const_iterator it=apDump::fmap.find(dname+"."+cname);
  if(it!=fmap.end()){ // found in registered components
    it->second(comp,this,lim);
  }
}


map<string,apDump::fdump> apDump::fmap;

/// dump object for debug output
apDump theDumpObj, *theDump=&theDumpObj; 

/*FILE *apDump::open_file(const string &fname, int *count, int group_ind, const char *mode){
  return fm.open_file(config->outdir+fname,count,group_ind, mode);
}*/


# ifndef NO_GNUDUMP

int apDump::DumpRegion(const Region_3 *reg,const string &fname, bool if_lim){
  if(comm->get_myrank())return 0;

  RegDumper<3> *dmp=reg->CreateDumper();
  if(!dmp)
    return -1;
  vector<VecContour<> > sides;
  dmp->Dump(sides, if_lim ? lim.ptr() : NULL);
  if(!if_lim && sides.size()==0)
    dmp->Dump(sides,lim.ptr());
  delete dmp;
  if(sides.size()==0)
    return 0;

  int count;
  if(format&GNUPLOT){
    string fullname=fname+".pol"; // name convention for gnuplot 
    FILE *f=fm.open_file(config->GetOutputDir()+fullname,&count);
    if(!f)
      return LOGERR(-1,fmt("emDump::DumpRegion: can't (re)open file '%s' for writing!\n",fname.c_str()),0); 
    ::DumpContours(f,sides,!count);
    fflush(f);
  }
  if(format&VTK){
    string fullname=fname+".vtp"; // name convention for VTK poly
    FILE *f=fm.open_file(config->GetOutputDir()+fullname,&count);
    if(!f)
      return LOGERR(-1,fmt("emDump::DumpRegion: can't (re)open file '%s' for writing!\n",fname.c_str()),0); 
    ::DumpContoursVTK(f,sides,!count);
    fflush(f);
  }
  return 1;
}

int apDump::DumpDirections(int n, const Vector_3 *origin, const Vector_3 *k, const string &fname, bool if_lim){
  if(comm->get_myrank())return 0;
  int count;
  FILE *f=fm.open_file(config->GetOutputDir()+fname,&count);
  if(!f)
    return LOGERR(-1,fmt("emDump::DumpDirections: can't (re)open file '%s' for writing!\n",fname.c_str()),0); 
  ::DumpDirections(f,n,origin,k,if_lim ? lim.ptr() : NULL,!count);
  fflush(f);
  return 1;
}

# else // NO_GNUDUMP

int apDump::DumpRegion(const Region_3 *reg,const string &fname, bool if_lim){
  return 0;
}

int apDump::DumpDirections(int n, const Vector_3 *origin, const Vector_3 *k, const string &fname, bool if_lim){
  return 0;
}

#define INSTANTIATE_GetRegionDumper(type) \
	template RegDumper<type::dimension> *GetRegionDumper< type >(const type *);

// Explicit template specification is needed for MSVC
INSTANTIATE_GetRegionDumper(Region<2>);
INSTANTIATE_GetRegionDumper(Circle);
INSTANTIATE_GetRegionDumper(Region<3>);
INSTANTIATE_GetRegionDumper(Box);
INSTANTIATE_GetRegionDumper(Polyhedron<Box::plane_it>);
INSTANTIATE_GetRegionDumper(Polyhedron<Plane_3*>);
INSTANTIATE_GetRegionDumper(Sphere);
INSTANTIATE_GetRegionDumper(Cylinder<Circle>);
INSTANTIATE_GetRegionDumper(Cone<Circle>);
INSTANTIATE_GetRegionDumper(ConfinedRegion<Region_3>);
INSTANTIATE_GetRegionDumper(ConfinedRegion<Polyhedron_3>);
INSTANTIATE_GetRegionDumper(ConfinedRegion<Cylinder<Circle> >);
INSTANTIATE_GetRegionDumper(ConfinedRegion<Cone<Circle> >);
INSTANTIATE_GetRegionDumper(ConfinedRegion<Sphere>);
INSTANTIATE_GetRegionDumper(StretchedRegion<Region<3> >);
INSTANTIATE_GetRegionDumper(Inverse<Region_3>);
INSTANTIATE_GetRegionDumper(Inverse<Polyhedron_3>);


# endif