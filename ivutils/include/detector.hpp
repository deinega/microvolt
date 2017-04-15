/*****************************************************************************
 *
 *   Kinetic Technologies Ltd. 2006       
 *
 *   Authors	: Ilya Valuev, Alexei Deinega,  KINTECH, Moscow, Russia
 *
 *   Project	: FDTD-II
 *
 *   $Revision: 1.3 $
 *   $Date: 2013/10/27 22:30:16 $
 *   @(#) $Header: /home/dev/Photon/ivutils/include/detector.hpp,v 1.3 2013/10/27 22:30:16 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/include/detector.hpp,v $
$Revision: 1.3 $
$Author: valuev $
$Date: 2013/10/27 22:30:16 $
*/
/******************************************************************************
 * $Log: detector.hpp,v $
 * Revision 1.3  2013/10/27 22:30:16  valuev
 * fixed support for COMPLEX
 *
 * Revision 1.2  2013/07/13 19:46:52  lesha
 * recorder can be not active, in this case detector doesn't record anything
 *
 * Revision 1.1  2013/04/23 01:05:21  lesha
 * *** empty log message ***
 *
 * Revision 1.2  2013/03/20 17:38:38  valuev
 * restructure with proper template instantiation
 *
 * Revision 1.1  2012/12/06 02:24:05  lesha
 * *** empty log message ***
 *
 * Revision 1.18  2012/11/20 05:38:40  lesha
 * dimensions sequence is changed in detectors
 *
 * Revision 1.17  2012/11/14 03:55:57  lesha
 * *** empty log message ***
 *
 * Revision 1.16  2012/11/04 04:28:04  lesha
 * put_request: action is made as a last argument
 *
 * Revision 1.15  2012/10/16 15:35:07  lesha
 * put_full_request is modified
 *
 * Revision 1.14  2012/10/07 19:51:42  lesha
 * contaoner_t is removed from detector
 *
 * Revision 1.13  2012/09/28 21:51:07  lesha
 * *** empty log message ***
 *
 * Revision 1.12  2012/09/27 19:15:12  lesha
 * documentation
 *
 * Revision 1.11  2012/09/26 03:00:59  lesha
 * *** empty log message ***
 *
 * Revision 1.10  2012/09/26 02:47:10  lesha
 * TextRecorder is documented
 *
 * Revision 1.9  2012/09/24 13:44:16  lesha
 * pencil, seqpack and table_function are excluded from transport
 *
 * Revision 1.8  2012/06/19 12:57:18  xander
 * const was added for dim_type
 *
 * Revision 1.7  2012/04/26 17:37:16  lesha
 * file unopening is treated
 *
 * Revision 1.6  2012/04/12 04:30:02  lesha
 * documentation
 *
 * Revision 1.5  2012/03/02 19:43:44  lesha
 * Detector::StartRecord message is added
 *
 * Revision 1.4  2012/02/22 04:21:23  lesha
 * MPI_VEC_TYPE is used
 *
 * Revision 1.3  2012/02/16 02:02:43  lesha
 * VTKGridDetectorSet is moved to photonic
 *
 * Revision 1.2  2012/02/15 23:25:17  lesha
 * *** empty log message ***
 *
 * Revision 1.47  2012/02/15 23:17:43  lesha
 * emtl.h is excluded from include files of detector.h
 *
 * Revision 1.46  2012/02/15 22:30:02  lesha
 * emtl.h is excluded from headers
 *
 * Revision 1.45  2012/02/15 20:22:11  lesha
 * valtype - vec_type
 *
*******************************************************************************/

# ifndef _DETECTOR_HPP
# define _DETECTOR_HPP

#include "detector.h"
#include "string_utils.h"

template<class vset_t>
void TextRecorder<vset_t>::Flush(){
  if(!if_record || inifile<1)return;
  fflush(f);
}

template<class vset_t>
int TextRecorder<vset_t>::StartRecord(int it){
  if(!if_record)return 0;
  if(inifile<1){
    f=fopen((config->GetRawDataDir()+fname).c_str(), it ? "a" : "w");
    if(!f)
      inifile=-1;
    else{
      if(!it){
        fprintf(f,"# ");
        for(size_t i=0;i<columns.size();i++)
          fprintf(f,"%d-%s\t",int(i+1),columns[i].c_str());
        fprintf(f,"\n");
      }
      inifile=1;
    }
    return inifile;
  }
  else
    return -1;
}

template<class vset_t>
int TextRecorder<vset_t>::NextRecord(void *filebuf){
  if(!if_record)return 0;
  if(inifile<1)return -1;
  
  size_t fnum=bytes_size/sizeof(vec_type);
  vec_type *ptr=(vec_type *)filebuf;
  typename vset_t::iterator it=vset->begin();
  for(size_t i=0;i<vnum;i++,++it){
    Vector_3 pos=*it;

    if(argtype&argx)fprintf(f,"%g\t",pos[0]);
    if(argtype&argy)fprintf(f,"%g\t",pos[1]);
    if(argtype&argz)fprintf(f,"%g\t",pos[2]);

    for(size_t j=0;j<fnum;j++){
      fprintf(f,format.c_str(),ptr[fnum*i+j]);
      fprintf(f,"\t");
    }
    fprintf(f,"\n");
  }
  return 1;
}

template<class vset_t>
int TextRecorder<vset_t>::EndRecord(){
  if(!if_record)return 0;
  if(inifile<0)return -1;
  inifile=-1;
  return fclose(f);
}

template<class vset_t>
int TextFluxRecorder<vset_t>::StartRecord(int it){
  if(!if_record)return 0;
  if(inifile<1){
    f=fopen((config->GetRawDataDir()+fname).c_str(), it ? "a" : "w");
    if(!f)
      inifile=-1;
    else{
      if(!it){
        fprintf(f,"# ");
        fprintf(f,"%s",extra.c_str());
        for(size_t i=0;i<columns.size();i++)
          fprintf(f,"%d-%s\t",int(i)+1,columns[i].c_str());
        fprintf(f,"\n");
        fflush(f);
      }
      inifile=1;
    }
    return inifile;
  }
  else
    return -1;
}

template<class vset_t>
int TextFluxRecorder<vset_t>::NextRecord(void *filebuf){
  if(!if_record)return 0;
  if(inifile<1)return -1;

  int dim=0;
  for(int di=0;di<3;di++){
    if(argtype&(argx<<di))
      dim++;
  }
  size_t jnum=rec_size/dim;
  vec_type *ptr=(vec_type *)filebuf;

  fprintf(f,"%s",extra.c_str());
  for(size_t j=0;j<jnum;j++){
    vec_type flux=0;
    typename vset_t::iterator it=vset->begin();
    for(size_t i=0;i<vnum;i++,++it){
      Vector_3 pos=*it;
      Vector_3 ds=it.ds();
      Vector_3 J;
      for(int di=0,dic=0;di<3;di++){
        if(!(argtype&(argx<<di)))
          continue;
        J[di]=ptr[rec_size*i+j*dim+dic];
        dic++;
//        fprintf(f,"%g\t",J[di]);
      }
//      fprintf(f,"\n");

      flux+=J*ds;
    }
    fprintf(f,"%g\t",flux);
  }
  fprintf(f,"\n");
  fflush(f);
  return 1;
}



template<class vset_t, class xfer_t>
int DetectorSet<vset_t,xfer_t>::find_major_rank() {
  int np=comm->get_nproc();
  int *ranks = new int[np]; // vecrors numbers on each processor
  memset(ranks, 0, np*sizeof(int));

  vec_it it=vsetp->begin();
  vec_it end=vsetp->end();

  for(; it!=end; ++it){
    Vector_3 pos=transform.ptr() ? (*transform.ptr())(*it) : *it;
    int rank=xfer.test_rank(pos);
    if(rank<0 && !rec_out){
      delete []ranks;
      return ::pmessage(vblWARN,-1,"Point (%g, %g, %g) of detector set '%s' is outside the calculated space.\n",pos[0],pos[1],pos[2],apComponent::name.c_str());
    }
    ranks[rank]++;
  }
  int max_sz=0;
  int rank=0;
  for (int i=0; i<np; i++) {
    if (ranks[i]>max_sz) {
      max_sz=ranks[i];
      rank=i;
    }
  }
  delete []ranks;
  if (max_sz==0) 
    return ::pmessage(vblWARN,-1,"Detector set '%s' is empty.\n", apComponent::name.c_str()); // empty vset
  return rank;
}

template<class vset_t, class xfer_t>
void DetectorSet<vset_t,xfer_t>::AddDefaultTimers(const vector<int> &p_id){
  if(p_id.size()!=2)return;
  int build=p_id[0];
  int calc=p_id[1];
  tid.det_pre=AddTimersField("det_pre", build);
  tid.det_rec=AddTimersField("det_rec", calc);
  xfer.InitComponent(name,dump,ut,makevec<int>(2,tid.det_pre,calc));
}

template<class vset_t, class xfer_t>
void DetectorSet<vset_t,xfer_t>::SetWorkFlag(int wf_){
  wf=wf_;
  R->SetWorkFlag(wf);
  if(wf&1){
    for(size_t i=0;i<xfer.get_full_data_size()*bufsize;i++)
      filebuf[i]=0;
  }
}

template<class vset_t, class xfer_t>
int DetectorSet<vset_t,xfer_t>::Init(){
  if(inirec>0)return -1;
  start(tid.det_pre);
  message(vblMESS1,0,"Building %d detectors for '%s'...\n", vsetp->size(), name.c_str());

  R->SetSize(vsetp->size());

  int major_rank=find_major_rank(); // perform find_major_rank even if local_write==true for vsetp checking
  if(major_rank<0)
    comm->set_exit();
  if(comm->check_exit())
    return LOGERR(-1,fmt("Can't find proper major rank for detector set '%s'",name.c_str()),LINFO);

  int wr_rank=local_write&1 ? 0 : major_rank;
  int my_rank=comm->get_myrank();
  R->SetRecordingFileNode(wr_rank);

  if (my_rank==wr_rank || local_write&2){

    bufnum=xfer.register_buffer(NULL);
    bufsize=0;

    int vsz=vsetp->size();
    int deg=int(floor(log10(vec_type(vsz))+.5))-1;
    if(deg<0)deg=0;
    int ufreq=int(pow(10.,deg)); // how often print a status message
    if(ufreq<100)
      ufreq=100;

//    int ufreq=100; // how often print a status message
//    for(;vsetp->size()>(size_t)ufreq*100;ufreq*=10){}

    vec_it it=vsetp->begin(), end=vsetp->end();

    for(int i=0; it!=end; ++it){
      Vector_3 pos=transform.ptr() ? (*transform.ptr())(*it) : *it;

      int rank=xfer.test_rank(pos);
      if(rank<0) // will be only if rec_out==1
        rank=wr_rank;
      R->record(local_write&2 ? rank : wr_rank);
      if(local_write&2){
        if(my_rank!=rank)
          continue;
      }

      int res=xfer.put_full_request(bufnum, bufsize, pos);

      if(res==-1)
        return ::pmessage(vblWARN,-1,"Point (%g, %g, %g) of detector set '%s' is outside the calculated space.\n",pos[0],pos[1],pos[2],apComponent::name.c_str());
      else if(res==-2)
        return ::pmessage(vblWARN,-1,"Not enough memory for detector set '%s'.\n",apComponent::name.c_str());
      i++;
//      if(i%ufreq==0)message(vblMESS1,0,"%d positions...\n",i);
    }

    if(filebuf) delete []filebuf;
    if(bufsize){
      size_t fds=xfer.get_full_data_size();
      filebuf = new char[fds*bufsize];
      xfer.substitute_buffer(bufnum, (typename xfer_t::val_t *)filebuf); // questionable construction
      log->add_memusage(name+", buf",fds*bufsize+R->data_size());
    }
  }
//  message(vblMESS1,0,"Done.\n");

  xfer.prepare_transfers();
  if_sent=false;
  inirec=1;

  log->add_memusage(name+", transfers",&xfer);
  log->add_mpiusage(name+", send",1,&xfer);
  log->add_mpiusage(name+", receive",0,&xfer);

  stop(tid.det_pre);
  return inirec;
}

template<class vset_t, class xfer_t>
void DetectorSet<vset_t,xfer_t>::Flush(){
  if(inirec<0)return;
  R->Flush();
}

template<class vset_t, class xfer_t>
int DetectorSet<vset_t,xfer_t>::StartRecord(int it){
  if(inirec<0)
    return ::pmessage(vblWARN,-1,"DetectorSet::StartRecord is called while inirec<0.\n");
  if(!R->active())return 0;
  return R->StartRecord(it);
}

template<class vset_t, class xfer_t>
int DetectorSet<vset_t,xfer_t>::StartNextRecord(){
  if(inirec<0)return -1;
  if(!R->active())return 0;
  if(if_sent)return 0;

  if(!(wf&1)){
    xfer.start_transfers();
    xfer.compute_local();
  }
  if_sent=true;
  return 0;
}

template<class vset_t, class xfer_t>
int DetectorSet<vset_t,xfer_t>::CompleteNextRecord(){
  if(inirec<0)return -1;
  if(!R->active())return 0;
  if(if_sent==false)return 0;

  if(!(wf&1))
    xfer.complete_transfers();
  if_sent=false;

  start(tid.det_rec);
  int result=R->NextRecord(filebuf);
  stop(tid.det_rec);
  return result;
}

template<class vset_t, class xfer_t>
int DetectorSet<vset_t,xfer_t>::EndRecord() {
  if(inirec<0)return -1;
  if(!R->active())return 0;
#ifdef USE_MPI
  //MPI_Barrier(comm->get_all());
  comm->barrier();
#endif
  return R->EndRecord();
}


# endif