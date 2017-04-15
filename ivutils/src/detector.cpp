/*****************************************************************************
 *
 *   Kinetic Technologies Ltd. 2006       
 *
 *   Authors	: Ilya Valuev, Alexei Deinega,  KINTECH, Moscow, Russia
 *
 *   Project	: FDTD-II
 *
 *   $Revision: 1.6 $
 *   $Date: 2013/10/18 19:28:58 $
 *   @(#) $Header: /home/dev/Photon/ivutils/src/detector.cpp,v 1.6 2013/10/18 19:28:58 lesha Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/src/detector.cpp,v $
$Revision: 1.6 $
$Author: lesha $
$Date: 2013/10/18 19:28:58 $
*/
/******************************************************************************
 * $Log: detector.cpp,v $
 * Revision 1.6  2013/10/18 19:28:58  lesha
 * documentation of new Fourier update
 *
 * Revision 1.5  2013/10/17 19:10:15  lesha
 * table is added to fourier transform recorded
 *
 * Revision 1.4  2013/09/15 04:50:29  lesha
 * ufreq>1 case is fixed
 *
 * Revision 1.3  2013/07/13 19:46:52  lesha
 * recorder can be not active, in this case detector doesn't record anything
 *
 * Revision 1.2  2013/04/20 03:39:06  lesha
 * *** empty log message ***
 *
 * Revision 1.1  2013/04/20 03:32:12  lesha
 * *** empty log message ***
 *
 * Revision 1.2  2013/03/20 17:38:00  valuev
 * restructure with extensible contour table
 *
 * Revision 1.1  2013/02/12 02:46:40  lesha
 * include and src
 *
 * Revision 1.34  2012/10/07 19:50:59  lesha
 * contaoner_t is removed from detector
 *
 * Revision 1.33  2012/02/16 02:04:21  lesha
 * VTKGridDetectorSet is moved to near_to_far
 *
 * Revision 1.32  2012/02/15 20:13:55  lesha
 * emDetectorSet: prefix em is excluded
 *
 * Revision 1.31  2011/04/14 05:57:15  lesha
 * cxfer_t is added
 *
 * Revision 1.30  2011/03/09 20:06:11  lesha
 * em part is moved out from detectors
 *
 * Revision 1.29  2010/10/24 23:31:48  valuev
 * added parallel VTK writer
 *
 * Revision 1.28  2010/06/06 08:58:30  lesha
 * emFieldsRecorder is made non-template
 *
 * Revision 1.27  2009/01/30 14:02:15  valuev
 * restructured as a library
 *
*******************************************************************************/

#include "detector.hpp"


void FileRecorder::Flush(){
  if(!if_record || inifile<1)return;
  SR.Flush();
}

int FileRecorder::StartRecord(int it){
  if(!if_record)return 0;
  if(inifile<1){
    SR.SetSeqLen(vnum);
    inifile=SR.OpenRecord(config->GetRawDataDir()+fname, it ? "a" : "w", it)>=0 ? 1 : -1;
    return inifile;
  }
  else
    return -1;
}

int FileRecorder::NextRecord(void *filebuf){
  if(!if_record)return 0;
  if(inifile<1)return -1;
  if(SR.AppendData(filebuf, vnum)<0)return -1;
  return SR.NextSlice();
}

int FileRecorder::EndRecord(){
  if(!if_record)return 0;
  if(inifile<0)return -1;
  inifile=-1;
  return SR.CloseRecord();
}


void emFourierRecorder::record(int node){
  if(ti<0){
//    nodes.start_record();
    nodes.clear();
    ti=0;
  }
//  nodes.next_record(node);
  nodes.push_back(node);
  if(comm->get_myrank()==node)bufsize++;
}

void emFourierRecorder::Flush(){
  if(fname.length() && fnum)
    FileRecord(fname);
}

// it!=0 case is not realized
int emFourierRecorder::StartRecord(int it){
//  if(ti==0)
//    nodes.end_record();
  if(bufsize && fnum){
    try{
      F=new vec_type[fnum*bufsize*2*rec_size];
      for(size_t i=0;i<fnum*bufsize*2*rec_size;i++)
        F[i]=0;
      if(exp_upd){
        exp_table_dt=new cvec_type[fnum];
        exp_table_t=new cvec_type[fnum];
        for(int fi=0;fi<fnum;fi++){
          vec_type f=fmin+fi*df;
          vec_type phase=ufreq*exp_dir*2*M_PI*f*dt;
          exp_table_dt[fi]=exp(cvec_type(0,phase));
          exp_table_t[fi]=1;
        }
      }
    }
    catch(...){
      ::pmessage(vblMESS1,0,"emFourierRecorder: not enough memory (%s are requested)\n",
        size_in_bytes(2*bytes_size*fnum*bufsize).c_str());
      comm->set_exit();
    }
  }
  if(comm->check_exit())return -1;
  return 1;
}

int emFourierRecorder::NextRecord(void *filebuf){
  if(!bufsize || !fnum)return 0;
  if(!(wf&1)){
    for(int fi=0;fi<fnum;fi++){
      cvec_type e;
      if(!exp_upd || ti%exp_upd==0){
        vec_type t=ti*dt, f=fmin+fi*df;
        vec_type phase=exp_dir*2*M_PI*f*t;
        e=exp(cvec_type(0,phase));
        if(exp_upd)
          exp_table_t[fi]=e;
      }
      else
        e=exp_table_t[fi];
      vec_type cs=e.real();
      vec_type sn=e.imag();
//      vec_type cs=cos(phase);
//      vec_type sn=sin(phase);
      for(size_t vi=0;vi<bufsize;vi++){
        for(size_t i=0;i<rec_size;i++){
          F[(fi*bufsize+vi)*2*rec_size+i]+=((vec_type *)filebuf)[vi*rec_size+i]*cs*ufreq;
          F[(fi*bufsize+vi)*2*rec_size+rec_size+i]+=((vec_type *)filebuf)[vi*rec_size+i]*sn*ufreq;
        }
      }
    }
    if(exp_upd){
      for(int fi=0;fi<fnum;fi++)
        exp_table_t[fi]*=exp_table_dt[fi];
    }
  }
  ti+=ufreq;
  return 1;
}

int emFourierRecorder::EndRecord(){
  if(fname.length() && fnum)
    FileRecord(fname);
  if(F)delete []F;
  if(exp_table_dt)delete []exp_table_dt;
  if(exp_table_t)delete []exp_table_t;
  F=NULL, exp_table_dt=NULL, exp_table_t=NULL;
  return 1;
}

void emFourierRecorder::FileRecord(const string &fname){
  int my_rank=comm->get_myrank();
  int nproc=comm->get_nproc();

  SeqRecord SR(2*bytes_size);
  vec_type *F2 = NULL; // for data transfers on writing processor, if not some data is collected on other processors
  int *tags = NULL; // used on writing processor for receiving parallel messages
  int tag=0; // used on non-writing processors for sending parallel messages

  if(my_rank==wr_rank){
    SR.SetSeqLen(vnum);
    SR.OpenRecord(config->GetRawDataDir()+fname,"w");

    if(bufsize!=vnum){ // not all data is on writing processor
      F2 = new vec_type[vnum*2*rec_size];
      for(size_t i=0;i<vnum*2*rec_size;i++)
        F2[i]=0;
      tags = new int[nproc];
      for(int i=0;i<nproc;i++)
        tags[i]=0;
    }
  }

  for(int fi=0;fi<fnum;fi++){
    if(my_rank==wr_rank){
      if(bufsize==vnum){ // all data is on writing processor
        SR.AppendData(F+fi*vnum*2*rec_size,vnum);
      }
#ifdef USE_MPI
      else{
        // проходит по пакеру noded и последовательно скидывает данные с доменов на F2 пишущего домена
        int count=0; // размер текущей порции точек, которая записывается с одного домена
        int displ=0; // текущий сдвиг в F2
        int my_displ=0; // текущий сдвиг в F
        int cur_rank=0; // текущий домен
//        for(int_pack::iterator it=nodes.begin(),e=nodes.end();;++it){
        for(vector<int>::iterator it=nodes.begin(),e=nodes.end();;++it){
          if(!(it!=e) || cur_rank!=*it){ // кончился пакер или осуществился переход на новый домен
            if(cur_rank==my_rank){
              for(size_t i=0;i<count*2*rec_size;i++)
                F2[i+displ*2*rec_size]=F[i+(fi*bufsize+my_displ)*2*rec_size];
              my_displ+=count;
            }
            else if(count){
              //MPI_Status status;
              //MPI_Recv(F2+displ,count*12,MPI_VEC_TYPE,cur_rank,tags[cur_rank],comm->get_all(),&status);
              comm->recv_vec_type(F2+displ*rec_size*2,count*rec_size*2,cur_rank,tags[cur_rank]);
              tags[cur_rank]++;
              if(tags[cur_rank]>=32767)
                tags[cur_rank]=0;
            }
            cur_rank=(it!=e)?*it:-1; // меняем текущий домен
            displ+=count; // увеличиваем сдвиг в F2 на количество записанных точек
            count=1;
          }
          else
            count++; // увеличиваем размер порции точек на текущем домене
          if(!(it!=e))break;
        }
        SR.AppendData(F2,vnum);
      }
#endif
      SR.NextSlice();
    }
#ifdef USE_MPI
    else if(bufsize){
      int count=0; // размер текущей порции точек, которая записывается с одного домена
      int displ=0; // текущий сдвиг в F
//      for(int_pack::iterator it=nodes.begin(),e=nodes.end();;++it){
      for(vector<int>::iterator it=nodes.begin(),e=nodes.end();;++it){
        if(it!=e && my_rank==*it)
          count++;
        else if(count){ // сменился домен с текущего на другой, делаем отправку порции данных
          //MPI_Send(F+fi*bufsize+displ,count*12,MPI_VEC_TYPE,wr_rank,tag,comm->get_all());
          comm->send_vec_type(F+(fi*bufsize+displ)*2*rec_size,count*12,wr_rank,tag);
          tag++;
          if(tag>=32767)
            tag=0;
          displ+=count; // увеличиваем сдвиг в F на количество отправленных точек
          count=0;
        }
        if(!(it!=e))break;
      }
    }
    //MPI_Barrier(comm->get_all());
    comm->barrier();
#endif
  }
  if(my_rank==wr_rank){
    if(bufsize!=vnum){
      delete []F2;
      delete []tags;
    }
    SR.CloseRecord();
  }
}

int emFourierRecorder::write_read(FILE *f, write_read_f *fun){
  if(bufsize && fnum){
    fun(&ti,sizeof(int),1,f);
    fun(F,2*bytes_size,fnum*bufsize,f);
  }
  return 1;
}

