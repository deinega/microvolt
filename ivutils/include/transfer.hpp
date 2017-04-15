/*****************************************************************************
 *
 *   Kinetic Technologies Ltd. 2006-2009       
 *
 *   Authors	: Ilya Valuev, Alexei Deinega,  KINTECH, Moscow, Russia
 *
 *   Project	: FDTD-II
 *
 *   $Revision: 1.10 $
 *   $Date: 2013/11/03 18:15:57 $
 *   @(#) $Header: /home/dev/Photon/ivutils/include/transfer.hpp,v 1.10 2013/11/03 18:15:57 lesha Exp $
 *
 *****************************************************************************/
/*
$Source: /home/dev/Photon/ivutils/include/transfer.hpp,v $
$Revision: 1.10 $
$Author: lesha $
$Date: 2013/11/03 18:15:57 $
*/
/******************************************************************************
 * $Log: transfer.hpp,v $
 * Revision 1.10  2013/11/03 18:15:57  lesha
 * *** empty log message ***
 *
 * Revision 1.9  2013/11/02 17:38:34  lesha
 * operator * for argument
 *
 * Revision 1.8  2013/10/14 07:53:02  lesha
 * InterpTransfer
 *
 * Revision 1.7  2013/10/03 12:45:15  valuev
 * edited comments
 *
 * Revision 1.6  2013/09/29 06:22:57  lesha
 * slight change
 *
 * Revision 1.5  2013/09/29 06:11:57  lesha
 * transfer documentation
 *
 * Revision 1.4  2013/09/29 01:33:53  lesha
 * documentation
 *
 * Revision 1.3  2013/09/28 21:21:01  lesha
 * transfer documentation
 *
 * Revision 1.2  2013/09/05 08:44:01  valuev
 * continued restructure for 64-bits (added index_t)
 *
*******************************************************************************/

#ifndef _TRANSFER_HPP
#define _TRANSFER_HPP

#include "transfer.h"

template<class arg_tt, class container_tt, class reg_tester_t, class interp_packer_t, class index_t>
InterpTransfer<arg_tt, container_tt, reg_tester_t, interp_packer_t, index_t>::InterpTransfer(container_t *cont_,interp_packer_t *packer_,interp_packer_t *packer_loc_):
cont(cont_), ip(packer_),iplocal(packer_loc_),srequests(NULL), rrequests(NULL), sstatuses(NULL), rstatuses(NULL), 
rec_state(0-(cont_==NULL)), commtbl((size_t)comm->get_nproc(),1), nsend(0), nrecv(0), mpidump(0){

# ifdef USE_MPI
  // creating specific communicator
  //MPI_Comm_dup(comm->get_all(), &mycomm);
  comm->dup_all(&mycomm);
  //MPI_Comm_rank(mycomm,&myrank);
  comm->comm_rank(mycomm,&myrank);
  //MPI_Type_contiguous( sizeof(irequest_t), MPI_BYTE, &ireq_mpit); 
  //MPI_Type_commit(&ireq_mpit);
  comm->create_type(sizeof(irequest_t),&ireq_mpit);
  //MPI_Type_contiguous( sizeof(intpair_t), MPI_BYTE, &intpair_mpit); 
  //MPI_Type_commit(&intpair_mpit);
  comm->create_type(sizeof(intpair_t),&intpair_mpit); 
# else
  myrank=0;
# endif

  np=comm->get_nproc();

  for(size_t i=0;i<commtbl.size();i++){
    commtbl[i]= new comm_t();
# ifdef USE_MPI
    commtbl[i]->mycomm=mycomm;
# endif
  }

  if(mpidump){
    char str[100];
    sprintf(str,"mpi%d_p%d.log",idcount,myrank);
    logfile=fopen(str,"wt");
  }
  else logfile=NULL;

  idcount++;    
  rankvec.resize((size_t)comm->get_nproc());
  // reserved output buffers
  for(int i=0;i<np;i++)
    register_buffer(NULL);    
}

template<class arg_tt, class container_tt, class reg_tester_t, class interp_packer_t, class index_t>
void InterpTransfer<arg_tt, container_tt, reg_tester_t, interp_packer_t, index_t>::AddDefaultTimers(const vector<int> &p_id){
  if(p_id.size()!=2)return;
  int build=p_id[0];
  int calc=p_id[1];
  tid.tr_prep=AddTimersField("tr_prep",build);
  tid.tr_lcl=AddTimersField("tr_lcl",calc);
  tid.tr_send=AddTimersField("tr_send",calc);
  tid.tr_recv=AddTimersField("tr_recv",calc);
  tid.tr_gthr=AddTimersField("tr_gthr",calc);
  tid.tr_sctr=AddTimersField("tr_sctr",calc);
}

template<class arg_tt, class container_tt, class reg_tester_t, class interp_packer_t, class index_t>
void InterpTransfer<arg_tt, container_tt, reg_tester_t, interp_packer_t, index_t>::init_mpi_table(){
  if(mpidump){
    fprintf(logfile,"initializing MPI transfers...\n");
    fflush(logfile);
  }
    
# if USE_MPI
  // counting the nonzero transfer sizes
  nsend=nrecv=0;
  for(size_t i=0;i<commtbl.size();i++){
    commtbl[i]->init();
    commtbl[i]->reset_buffer(1);
    if(commtbl[i]->insz!=0)nrecv++;
    if(commtbl[i]->outsz!=0)nsend++; 
  }

  if(srequests)free_mpi_requests();

  size_t span=comm->create_nblck_array(nsend,&srequests,&sstatuses);
  comm->create_nblck_array(nrecv,&rrequests,&rstatuses);

  // making requests for nonzero size transfers
  nsend=nrecv=0;
  for(size_t i=0;i<commtbl.size();i++){
    request_t *s=NULL, *r=NULL; // array for MPI_Request
    if(commtbl[i]->insz!=0){
      r=(char *)rrequests+nrecv*span;
      nrecv++;
    }
    if(commtbl[i]->outsz!=0){
      s=(char *)srequests+nsend*span;
      nsend++; 
    }
    commtbl[i]->make_requests(reqpair_t(s,r,(int)i),comm.ptr());
    if(mpidump){
      if(commtbl[i]->outsz!=0)
        fprintf(logfile,"send to   %lu: %lu values, request %ld\n",i,commtbl[i]->outsz,*(long *)s);
      if(commtbl[i]->insz!=0)
        fprintf(logfile,"recv from %lu: %lu values, request %ld\n",i,commtbl[i]->insz,*(long *)r);
    }      
  }
  if(mpidump){
    fprintf(logfile,"\n");
    fflush(logfile);
  }
# endif
}

template<class arg_tt, class container_tt, class reg_tester_t, class interp_packer_t, class index_t>
void InterpTransfer<arg_tt, container_tt, reg_tester_t, interp_packer_t, index_t>::free_mpi_requests(){
# ifdef USE_MPI
  comm->free_nblck_array(nsend,&srequests,&sstatuses);
  comm->free_nblck_array(nrecv,&rrequests,&rstatuses);
# endif
}

template<class arg_tt, class container_tt, class reg_tester_t, class interp_packer_t, class index_t>
int InterpTransfer<arg_tt, container_tt, reg_tester_t, interp_packer_t, index_t>::put_vrequest(int ibuf, int ind, 
vector<arg_t> &argsum, coef_t c_sum, int action/*, emDumpInfo *from, emDumpInfo *to*/){
  if(rec_state!=1){ // adjust record state
    ip.start_record();
    iplocal.start_record();
    rec_state=1;
  }
  interp_form_t coeff;
  ibuf*=2;
  int mibuf=ibuf+action;
  for(size_t k=0;k<argsum.size();k++){
    // determining the rank of the requested position
    int rank=cont->test_rank(argsum[k].pos, container_t::MODIFY);
    if(rank<0) // no action for undefined point
      continue;
     
    // local or non-local?
    if(rank==myrank){ // local, put directly to the local hash
      vector<arg_t> nonlocal;
      if (coeff.valid()) {
        interp_form_t coeff_k=cont->create_interpolation(argsum[k],0,&nonlocal);
        coeff=linear_mix(coef_t(1.),coeff,c_sum,coeff_k);
      }
      else {
        coeff=cont->create_interpolation(argsum[k],1,&nonlocal);
        if (c_sum!=1.)
          coeff*=c_sum; // insert c_summ into create_interpolation
      }
        
      // nonlocal, level-1 requests
      for(size_t i=0;i<nonlocal.size();i++){
        int rank1=cont->test_rank(nonlocal[i].pos, container_t::MODIFY);
        if (rank1<0)
          continue;
        // coding the action in ibuf
        rankvec[rank1].push_back(irequest_t(myrank,ibuf,ind,1,nonlocal[i]*c_sum));
      }
#if 0
// THIS IS OLD DUMPER, SHOULD BE MODIFIED 
      if(from && to && coeff.valid()){// currently only local transfers are dumped
        if(from){
          from->pos=argsum[k].pos;
          from->fdirv=argsum[k].fdir; 
        }
        dumper_t mydump(this,from,to);
        mydump.dump_on_init(ind,coeff);
      }
#endif
    }
    else{ // non-local, save for negotiating (basic level-0 request)
      rankvec[rank].push_back(irequest_t(myrank,ibuf,ind,0,argsum[k]*c_sum));
    }
  }
  // recording the local part
  if(coeff.valid() || action){ // unvalid coeff is recorded to zero the buffer element (prepare for nonlocal sum)
    iplocal.next_record(make_pair(coeff,make_pair(mibuf,ind)));
  }
  return 1;
}

template<class arg_tt, class container_tt, class reg_tester_t, class interp_packer_t, class index_t>
int InterpTransfer<arg_tt, container_tt, reg_tester_t, interp_packer_t, index_t>::fill_buffers(iterator it, iterator e){
  int nb=(int)buff.size();
  for(;it!=e;++it){
    typename base_t::value_t iv=*it;
    int ibuf=iv.second.first;
    if (buff[ibuf/2]==NULL)
      continue; // some buffers can be turned off 
    int ind=iv.second.second;
    // computing the value
    val_t v=cont->get_interp_value(&iv.first);  // iv.first is interp_form_t
//    reval(v,v);
    if(ibuf%2)
      buff[ibuf/2][ind]=v; // rerecord
    else 
      buff[ibuf/2][ind]+=v; // add
  }
  return 1;
}

#if 0
template<class arg_tt, class container_tt, class reg_tester_t, class interp_packer_t>
int InterpTransfer<arg_tt, container_tt, reg_tester_t, interp_packer_t>::fill_buffers(iterator it, iterator e, int *shifts){
  int nb=(int)buff.size();
  for(;it!=e;++it){
    // computing the value
    typename base_t::value_t iv=*it;
    val_t val=cont->get_interp_value(&iv.first);

    int ibuf=iv.second.first;
    int ind=iv.second.second;
    ind=shifts[ind];
    if(ibuf%2)
      buff[ibuf/2][ind]=val;
    else 
      buff[ibuf/2][ind]+=val;
  }
  return 1;
}
#endif

template<class arg_tt, class container_tt, class reg_tester_t, class interp_packer_t, class index_t>
int InterpTransfer<arg_tt, container_tt, reg_tester_t, interp_packer_t, index_t>::transfer_requests(){
  // communicating sizes
  for(int i=0;i<np;i++)
    svecsz[i]=(int)rankvec[i].size(); // number of interpolation requests for each processor
# ifdef USE_MPI
  // creating specific communicator
  //MPI_Comm mycomm1;
  void *mycomm1;
  //MPI_Comm_dup(comm->get_all(), &mycomm1);
  comm->dup_all(&mycomm1);
  //MPI_Alltoall(svecsz,1,MPI_INT,rvecsz,1, MPI_INT,mycomm1);
  comm->int_alltoall(svecsz,1,rvecsz,1,mycomm1);
# else
  rvecsz[0]=svecsz[0];
# endif
  // preparing request send-recv buffers
  int ssz=0, rsz=0;
  for(int i=0;i<np;i++){
    rdispl[i]=rsz;
    sdispl[i]=ssz;
    rsz+=rvecsz[i];
    ssz+=svecsz[i];
  }
  irequest_t *rbuf= rsz? new irequest_t[rsz] : NULL;
  irequest_t *sbuf= ssz? new irequest_t[ssz] : NULL;
  // copying interpolation requests
  for(int i=0;i<np;i++){
    for(int j=0;j<svecsz[i];j++)
      sbuf[sdispl[i]+j]=rankvec[i][j];
    rankvec[i].clear();
  }
  // communicating interpolation requests
# ifdef USE_MPI
  //MPI_Alltoallv(sbuf, svecsz, sdispl, ireq_mpit, rbuf, rvecsz, rdispl, ireq_mpit, mycomm1);
  comm->alltoallv(sbuf, svecsz ,sdispl, ireq_mpit, rbuf, rvecsz, rdispl, ireq_mpit, mycomm1);
  //MPI_Comm_free(&mycomm1);
  comm->remove_comm(mycomm1);
# else
  for(int i=0;i<svecsz[0];i++)
    rbuf[i]=sbuf[i];
# endif
  // unpacking received requests and decoding them
  for(int i=0;i<np;i++){
    for(int j=0;j<rvecsz[i];j++){
      vector<arg_t> nonlocal;
      irequest_t &req=rbuf[rdispl[i]+j];
      arg_t iarg(req.arg);
      interp_form_t coeff=cont->create_interpolation(iarg,0,(req.level>0 ? NULL:&nonlocal)); 
      if(req.dest==myrank) // this rank, interpolated value should be recorded directly to destination
        iplocal.next_record(make_pair(coeff,make_pair(req.ibuf,req.ind)));
      else{ // other rank, put in send buffer for the given rank
        commtbl[req.dest]->destpack.next_record(make_pair(req.ibuf,req.ind));
        // other rank, interpolated value should be recorded to input buffer
        ip.next_record(make_pair(coeff,make_pair(2*req.dest+1,(int)commtbl[req.dest]->outsz)));
        commtbl[req.dest]->outsz++;
      }
      // nonlocal, level 1 secondary request
      for(size_t k=0;k<nonlocal.size();k++){
        int rank1=cont->test_rank(nonlocal[k].pos, container_t::MODIFY);
        if (rank1<0) 
          continue;
        rankvec[rank1].push_back(irequest_t(req.dest,req.ibuf,req.ind,1,arg_t(nonlocal[k].c,nonlocal[k].pos,nonlocal[k].ftype,nonlocal[k].fdir)));
      }
    }
  }
  if(rbuf) delete [] rbuf;
  if(sbuf) delete [] sbuf;
  return 1;
}

template<class arg_tt, class container_tt, class reg_tester_t, class interp_packer_t, class index_t>
int InterpTransfer<arg_tt, container_tt, reg_tester_t, interp_packer_t, index_t>::prepare_transfers(){
  start(tid.tr_prep);
  rvecsz=new int[np];
  svecsz=new int[np];
  rdispl=new int[np];
  sdispl=new int[np];
  // phase I: communicate primary interpolation requests
  transfer_requests();
  // phase II: communicate secondary interpolation requests
  transfer_requests();
  //rankvec.clear(); // why not?
  complete_record();
  for(int i=0;i<np;i++){
    commtbl[i]->destpack.end_record();
  }
  // finished interpolation collection

  // phase III: collect destination data, organize receives
  // communicating output buffer sizes
  for(int i=0;i<np;i++)
    svecsz[i]=(int)commtbl[i]->outsz;
# ifdef USE_MPI
  // creating specific communicator
  //MPI_Comm mycomm1;
  void *mycomm1;
  //MPI_Comm_dup(comm->get_all(), &mycomm1);
  comm->dup_all(&mycomm1);
  //MPI_Alltoall(svecsz,1,MPI_INT,rvecsz,1, MPI_INT,mycomm1);
  comm->int_alltoall(svecsz,1,rvecsz,1,mycomm1);
# else
  rvecsz[0]=svecsz[0];
# endif
  // preparing request send-recv buffers
  int ssz=0, rsz=0;
  for(int i=0;i<np;i++){
    rdispl[i]=rsz;
    sdispl[i]=ssz;
    rsz+=rvecsz[i];
    ssz+=svecsz[i];
  }
  intpair_t *rbuf= rsz? new intpair_t[rsz] : NULL;
  intpair_t *sbuf= ssz? new intpair_t[ssz] : NULL;
  // unpacking and copying buffer destinations
  for(int i=0;i<np;i++){
    typename comm_t::packer_t::iterator it=commtbl[i]->destpack.begin(), e=commtbl[i]->destpack.end();
    for(int j=0;it!=e;++it,++j){
      sbuf[sdispl[i]+j]=*it;
    }
    commtbl[i]->destpack.clear();
  }
  // communicating interpolation requests
# ifdef USE_MPI
  //MPI_Alltoallv(sbuf,svecsz,sdispl, intpair_mpit, rbuf, rvecsz, rdispl, intpair_mpit,mycomm1);
  comm->alltoallv(sbuf,svecsz,sdispl, intpair_mpit, rbuf, rvecsz, rdispl, intpair_mpit,mycomm1);
  //MPI_Comm_free(&mycomm1);
  comm->remove_comm(mycomm1);
# else
  for(int i=0;i<svecsz[0];i++)
    rbuf[i]=sbuf[i];
# endif
  // packing buffer destinations
  for(int i=0;i<np;i++){
    commtbl[i]->destpack.start_record();
    for(int j=0;j<rvecsz[i];j++){
      commtbl[i]->destpack.next_record(rbuf[rdispl[i]+j]);        
    }
    commtbl[i]->insz=rvecsz[i];
    // finish record and allocate buffers
  }
  if(rbuf) delete [] rbuf;
  if(sbuf) delete [] sbuf;
   
  init_mpi_table();
  // binding output buffers
  for(int i=0;i<np;i++){
    substitute_buffer(i,commtbl[i]->outbuf);
  }
  delete [] rvecsz;
  delete [] svecsz;
  delete [] rdispl;
  delete [] sdispl;
  stop(tid.tr_prep);
  return 1;
}

template<class arg_tt, class container_tt, class reg_tester_t, class interp_packer_t, class index_t>
void InterpTransfer<arg_tt, container_tt, reg_tester_t, interp_packer_t, index_t>::start_transfers(){
  start(tid.tr_gthr);
  // prepare output buffers
  fill_buffers(begin(),end());
  stop(tid.tr_gthr);
  // start non-blocking transfers
  if(mpidump){
    fprintf(logfile,"starting transfers...");
    fflush(logfile);
  }
# ifdef USE_MPI
//#if 0 // test with turned off MPI (for balanced domain decomposition)
  start(tid.tr_send);
  if(nsend)
    //MPI_Startall(nsend,srequests);
    comm->startall(nsend,srequests);
  if(nrecv)
    //MPI_Startall(nrecv,rrequests);
    comm->startall(nrecv,rrequests);
  stop(tid.tr_send);
# endif
  if(mpidump){
    fprintf(logfile,"ok\n");
    fflush(logfile);
  }
}

template<class arg_tt, class container_tt, class reg_tester_t, class interp_packer_t, class index_t>
//void InterpTransfer<container_tt, reg_tester_t, interp_packer_t, index_t>::complete_transfers(int *shifts){
void InterpTransfer<arg_tt, container_tt, reg_tester_t, interp_packer_t, index_t>::complete_transfers(){
  if(mpidump){
    fprintf(logfile,"completing transfers...");
    fflush(logfile);
  }
# ifdef USE_MPI
//#if 0 // test with turned off MPI (for balanced domain decomposition)
//  int flag=0;
    
  //MPI_Testall((int)commtbl.size(),rrequests,&flag,rstatuses);

  if(nsend){
    //MPI_Testall(nsend,srequests,&flag,sstatuses);
    if(mpidump){
      fprintf(logfile,"tst_send...");
      fflush(logfile);
    }
    start(tid.tr_recv);
    //MPI_Waitall(nsend,srequests,sstatuses);
    comm->waitall(nsend,srequests,sstatuses);
    stop(tid.tr_recv);
  }
  if(nrecv){
    //MPI_Testall(nrecv,rrequests,&flag,rstatuses);
    if(mpidump){
      fprintf(logfile,"tst_recv...");
      fflush(logfile);
    }
    start(tid.tr_recv);
    //MPI_Waitall(nrecv,rrequests,rstatuses);
    comm->waitall(nrecv,rrequests,rstatuses);
    stop(tid.tr_recv);
  }
  # endif
  if(mpidump){
    fprintf(logfile,"ok\n");
    fflush(logfile);
  }
  start(tid.tr_sctr);
  for(size_t i=0;i<commtbl.size();i++){
//    if (shifts)
//      commtbl[i]->transfer_dest(buff, shifts);
//    else
      commtbl[i]->transfer_dest(buff);
  }
  stop(tid.tr_sctr);
}

template<class arg_tt, class container_tt, class reg_tester_t, class interp_packer_t, class index_t>
int InterpTransfer<arg_tt, container_tt, reg_tester_t, interp_packer_t, index_t>::idcount=0;

#endif
