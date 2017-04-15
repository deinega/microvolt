/// \file \brief Template function definitions for  grid.h

#include "grid.h"

template <class value_tt, size_t il, class index_t>
index_t UniformGrid<value_tt, il, index_t>::iterator::operator-(const typename UniformGrid<value_tt, il, index_t>::iterator& other) const {
  index_t sz[3], df[3];
  for(int i=0;i<3;i++){
    sz[i]=parent->ien[i]-parent->ist[i]+1;
    df[i]=ind[i]-other.ind[i];
  }
  return (df[0]*sz[1]+df[1])*sz[2]+df[2];
}

template <class value_tt, size_t il, class index_t>
typename UniformGrid<value_tt, il, index_t>::iterator& UniformGrid<value_tt, il, index_t>::iterator::operator++(){ //prefix
  int ci=2;
  for(;ci>=0;ci--){
    ind[ci]++;
    if(ind[ci]>parent->ien[ci]){
      ind[ci]=parent->ist[ci];
    }
    else break;
  }
  if(ci<0){
    ie=1; // end of sequence
  }
  return *this;
}

template <class value_tt, size_t il, class index_t>
typename UniformGrid<value_tt, il, index_t>::iterator& UniformGrid<value_tt, il, index_t>::iterator::operator+=(index_t incr){
  while(incr>0 && !ie){
    int ci=2;
    for(;ci>=0;ci--){
      if(ci==2){ // only first index lowers increment
        ind[2]+=incr;
        if(ind[2]>parent->ien[2]){
          incr=ind[2]-parent->ien[2];
          ind[2]=parent->ist[2];
        }
        else{
          incr=0;
          break;
        }
      }
      else{
        ind[ci]++;
        if(ind[ci]>parent->ien[ci])ind[ci]=parent->ist[ci];
        else break;
      }
    }
    if(ci<0){
      ie=1;  
    }
  }
  return *this;
}

template <class value_tt, size_t il, class index_t>
void UniformGrid<value_tt, il, index_t>::init(const Vector_3 &v1, const Vector_3&v2, const int dir, const index_t *sz, const index_t *start, const index_t *end){
  for(int i=0;i<3;i++){
    if(v2[i]<v1[i])
      drc|=1<<i;
  }
  b.init(v1,v2);
  Vector_3 bsz=b.GetSize();
  for(int i=0;i<3;i++){
    ngr[i]=sz[i];
    if(ngr[i]<2){
      ngr[i]=1;
      dx[i]=bsz[i];
      pref[i]=(v1[i]+v2[i])/2.;
    }
    else{
      dx[i]=bsz[i]/(ngr[i]-1);
      pref[i]=b.get_p1()[i];
    }
    nst[i]=0;
    nen[i]=ngr[i]-1;
  }

  SetIteratorRange(start,end); 
  SetInterpolationRange(start,end);
  SetMemoryRange(start,end);
  int dir0=abs(dir)-1;
  if(dir&&!bsz[dir0]){
    vec_type s=1;
    for(int i=1;i<3;i++){
      int dir1=(dir0+i)%3;
      s*=bsz[dir1]/(sz[dir1]-(sz[dir1]>1));
    }
    ds[dir0]=(dir>0?1:-1)*s;
  }
}

template <class value_tt, size_t il, class index_t>
void UniformGrid<value_tt, il, index_t>::set_range(index_t *st, index_t *en, index_t *gr, const index_t *start, const index_t *end){
  for(int i=0;i<3;i++){
    if(start) st[i]=start[i];
    else st[i]=nst[i];
    if(end)en[i]=end[i];
    else en[i]=nen[i];
    gr[i]=en[i]-st[i]+1;
  }
}

template <class value_tt, size_t il, class index_t>
void UniformGrid<value_tt, il, index_t>::get_1d_coeff(const Vector_3 &place, index_t *ind, vec_type *c, bool nonlocal) const{
  Vector_3 p=ch_drc(place)-b.get_p1();
  for(int i=0;i<3;i++){
    ind[i]=(index_t)floor(acdiv(p[i],dx[i]));
    if(nonlocal)
      c[i]=1.-fabs(acdiv(p[i],dx[i])-(vec_type)ind[i]);
    else{
      if(ind[i]<nst[i]-1)ind[i]=nst[i]-1;
      if(ind[i]>nen[i])ind[i]=nen[i];

      if(ind[i]<nst[i])
        c[i]=0; 
      else if(ind[i]>=nen[i])
        c[i]=1;
      else
        c[i]=1.-fabs(acdiv(p[i],dx[i])-(vec_type)ind[i]);
    }
  }
}

template <class value_tt, size_t il, class index_t>
int UniformGrid<value_tt, il, index_t>::GetCoeff(const Vector_3 &place, index_t *int_ind, vec_type *values, int force_external, vector<Vector_3> *nonlocal, bool all_nonlocal) const{
  if(!force_external && !b.TestPoint(place))return -1; // not in the region
  index_t ind[3];
  vec_type c[3];
  get_1d_coeff(place,ind,c,nonlocal!=NULL);
/*  for(int i=0;i<3;i++){
    ind[i]=(int)floor(acdiv(p[i],dx[i]));
    if(nonlocal)
      c[i]=1.-fabs(acdiv(p[i],dx[i])-(vec_type)ind[i]);
    else{
      if(ind[i]<nst[i]-1)ind[i]=nst[i]-1;
      if(ind[i]>nen[i])ind[i]=nen[i];

      if(ind[i]<nst[i])
        c[i]=0; 
      else if(ind[i]>=nen[i])
        c[i]=1;
      else
        c[i]=1.-fabs(acdiv(p[i],dx[i])-(vec_type)ind[i]);
    }
  }*/
  index_t cv[4]={0,0,0,0};
  vec_type csum=0; // sum of all coefficients, should be equal to one
  int l=0, inloc=0;
  while(!cv[3]){
    // calculating the coefficient
    vec_type ct=1;
    for(int i=0;i<3;i++){
      if(!cv[i])ct*=c[i];
      else ct*=1.-c[i];
    }
    // filling
    if(fabs(ct)>VEC_ZERO){
      // check for regulary supported index
      if((all_nonlocal==false) && check_interpolation_ind(ind[0]+cv[0],ind[1]+cv[1],ind[2]+cv[2]))
        int_ind[l]=pack_ind(ind[0]+cv[0],ind[1]+cv[1],ind[2]+cv[2]);
      else if(nonlocal){
        inloc++;
        int_ind[l]=-inloc;
        nonlocal->push_back(Position(ind[0]+cv[0],ind[1]+cv[1],ind[2]+cv[2]));
      }
      else
        throw "GetCoeff"; // impossible
      values[l++]=ct;
      csum+=ct;
    }
    // index increment
    for(int ci=0;ci<4;ci++){
      cv[ci]++;
      if(cv[ci]>=2){
        cv[ci]=0;
      }
      else break;
    }
  }
//  if(l>0 && fabs(csum-1.)>1e-7) // recalibrating
//    for(i=0;i<l;i++)values[i]/=csum; // TEMPORARY COMMENT
  return l;
}

template <class value_tt, size_t il, class index_t>
int UniformGrid<value_tt, il, index_t>::test_local(const Vector_3 &place) const{
  Vector_3 p=place-b.get_p1();
  for(int i=0;i<3;i++){
    vec_type val=acdiv(p[i],dx[i]);
    if(val<nst[i] || val>nen[i]) 
      return 0;
  }
  return 1;
}

template <class value_tt, size_t il, class index_t>
typename UniformGrid<value_tt, il, index_t>::value_t UniformGrid<value_tt, il, index_t>::Interpolate(const Vector_3 &place, int force_external) const {
  vec_type values[8];
  value_tt sum=value_tt();
  index_t ind[8];
  int ni=GetCoeff(place,ind,values,force_external);
  for(int i=0;i<ni;i++)sum+=ptr[ind[i]]*values[i];
  return (value_tt)sum;
}

template <class value_tt, size_t il, class index_t>
int UniformGrid<value_tt, il, index_t>::AddValue(const Vector_3 &place, value_tt value, int force_external){
  vec_type values[8], sum=0;
  index_t ind[8];
  int ni=GetCoeff(place,ind,values,force_external);
  for(int i=0;i<ni;i++)ptr[ind[i]]+=value*values[i];
  return ni;
}

template <class value_tt, class index_t>
int InterpBox<value_tt, index_t>::GetCoeff(const Vector_3 &place, index_t *int_ind, vec_type *values, int force_external) const{
  if(!force_external && !Box::TestPoint(place))return -1; // not in the region
  Vector_3 p=place-Box::get_p1();
  int i, l=0; index_t ind[3];
  vec_type c[3], t[3];
  for(i=0;i<3;i++){
    c[i]= this->sz[i] > VEC_ZERO ?  1.-p[i]/sz[i] : 1.;
    if(t[i]<0)c[i]=0.;
    else if(c[i]>1.)c[i]=1.;
  }
  index_t cv[4]={0,0,0,0};
  while(!cv[3]){
    // calculating the coefficient
    vec_type ct=1.;
    for(i=0;i<3;i++){
      if(!cv[i])ct*=c[i];
      else ct*=1.-c[i];
    }
    // filling
    if(fabs(ct)>VEC_ZERO){
      int_ind[l]=pack_ind(cv[0],cv[1],cv[2]);
      values[l++]=ct;
    }
    // index increment
    for(int ci=0;ci<4;ci++){
      cv[ci]++;
      if(cv[ci]>=2){
        cv[ci]=0;
      }
      else break;
    }
  }
  return l;
}

template <class value_tt, size_t ili, class index_t>
void NonUniformGrid<value_tt, ili, index_t>::get_1d_coeff(const Vector_3 &p, index_t *ind, vec_type *c, bool nonlocal) const{
  for(int i=0;i<3;i++){

    typename vector<value_tt>::const_iterator lb=lower_bound(x[i].begin(),x[i].end(),p[i]);
    ind[i]=lb-x[i].begin()-1;

    if(ind[i]<nst[i])
      c[i]=0; 
    else if(ind[i]>=nen[i])
      c[i]=1;
    else if(p[i]==*lb){
      ind[i]++;
      c[i]=1;
    }
    else{
      value_tt neigh[2]={x[i][ind[i]],x[i][ind[i]+1]};
      value_tt dneigh=neigh[1]-neigh[0];
      c[i]=acdiv(neigh[1]-p[i],dneigh);
    }
  }
}
