# include "seqpack.h"
# include "conio.h"

//# define INT
//# define ITER
//# define DATA
//# define PAIR
//# define ITERDATA
# define  GROUP

//template <>
//int get_group_count<int, equal_to<int>, vector_set<int, equal_to<int> >(group_pack<<int, equal_to<int>, vector_set<int, equal_to<int> >>::iterator &it);

typedef double valtype;
typedef pair<valtype, valtype> valtype2;
typedef pair<valtype, valtype2 > valtype3;
typedef pair<valtype2, valtype2> valtype4;

typedef pair_pack<valtype, valtype> valtype2_pack;
//typedef pair_pack<valtype, valtype, data_unpacked<valtype> > valtype11_pack;
//typedef pair_pack<valtype2, valtype2, valtype11_pack, valtype2_pack> valtype13_pack;
typedef pair_pack<valtype2, valtype2, valtype2_pack, valtype2_pack> valtype13_pack;


/// function to test integer sequence packing
int test_seqpack(int len,  double seqprob, int seqlenav, int spec_val, double specprob);

/// reads the input to int_pack from the string in the form of a list: 1-10,3,4,3-65
/// the integers must be positive inside a range [imin,imax]
/// the ranges num1-num2 are by default with increment 1 (-1)
/// delimiters (,-) are specified in the delims string
/// special value "all" means the whole range [imin,imax]
/// @return:  -1 wrong arguments,
///           -2 read (format) error
///           -3 outside the range
///            1 OK
int add_ranges(int_pack &packer, const string &spec, int imin=0, int imax=1000, bool finalize=true,const string &delims=",-");


void main() {
  int i=0;
  int res;
# ifdef INT
  test_seqpack(10u0, 0.5, 5, 1, 0.5);
  getch();
# endif
# ifdef ITER
  {
  int N=100;
  int *ptr(new int[N]);
  for (i=0; i<N; i++)
    ptr[i]=i;

  int len=6;
  int *unpacked(new int[len]);
  unpacked[0]=0;
  unpacked[1]=1;
  unpacked[2]=1;
  unpacked[3]=1;
  unpacked[4]=1;
  unpacked[5]=2;

  iter_pack<int *> IP(ptr);

  IP.start_record();
  for (i=0; i<len; i++) {
    IP.next_record(unpacked[i]);
  }
  IP.end_record();

  printf("%d records, packed size is %d, effective packed size is %d\n",len,IP.size(),IP.packed_size());
  printf("Comparing:\n");
  iter_pack<int *>::iterator it=IP.begin(), end=IP.end();
  int cur_ind(0);
  for (i=0; i<len; i++) {
    if(!(it!=end)){
      printf("Unexpected end of sequence at %d\n",i);
      res=0;
      break;
    }
    cur_ind+=unpacked[i];
    int vu=ptr[cur_ind];
    int vp=*(*it);
    if(vu!=vp){
      printf("! %d: %d  %d\n",i,vu,vp);
      res=0;
    }
    it++;
  }
  if(!res){
    printf("Sequences differ.\n");
  }
  else{
    printf("Sequences match.\n");
  }
  getch();
  delete []ptr;
  delete []unpacked;
  }
# endif
# ifdef DATA
  {
  int N=100;
  int *ptr(new int[N]);
  for (i=0; i<N; i++)
    ptr[i]=i;

  int len=6;
  int *unpacked(new int[len]);
  unpacked[0]=0;
  unpacked[1]=1;
  unpacked[2]=1;
  unpacked[3]=1;
  unpacked[4]=1;
  unpacked[5]=2;

  data_pack<int, equal_to<int>, vector_set<int>> D;
  D.start_record();
  for (i=0; i<len-2; i++) {
    D.next_record(unpacked[i]);
  }
  D.next_group(unpacked+4, 1);
  D.next_group(unpacked+5);
  D.end_record();

  printf("%d records, packed size is %d, effective packed size is %d\n",len,D.size(),D.packed_size());
  printf("Comparing:\n");
  data_pack<int, equal_to<int>, vector_set<int> >::iterator it=D.begin(), end=D.end();
  for (i=0; i<len; i++) {
    if(!(it!=end)){
      printf("Unexpected end of sequence at %d\n",i);
      res=0;
      break;
    }
    int vu=ptr[unpacked[i]];
    int vp=*it;
    get_group_count(it);
    if(vu!=vp){
      printf("! %d: %d  %d\n",i,vu,vp);
      res=0;
    }
    it++;
  }
  if(!res){
    printf("Sequences differ.\n");
  }
  else{
    printf("Sequences match.\n");
  }
  getch();
  delete []ptr;
  delete []unpacked;

  }
# endif
# ifdef PAIR 
  {

  vector<int> vct;
  vct.push_back(1);
  vct.push_back(2);
  for (vector<int>::iterator iv=vct.begin(); iv!=vct.end(); iv++) {
    (*iv)++;
  }

  int N=100;
  int *ptr(new int[N]);
  for (i=0; i<N; i++)
    ptr[i]=i;

  int len=6;
  int *unpacked(new int[len]);
  unpacked[0]=0;
  unpacked[1]=1;
  unpacked[2]=1;
  unpacked[3]=1;
  unpacked[4]=1;
  unpacked[5]=2;

  pair_pack<int, int> P;
  
  P.start_record();
  for (i=0; i<len-2; i++) {
    P.next_record(unpacked[i], unpacked[i]);
  }
  pair<int, int> gr[2];
  gr[0].first=gr[0].second=unpacked[4];
  gr[1].first=gr[1].second=unpacked[5];

  P.next_group(gr,2);
  P.end_record();

  printf("%d records, packed size is %d, effective packed size is %d\n",len,P.size(),P.packed_size());
  printf("Comparing:\n");
  pair_pack<int, int>::iterator it=P.begin(), end=P.end();
  for (i=0; i<len; i++) {
    if(!(it!=end)){
      printf("Unexpected end of sequence at %d\n",i);
      res=0;
      break;
    }
    int vu=ptr[unpacked[i]];
    pair<int, int> vp=*it;
    if((vu!=vp.first)||(vu!=vp.second)){
      printf("! %d: %d  %d\n",i,vu,vp.first);
      res=0;
    }
    it++;
  }
  if(!res){
    printf("Sequences differ.\n");
  }
  else{
    printf("Sequences match.\n");
  }
  getch();
  delete []ptr;
  delete []unpacked;

  }
# endif
# ifdef ITERDATA
  {
  int N=100;
  int *ptr(new int[N]);
  for (i=0; i<N; i++)
    ptr[i]=i;

  int len=6;
  int *unpacked(new int[len]);
  unpacked[0]=0;
  unpacked[1]=1;
  unpacked[2]=1;
  unpacked[3]=1;
  unpacked[4]=1;
  unpacked[5]=2;

  iter_data_pack<int *, int> ID(ptr);

  ID.start_record();
  for (i=0; i<len; i++) {
    ID.next_record(unpacked[i], unpacked[i]);
  }
  ID.end_record();

  printf("%d records, packed size is %d, effective packed size is %d\n",len,ID.size(),ID.packed_size());
  printf("Comparing:\n");
  iter_data_pack<int *, int>::iterator it=ID.begin(), end=ID.end();
  int cur_ind(0);
  for (i=0; i<len; i++) {
    if(!(it!=end)){
      printf("Unexpected end of sequence at %d\n",i);
      res=0;
      break;
    }
    cur_ind+=unpacked[i];
    int vu1=ptr[cur_ind];
    int vu2=ptr[unpacked[i]];
    pair<int *, int> vp=*it;
    if((vu1!=*(vp.first))||(vu2!=vp.second)){
      printf("! %d: %d  %d\n",i,vu1,*(vp.first));
      res=0;
    }
    it++;
  }
  if(!res){
    printf("Sequences differ.\n");
  }
  else{
    printf("Sequences match.\n");
  }
  getch();
  delete []ptr;
  delete []unpacked;

  }
# endif

# ifdef GROUP
  {
  int N=100;
  /// must return 1
  i=get_group_count(&N);
  int *ptr(new int[N]);
  for (i=0; i<N; i++)
    ptr[i]=i;

  int len=6;
  int *unpacked(new int[len]);
  unpacked[0]=0;
  unpacked[1]=1;
  unpacked[2]=1;
  unpacked[3]=1;
  unpacked[4]=1;
  unpacked[5]=2;

  group_pack<int, equal_to<int>, vector_set<int, equal_to<int> >, subgroup_test<equal_to<int> > > G;
//template<class T, class comp_pr=equal_to<T>, class set_t=vector_set<T, comp_pr> , class group_pr=subgroup_test<typename set_t::key_compare> >
  G.start_record();
  G.next_record(unpacked[0]);
  G.next_group(&unpacked[1], 2);
  G.next_group(&unpacked[3], 2);
  G.next_group(&unpacked[5], 1);
  G.end_record();

//  group_pack<int>::iterator it=G.begin();
//  *it;



  printf("%d records, packed size is %d, effective packed size is %d\n",len,G.size(),G.packed_size());
  printf("Comparing:\n");
  group_pack<int>::iterator it=G.begin(), end=G.end();
  for (i=0; i<len; i++) {
    if(!(it!=end)){
      printf("Unexpected end of sequence at %d\n",i);
      res=0;
      break;
    }
    int vu=ptr[unpacked[i]];
    int vp=*it;
    get_group_count(it);
    if(vu!=vp){
      printf("! %d: %d  %d\n",i,vu,vp);
      res=0;
    }
    it++;
  }
  if(!res){
    printf("Sequences differ.\n");
  }
  else{
    printf("Sequences match.\n");
  }
  getch();
  delete []ptr;
  delete []unpacked;

  }
# endif

}




int test_seqpack(int len,  double seqprob, int seqlenav, int spec_val, double specprob){
  vector<int> unpacked(len);
  int_pack packed(spec_val);

  // generating the random sequence
  int i=0;
  // starting value
  int vcur=12, interval, seq=0;
  packed.start_record();
  do{
    if(i!=0){
      if(seq<=0){
        // sequence or not
        seq=(((double)rand())/RAND_MAX < seqprob ? 1 : 0);
        if(seq){
          // the length of the sequence
          seq=(int)(seqlenav*(double)rand()/RAND_MAX);
        }
        // special interval or not
        int spec=(((double)rand())/RAND_MAX < specprob ? 1 : 0);
        if(!spec){
          interval=(int)(100.*((double)rand())/RAND_MAX);
        }
        else interval=spec_val;
      }
      else seq--;
      int v=vcur+interval;
      if(v<0)break; // end of integer range, can not be packed further
      vcur=v;
    }
    unpacked[i]=vcur;
    if(i==16329){
      i++;
      i--;
    }
    packed.next_record(vcur);
    i++;
  }while(i<len);
  packed.end_record();

  printf("%d records, packed size is %d, effective packed size is %d\n",i,(int)packed.size(),int(packed.packed_size()/sizeof(int)));
  printf("Comparing:\n");
  int_pack::iterator it=packed.begin(), end=packed.end();
  int res=1;
  for(i=0;i<(int)unpacked.size();i++){
    if(!(it!=end)){
      printf("Unexpected end of sequence at %d\n",i);
      res=0;
      break;
    }
    int vu=unpacked[i];
    int vp=*it;
    //printf("%d: %d  %d\n",i,vu,vp);
    if(vu!=vp){
      printf("! %d: %d  %d\n",i,vu,vp);
      res=0;
    }
    if(i==16328){
      i++;
      i--;
    }
    it++;
  }
  if(!res){
    printf("Sequences differ.\n");
  }
  else{
    printf("Sequences match.\n");
  }
  return res;
}

int add_ranges(int_pack &packer, const string &spec, int imin, int imax, bool finalize,const string &delims){
  int res=1;
  if(spec=="all"){
    for(int i=imin;i<=imax;i++)
      packer.next_record(i);
  }
  else{
    if(delims.size()<2)
      return -1;
    char comma[]=",", range[]="-";
    comma[0]=delims[0];
    range[0]=delims[1];
    char *buff=new char[spec.size()+1];
    strncpy(buff,spec.c_str(),spec.size()+1);
    
    
    char *str1=strtok(buff,comma);
    while(str1){
      int n1, n2;
      if(sscanf(str1,"%d",&n1)!=1){
        res= -2;
        break;
      }
      char *p=strstr(str1,range);
      if(p){
        *p=0;
        if(sscanf(p+1,"%d",&n2)!=1){
          res=-2;
          break;
        }
        *p=range[0];
      }
      else
        n2=n1;
      if(n1<0 || n2<0){
        res= -3;
        break;
      }
      int dn= n1==n2? 0 : n1<n2 ? 1 : -1; // this is to be replaced by reading optional step
      n2+=dn;
      do{
        if(imin<=n1 && n1<=imax) // range control
          packer.next_record(n1);
        n1+=dn;
      }while(n1!=n2);
      str1=strtok(NULL,comma);
    }
    delete [] buff;
  }
  if(res>0 && finalize)
    packer.end_record();
  return res; 
}
