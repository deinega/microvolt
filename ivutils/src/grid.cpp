/// \file \brief Non-template function definitions for  grid.h 
#include "grid.h"
#include "linsolv.hpp"

Vector_3 get_iterator_dx(const UniformGrid<Vector_3, 1>::iterator &it){
  return it.dx();
}

vec_type volume_integral_grid_t::operator()(const Vector_3 &p1, const Vector_3 &p2) const{
  vec_type vol=Box(p1,p2).Volume();
  if(!grid.GetPtr())
    return vol;
  Vector_3 dx=grid.get_dx();
  int n[3]; // discretization of the volume
  Vector_3 st; // start space point
  for(int i=0;i<3;i++){
    n[i]=int(floor(acdiv(p2[i]-p1[i],dx[i])))+1;
    st[i]=p1[i]+((p2[i]-p1[i])-(n[i]-1)*dx[i])/2;
  }
  vec_type sum=0;
  Vector_3 pos;
  int ic[3];
  for(ic[0]=0;ic[0]<n[0];ic[0]++){
    pos[0]=st[0]+ic[0]*dx[0];
    for(ic[1]=0;ic[1]<n[1];ic[1]++){
      pos[1]=st[1]+ic[1]*dx[1];
      for(ic[2]=0;ic[2]<n[2];ic[2]++){
        pos[2]=st[2]+ic[2]*dx[2];
        sum+=grid.Interpolate(pos);
      }
    }
  }
  return vol*sum/(n[0]*n[1]*n[2]);
}

class mesh_generation_t{
  vec_type q;
  int n;
  vec_type dx1; 
  vec_type dx2;
  vec_type dist;

  int find_q(vec_type *y);
  int corr_q(vec_type *y);

public:

  mesh_generation_t(vec_type dx1_, vec_type dx2_, vec_type dist_):dx1(dx1_),dx2(dx2_),dist(dist_){
    q = dx1<dx2 ? 2 : .5;
  }

  int mesh_generation(vector<vec_type> &dx);
};

/* 
sum_1^n dx1*q^i = dist
dx1*q^{n+1} = dx2

( if dx1 <-> dx2 then q -> 1/q )

sum_1^n q^i = q(1-q^n)/(1-q) = dist/dx1
q^{n+1} = dx2/dx1

n+1 = log_q(dx2/dx1)

q(1-q^{log_q(dx2/dx1)-1})/(1-q) = dist/dx1
*/

int mesh_generation_t::find_q(vec_type *y){
  vec_type dx2_dx1=dx2/dx1;
  vec_type dist_dx1=dist/dx1;

  vec_type log_q = log(dx2_dx1)/log(q);
  vec_type pow_q=pow(q,log_q-1);
  vec_type sum=q*(1-pow_q)/(1-q);
  *y=sum-dist_dx1;

  return 1;
}

int mesh_generation_t::corr_q(vec_type *y){
  vec_type dist_dx1=dist/dx1;
  vec_type pow_q=pow(q,n);
  vec_type sum=q*(1-pow_q)/(1-q);
  *y=sum-dist_dx1;

  return 1;
}

int mesh_generation_t::mesh_generation(vector<vec_type> &dx){

  dx.clear();

  newton_raphson<vec_type,linear_solver<vec_type> > nr(1);
  if(nr(&q,this,&mesh_generation_t::find_q,10)<0)
    return -1;

  vec_type dx2_dx1=dx2/dx1;
  vec_type nval = log(dx2_dx1)/log(q)-1;
  n=int(floor(nval+.5));

  if(nr(&q,this,&mesh_generation_t::corr_q,10)<0)
    return -1;

  vec_type qq=q;
  vec_type d=0;
  for(int i=0;i<n;i++,qq*=q){
    vec_type step=dx1*qq;
    d+=step;
    dx.push_back(step);
  }

  return 1;
}

int generate_1d_mesh(Vector_3 rough, vector<Vector_3> dence, vector<vec_type> &mesh){

  if(rough[1]-rough[0]<=rough[2] || rough[2]<=0)
    return -1; // incorrect mesh definition

  vector<int> inter; // approximate size of intermediate region in rough mesh steps

  for(size_t i=0;i<dence.size();i++){ // cycle along dence meshes inserted to the rough mesh
    if(dence[i][1]-dence[i][0]<dence[i][2] || dence[i][2]<=0 || dence[i][2]>=rough[2])
      return -1; // incorrect mesh definition
    inter.push_back(int(ceil(acdiv(rough[2],dence[i][2]))));
  }

  // high resolution meshes should be at the distance at least 2*(inter+1) from each other
  // otherwise current algorithm will not work
  for(size_t i=0;i<dence.size();i++){
    if(dence[i][0]<rough[0] || dence[i][1]>rough[1])
      return -1; // incorrect mesh definition
    for(size_t j=i+1;j<dence.size();j++){
      vec_type dist=rough[2]*(inter[i]+inter[j]+2); // minimal distance between two high resolution meshes
      if(dence[i][0]<dence[j][0]){
        if(dence[i][1]+dist>dence[j][0])
          return -1; // two high resolution meshes are too close to each other
      }
      else{
        if(dence[i][0]-dist<dence[j][1])
          return -1; // two high resolution meshes are too close to each other
      }
    }
  }

  int rough_n=int(ceil(acdiv(rough[1]-rough[0],rough[2]))); // mesh size of rough mesh
  rough[2]=(rough[1]-rough[0])/rough_n; // correcting step of rough mesh
  
  // nodes from rough mesh that will not be included to generated nonuniform mesh
  vector<pair<int,int> > excl;
  
  for(size_t i=0;i<dence.size();i++){ // cycle along dence meshes inserted to the rough mesh

    int dence_n=int(ceil(acdiv(dence[i][1]-dence[i][0],dence[i][2]))); // mesh size of high resolution mesh
    dence[i][2]=(dence[i][1]-dence[i][0])/dence_n; // correcting step of high resolution mesh

    for(int j=0;j<=dence_n;j++)
      mesh.push_back(dence[i][0]+j*dence[i][2]); // generated mesh will have internal nodes of high resolution mesh

    // left and right rough mesh node index.
    // between these indices we don't use rough mesh
    int ex[2]; 

    // generation of intermediate steps between rough mesh and inserted dence mesh
    // we consider first space at the left to the dence mesh,
    // and then space at the right
    for(int lr=0;lr<2;lr++){
      int sign = lr ? 1 : -1;
      vec_type x=dence[i][lr]; // x1 (left) or x2 (right)
      vec_type div=acdiv(x,rough[2]); // coordinate of x at the rough mesh
      ex[lr] = int(lr ? ceil(div) : floor(div)); // shift to the neighbour node of the rough mesh
      ex[lr]+=sign*inter[i]; // shift on inter rough mesh nodes
      div = lr ? ceil(div)-div : div-floor(div);
      vec_type dist=(inter[i]+div)*rough[2]; // distance of intermediate region between dence and rough meshes

      // we cover this distance by noniniform mesh step
      mesh_generation_t mg(dence[i][2],rough[2],dist);
      vector<vec_type> dx;
      if(mg.mesh_generation(dx)<0) // generation
        return -1; // not successful

      vec_type dxc=0;

      vec_type d=0; // position of node of intermediate region with nonuniform step
      // check if intermediate region is within rough mesh
      for(size_t j=0;j<dx.size();j++){
        d+=dx[j];
        vec_type pos=x+sign*d;
        if(!lr && pos<rough[0]){ // some part of intermediate region is at the left to the first node of the rough mesh
          dist=dence[i][0]-rough[0];
          ex[0]=0;
          dxc=dx[j];
          break;
        }
        else if(lr && pos >rough[1]){ // some part of intermediate region is at the right to the last node of the rough mesh
          dist=-(dence[i][1]-rough[1]);
          ex[1]=rough_n;
          dxc=dx[j];
          break;
        }
      }
      // some part of intermediate region is outside rough mesh
      if(dxc){ // new nonuniform mesh generation for shortened intermediate region
        mesh_generation_t mg(dence[i][2],dxc,dist);
        if(mg.mesh_generation(dx)<0)
          return -1;
      }

      d=0;
      for(size_t j=0;j<dx.size();j++){
        d+=dx[j];
        vec_type pos=x+sign*d;
        mesh.push_back(pos); // adding nonuniform part to generated mesh
      }
//      accomp(d,dist); // just check
    }
    excl.push_back(make_pair(ex[0],ex[1]));
  }
  for(int i=0;i<=rough_n;i++){ // listing nodes of rough mesh
    int ex=0;
    for(size_t ie=0;ie<excl.size();ie++){
      if(i>=excl[ie].first && i<=excl[ie].second){
        ex=1;
        break;
      }
    }
    if(ex)
      continue; // already covered by high resolution mesh or nonuniform mesh

    mesh.push_back(rough[0]+i*rough[2]); // adding nodes of rough meesh
  }
  sort(mesh.begin(),mesh.end());

  if (mesh.size() > 1) {
    // Remove coinciding points
    vector<vec_type>::iterator it1 = mesh.begin();
    for(vector<vec_type>::iterator it2 = it1+1; it2 != mesh.end(); it2++) {
      if (!accomp(*it1, *it2))
        *(++it1) = *it2;
    }
    mesh.erase(it1+1, mesh.end());
  }

  return (int)mesh.size();
}
