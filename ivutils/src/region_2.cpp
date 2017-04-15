#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#define _USE_MATH_DEFINES
#endif

#include "region_2.hpp"

Polygon_2 *GetRegularPolygon(const Vector_2 &center, vec_type d, int n){
  Polygon_2 *poly = new Polygon_2();
  vec_type dfi=2*M_PI/n, fi1=0;//M_PI/2;
  for (int i=0; i<n; i++)
    poly->add(Vector_2(center[0]+d*cos(fi1+i*dfi), center[1]+d*sin(fi1+i*dfi)));
  return poly;
}
