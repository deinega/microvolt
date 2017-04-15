#include "linsolv.h"

/*
System of linear equations:

 x - y = 1
 x + y = 2

 solution: x = 1.5, y = 0.5
*/

#define PARDISO

int main(){

#ifndef PARDISO
  linear_solver<double> solv;
#else
  pardiso_solver<double> solv;
  solv.set_sparce(2); // set assumed maximal elements number at each row
#endif
  // 2 - system size, 1|2 - two bits for managing RHS column and matrix
  solv.init(2,1|2);
  solv.start_record();
  // assign matrix elements values
  solv.set_m(0,0,1);
  solv.set_m(0,1,1);
  solv.set_m(1,0,-1);
  solv.set_m(1,1,1);
  // assign RHS column value
  solv.v(0)=1;
  solv.v(1)=2;
  solv.end_record();
  double x[2];
  // solution is recorded in x
  int res=solv(x); // res is 1 if everything ok

  return 0;
}
