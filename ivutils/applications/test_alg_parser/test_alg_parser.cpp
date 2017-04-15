#include "alg_parser.h"

int main(){

  AlgebraicParser AP;
  double val; // here calculated value will be recorded
  int res=AP("sqrt(5+cos(3*pi))", &val); // if res=1 expression is correct
  res=AP("1e-5", &val); // if res=1 expression is correct

  return 0;
}
