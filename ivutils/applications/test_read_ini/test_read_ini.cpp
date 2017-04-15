#include "read_ini.h"

// example how to read ini file ini.dat

int main(){

  Ini I;
  if(I.read("..\\ini.dat")<0)
    return -1; // ini file is not found or incorrect

  double center[3], radius, double_size;
  I.parse_key("center", "ddd", NULL, center, NULL);
  I.parse_key("radius", "d", NULL, &radius, NULL);
  I.parse_key("double_size", "d", NULL, &double_size, NULL);
  int num;
  I.parse_key("num", "i", &num, NULL, NULL);
  string fname;
  I.parse_key("file_name", "s_s", NULL, NULL, &fname);

  return 0;
}
