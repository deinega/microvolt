#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#define _USE_MATH_DEFINES
#endif

#include <cmath>
#include <limits>
#include <stdlib.h>
#include "alg_parser.h"
#include "string_utils.h"

AlgebraicParser::AlgebraicParser() {
  // fill maps with constants and functions
  // one can extend given list of constants or functions
  constants["pi"]=M_PI;
  constants["inf"]=numeric_limits<double>::max();

  functions["sqrt"]=sqrt, functions["sin"]=sin, functions["cos"]=cos, functions["tan"]=tan;
  functions["asin"]=asin, functions["acos"]=acos, functions["atan"]=atan, functions["abs"]=fabs;
  functions["floor"]=floor, functions["ceil"]=ceil;
    
  bifunctions["pow"]=pow;
};

int AlgebraicParser::find_symbol(const string &arg, char symb) {
  int br=0; //����������� ������
  int len=(int)arg.length();
  for (int i=len-1; i>=0; i--){
    if (arg[i]==')')
      br++;
    else if (arg[i]=='(') {
      br--;
      if (br<0)
        return -2;
    }
    else if(!br && arg[i]==symb ){ // !br ��������, ��� ������ ������ ��� ������
      if(symb=='-' || symb=='+'){// checking for exponential format
        if(i>=2 && i<len-1 && isdigit(arg[i+1]) && isdigit(arg[i-2]) && tolower(arg[i-1])=='e') // this is exp format
          continue;
      }
      return i;
    }
  }
  if (br) return -2;
  return -1;
}

int AlgebraicParser::atof_map(const string &str, double *val) {
  *val=atof(str.c_str());
  if (*val==0 && str!="0") {
    if(str.length()==0)
      return -1; // empty string
    if(str[str.length()-1]==')') { // function
      int i;
      for (i=0; i<(int)str.length(); i++) {
        if (str[i]=='(') break;
      }
      if (i==(int)str.length()) return -1;
      string func=SUBSTR(str, 0, i);
      if (zerofunctions.find(func)!=zerofunctions.end()) {
        if (i!=(int)str.length()-2) return -1;
        *val=zerofunctions[func]();
      }
      else if (functions.find(func)!=functions.end()) {
        double val2;
        if (operator()(SUBSTR(str, i+1, (int)str.length()-(i+1)-1), &val2)<1) return -1;
        *val=functions[func](val2);
      }
      else if (bifunctions.find(func)!=bifunctions.end()) {
        int j;
        int br=0;
        for (j=i+1; j<(int)str.length(); j++) {
          if (str[j]=='(')
            br++;
          else if (str[j]==')') {
            if (br==0)
              return -1;
            br--;
          }
          else if (str[j]==',' && br==0) break;
        }
        if (j>=(int)str.length()-2 || j==i+1) return -1;
        double val2, val3;
        if (operator()(SUBSTR(str, i+1, j-(i+1)), &val2)<1) return -1;
        if (operator()(SUBSTR(str, j+1, (int)str.length()-(j+1)-1), &val3)<1) return -1;
        *val=bifunctions[func](val2, val3);
      }
      else if (trifunctions.find(func)!=trifunctions.end()) {
        int j, j1=-2;
        int br=0;
        for (j=i+1; j<(int)str.length(); j++) {
          if (str[j]=='(')
            br++;
          else if (str[j]==')') {
            if (br==0)
              return -1;
            br--;
          }
          else if (str[j]==',' && br==0) {
            if (j1<0)
              j1=j;
            else
              break;
          }
        }
        if (j>=(int)str.length()-2 || j1==i+1 || j-j1<2) return -1;
        double val1,val2, val3;
        if (operator()(SUBSTR(str, i+1, j1-(i+1)), &val1)<1) return -1;
        if (operator()(SUBSTR(str, j1+1, j-(j1+1)), &val2)<1) return -1;
        if (operator()(SUBSTR(str, j+1, (int)str.length()-(j+1)-1), &val3)<1) return -1;
        *val=trifunctions[func](val1,val2, val3);
      }
      else
        return -1;
    }
    else { // constant
      if (constants.find(str)==constants.end()) return -1;
      *val=constants[str];
    }
  }
  return 1;
}

int AlgebraicParser::operator()(const string &str, double *val){

 // coding of arythmetical operations
  const int extOP_PLUS=1;
  const int extOP_MINUS=-1;
  const int extOP_MUL=2;
  const int extOP_DIV=-2;

  /* ������� ��������� ������ �� ���������, ��������� + � -, 
  ��������� ��� �������� �������� ������� �����������.
  ����� ��������� ���������, ��������� * � /.
  ������ ������� ��������� ������ �� ��� ���������.
  ���������� �������� �������������� � ������� ������������ ������ ���� �� �������.
  */

  if (str.length()==0)
    return -2;

  int direct=0; // 1 - ������ ��������� ���� �����, 0 - ������ ��������� ���� ���������

  //���� � ����� ������ ������ + ��� - ��� ������ ()
  int op=extOP_MINUS; // �������� �������������� �������� (��������� �������� ����������� � ������� define)
  int c1=find_symbol(str,'-');
  int c2=find_symbol(str,'+');
  if (c1<-1 || c2<-1) return -2; // ������ � ������ ����������� �����������

  if(c1<0 && c2<0){ // none of +/- found
    //���� � ����� ������ ������ * ��� / ��� ������ ()
    op=extOP_MUL;
    c1=find_symbol(str,'*');
    c2=find_symbol(str,'/');
    if (c1<-1 || c2<-1) return -2; // ������ � ������ ����������� �����������
    direct=1; // ������������, ��� ������ ��������� ���� �����
  }

  if(c2>c1){ // ������ � ����� ����������� ��������, ��������������� ���������� �������� (+ ��� / )
    op=-op;
    c1=c2;
  }
  double val1=0, val2=0; // � ��� ���������� ����� �������� �������� ��� ������� � ������� ����������

  if (c1<0 && c2<0) { // �� ������ ������� + - * / � ������ �� ���������
    if(str[0]=='(' && str[str.length()-1]==')') { // ������ ���� ��������� � ������� (...)
      string token1=SUBSTR(str, 1, (int)str.length()-2); // �������, ��� ������ ������
      return operator()(token1, val); // ������� ��� � ����������
    }
    // ������ ���� ������ �����. � ���� ������ val1 ����� ��������, val2 ����� ���� ������,
    // � �������� ����� ���������� (� ���������� �� ������� 1*val2=val2)
    val1=1;
  }
  else {
    string token1=SUBSTR(str, 0, c1);
    trim(token1);
    if(token1.length()==0){ // ����� ��������� ����� �������� ����� ����
      if(op==extOP_MINUS)val1=0.; // ������� ����� (��������, -5)
      else return -1; // ������
    }
    else{
      if(operator()(token1,&val1)<0) // ������� val1
        return -1;
    }
  }

//  string token2=str.substr(c1+1, str.length()-(c1+1));
  string token2=my_substr(str, c1+1, (int)str.length()-(c1+1));
  if (token2.length()==0)
    return -1; // ������
  trim(token2);
  if(token2[0]=='(' && token2[token2.length()-1]==')'){ // ������ ��������� ���� ��������� � �������, � �� �����
    direct=0;
  }
  if(direct){ // ������ ��������� ���� ��������������� �����
    if(atof_map(token2, &val2)<1)return -1;
  }
  else{
    if(operator()(token2,&val2)<0)
      return -1;
  }
  // �� ��������� ��� ��������� val1 � val2 � ���������� �������������� ��������  
  switch(op){
    case extOP_PLUS:
      *val=val1+val2;
    break;
    case extOP_MINUS:
      *val=val1-val2;
    break;
    case extOP_DIV:
      *val=val1/val2;
    break;
    case extOP_MUL:
      *val=val1*val2;
    break;
  }
  return 1;
}
