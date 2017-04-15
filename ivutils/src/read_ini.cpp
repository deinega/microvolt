#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "string_utils.h"
#include "read_ini.h"
#include "alg_parser.h"
#include "logexc.h"

/// some useful math functions

/// random number in [0,1]
double random_0_1() {
  return double(rand())/double(RAND_MAX);
}

//double my_time() {
//  return double(time(NULL));
//}

double my_round(double arg) {
  return floor(arg+0.5);
}

//double condition(double arg1, double arg2, double arg3) {
//  return arg1 ? arg2 : arg3;
//}

int Ini::read(const char *filename, const int n_ini_str, const string *ini_str){

//  stdfile_logger log("",1,stdout,stderr,vblALLBAD,vblFATAL,1);

  FILE *f=NULL;
  if((f=fopen(filename, "r"))==NULL) 
    return LOGERR(-1,fmt("Input file is not found\n"),0);
  string key, val;
  map<string, bool> ini_strings; // line and flag if this line is already found in ini file
  mapstring.clear();
  vector<int> ifs; // vector for nested conditions

  for (int i=0; i<n_ini_str; i++) { // read ini_str
    if (read_line(ini_str[i], key, val)) { // extract key and value
      ini_strings[key]=false;
      remove_comments(val);
      mapstring[key]=val;
    }
  }

  for (int strnum=0; !feof(f); strnum++) { // read ini file
    char buf[CH_MAX];
    if (fgets(buf, CH_MAX, f)==NULL) // end of the file
      break;
    int len=(int)strlen(buf);
    if (len==CH_MAX)
      return LOGERR(-1,fmt("String %d in input file is too large\n", strnum),0);
    if (len==0)
      continue;
    if (buf[len-1]=='\n')
      len--;
    string sbuf(buf, len);
    remove_comments(sbuf);

    int if_res=if_test(sbuf, ifs);
    if (if_res<0)
      return LOGERR(-1,fmt("Incorrect 'if-else' using (string %d)\n", strnum),0);
    if (if_res==0) // this is one of the condition operator or condition is not satisfied
      continue;

    if (read_line(sbuf, key, val)) { // extract key and value
      if (ini_strings.find(key)!=ini_strings.end()) {
        if (ini_strings[key])
          return LOGERR(-1,fmt("Key '%s' is assigned twice\n", key.c_str()),0);
        else
          ini_strings[key]=true;
      }
      else if (mapstring.find(key)!=mapstring.end())
        return LOGERR(-1,fmt("Key '%s' is assigned twice\n", key.c_str()),0);
      else
        mapstring[key]=val; // adding key and value to the map
    }
  }
  if (ifs.size())
    return LOGERR(-1,fmt("Incorrect 'if-else' using: 'if' number is larger then 'end' one\n"),0);
  fclose(f);


  for (map_it it=mapstring.begin(), end=mapstring.end(); it!=end; it++) {
    int level=0;
    if(process_key_dollars((*it).first, (*it).second, level)<0) // substituting references $ values
      return -1;
  }

  return 1;
}

int Ini::read_line(const string &s, string &strkey, string &strval) {
  char buff[CH_MAX+1]; // character buffer to record
  int cnt=0; // current position in buffer
  bool if_value=false; // false if key is recorded; true if value is recorded

  for(size_t i=0;i<s.length();i++){

    if (s[i]==':' && !if_value) {
      if (cnt!=0) { // record key
        buff[cnt]=0;
        strkey=string(buff);
        trim(strkey);
        if (!strkey.length())
          return 0; // key is empty
      }
      else // key is empty
        return 0;
      if_value=true;
      cnt=0;
    }
    else if (s[i]==';') // comment
      break;
    else if (cnt>=CH_MAX) 
      return 0; // to long key or value string
    else
      buff[cnt++]=s[i]; // filling the buffer
  }

  if (cnt!=0 && if_value) { // record value
    buff[cnt]=0;
    strval=string(buff);
    trim(strval);
    return 1;
  }

  return 0; // value is not defined
}

int Ini::process_value_dollars(const string &val, string &result, int &level) {
  const char *value=val.c_str();
  char new_value[CH_MAX];
  char ref[CH_MAX];
  int len((int)strlen(value)), nv(0), r(0);
  bool if_dollar(false), if_arithmetic_dollar(false), if_braces(false);

  for (int i=0; i<=len; i++) {
    if (value[i]=='$' || value[i]=='&') {
      if (if_dollar) 
        return -1;
      if_dollar=true;
      if_arithmetic_dollar=(value[i]=='&');
    }
    else if (!if_dollar) {
      if (i==len)
        new_value[nv]=0;
      else {
        if (nv>=CH_MAX) 
          return -2;
        new_value[nv++]=value[i];
      }
    }
    else if (value[i]=='{') {
      if (value[i-1]=='$' || value[i]=='&')
        if_braces=true;
      else
        return -1;
    }
    else if (i==len || value[i]=='}' || value[i]==' ' || value[i]==char(9) || value[i]==',' || value[i]=='(' || value[i]==')' || value[i]=='+' || value[i]=='-' || value[i]=='*' || value[i]=='/' || value[i]=='\\') {
      if (if_braces) {
        if (i==len || value[i]!='}')
          return -1;
        else
          if_braces=false;
      }
      else if (value[i]=='}')
        return -1;
      ref[r]=0;
      level++;
      string dollar_value;
      if(process_key_dollars(string(ref), dollar_value, level)<0)
        return -1;
      if (if_arithmetic_dollar)
        dollar_value=string("(")+dollar_value+string(")");
      level--;
      if (!(i==len || value[i]=='}')) {
        char symb[2]={value[i],0};
        dollar_value+=string(symb);
      }
      int dl=(int)dollar_value.length();
      if (dl>CH_MAX-nv-1) 
        return -2;
      new_value[nv]=0;
      strcat(new_value, dollar_value.c_str());
      nv+=dl;
      r=0;
      if_dollar=false;
    }
    else {
      if (r>=CH_MAX) 
        return -3;
      ref[r++]=value[i];
    }
  }

  result=string(new_value);
  return 1;
}

int Ini::process_key_dollars(const string &key, string &result, int &level) {
  if (level>MAX_LEVEL_REFERENCES) 
    return LOGERR(-1,fmt("Too deep indirect reference, possible references loop (check for example key %s)\n", key.c_str()),0);
  if (mapstring.find(key)==mapstring.end())
    return LOGERR(-1,fmt("Reference %s is absent\n", key.c_str()),0);
  int res=process_value_dollars(mapstring[key], result, level);
  if (res==-1)
    return LOGERR(-1,fmt("Incorrect reference in %s\n", key.c_str()),0);
  else if (res==-2)
    return LOGERR(-1,fmt("Too large definition of %s\n", key.c_str()),0);
  else if (res==-3)
    return LOGERR(-1,fmt("Too large reference in %s\n", key.c_str()),0);
  return res;
}


void Ini::remove_comments(string &str) {
  for (size_t i=0; i<str.length(); i++) {
    if (str[i]==';') {
      str=SUBSTR(str, 0, (int)i); // remove everything after comment
      return;
    }
  }
}

int Ini::parse_next_substr(const string &str, string &token, int &i, int start, int end) {

  const int MAX_INI_WORD=100;
  const int MAX_INI_SYMBOL=10;

  char word[MAX_INI_WORD];
  char symb[MAX_INI_SYMBOL];
  // current position in buffer
  int iw=0;
  int is=0;
  // status for buffer (write or not wtite)
  bool rw=false;
  bool rs=false;
  int state=0; // current status (0 - space, 1 - punctiation, 2 - word)
  int pr_state=0; // previous status

  if (end==-1)
    end=(int)str.length();

  for (i=start; i<end; i++) {
    if (str[i]==' ' || str[i]==char(9)) 
      state=0;
    else if (str[i]==',' || str[i]==':' || str[i]==';')
      state=1;
    else
      state=2;

    if (state==0) { // space
      if (pr_state==1) {
        rs=true;
        break;
      }
      else if (pr_state==2) {
        rw=true;
        break;
      }
    }
    else if (state==1) { // punctuation
      if (pr_state==2) {
        rw=true;
        break;
      }
      else {
        if (is==MAX_INI_SYMBOL-1) // fill buffer for punctiation
          return -1;
        symb[is]=str[i];
        is++;
      }
    }
    else { // word
      if (pr_state==1) {
        rs=true;
        break;
      }
      else {
        if (iw==MAX_INI_WORD-1) // fill buffer for word
          return -1;
        word[iw]=str[i];
        iw++;
      }
    }
    pr_state=state;
  } // for
  if (rs || (i==end && state==1)) {
    symb[is]=0;
    token=string(symb); // next token is punctiation
    return 1;
  }
  else if (rw || (i==end && state==2)) {
    word[iw]=0;
    token=string(word); // next token is word
    return 2;
  }
  return 0;
}

int Ini::test_condition(const string &str) {
  const char* cstr=str.c_str();
  if (strcmp(cstr, "if")==0)
    return 1;
  else if (strcmp(cstr, "elseif")==0)
    return 2;
  else if (strcmp(cstr, "else")==0)
    return 3;
  else if (strcmp(cstr, "end")==0)
    return 4;
  else return 0;
}

int Ini::parse_equal_evaluate(const string &str) {
  int action; // 1 - equality, 0 - inequality

  int c=int(str.find(string("==")));
  if (c==0 || c>=int(str.length())-2)
    return -1;
  else if (c>0)
    action=1;
  else {
    c=int(str.find(string("!=")));
    if (c==0 || c>=int(str.length())-2)
      return -1;
    else if (c>0)
      action=0;
    else
      return -1;
  }
  string token1=SUBSTR(str,0, c);
  string token2=SUBSTR(str,c+1+1, (int)str.length()-(c+1+1));
  trim(token1);
  trim(token2);
  if (token1.length()==0 || token2.length()==0)
    return -1;
  if (int(token1.find(string("==")))>=0)
    return -1;
  if (int(token2.find(string("==")))>=0)
    return -1;
  if (int(token1.find(string("!=")))>=0)
    return -1;
  if (int(token2.find(string("!=")))>=0)
    return -1;

  int level=0;
  string rtoken1, rtoken2;
  int res;
  res=process_value_dollars(token1, rtoken1, level);
  if (res==1)
    res=process_value_dollars(token2, rtoken2, level);
  if (res==-1)
    return LOGERR(-1,fmt("Incorrect reference in 'if-else'\n"),0);
  else if (res==-2)
    return LOGERR(-1,fmt("Too large condition in 'if-else'\n"),0);
  else if (res==-3)
    return LOGERR(-1,fmt("Too large reference in 'if-else'\n"),0);

  return (rtoken1.compare(rtoken2)==0) ? action : !action;
}

int Ini::if_test(const string &str, vector<int> &ifs) {
  string token;
  int i;
  int next=parse_next_substr(str, token, i);
  if(next<0)return -1;
  int key_word=(next==2) ? test_condition(token) : 0;
  int n=(int)ifs.size();

  if (key_word==0) { // no condition operator
    for (int j=0; j<n; j++) {
      if ((ifs[j]&1)==0) // check all nested condition operators
        return 0; // some of condition is not satisfied
    }
    return 1; // all possible condition are satisfied
  }
  else if (key_word==4) { // end
    if (n==0)
      return -1; // using end without using if before
    ifs.pop_back(); // exit from current condition 
  }
  else if (key_word==1) { // if
    string rest=SUBSTR(str,i,(int)str.length()-i); // reading condition argument (like a==b)
    int val=parse_equal_evaluate(rest);
    if(val<0)
      return -1; // condition argument is wrong
    ifs.push_back(val | (val<<1));
  }
  else { // elseif or else
    if (n==0)
      return -1; // using elseif or else without using if
    int &cur=ifs[n-1]; // result of previous if or elseif
    if (cur&4)
      return -1; // using elseif or else after else
    if (key_word==3) // else
      cur|=4;

    if (cur&2) {
      cur&=~1; // prevous condition was satisfied, this condition is not satisfied
    }
    else { // previous condition was not satisfied
      if (key_word==3) // else
        cur|=3;
      else {
        string rest=SUBSTR(str, i, (int)str.length()-i); // elseif
        int val=parse_equal_evaluate(rest);
        if(val<0)
          return -1; // condition argument is wrong
        cur|=val | (val<<1);
      }
    }
  }
  return 0;
}

int Ini::parse_string(const string &value, const char *format, int *ints, double *doubles, string *strings, int *optionals) {

  AlgebraicParser alg; // to parse algbraic expressions

  alg.functions["round"]=my_round;
  alg.zerofunctions["random"]=random_0_1;
//  alg.zerofunctions["time"]=my_time;
//  alg.trifunctions["condition"]=condition;

  int i=0, len=value.length(), frmlen=strlen(format), cnt=0, num=0, result=0;
  int is_equal=0; // a (0) || a = (1) || a = b (2)
  bool is_comma=false; // if after comma
  bool is_wait=false; // if waiting for the value
  bool is_optionals=false;
  char ch[CH_MAX+1]; // for current word recording

  if (optionals)
    *optionals=0;

  for (; i<=len; i++) {
    if (value[i]=='=') {
      if (is_equal) return -3; // a = = (1) || a = b = (2)
      if (is_comma) return -6; // a, =
      cnt=0;
      is_wait=false;
      is_equal=1; // first '='
    }
    else if (value[i]==',') {
      if (is_comma) return -4; // a, ,
      if (is_equal==1) return -5; // a =,
      is_comma=true;
      is_wait=false;
      if (cnt>=CH_MAX) return num+1; // number of too large substring
      else ch[cnt++]=value[i];
    }
    else if (value[i]==' ' || value[i]==char(9)) {
      if (is_comma) continue;
      if (is_equal==1) continue;
      else is_wait=true; // now waiting for any significant letter and record previous word
    }
    else {
      if ((i==len)||(is_wait)) {
        if ((i==len)&&((is_equal==1)||(is_comma))) return -1; // a = || a,
        // record
        if (cnt) {
          ch[cnt]=0;


          if (optionals) {
            if (*optionals)
              (*optionals)++;
          }

          if (format[num]=='_') {
            num++;
            is_optionals=true;
            if (optionals)
              *optionals=1;
          }

          if (format[num]=='i') {
//            *ints=atoi(ch);
//            result=((*ints==0)&&(*ch!='0')) ? 0 : 1;
            double aux;
            result=alg(string(ch), &aux);
            *ints=(int)aux;
            ints++;
          }
          else if (format[num]=='d') {
//            *doubles=atof(ch);
//            result=((*doubles==0)&&(*ch!='0')) ? 0 : 1;
            result=alg(string(ch), doubles);
            doubles++;
          }
          else if (format[num]=='s') {
            *strings=string(ch);
            result=1;
            strings++;
          }
          num++;
          if (result<1) return num; // number of incorrect substring
          if (num==frmlen) break; // MAY BE check for superfluos letters
          cnt=0;
        }
        is_wait=false;
        if (is_equal==2) is_equal=0;
      }
      else {
        if (is_comma) is_comma=false; // a, b
        if (is_equal==1) is_equal++; // a = b
      }
      if (cnt>=CH_MAX) 
        return num+1; // number of too large substring
      else ch[cnt++]=value[i];
    }
  }
  if (num<frmlen) { // substring number are less than require. Possible error with optional parameters
    if (!is_optionals && format[num]!='_') 
      return -1;
//    if (!optionals) return -1;
//    else if(!(*optionals) && format[num]!='_') return -1;
  }
  if (i<len)
    return -2;
  return 0;
/*
  int i(0), len(strlen(value)), frmlen(strlen(format)), cnt(0), num(0), result;
  char ch[CH_MAX+1];

  for (; i<=len; i++) {
    if ((value[i]==' ')||(i==len)) {
      if (cnt) {
        ch[cnt]=0;
        if (format[num]=='i') {
          *ints=atoi(ch);
          result=((*ints==0)&&(*ch!='0')) ? 0 : 1;
          ints++;
        }
        else if (format[num]=='d') {
          *doubles=atof(ch);
          result=((*doubles==0)&&(*ch!='0')) ? 0 : 1;
          doubles++;
        }
        else if (format[num]=='s') {
          *strings=string(ch);
          result=1;
          strings++;
        }
        num++;
        if (!result) return num; // number of incorrect substring
        if (num==frmlen) break;
        cnt=0;
      }
    }
    else {
      if (cnt>=CH_MAX) return num+1; // number of too large substring
      else ch[cnt++]=value[i];
    }
  }
  if (num<frmlen) return -1; // substring number are less than require
  return 0;
*/
}

int Ini::parse_key(const string &key, const char *format, int *ints, double *doubles, string *strings, int *optionals, bool force) {
  if(!(mapstring.find(key)!=mapstring.end())) {
    if(force)return LOGERR(-1,fmt("Description '%s' is absent\n", key.c_str()),0);
    return 0;
  }
  int result=parse_string(mapstring[key], format, ints, doubles, strings, optionals);
  if (result>0) return LOGERR(-1,fmt("Incorrect parameter %d definition in %s\n", result, key.c_str()),0);
  else if (result==-1) return LOGERR(-1,fmt("Too few parameters in %s\n", key.c_str()),0);
  else if (result==-2) return LOGERR(-1,fmt("Too many parameters in %s\n", key.c_str()),0);
  else if (result<-2) return LOGERR(-1,fmt("Desription '%s' is incorrect\n", key.c_str()),0);
  return 1;
}

string Ini::get_ini_line(const string &key){
  map_t::iterator it=mapstring.find(key);
  if(it!=mapstring.end()){
    return key+": "+it->second;
  }
  return "";
}
