#ifndef READ_INI
#define READ_INI

/// @file read_ini.h \brief Classes to read special ini files

#include <vector>
#include <map>
#include <string>

using namespace std;

/* 
This class allows to parse simple ini files in the following format:
key: value1 value2 value3

One can express values using algebraic expression (f.e. 1+1)
Symbol '$' allows to substitute value of some defined key (f.e. 1+$key)
Condition operator if-elseif-else-end can be used to skip some parts of ini file
Comments are started with semicolon ;

Example of ini file can be found in test project

Keys and values of processed ini file are saved to map using function read,
and then they can be extracted using function parse_key
*/
class Ini{

  static const int MAX_LEVEL_REFERENCES=100; // maximal level on nested references $
  static const int CH_MAX=4000; // maximal string length for key or values

  typedef map<string, string, less<string> > map_t;
  typedef map<string, string, less<string> >::iterator map_it;

  map_t mapstring; // map to store keys and values

  // read one line in the format "key: value1 value2" and record key and value to the mapstring
  int read_line(const string &s, string &strkey, string &strvalue);

  // takes string val and substitue all strings after $ to values of corresponding keys
  // level is recursion level of the function:
  // if function is called from the same function the level is increasing
  int process_value_dollars(const string &val, string &result, int &level);

  // call process_value_dollars for value of the given key
  int process_key_dollars(const string &key, string &result, int &level);

  // remove part of string starting from ";"
  void remove_comments(string &str);

  // auxiliary funtion used in if_test
  // analyze string str from start to end (if end=-1 to str.length()
  // and record token which is word (in this case function returns 2) 
  // or combination of symbols .:; (in this case function returns 1)\
  // if there is no token returns 0
  // if token is too long returns -1
  // current position after reading token is recorded to i
  int parse_next_substr(const string &str, string &token, int &i, int start=0, int end=-1);

  // test if str is condition operator and return the code of condition:
  // 1 - if, 2- elseif, 3 - else, 4 - end (otherwise returns 0)
  int test_condition(const string &str);

  // check the condition (str is smth like 1==1 or 1!=1)
  // returns if condition is true, 0 if condition is false, and -1 if something wrong
  int parse_equal_evaluate(const string &str);

  // read the string, change ifs (vector of nested conditions)
  // returns 1 if all nested condition are satisfied
  // 0 if some nested condition is not satisfied or this is condition line
  // -1 if using condition operators is incorrect
  int if_test(const string &str, vector<int> &ifs);

//public: // for compatibility with experiments project: PLEASE CHECK!
// WHICH EXPERIMENTS PROJECT DO YOU MEAN ???

public:

  // parse ini file with filename and record the result to mapstring
  // some of parameters of ini file can be redefined using n_ini_str and ini_str
  // in this case ini_str is array of n_in_str strings of the ini format:
  // key: value1 valu2 value3
  // using this option can be useful to change ini file by comand line arguments
  int read(const char *filename, const int n_ini_str=0, const string *ini_str=NULL);

  // see comments to parse_key
  static int parse_string(const string &value, const char *format, int *ints, double *doubles, string *strings, int *optionals=NULL);

  /* extract values for given key in integer, double or string format
   format is string of characters i (for integer), d (for double), s (for string)
   and _ (separate ordinary parameters from optional parameters)

   for example you have:
   key: 1.5 0 2 name1 name2 33
   
   int ints[3];
   dobule prm;
   string str[2];
   int opt;
   parse_key("key", "diiss_ii",ints,&prm,str,&opt);
   
   result will be int[3]={0,2,33}, prm=1.5, str[2]={"name1","name2"}, opt=1 (one optional parameters is used)

   see test example for better understanding

   returns 1 if everything ok or -1 otherwise
   if key is absent then returns -1 and shows error message (if forse=true) or just returns 0 (otherwise)
   */
  int parse_key(const string &key, const char *format, int *ints, double *doubles, string *strings, int *optionals=NULL, bool force=true);

  // return string in format key: value
  string get_ini_line(const string &key);

/*  // for compatibility with experiments project: PLEASE CHECK!
// WHICH EXPERIMENTS PROJECT DO YOU MEAN ???
  static int catoi(const char *str, int *val){
    *val=atoi(str);
    return -(*val==0 && (strcmp(str,"0")));
  }*/
};

#endif
