# ifndef _ALG_PARSER_H
# define _ALG_PARSER_H

/// @file alg_parser.h \brief Interface for parsing algebraic expressions

#include <string>
#include <map>

using namespace std;

/// This class allows to parse and calculate algebraic expression 
/// with + - * /, brackets () and some constants and functions.
/// constants and functions should be recorded in corresponding maps (see below).
/// f. e. "sqrt(5+cos(3*pi))"
/// Main function of the class is operator() which calculates given algebraic expression.
/// Maps for constants and functions are filled in constructor 
/// with some common contans and functions (as pi or cos).
/// User can add some other funtions to these maps, f. e. functions["sqrt"]=sqrt
class AlgebraicParser{

  // finds a position of the first character symb in the string arg (starting from the end of the string).
  // if symbol is not found returns -1.
  // characters inside brackets () are ingnored
  // if brackets are used incorrectly (f.e. "5)(-38*4("), returns -2.
  // skips '+' or '-' symbols in case they are part of the exponential number format 
  int find_symbol(const string &arg, char symb);

  // str is given string, calculated value is recorded to val,
  // if str is correct returns 1, otherwise returns -1.
  // calls C function atof, or calculates constant or function value using content of corresponding maps.
  // this function is used in operator()
  int atof_map(const string &str, double *val);

public:

  typedef double (*zerofun)();
  typedef double (*fun)(double);
  typedef double (*bifun)(double, double);
  typedef double (*trifun)(double, double, double);

  // maps of names for constants and functions with 0, 1, 2 and 3 arguments
  // this names will be recognized during estimation of algebraic expression
  map<string, double> constants; 
  map<string, zerofun> zerofunctions;
  map<string, fun> functions;
  map<string, bifun> bifunctions;
  map<string, trifun> trifunctions;

  AlgebraicParser();

  /** 
  calculate the algebraic expression str and record the result in val.
  if str is correct returns 1, otherwise returns -1.  
  str can have + - * /, brackets () 
  and some constants and functions which are recorded in corresponding maps (see above).
  f.e., str="sqrt(5+cos(3*pi))".
  T could be float or double. */
  int operator()(const string &str, double *val);
};

# endif
