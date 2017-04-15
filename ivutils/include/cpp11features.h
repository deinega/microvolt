#ifndef CPP11_FEATURES_H
#define CPP11_FEATURES_H

/// @file cpp11features.h \brief Definitions of macros for C++11 features available

/** \if EN
  The header file contains macros which are enabled if the respective features of C++11 standard
  (ISO/IEC 14882:2011) are available. The definitions are based on the following pages:
    Table of status of C++11 language features in compilers:
      http://wiki.apache.org/stdcxx/C++0xCompilerSupport
    List of compiler version macros
      http://sourceforge.net/p/predef/wiki/Compilers/
 
  Supported compilers:
    Clang 3.0,
    GCC 4.0,
    Intel C++ Compiler 9.0,
    MS Visual C++ 8.0 (Visual Studio 2005),
  and their newer versions

  Note that Clang, GCC and Intel C++ Compiler need special command line options in order to enable
  C++11 features (-std=c++11 or -std=c++0x, see compiler documentation). MSVC always enables all
  implemented C++11 features.

  In order to disable all HAS_xxx macros listed below and turn off usage of C++11 features in this
  library one can globally define macro DISABLE_CPP11 (it doesn't turn off C++11 support, but C++11
  features won't be used in the library, which may result in a certain increase of calculation time
  and/or memory usage).

  The list is significantly shortened and simplified in order to support features which are
  relevant to this library only.

  The following macros are defined if the respective features are available (the section numbers
  below refer to the respective sections of ISO/IEC 14882:2011(E) International Standard:
  Information technology � Programming languages � C++):

    HAS_ARRAY - defined if header <array> (which contains std::array class) is available
      See sections 23.3.1 and 23.3.2

    HAS_DEFAULTED - defined if defaulted and deleted class members are available
      See sections 8.4.2 and 8.4.3 or http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2007/n2346.htm

    HAS_FRIENDEXT - defined if extended friend declarations are supported
      (e.g., friend T; where T is a template parameter).
      See section 11.3 or http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2005/n1791.pdf

    HAS_INITLIST - defined if initializer lists for arrays and structs can be used in expressions
      See sections 8.5.4 and 18.9 or http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2008/n2672.htm

    HAS_LAMBDAS - defined if lambda expressions and std::function are supported.
      See sections 5.1.2 and 20.8.11. The macro may contain one of the following values depending
      on the feature version implemented:
        09 - lambdas v.0.9: stateless lambdas, no conversion of a closure object to an object of an
          "ordinary" type. See http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2008/n2550.pdf
        10 - lambdas v.1.0: mutable lambdas were added.
          See http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2008/n2658.pdf
        11 - lambdas v.1.1: std::function<> template class was added, a number of corner points was
          clarified. See http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2009/n2927.pdf
      Such detalization may appear to be redundant, and further revisions of this may shorten it to
      checking for the final (1.1) version. If so, the version number in the macro definition will
      be omitted.

    HAS_RVALUE_REFS - defined if rvalue reference type (Type&&) is defined. This doesn't imply that
      member functions reference qualifiers (for definitions like T::f()&) are defined as well (if
      these will be needed, then a separate macro should be added to this file).
      See section 8.3.2 or http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2006/n2118.html
      There were several revisions of this feature which clarified some points. These points may or
      may not be relevant to our code. For now, the macro is defined if any implementation of
      rvalue references is available, but this may be changed in a way similar to lambdas.

    HAS_SMART_PTR - defined if std::shared_ptr, std::weak_ptr and std::unique_ptr are defined.
      See section 20.7.

    HAS_STDFUNCTION - defined if std::function<> is supported and lambda expressions can be cast
      to objects of this class. The macro is always defined if HAS_LAMBDAS >= 11, is always
      undefined if HAS_LAMBDAS is undefined and may be defined otherwise (some compilers don't
      support lambdas v.1.1 completely, but support std::function<>).

    HAS_VARTEMPLATES - defined if templates with arbitrary number of parameters are available.
      See section 14.5.3 or http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2007/n2242.pdf and
      http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2008/n2555.pdf
      
  The following macros help writing code which uses C++11 features if they are available:

    SHARED_OR_AUTO(type) - is replaced with std::shared_ptr<type> if it's available, and with
      std::auto_ptr<type> otherwise. Note that std::auto_ptr has no shared ownership and can't
      be used in containers.

    UNIQUE_OR_AUTO(type) - is replaced with std::unique_ptr<type> if it's available, and with
      std::auto_ptr<type> otherwise. Note that std::auto_ptr can't be used in containers.

    SHARED_OR_RAWPTR(type) - is replaced with std::shared_ptr<type> if it's available, and with
      type* otherwise. Make sure that the pointed object will be deleted exactly once if it is raw.

    UNIQUE_OR_RAWPTR(type) - is replaced with std::unique_ptr<type> if it's available, and with
      type* otherwise. Make sure that the pointed object will be deleted exactly once if it is raw.

    MAKE_SHARED_OR_AUTO(type, args) - is replaced with (std::make_shared<type> args), if
      std::shared_ptr is supported, and with (std::auto_ptr<type>(new type args)) otherwise. Note
      that the arguments for the constructor should be specified in additional braces, e. g.:
      MAKE_SHARED_OR_AUTO(Type, (arg1, arg2, arg3));
      std::shared_ptr can be initialized with a raw pointer too, but the created object consumes
      slightly more memory since the memory is allocated for the control block and the pointed
      object separately.

    MAKE_SHARED_OR_RAW(type, args) - are replaced with (std::make_shared<type> args), if
      std::shared_ptr is supported, and with (new type args) otherwise. Note that the arguments for
      the constructor should be specified in additional braces, e. g.:
      MAKE_SHARED_OR_AUTO(Type, (arg1, arg2, arg3));

    DELETED - is replaced with =delete if deleted functions are available and is empty otherwise.

\endif
\if RU
  � ������ ����� ������������ �������, ��������������� ��������� ������������, ��������� �
  ��������� C++11 (ISO/IEC 14882:2011). ����������� �������� �� ����������, �������������� ��
  ��������� ������:
    ������� ������� ������������ C++11 � ��������� ������������:
      http://wiki.apache.org/stdcxx/C++0xCompilerSupport
    ������ ��������, ������������ ������ �����������
      http://sourceforge.net/p/predef/wiki/Compilers/
    
  �������������� �����������
    Clang 3.0,
    GCC 4.0,
    Intel C++ compiler 9.0,
    MS Visual C++ 8.0 (Visual Studio 2005),
  � ����� ����� ������

  ��� ��������� ��������� ������������ C++11 � Clang, GCC � ����������� Intel ���������� �������
  ��������������� �������� � ��������� ������ (-std=c++11 ��� -std=c++0x, ��. ��������
  �����������). � MSVC ��������� C++11 �������� ������.

  ����� ��������� ��� ����������� ����� ������� HAS_xxx �, ��������������, ������ �������������
  ������������ C++11 � ������ ����������, ���� ��������� ���������� ������ DISABLE_CPP11
  (��������� C++11 ������������ ��� ���� �� �����������, �� ����������� C++11 ��������������
  � ���������� �� �����; ��� ����� �������� � ���������� ���������� ������� ����� �/��� ������
  ������������ ������).

  ����� ����������� ����������� �������� � ������� � ����� ��������� ������������� ��� �����������,
  ������� ������������ �������� ��� ������ ����������

  ��������� ������� ����������, ���� ��������������� ����������� �������������� ������������
  (������ �� ������ �������� ��������� � ��������� ISO/IEC 14882:2011(E) (International Standard:
  Information technology � Programming languages � C++):

    HAS_ARRAY - ���������, ���� ���������� ������������ ���� <array> (���������� ����� std::array).
      ��. ������� 23.3.1 � 23.3.2

    HAS_DEFAULTED - ���������, ���� �������������� ��������� ������ � ������ �� ���������
      (������������� =delete � =default).
      ��. ������� 8.4.2 � 8.4.3 ��� http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2007/n2346.htm

    HAS_FRIENDEXT - ���������, ���� �������������� ����������� ����������� ������
      (��������, friend T; ��� T - �������� �������).
      ��. ������ 11.3 ��� http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2005/n1791.pdf

    HAS_INITLIST - ���������, ���� ������ ������������� �������� � �������� ����� �������������� �
      ����������.
      ��. ������� 8.5.4 � 18.9 ��� http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2008/n2672.htm

    HAS_LAMBDAS - ���������, ���� �������������� ������-��������� � ����� std::function.
      ��. ������� 5.1.2 � 20.8.11. ������ ������ ����� ��������� ���� �� ��������� �������� �
      ����������� �� ����, ����� ������ ������ ���������������� �����������:
        09 - ������ v.0.9: ������ ��� ���������, ��� ����������� �������������� ������� ���������
          � "�������" ���. ��. http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2008/n2550.pdf
        10 - ������ v.1.0: ��������� ����������� �������� ���������� ����� (����� � ����������).
          ��. http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2008/n2658.pdf
        11 - ������ v.1.1: �������� ��������� ����� std::function<>, ��������� ��������� � ����
          �������� ��������. ��. http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2009/n2927.pdf
      ��������� ����������� ����� ��������� �� �������� ��������, � ������� ������� ������� �����
      �������� ���������� �� �������� ������� ��������� �������� (1.1) ������ �����. � ���� ������
      ����� ������ � ����������� ������� ������� ����� ������.

    HAS_RVALUE_REFS - ���������, ���� �������������� ��� "������ �� rvalue" (Type&&). ��� ����
      ������������� ������ �� this (T::f()&) ����� �� �������������� (���� ������� ����������� �
      �������, �� ��� ��� ������� ���������� ��������� ������).
      ��. ������ 8.3.2 ��� http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2006/n2118.html
      ���� �������� ��������� ������� ���� ����������������, ����������� ��� ��������. ��� �����
      ���� ����� ��� ���� ������ ����������, ����� � �� ����. �� ������ ������ ��������� ������
      ��������� ��� ������� ����� ���������� ������ �� rvalue, ��, ��������, ����������� ��������
      ������, ��� ��� ���� ������� ��� �����.

    HAS_SMART_PTR - ���������, ���� �������������� std::shared_ptr, std::weak_ptr � std::unique_ptr.
      ��. ������ 20.7.

    HAS_STDFUNCTION - ���������, ���� �������������� ��������� ����� std::function<> � ��������
      �������������� ������-��������� � ������� ����� ������. ������ ������ ������ ���������, ����
      HAS_LAMBDAS >= 11, ������ �� ���������, ���� �� ��������� HAS_LAMBDAS � ����� ���� ���������,
      ���� HAS_LAMBDAS ��������� � HAS_LAMBDAS < 11 (��������� ����������� �� ������������ ������
      v.1.1 ���������, �� ������������ std::function<>).

    HAS_VARTEMPLATES - ���������, ���� �������������� ������� � ������������ ������ ����������.
      ��. ������ 14.5.3 ��� http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2007/n2242.pdf �
      http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2008/n2555.pdf
      
  ��������� ������� ����� ������������ ��� ��������� ����, ������� ���������� ����������� C++11
  ��� ������� �������.

    SHARED_OR_AUTO(type) - ���������� �� std::shared_ptr<type>, ���� �� ��������������, � ��
      std::auto_ptr<type> � ��������� ������. �������� ��������, ��� std::auto_ptr �� ������������
      ����������� �������� � �� ����� �������������� � �����������.

    UNIQUE_OR_AUTO(type) - ���������� �� std::unique_ptr<type>, ���� �� ��������������, � ��
      std::auto_ptr<type> � ��������� ������. �������� ��������, ��� std::auto_ptr �� �����
      �������������� � �����������.

    SHARED_OR_RAWPTR(type) - ���������� �� std::shared_ptr<type>, ���� �� ��������������, � ��
      type* � ��������� ������. �� �������� ���������� �������� ������� ����� ���� ��� � ������,
      ���� �� ���� ��������� ������� ���������.

    UNIQUE_OR_RAWPTR(type) - ���������� �� std::unique_ptr<type>, ���� �� ��������������, � ��
      type* � ��������� ������. �� �������� ���������� �������� ������� ����� ���� ��� � ������,
      ���� �� ���� ��������� ������� ���������.

    MAKE_SHARED_OR_AUTO(type, args) - ���������� �� (std::make_shared<type> args), ����
      std::shared_ptr ��������������, � �� (std::auto_ptr<type>(new type args)) � ��������� ������.
      �������� ��������, ��� ��������� ������������ ������ ����������� � �������������� �������,
      ��������:
      MAKE_SHARED_OR_AUTO(Type, (arg1, arg2, arg3));
      std::shared_ptr ���� ����� ���� ��������������� ����������, ���������� ����� new, ������
      ����������� ����� ������� ������ ��������� �������� ������� ������ ����� (��� ��� ������ ���
      ������������ ����� � ������������ ������� ���������� ��������).

    MAKE_SHARED_OR_RAW(type, args) - ���������� �� (std::make_shared<type> args), ����
      std::shared_ptr ��������������, � �� (new type args) � ��������� ������. �������� ��������,
      ��� ��������� ������������ ������ ����������� � �������������� �������, ��������:
      MAKE_SHARED_OR_AUTO(Type, (arg1, arg2, arg3));

    DELETED - ���������� �� =delete, ���� ��������� ������ ��������������, � �� ������ ������ �
      ��������� ������.
\endif
*/

#ifndef DISABLE_CPP11

#if defined(__clang__)
#  define HAS_FRIENDEXT
#  define HAS_SMART_PTR
#  if __has_extension(cxx_rvalue_references)
#    define HAS_RVALUE_REFS
#  endif
#  if __has_include(<array>)
#    define HAS_ARRAY
#  endif
#  if __has_extension(cxx_generalized_initializers)
#    define HAS_INITLIST
#  endif
#  if __has_extension(cxx_defaulted_functions)
#    define HAS_DEFAULTED
#  endif
#  if __has_extension(cxx_variadic_templates)
#    define HAS_VARTEMPLATES
#  endif
#  if __has_extension(cxx_lambdas)
#    define HAS_LAMBDAS 11
#  endif
#elif defined(__INTEL_COMPILER) && ((defined(__GNUC__) && defined(__GXX_EXPERIMENTAL_CXX0X__)) \
  || (defined(_MSC_VER) && defined(__INTEL_CXX11_MODE__)))
// On Windows environment when using Intel C++ compiler with Visual Studio 2010* or 2012*,
// the C++11 features supported by Visual C++ 2010/2012 are enabled by default.
// Use "/Qstd=c++11" to turn on the support for all other cases.
// If this command line option is specified, then macro __INTEL_CXX11_MODE__
// is likely to be defined, but it isn't documented, so we may be wrong here
#  if _MSC_VER+0 >= 1600 || (__GNUC__+0)*100 + __GNUC_MINOR___+0 >= 403
// Presence of <array> depends on library
#    define HAS_ARRAY
#  endif
#  if __INTEL_COMPILER >= 1300
#    define HAS_INITLIST
// Intel C++ Compiler 13 doesn't contain lambdas v.1.1, but supports conversion of lambdas
// to std::function
#    define HAS_STDFUNCTION
#  endif
#  if __INTEL_COMPILER >= 1201
#    define HAS_VARTEMPLATES
#  endif
#  if __INTEL_COMPILER >= 1200
#    define HAS_DEFAULTED
#    define HAS_LAMBDAS 10
#  elif __INTEL_COMPILER >= 1100
#    define HAS_LAMBDAS 9
#  endif
#  if __INTEL_COMPILER >= 1101
#    define HAS_RVALUE_REFS
#    define HAS_SMART_PTR
#  endif
#  if __INTEL_COMPILER >= 1100
#    define HAS_FRIENDEXT
#  endif
#elif defined(_MSC_VER)
// Microsoft added some new C++11 features in MSVC Nov.2012 CTP
// (customer technology preview), but it's a bit buggy. If we
// encounter too much problems, then it may be better to replace
// the condition below with _MSC_VER >= 1800
#  if _MSC_FULL_VER >= 170051025
#    define HAS_INITLIST
#    define HAS_VARTEMPLATES
#  endif
#  if _MSC_VER >= 1700
#    define HAS_LAMBDAS 11
#  elif _MSC_VER >= 1600
#    define HAS_LAMBDAS 10
#  endif
#  if _MSC_VER >= 1600
#    define HAS_ARRAY
#    define HAS_FRIENDEXT
#    define HAS_RVALUE_REFS
#    define HAS_SMART_PTR
// Though MSVC10 has no complete implementation of lambdas v.1.1, it has std::function class
// and lambdas can be converted to it
#    define HAS_STDFUNCTION
#  endif
// MSVC didn't implement defaulted and deleted functions yet
//#  if _MSC_VER >= 9999
//#    define HAS_DEFAULTED
//#  endif
#elif defined(__GNUC__) && defined(__GXX_EXPERIMENTAL_CXX0X__)
// assume that we use GNU stdc++ library
#  if __GNUC__*100 + __GNUC_MINOR__ >= 407
#    define HAS_FRIENDEXT
#  endif
#  if __GNUC__*100 + __GNUC_MINOR__ >= 405
#    define HAS_LAMBDAS 11
#  endif
#  if __GNUC__*100 + __GNUC_MINOR__ >= 404
#    define HAS_INITLIST
#    define HAS_DEFAULTED
#  endif
#  if __GNUC__*100 + __GNUC_MINOR__ >= 403
#    define HAS_ARRAY
#    define HAS_RVALUE_REFS
#    define HAS_SMART_PTR
#    define HAS_VARTEMPLATES
#  endif
#endif

#if HAS_LAMBDAS >= 11
#define HAS_STDFUNCTION
#endif

#endif

// We don't want to reinvent smart pointers, but if there is no C++11
// smart pointer support, then the user should decide himself which
// tools to use (raw pointers, std::auto_ptr, boost or other libs),
// since none are completely compatible with C++11.
// The macros below provide automatic choice between
// std::shared_ptr/std::unique_ptr and std::auto/raw pointer types
// depending on presence of C++11 smart pointer support.
// Though these definitions save the user from using bulky and
// duplicating #ifdef switching in definitions, they don't provide
// memory management required for raw pointers and, of course, they
// don't make any pointers safe and compatible with each other.
#ifdef HAS_SMART_PTR
#define SHARED_OR_AUTO(type) std::shared_ptr<type>
#define UNIQUE_OR_AUTO(type) std::unique_ptr<type>
#define SHARED_OR_RAWPTR(type) std::shared_ptr<type>
#define UNIQUE_OR_RAWPTR(type) std::unique_ptr<type>

// Note that all compilers implemented rvalue references before variadic templates (or
// simultaneously with them), so a separate check is not needed
#define MAKE_SHARED_OR_AUTO(Type, Args) (std::make_shared<Type> Args)
#define MAKE_SHARED_OR_RAW MAKE_SHARED_OR_AUTO
#else
#define SHARED_OR_AUTO(type) std::auto_ptr<type>
#define UNIQUE_OR_AUTO(type) std::auto_ptr<type>
#define SHARED_OR_RAWPTR(type) type*
#define UNIQUE_OR_RAWPTR(type) type*

#define MAKE_SHARED_OR_AUTO(type, args) (std::auto_ptr<type>(new type args))
#define MAKE_SHARED_OR_RAW(type, args) (new type args)
#endif

#ifdef HAS_DEFAULTED
#define DELETED = delete
#else
#define DELETED
#endif

#endif
