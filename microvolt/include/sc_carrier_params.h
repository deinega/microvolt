#ifndef SC_CARRIER_PARAMS
#define SC_CARRIER_PARAMS

#pragma warning(disable:4996)

#include <algorithm>
#include <functional>

#include "cpp11features.h"

#ifdef HAS_ARRAY
#include <array>
#else
#include <stdexcept>
#endif

// The header is needed just for the next typedef
#include "vector_3.h"
typedef vec_type valtype; // double or float (if macro SINGLE_PRECISION is defined)

// Flags for output data
enum scoutTYPE{
  outPsi=0x1,    // Potential
  outFermi=0x2,  // Quasi-Fermi levels
  outBands=0x4, // Condunction and valence bands
  outConc=0x8,   // Concentrations of charge carriers
  outDop=0x10,    // Concentrations of dopants
  outR_sc=0x20,  // Recombination source term
  outG_sc=0x40,  // Electron-hole generation term
  outJ=0x80,    // Components of the current density vector (for electrons or for holes)
  outJtot=0x100, // Components of the current density vector (for electrons + holes)
  outE_sc=0x200, // Electric field
  outResidue=0x400, // Residue for disretized equation
  outExcitons=0x800, // Concentrations of singlet and triplet excitons
};

// enumeration to denote different type of variables
enum{ // change order
  rpConc=0x1, // concentration
  rpFermi=0x2, // quasi-fermi potential
  rpSB=0x4, // slotboom variable
  rpC = rpConc|rpFermi|rpSB,

  rpPsi=0x8, // potential
  rpSBpsi=0x10, // exp(psi*q/kT)
  rpP = rpPsi|rpSBpsi,

  rpJ=0x20, // curent

  rpE=0x40, // electrones
  rpH=0x80, // holes
  rpEH = rpE|rpH,

  rpExcitons=0x100, // excitons
  rpS=0x200, // singlet excitons
  rpT=0x400, // triplet excitons
  rpST = rpS|rpT,

  rpAll=0xfffff
};

// index 0 is for electrons, 1 is for holes
#ifdef HAS_ARRAY
struct CarrierParams : public std::array<valtype, 2> {
#else
// Declare CarrierParams similarly to std::array from C++11
// if we have no <array> header
struct CarrierParams {
  typedef valtype value_type;
  typedef value_type* pointer;
  typedef const value_type* const_pointer;
  typedef value_type& reference;
  typedef const value_type& const_reference;
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;

private:
  template <class ptrType, class refType>
  class iterator_base : public std::iterator<random_access_iterator_tag,
    value_type, difference_type, ptrType, refType> {
  private:
    typedef iterator_base<ptrType, refType> iterType;
    friend struct CarrierParams;
    ptrType pos;

    explicit iterator_base(ptrType pos_) : pos(pos_) { }
  public:
    iterator_base() : pos(0) { }

    refType operator*() const { return *pos; }
    ptrType operator->() const { return &**this; }
    iterType& operator++() { ++pos; return *this; }
    iterType operator++(int) { iterType tmp(*this); ++*this; return tmp; }
    iterType& operator--() { --pos; return *this; }
    iterType operator--(int) { iterType tmp(*this); --*this; return tmp; }
    iterType& operator+=(difference_type offset) { pos += offset; return *this; }
    iterType operator+(difference_type offset) const { iterType tmp(*this); return tmp += offset; }
    iterType& operator-=(difference_type offset) { pos -= offset; return *this; }
    iterType operator-(difference_type offset) const { iterType tmp(*this); return tmp -= offset; }
    difference_type operator-(const iterType& other) const { return pos - other.pos; }
    refType operator[](difference_type offset) const { return *(*this + offset); }
    bool operator==(const iterType& other) const { return pos == other.pos; }
    bool operator!=(const iterType& other) const { return !(*this == other); }
    bool operator<(const iterType& other) const { return pos < other.pos; }
    bool operator>(const iterType& other) const { return other < *this; }
    bool operator<=(const iterType& other) const { return !(other < *this); }
    bool operator>=(const iterType& other) const { return !(*this < other); }
  };

public:
  typedef iterator_base<pointer, reference> iterator;
  typedef iterator_base<const_pointer, const_reference> const_iterator;
  typedef std::reverse_iterator<iterator> reverse_iterator;
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

  void fill(const value_type& u) { std::fill_n(begin(), size(), u); }
  void swap(CarrierParams& other) { std::swap_ranges(begin(), end(), other.begin()); }

  // iterators:
  iterator begin() { return iterator(elems); }
  const_iterator begin() const { return cbegin(); }
  iterator end() { return begin() + size(); }
  const_iterator end() const { return cend(); }
  reverse_iterator rbegin() { return reverse_iterator(end()); }
  const_reverse_iterator rbegin() const { return crbegin(); }
  reverse_iterator rend() { return reverse_iterator(begin()); }
  const_reverse_iterator rend() const { return crend(); }
  const_iterator cbegin() const { return const_iterator(elems); }
  const_iterator cend() const { return cbegin() + size(); }
  const_reverse_iterator crbegin() const { return const_reverse_iterator(cend()); }
  const_reverse_iterator crend() const { return const_reverse_iterator(cbegin()); }

  // capacity:
  size_type size() const { return N; }
  size_type max_size() const { return N; }
  bool empty() const { return false; }

  // element access:
  reference operator[](size_type n) { return elems[n]; }
  const_reference operator[](size_type n) const { return elems[n]; }

  const_reference at(size_type n) const {
    if(n >= size())
      throw std::out_of_range("Invalid charge carrier index");
    return (*this)[n];
  }
  reference at(size_type n) {
    if(n >= size())
      throw std::out_of_range("Invalid charge carrier index");
    return (*this)[n];
  }

  reference front() { return (*this)[0]; }
  const_reference front() const { return (*this)[0]; }
  reference back() { return (*this)[size()-1]; }
  const_reference back() const { return (*this)[size()-1]; }
  pointer data() { return elems; }
  const_pointer data() const { return elems; }

  bool operator==(const CarrierParams& other) const {
    return std::equal(cbegin(), cend(), other.cbegin());
  }
  bool operator!=(const CarrierParams& other) const { return !(*this == other); }
  bool operator<(const CarrierParams& other) const {
    return std::lexicographical_compare(cbegin(), cend(), other.cbegin(), other.cend());
  }
  bool operator>(const CarrierParams& other) const { return other < *this; }
  bool operator<=(const CarrierParams& other) const { return !(other < *this); }
  bool operator>=(const CarrierParams& other) const { return !(*this < other); }

private:
  static const size_type N = 2;
  value_type elems[N];

public:
#endif
  explicit CarrierParams(value_type val = 0) { fill(val); }
  
  CarrierParams(value_type n, value_type p) {
    (*this)[0] = n;
    (*this)[1] = p;
  }

#ifdef HAS_ARRAY
  CarrierParams(const std::array<value_type, 2>& other)
    : std::array<value_type, 2>(other) { }
#endif
};

#ifndef HAS_ARRAY
inline void swap(CarrierParams& x, CarrierParams& y) { return x.swap(y); }
#endif

// Arithmetic operations and functions on CarrierParams
// We could use valarray as an underlying type, but it's allocated dynamically,
// which is inefficient for arrays with few elements like this one
inline CarrierParams& operator+=(CarrierParams& x, const CarrierParams& y) {
  std::transform(x.begin(), x.end(), y.begin(), x.begin(), std::plus<CarrierParams::value_type>());
  return x;
}

inline CarrierParams& operator-=(CarrierParams& x, const CarrierParams& y) {
  std::transform(x.begin(), x.end(), y.begin(), x.begin(), std::minus<CarrierParams::value_type>());
  return x;
}

inline CarrierParams& operator*=(CarrierParams& x, const CarrierParams& y) {
  std::transform(x.begin(), x.end(), y.begin(), x.begin(), std::multiplies<CarrierParams::value_type>());
  return x;
}

inline CarrierParams& operator/=(CarrierParams& x, const CarrierParams& y) {
  std::transform(x.begin(), x.end(), y.begin(), x.begin(), std::divides<CarrierParams::value_type>());
  return x;
}

inline CarrierParams& operator+=(CarrierParams& x, CarrierParams::value_type y) {
  std::transform(x.begin(), x.end(), x.begin(), std::bind2nd(std::plus<CarrierParams::value_type>(), y));
  return x;
}

inline CarrierParams& operator-=(CarrierParams& x, CarrierParams::value_type y) {
  std::transform(x.begin(), x.end(), x.begin(), std::bind2nd(std::minus<CarrierParams::value_type>(), y));
  return x;
}

inline CarrierParams& operator*=(CarrierParams& x, CarrierParams::value_type y) {
  std::transform(x.begin(), x.end(), x.begin(), std::bind2nd(std::multiplies<CarrierParams::value_type>(), y));
  return x;
}

inline CarrierParams& operator/=(CarrierParams& x, CarrierParams::value_type y) {
  std::transform(x.begin(), x.end(), x.begin(), std::bind2nd(std::divides<CarrierParams::value_type>(), y));
  return x;
}

inline CarrierParams operator+(CarrierParams x) {
  return x;
}

inline CarrierParams operator-(CarrierParams x) {
  std::transform(x.begin(), x.end(), x.begin(), std::negate<CarrierParams::value_type>());
  return x;
}

inline CarrierParams operator+(CarrierParams x, const CarrierParams& y) { return x += y; }
inline CarrierParams operator-(CarrierParams x, const CarrierParams& y) { return x -= y; }
inline CarrierParams operator*(CarrierParams x, const CarrierParams& y) { return x *= y; }
inline CarrierParams operator/(CarrierParams x, const CarrierParams& y) { return x /= y; }
inline CarrierParams operator+(CarrierParams x, CarrierParams::value_type y) { return x += y; }
inline CarrierParams operator-(CarrierParams x, CarrierParams::value_type y) { return x -= y; }
inline CarrierParams operator*(CarrierParams x, CarrierParams::value_type y) { return x *= y; }
inline CarrierParams operator/(CarrierParams x, CarrierParams::value_type y) { return x /= y; }
inline CarrierParams operator+(CarrierParams::value_type x, CarrierParams y) { return y += x; }
inline CarrierParams operator-(CarrierParams::value_type x, CarrierParams y) { return x + (-y); }
inline CarrierParams operator*(CarrierParams::value_type x, CarrierParams y) { return y *= x; }

inline CarrierParams operator/(CarrierParams::value_type x, CarrierParams y) {
  std::transform(y.begin(), y.end(), y.begin(), std::bind1st(std::divides<CarrierParams::value_type>(), x));
  return y;
}

inline CarrierParams exp(CarrierParams x) {
  std::transform(x.begin(), x.end(), x.begin(), std::ptr_fun<CarrierParams::value_type, CarrierParams::value_type>(std::exp));
  return x;
}

inline CarrierParams log(CarrierParams x) {
  std::transform(x.begin(), x.end(), x.begin(), std::ptr_fun<CarrierParams::value_type, CarrierParams::value_type>(std::log));
  return x;
}

inline CarrierParams sqrt(CarrierParams x) {
  std::transform(x.begin(), x.end(), x.begin(), std::ptr_fun<CarrierParams::value_type, CarrierParams::value_type>(std::sqrt));
  return x;
}

template <class Type>
inline CarrierParams pow(CarrierParams x, Type n) {
  std::transform(x.begin(), x.end(), x.begin(), std::bind2nd(
    std::ptr_fun<CarrierParams::value_type, Type, CarrierParams::value_type>(std::pow), n));
  return x;
}


// Functions relevant to CarrierParams only
// Apply negation to electrons only
inline CarrierParams negate_first(CarrierParams x) {
  x[0] = -x[0];
  return x;
}

inline CarrierParams negate_first(CarrierParams::value_type x) { return negate_first(CarrierParams(x)); }

// Returns {exp(-x[0]), exp(x[1])}
inline CarrierParams altexp(CarrierParams x) { return exp(negate_first(x)); }

// Returns {-log(x[0]), log(x[1])}
inline CarrierParams altlog(CarrierParams x) { return negate_first(log(x)); }

#endif
