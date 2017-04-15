#ifndef SC_CONTEXT_H
#define SC_CONTEXT_H

/// @file sc_context.h Context for parameters which depend on calculated values.
/// This is used for extension of Microvolt for organic solar cells.

#include <map>
#include <set>

#include "cpp11features.h"

// The header is needed just for the next typedef
#include "vector_3.h"
typedef vec_type valtype; // double or float (if macro SINGLE_PRECISION is defined)


// Context for parameters which depend on calculated values
// (like electric field or carrier concentration)
class scPointContext {
public:
  // Identifiers for variables needed in the context
  // In general, these should depend on the equations involved
  // thus this list may be removed from here in future
  enum VarID {
    varPsi,      // Potential
    varEField,   // Electric field (absolute value)
    varFermiN,   // Quasi-Fermi level (inorganic) for electrons, eV
    varFermiP,   // Quasi-Fermi level (inorganic) for holes, eV
    varCondBand, // Condunction band, eV
    varValBand,  // Valence band, eV
    varConcN,    // Concentration of electrons, 1/m3
    varConcP,    // Concentration of holes, 1/m3
  };

  void set(VarID id, valtype value) { vars[id] = value; }

  valtype get(VarID id) const {
#if defined(_MSC_VER) && _MSC_VER < 1600
    // There is no map::at in C++03. MSVC supports it starting from VS2010.
    // Maybe a condition for some other compilers is also needed
    return 
      contains(id) ? 
      vars.find(id)->second : 
      throw std::out_of_range("Invalid variable id");
#else
    return vars.at(id);
#endif
  }

  bool contains(VarID id) const { return vars.count(id) > 0; }

  std::set<VarID> getVarList() const {
    std::set<VarID> result;
    for(std::map<VarID, valtype>::const_iterator it = vars.begin(), e = vars.end(); it != e; it++)
      result.insert(it->first);

    return result;
  }

  bool empty() const { return vars.empty(); }

  // The comparison may be slow, maybe it's better to generate a hash
  bool operator==(const scPointContext& other) const { return vars == other.vars; }
  bool operator!=(const scPointContext& other) const { return !(*this == other); }

private:
  // std::map seems to be too slow since it needs to be stored in heap
  // It also isn't capable of storing different types of variables
  // Wouldn't an explicit list of fields with a bit field (which specifies the needed variables) be better?
  std::map<VarID, valtype> vars;
};

#endif
