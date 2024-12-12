#ifndef DIGAONTPOLY_H
#define DIGAONTPOLY_H

#include <vector>
#include <memory>
#include "module_base/macros.h"   // GetRealType
#include "module_hamilt_general/hamilt.h"
#include "module_base/global_variable.h"
#include "module_basis/module_ao/parallel_orbitals.h"

namespace hsolver
{

template <typename T>
class DiagoNTPoly
{
  private:
    using Real = typename GetTypeReal<T>::type;
    static std::vector<double> mu_buffer;

  public:
    DiagoNTPoly(const Parallel_Orbitals* ParaV_in);
    void diag(hamilt::Hamilt<T>* phm_in, psi::Psi<T>& psi, Real* eigenvalue_in);
    const Parallel_Orbitals* ParaV;
    std::vector<T*> DM;
    std::vector<T*> EDM;
    double energy;
    double chemical_potential;
    ~DiagoNTPoly();
};
} // namespace hsolver

#endif
