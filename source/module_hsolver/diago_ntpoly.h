#ifndef DIGAONTPOLY_H
#define DIGAONTPOLY_H

#include <vector>
#include <memory>
#include "module_base/macros.h"   // GetRealType
#include "module_hamilt_general/hamilt.h"
#include "module_base/global_variable.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_ntpoly/ntpoly_solver.h"

namespace hsolver
{

template <typename T>
class DiagoNTPoly : public DiagH<T>
{
  private:
    using Real = typename GetTypeReal<T>::type;
    static std::vector<double> mu_buffer;

  public:
    DiagoNTPoly(const Parallel_Orbitals* ParaV_in);
    void diag(hamilt::Hamilt<T>* phm_in, psi::Psi<T>& psi, Real* eigenvalue_in) override;
    const Parallel_Orbitals* ParaV;
    std::vector<T*> DM;
    std::vector<T*> EDM;
    double energy;
    double chemical_potential;
    std::unique_ptr<NTPOLY::NTPOLY_Solver> ps; // pointer to NTPOLY_Solver
    ~DiagoNTPoly();
};
} // namespace hsolver

#endif
