#ifndef DIAGOELSI_H
#define DIAGOELSI_H

#include <map>
#include "diagh.h"
#include "module_orbital/parallel_orbitals.h"
#include <elsi.h>

namespace hsolver
{
extern bool is_new_e_iteration;
extern elsi_handle KS_handle;
extern bool is_KS_handle_inited;

class DiagoElsi : public DiagH<double>
{

  public:
  // Elsi could be called by different functions, and each call needs a 
  // independent elsi handle because the parameters may not be the same.
  // The ELSI_HANDLE_POOL is used to save all elsi handles and the key 
  // is the name of the function who calls elsi.
    static std::map<std::string, elsi_handle> ELSI_HANDLE_POOL;
    
  // The ELSI_NEW_ITER_POOL is used to save the states of interations
  // for elsi handles because there may be more than one elsi handle is
  // used when the program is running.
  // The key should be the same key of the handle in ELSI_HANDLE_POOL
    static std::map<std::string, int> ELSI_ITERATION_POOL;
    
    DiagoElsi(const std::string elsi_handle_key, const MPI_Comm COMM_DIAG, const int desc[]);
    ~DiagoElsi();
    void diag(hamilt::Hamilt<double>* phm_in, psi::Psi<double>& psi, double* eigenvalue_in) override;
    void diag(hamilt::Hamilt<double>* phm_in, psi::Psi<std::complex<double>>& psi, double* eigenvalue_in) override;
    
    // finalize an elsi handle
    static void finalize(const std::string elsi_handle_key);
  private:
    std::string elsi_handle_key_;
    elsi_handle h_elsi_;
};

} // namespace hsolver

#endif
