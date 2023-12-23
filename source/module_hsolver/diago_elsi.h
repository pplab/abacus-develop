#ifndef DIAGOELSI_H
#define DIAGOELSI_H

#include "diagh.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include <elsi.h>

namespace hsolver
{
    extern bool is_new_e_iteration;
    extern elsi_handle KS_handle;
    extern bool is_KS_handle_inited;
    
    template<typename T>
    class DiagoElsi : public DiagH<T>
    {
    private:
        using Real = typename GetTypeReal<T>::type;

    public:
        void diag(hamilt::Hamilt<T>* phm_in, psi::Psi<T>& psi, Real* eigenvalue_in) override;        
    };
} // namespace hsolver

#endif
