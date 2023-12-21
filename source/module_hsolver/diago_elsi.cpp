#include "diago_elsi.h"

#include "module_base/global_variable.h"
#include "module_base/lapack_connector.h"
#include "module_base/timer.h"
#include "module_base/tool_quit.h"
#include "module_hamilt_general/operator.h"
extern "C"
{
#include "module_base/blacs_connector.h"
#include "module_base/scalapack_connector.h"
}

// for debug
#include <sstream>
#include <fstream>
#include <string>
// end debug

typedef hamilt::MatrixBlock<double> matd;
typedef hamilt::MatrixBlock<std::complex<double>> matcd;

namespace hsolver
{

bool is_new_e_iteration;
elsi_handle ELSI_handle;
bool is_ELSI_handle_inited=false;

inline void set_ELSI_handle(const MPI_Comm COMM_DIAG, const int desc[])
{
    if(!is_ELSI_handle_inited)
    {
        int blacs_ctxt=desc[1];
        int n_basis=desc[2];
        int nblk=desc[4];
        int solver=1; // 1: ELPA
        int par_mode=1; // 1: multi-processes
        int mat_format=0; // 0: blacs dense matrix format
        const int n_electron = static_cast<int> (GlobalV::nelec); // number of electrons in the system. 
        // It is not used when call ELPA kernel, but is needed to init elsi handle.
        const int n_state = GlobalV::NBANDS;
        c_elsi_init(&ELSI_handle, solver, par_mode, mat_format, n_basis, n_electron, n_state);
        c_elsi_set_mpi(ELSI_handle, COMM_DIAG);
        c_elsi_set_blacs(ELSI_handle, blacs_ctxt, nblk);
        is_ELSI_handle_inited=true;
    }
}

void DiagoElsi::diag(hamilt::Hamilt<double> *phm_in, psi::Psi<std::complex<double>> &psi, double *eigenvalue_in)
{
    ModuleBase::TITLE("DiagoElsi", "diag");
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"DiagElsi for complex started");
#ifdef __MPI
    matcd h_mat, s_mat;
    phm_in->matrix(h_mat, s_mat);

    if(is_ELSI_handle_inited)
    {
        c_elsi_reinit(ELSI_handle); 
        // complex type is run for k point and the S matrix can not be reused within electron steps
        // because it is used for different k point in every single step.
    }
    else
    {
        set_ELSI_handle(MPI_COMM_WORLD, (const int*)h_mat.desc);
    }

    std::vector<double> eigen(GlobalV::NLOCAL, 0.0);

    auto h=reinterpret_cast<double _Complex*>(h_mat.p);
    auto s=reinterpret_cast<double _Complex*>(s_mat.p);
    double* eval=eigen.data();
    auto evec=reinterpret_cast<double _Complex*>(psi.get_pointer());

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"start solving...");
    c_elsi_ev_complex(h_elsi_, h, s, eval, evec);
   
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"c_elsi_ev_complex ended, eval[0]", eval[0]);

    const int inc=1;
    BlasConnector::copy(GlobalV::NBANDS, eigen.data(), inc, eigenvalue_in, inc);
#else
    ModuleBase::WARNING_QUIT("DiagoElsi", "DiagoElsi only can be used with macro __MPI");
#endif
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"DiagElsi for complex ended");
}

void DiagoElsi::diag(hamilt::Hamilt<double> *phm_in, psi::Psi<double> &psi, double *eigenvalue_in)
{
    ModuleBase::TITLE("DiagoElsi", "diag");
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"DiagElsi for double started");
#ifdef __MPI
    matd h_mat, s_mat;
    phm_in->matrix(h_mat, s_mat);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"h and s are inited");

    if(is_ELSI_handle_inited)
    {
        if(is_new_e_iteration) c_elsi_reinit(ELSI_handle);
    }
    else {
        set_ELSI_handle(MPI_COMM_WORLD, (const int*)h_mat.desc);
    }
    
    std::vector<double> eigen(GlobalV::NLOCAL, 0.0);
    double* eval=eigen.data();

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"start solving...");
    c_elsi_ev_real(ELSI_handle, h_mat.p, s_mat.p, eval, psi.get_pointer());
    
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "K-S equation was solved by elsi");
    const int inc=1;
    BlasConnector::copy(GlobalV::NBANDS, eigen.data(), inc, eigenvalue_in, inc);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "eigenvalues were copied to ekb");
#else
    ModuleBase::WARNING_QUIT("DiagoElsi", "DiagoElsi only can be used with macro __MPI");
#endif
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"DiagElsi for double ended");
}

} // namespace hsolver
