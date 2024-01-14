#include "diago_elsi.h"
#include "module_base/global_variable.h"
#include "module_base/lapack_connector.h"
#include "module_base/timer.h"
#include "module_base/tool_quit.h"
extern "C"
{
#include "module_base/blacs_connector.h"
#include "module_base/scalapack_connector.h"
}
#include "unistd.h"


typedef std::complex<double> double_complex;
typedef hamilt::MatrixBlock<double> matd;
typedef hamilt::MatrixBlock<double_complex> matcd;

using std::vector;

namespace hsolver
{
    bool is_new_e_iteration=true;
    elsi_handle ELSI_handle;
    bool is_ELSI_handle_inited=false;

    inline void set_ELSI_handle(const MPI_Comm COMM_DIAG, const int desc[])
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

    template<>
    void DiagoElsi<double_complex>::diag(hamilt::Hamilt<double_complex>* phm_in, psi::Psi<double_complex>& psi, Real* eigenvalue_in)
    {
        ModuleBase::TITLE("DiagoElsi", "diag");
        ModuleBase::timer::tick("DiagoElsi", "diag_dc");
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"DiagElsi for complex started");
    #ifdef __MPI
        const int inc=1;
        matcd h_mat, s_mat;
        phm_in->matrix(h_mat, s_mat);
        if(is_ELSI_handle_inited)
        {
            c_elsi_reinit(ELSI_handle);
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"ELSI_handle is reinited");
            // complex type is run for k point and the S matrix can not be reused within electron steps
            // because it is used for different k point in every single step.
        }
        else
        {
            set_ELSI_handle(MPI_COMM_WORLD, (const int*)h_mat.desc);
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"ELSI_handle is set");
        }

        vector<double> eigen(GlobalV::NLOCAL, 0.0);
        double* eval=eigen.data();
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"eigen is allocated, size ", GlobalV::NLOCAL);

        auto h=reinterpret_cast<double _Complex*>(h_mat.p);
        auto s=reinterpret_cast<double _Complex*>(s_mat.p);
        auto evec=reinterpret_cast<double _Complex*>(psi.get_pointer());

        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"start solving...");
        c_elsi_ev_complex(ELSI_handle, h, s, eval, evec);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"c_elsi_ev_complex ended, eval[0]", eval[0]);

        BlasConnector::copy(GlobalV::NBANDS, eval, inc, eigenvalue_in, inc);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"eval is copied to eigenvalue_in");
    #else
        ModuleBase::WARNING_QUIT("DiagoElsi", "DiagoElsi only can be used with macro __MPI");
    #endif
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"DiagElsi for complex ended");
        ModuleBase::timer::tick("DiagoElsi", "diag_dc");
    }

    template<>
    void DiagoElsi<double>::diag(hamilt::Hamilt<double>* phm_in, psi::Psi<double>& psi, Real* eigenvalue_in)
    {
        ModuleBase::TITLE("DiagoElsi", "diag");
        ModuleBase::timer::tick("DiagoElsi", "diag_d");
    #ifdef __MPI
        const int inc=1;
        matd h_mat, s_mat;
        phm_in->matrix(h_mat, s_mat);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"h and s are inited");

        if(is_ELSI_handle_inited)
        {
            if(is_new_e_iteration)
            {
                c_elsi_reinit(ELSI_handle);
                is_new_e_iteration=false;
                ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"ELSI_handle is reinited");
            }
        }
        else
        {
            set_ELSI_handle(MPI_COMM_WORLD, (const int*)h_mat.desc);
            is_new_e_iteration=false;
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"ELSI_handle is set");
        }
        
        vector<double> eigen(GlobalV::NLOCAL, 0.0);
        double* eval=eigen.data();

        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"start solving...");
        c_elsi_ev_real(ELSI_handle, h_mat.p, s_mat.p, eval, psi.get_pointer());
        
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "K-S equation was solved by elsi");
        BlasConnector::copy(GlobalV::NBANDS, eigen.data(), inc, eigenvalue_in, inc);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "eigenvalues were copied to ekb");
    #else
        ModuleBase::WARNING_QUIT("DiagoElsi", "DiagoElsi only can be used with macro __MPI");
    #endif
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"DiagElsi for double ended");
        ModuleBase::timer::tick("DiagoElsi", "diag_d");
    }
} // namespace hsolver
