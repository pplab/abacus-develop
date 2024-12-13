#ifdef __NTPOLY
#include <mpi.h>
#include <complex>
#include <memory>
#include "module_parameter/parameter.h"
#include "module_base/global_variable.h"
#include "module_base/tool_quit.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_ntpoly/simple_ntpoly.h"
#include "diago_ntpoly.h"

typedef hamilt::MatrixBlock<double> matd;
typedef hamilt::MatrixBlock<std::complex<double>> matcd;

namespace hsolver
{
template <typename T>

DiagoNTPoly<T>::DiagoNTPoly(const Parallel_Orbitals* ParaV_in)
{
    const int nspin = PARAM.inp.nspin==2 ? 2:1;

    this->ParaV = ParaV_in;

    this->DM.resize(nspin);
    this->EDM.resize(nspin);
    for (int i = 0; i < nspin; i++)
    {
        this->DM[i] = new T[ParaV->nrow * ParaV->ncol];
        this->EDM[i] = new T[ParaV->nrow * ParaV->ncol];
    }

}

template <typename T>
DiagoNTPoly<T>::~DiagoNTPoly()
{
    const int nspin = PARAM.inp.nspin==2 ? 2:1;
    for (int i = 0; i < nspin; i++)
    {
        delete[] this->DM[i];
        delete[] this->EDM[i];
    }

}

template <>
void DiagoNTPoly<double>::diag(hamilt::Hamilt<double>* phm_in, psi::Psi<double>& psi, double* eigenvalue_in)
{
    ModuleBase::TITLE("DiagoNTPoly", "diag");
    matd h_mat, s_mat;
    phm_in->matrix(h_mat, s_mat);
    int ik = psi.get_current_k();
    if(ntpoly::for_debug) ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "enter DiagoNTPoly<double>::diag, ik", ik);
    const int nelec=PARAM.inp.nelec;
    const int nspin = PARAM.inp.nspin==2 ? 2:1;
    const double converge_density = 1e-10;
    const double converge_overlap = 1e-10;
    const double threshold = 1e-10;
    ntpoly::simple_ntpoly(MPI_COMM_WORLD, h_mat.desc,
                      h_mat.row, h_mat.col,
                      converge_density, converge_overlap, threshold, 
                      nelec, nspin, h_mat.p, s_mat.p,
                      DM[ik], EDM[ik],
                      this->energy, this->chemical_potential);
    if(ntpoly::for_debug) ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "solved simple_ntpoly");
}

template <>
void DiagoNTPoly<std::complex<double>>::diag(hamilt::Hamilt<std::complex<double>>* phm_in,
                                            psi::Psi<std::complex<double>>& psi,
                                            double* eigenvalue_in)
{
    ModuleBase::TITLE("DiagoNTPoly", "diag");
    ModuleBase::WARNING_QUIT("DiagoNTPoly", "NTPoly is not completed for multi-k case");
}

template class DiagoNTPoly<double>;
template class DiagoNTPoly<std::complex<double> >;

} // namespace hsolver
#endif