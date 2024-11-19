#ifdef __NTPOLY
#include <mpi.h>
#include <complex>
#include <memory>
#include "diago_ntpoly.h"
#include "module_base/global_variable.h"
#include "module_base/tool_quit.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_ntpoly/simple_ntpoly.h"

typedef hamilt::MatrixBlock<double> matd;
typedef hamilt::MatrixBlock<std::complex<double>> matcd;

namespace hsolver
{
template <typename T>

template <typename T>
DiagoNTPoly<T>::DiagoNTPoly(const Parallel_Orbitals* ParaV_in)
{
    int nspin = GlobalV::NSPIN;
    if (GlobalV::NSPIN == 4)
    {
        nspin = 1;
    }

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
    int nspin = GlobalV::NSPIN;
    if (GlobalV::NSPIN == 4)
    {
        nspin = 1;
    }
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
    std::vector<double> eigen(GlobalV::NLOCAL, 0.0);
    int ik = psi.get_current_k();
    NTPoly::simple_ntpoly(this->ParaV->blacs_ctxt,
                      this->ParaV->nb,
                      this->ParaV->nrow,
                      this->ParaV->ncol,
                      h_mat.p,
                      s_mat.p,
                      DM[ik],
                      EDM[ik]);
    this->energy = this->ps->get_totalFreeEnergy();
    this->chemical_potential = this->ps->get_totalEnergyH();
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