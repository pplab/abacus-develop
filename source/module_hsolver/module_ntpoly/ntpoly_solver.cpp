#ifdef __NTPOLY
#include "module_base/parallel_global.h"
#include "ntpoly_solver.h"

#include <mpi.h>
#include <cstring>
#include <vector>

#include "module_base/global_variable.h"

extern MPI_Comm DIAG_WORLD;
extern MPI_Comm GRID_WORLD;
namespace ntpoly
{

void NTPOLY_Solver::NTPOLY_Solver(int desc[9], int nrow_in, int ncol_in, int nelec_in, int nspin_in, int nblk_in,
        double converge_density_in, double converge_overlap_in, double threshold_in)
        {
            std::memcpy(desc, desc_in, sizeof(int)*9);
            this->nrow=nrow_in;
            this->ncol=ncol_in;
            this->nelec=nelec_in;
            this->nspin=nspin_in;
            this->nblk=nblk_in;
            this->converge_density=converge_density_in;
            this->converge_overlap=converge_overlap_in;
            this->threshold=threshold_in;
        }

int NTPOLY_Solver::solve(const double* h,
                         const double* s,
                         double*& DM,
                         double*& EDM,
                         double& energy,
                         double& chemical_potential)
{
    int ret=simple_ntpoly(GRID_WORLD, this->desc, 
                this->nrow, this->ncol, 
                this->converge_density, this->converge_overlap, this->threshold, 
                this->nelec, this->nspin, H, S, DM, EDM, 
                energy, chemical_potential);
    return ret;
}

} // namespace ntpoly
#endif