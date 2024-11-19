#ifndef NTPOLY_Solver_H
#define NTPOLY_Solver_H

#include <vector>
#include "dist_bcd_matrix.h"
#include "DensityMatrixSolvers.h"
#include "Logging.h"
#include "PSMatrix.h"
#include "Permutation.h"
#include "ProcessGrid.h"
#include "SolverParameters.h"
#include "SquareRootSolvers.h"

namespace ntpoly
{
class NTPOLY_Solver
{
  public:
    void solve(const double* h,
               const double* s,
               double*& DM,
               double*& EDM,
               double& energy,
               double& chemical_potential);

};
} // namespace ntpoly
#endif // NTPOLY_Solver_H