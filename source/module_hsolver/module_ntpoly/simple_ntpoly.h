#ifndef SIMPLE_NTPOLY_H
#define SIMPLE_NTPOLY_H
#include <mpi.h>

// NTPoly headers
#include <PSMatrix.h>
namespace ntpoly
{
   const bool for_debug=true;
   static bool require_init_NTPOLY=true;
/**
 * Main function for performing NTPoly calculations.
 *
 * @param comm_2D The MPI communicator for the 2D grid.
 * @param desc The descriptor array for the BCD matrix.
 * @param nrow The number of rows in the local matrix.
 * @param ncol The number of columns in the local matrix.
 * @param converge_density Convergence threshold for density matrix.
 * @param converge_overlap Convergence threshold for overlap matrix.
 * @param threshold Threshold for matrix element values.
 * @param nelec Number of electrons.
 * @param H Hamiltonian matrix.
 * @param S Overlap matrix.
 * @param DM Density matrix.
 * @param EDM Energy density matrix.
 * @param energy the total energy of the system.
 * @param chemical_potential The chemical potential.
 * @return Returns 0 if successful, or an error code if an error occurs.
 */
int simple_ntpoly(const MPI_Comm comm_2D, const int desc[], 
    const int nrow, const int ncol, 
    const double converge_density, const double converge_overlap, const double threshold, 
    const int nelec, const int nspin, const double H[], const double S[], 
    double DM[], double EDM[], 
    double& energy, double& chemical_potential);

/**
 * Constructs a PSMatrix from a Block Cyclic Distributed (BCD) matrix.
 *
 * @param PSM The PSMatrix to be constructed.
 * @param comm_2D The MPI communicator for the 2D grid.
 * @param desc The descriptor array for the BCD matrix.
 * @param nrow The number of rows in the local matrix.
 * @param ncol The number of columns in the local matrix.
 * @param M The BCD matrix from which the PSMatrix will be constructed.
 * @param threshold The threshold below which values are considered zero.
 * @return Returns 0 if successful, or an error code if an error occurs.
 */
int constructPSMatrixFromBCD(NTPoly::Matrix_ps& PSM, 
    const MPI_Comm comm_2D, const int desc[],
    const int nrow, const int ncol, const double M[], const double threshold);

/**
 * Reads all non-zero values from a Block Cyclic Distributed (BCD) matrix into a TripletList.
 * 
 * @param tripletList The TripletList to which the non-zero values will be appended.
 * @param comm_2D The MPI communicator for the 2D grid.
 * @param desc The descriptor array for the BCD matrix.
 * @param nrow The number of rows in the local matrix.
 * @param ncol The number of columns in the local matrix.
 * @param M The BCD matrix from which non-zero values will be read.
 * @param threshold The threshold below which values are considered zero.
 * @return Returns 0 if successful, or an error code if an error occurs.
 */
int readTripletListFromBCD(NTPoly::TripletList_r& tripletList, 
    const MPI_Comm comm_2D, const int desc[],
    const int nrow, const int ncol, const double M[], const double threshold);

/**
 * Constructs a Block Cyclic Distributed (BCD) matrix from a PSMatrix.
 *
 * @param PSM The PSMatrix to be converted.
 * @param comm_2D The MPI communicator for the 2D grid.
 * @param desc The descriptor array for the BCD matrix.
 * @param nrow The number of rows in the local matrix.
 * @param ncol The number of columns in the local matrix.
 * @param M The BCD matrix to be filled.
 * @return Returns 0 if successful, or an error code if an error occurs.
 */
int constructBCDFromPSMatrix(NTPoly::Matrix_ps& PSM, 
    const MPI_Comm comm_2D, const int desc[],
    const int nrow, const int ncol, double M[]);

/**
 * Calculates the global index corresponding to a given local index 
 * in a block-cyclic distribution.
 *
 * @param localIndex The local index within the current process.
 * @param nblk The number of elements in each block.
 * @param nprocs The total number of process rows or columns.
 * @param myproc The rank of the current process row or column.
 * @return The global index corresponding to the given local index.
 */
static inline int globalIndex(const int localIndex, const int nblk, const int nprocs, const int myproc)
{
    int iblock, gIndex;
    iblock = localIndex / nblk;
    gIndex = (iblock * nprocs + myproc) * nblk + localIndex % nblk;
    return gIndex;
}

/**
 * Calculates the local index and process number corresponding to a given global index 
 * in a block-cyclic distribution.
 *
 * @param globalIndex The global index within the entire distributed matrix.
 * @param nblk The number of elements in each block.
 * @param nprocs The total number of process rows or columns.
 * @param lcoalProc The rank of the process row or column to which the local index belongs.
 * @return The local index corresponding to the given global index.
 */
static inline int localIndex(const int globalIndex, const int nblk, const int nprocs, int& lcoalProc)
{
    lcoalProc = int((globalIndex % (nblk * nprocs)) / nblk);
    return int(globalIndex / (nblk * nprocs)) * nblk + globalIndex % nblk;
}
}
#endif // SIMPLE_NTPOLY_H