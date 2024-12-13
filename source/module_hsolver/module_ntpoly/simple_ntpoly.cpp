#ifdef __NTPOLY
#include <iostream>
#include <vector>
#include <ProcessGrid.h>
#include <PSMatrix.h>
#include <TripletList.h>
#include <Triplet.h>
#include <Permutation.h>
#include <SolverParameters.h>
#include <SquareRootSolvers.h>
#include <DensityMatrixSolvers.h>
#include "Cblacs.h"
#include "simple_ntpoly.h"
#include "module_base/global_function.h"

namespace ntpoly
{
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
        double& energy, double& chemical_potential)
    {
        const int nFull=desc[2];
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "enter simple_ntpoly, nFull", nFull);
        // init default process grid
        int process_slice=1;
        NTPoly::ConstructGlobalProcessGrid(comm_2D, process_slice);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "GlobalProcessGrid is constructed");

        // init PSMatrices of Hamiltonian, Overlap, ISQOverlap, Density and EnergyDensity
        NTPoly::Matrix_ps Hamiltonian(nFull);
        NTPoly::Matrix_ps Overlap(nFull);  
        NTPoly::Matrix_ps ISQOverlap(nFull);
        NTPoly::Matrix_ps Density(nFull);
        NTPoly::Matrix_ps EnergyDensity(nFull);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "All PSMatrices are allocated, ActualDimension is", Hamiltonian.GetActualDimension());

        // convert H and S from BCD matrix to PSMatrix
        constructPSMatrixFromBCD(Hamiltonian, comm_2D, desc, nrow, ncol, H, threshold);
        constructPSMatrixFromBCD(Overlap, comm_2D, desc, nrow, ncol, S, threshold);        
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "H and S are converted to PSMatrix");

        // set purmutation
        NTPoly::Permutation permutation(nFull);
        permutation.SetRandomPermutation();

        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "permutation is done");

        // set solver parameters
        NTPoly::SolverParameters solver_parameters;
        solver_parameters.SetConvergeDiff(converge_overlap);
        solver_parameters.SetLoadBalance(permutation);
        solver_parameters.SetThreshold(threshold);
        solver_parameters.SetVerbosity(true);
        
        // InverseSquareRoot(Overlap, ISQOverlap, solver_parameters)
        NTPoly::SquareRootSolvers::InverseSquareRoot(Overlap, ISQOverlap, solver_parameters);

        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "ISQOverlap is done");

        // Solve the Density Matrix.
        // Change the solver variable for computing the density matrix.
        solver_parameters.SetConvergeDiff(converge_density);
        const double spin_degeneracy=static_cast<double>(nspin);
        const double trace=nelec/spin_degeneracy;
        NTPoly::DensityMatrixSolvers::TRS2(Hamiltonian, ISQOverlap, trace, 
                        Density, energy, chemical_potential, solver_parameters);
        Density.Scale(spin_degeneracy);
        
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Density Matrix is done");
        // convert DM from the PSMatrix to a BCD matrix
        constructBCDFromPSMatrix(Density, comm_2D, desc, nrow, ncol, DM);

        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Density Matrix is converted to BCD format");

        // Solve the Energy Density Matrix
        NTPoly::DensityMatrixSolvers::EnergyDensityMatrix(Hamiltonian, Density, EnergyDensity, threshold);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "EnergyDensity Matrix is done");

        constructBCDFromPSMatrix(EnergyDensity, comm_2D, desc, nrow, ncol, EDM);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "EnergyDensity Matrix is converted to BCD format");
        return 0;
    }

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
        const int nrow, const int ncol, const double M[], const double threshold)
    {
        // init PSMatrix
        const int nFull=desc[2];    
        //PSM.Resize(nFull);
        //ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "PSM is resized to", PSM.GetActualDimension());

        // read all non-zero values from BCD matrix into a tripletlist
        NTPoly::TripletList_r tripletList;
        if(comm_2D != MPI_COMM_NULL)
        {
            readTripletListFromBCD(tripletList, comm_2D, desc, nrow, ncol, M, threshold);
        }
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "the BCD Matrix is converted to tripletList, non-zero elements are:", tripletList.GetSize());
        // fill PSMatrix from tripletlist
        PSM.FillFromTripletList(tripletList);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "the PSMatrix is filled from tripletList, size is", PSM.GetSize());
        return 0;
    }

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
        const int nrow, const int ncol, const double M[], const double threshold)
    {
        int blacs_context=desc[1];
        const int nblk=desc[4];
        int nprow, npcol, myprow, mypcol;
        Cblacs_gridinfo(blacs_context, &nprow, &npcol, &myprow, &mypcol);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "enter readTripletListFromBCD, initial tripletList size is", tripletList.GetSize());
        // read all non-zero values from BCD matrix into a tripletlist
        NTPoly::Triplet_r tmp_t;
        for(int i=0; i<ncol; ++i)
        {
            tmp_t.index_column=globalIndex(i, nblk, npcol, mypcol);
            for(int j=0; j<nrow; ++j)
            {
                const int idx=i*nrow+j;
                const double val=M[idx];
                if(std::abs(val)<threshold) continue;
                tmp_t.index_row=globalIndex(j, nblk, nprow, myprow);
                tmp_t.point_value=val;
                tripletList.Append(tmp_t);
            }
        }
        return 0;
    }

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
        const int nrow, const int ncol, double M[])
    {
        int blacs_context=desc[1];
        const int nblk=desc[4];
        int nprow, npcol, myprow, mypcol;
        Cblacs_gridinfo(blacs_context, &nprow, &npcol, &myprow, &mypcol);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "enter constructBCDFromPSMatrix, nblk is", nblk);
        // gather matrix elements of current process to a tripletlist from the PSMatrix
        // and then fill the BCD matrix
        NTPoly::TripletList_r tripletList;
        for(int i=0; i<nrow; i+=nblk)
        {
            const int start_row=i;
            const int end_row=std::min(i+nblk, nrow);
            for(int j=0; j<ncol; j+=nblk)
            {
                const int start_col=j;
                const int end_col=std::min(j+nblk, ncol);
                PSM.GetMatrixBlock(tripletList, start_row, end_row, start_col, end_col);

                // fill the BCD matrix
                for(int k=0; k<tripletList.GetSize(); ++k)
                {
                    const NTPoly::Triplet_r tmp_t=tripletList.GetTripletAt(k);
                    const int gRow=tmp_t.index_row;
                    const int gCol=tmp_t.index_column;
                    const int lRow=localIndex(gRow, nblk, nprow, myprow);
                    const int lCol=localIndex(gCol, nblk, npcol, mypcol);
                    const int idx=lRow+lCol*nrow;
                    M[idx]=tmp_t.point_value;
                }
            }
        }
        return 0;
    }
}
#endif // NTPOLY