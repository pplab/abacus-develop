#ifdef __NTPOLY
#include <iostream>
#include <vector>
#include <string>
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
//#include "module_base/global_function.h"
#include "utils.hpp"

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
        if(for_debug) 
        {
            //ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "enter simple_ntpoly, nFull", nFull);
            outlog("enter simple_ntpoly, nFull", nFull);
            int myid;
            MPI_Comm_rank(comm_2D, &myid);
            if(myid == 0)
            {
                saveParametersToFile("parameters.dat", nFull, nelec, nspin, converge_density, converge_overlap, threshold);
            }            
            // int saveBCDMatrixToFile(const MPI_Comm comm, const int* desc, const int nrow, const int ncol, double* matrix, const std::string& filename)
            saveBCDMatrixToFile(comm_2D, desc, nrow, ncol, H, "H.dat");
            saveBCDMatrixToFile(comm_2D, desc, nrow, ncol, S, "S.dat");
        }
        // init default process grid
        int process_slice=1;
        if(require_init_NTPOLY)
        {
            NTPoly::ConstructGlobalProcessGrid(comm_2D, process_slice);
            require_init_NTPOLY=false;
        }        
        if(for_debug) //ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "GlobalProcessGrid is constructed");
        {
            outlog("GlobalProcessGrid is constructed");
        }

        // init PSMatrices of Hamiltonian, Overlap, ISQOverlap, Density and EnergyDensity
        NTPoly::Matrix_ps Hamiltonian(nFull);
        NTPoly::Matrix_ps Overlap(nFull);  
        NTPoly::Matrix_ps ISQOverlap(nFull);
        NTPoly::Matrix_ps Density(nFull);
        NTPoly::Matrix_ps EnergyDensity(nFull);
        if(for_debug) 
        {
            // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, 
            //     "All PSMatrices are allocated, ActualDimension is", 
            //     Hamiltonian.GetActualDimension());
            // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, 
            //     "LogicalDimension is", 
            //     Hamiltonian.GetLogicalDimension());
            outlog("All PSMatrices are allocated, ActualDimension is", Hamiltonian.GetActualDimension());
            outlog("LogicalDimension is", Hamiltonian.GetLogicalDimension());
        }

        // convert H and S from BCD matrix to PSMatrix
        constructPSMatrixFromBCD(Hamiltonian, comm_2D, desc, nrow, ncol, H, threshold);
        constructPSMatrixFromBCD(Overlap, comm_2D, desc, nrow, ncol, S, threshold);        
        //ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "H and S are converted to PSMatrix");
        if(for_debug)
        {
            outlog("H and S are converted to PSMatrix");
        }

        // set purmutation
        NTPoly::Permutation permutation(Hamiltonian.GetLogicalDimension());
        permutation.SetRandomPermutation();

        if(for_debug) //ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "permutation is done");
        {
            outlog("permutation is done");
        }

        // set solver parameters 
        NTPoly::SolverParameters solver_parameters;
        solver_parameters.SetConvergeDiff(converge_overlap);
        solver_parameters.SetLoadBalance(permutation);
        solver_parameters.SetThreshold(threshold);
        solver_parameters.SetVerbosity(true);
        
        // InverseSquareRoot(Overlap, ISQOverlap, solver_parameters)
        NTPoly::SquareRootSolvers::InverseSquareRoot(Overlap, ISQOverlap, solver_parameters);

        if(for_debug) //ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "ISQOverlap is done");
        {
            outlog("ISQOverlap is done");
        }

        // Solve the Density Matrix.
        // Change the solver variable for computing the density matrix.
        solver_parameters.SetConvergeDiff(converge_density);
        const double spin_degeneracy = nspin==1? 2.0: 1.0;
        const double trace=nelec/spin_degeneracy;
        NTPoly::DensityMatrixSolvers::TRS2(Hamiltonian, ISQOverlap, trace, 
                        Density, energy, chemical_potential, solver_parameters);
        Density.Scale(spin_degeneracy);
		energy *= spin_degeneracy;
        
        if(for_debug) //ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Density Matrix is done");
        {
            outlog("Density Matrix is done");
        }
        // convert DM from the PSMatrix to a BCD matrix
        constructBCDFromPSMatrix(Density, comm_2D, desc, nrow, ncol, DM);

        if(for_debug)
        {
            //ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Density Matrix is converted to BCD format");
            outlog("Density Matrix is converted to BCD format");
            saveBCDMatrixToFile(comm_2D, desc, nrow, ncol, DM, "DM.dat");
        } 

        // Solve the Energy Density Matrix
        NTPoly::DensityMatrixSolvers::EnergyDensityMatrix(Hamiltonian, Density, EnergyDensity, threshold);
        if(for_debug) //ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "EnergyDensity Matrix is done");
        {
            outlog("EnergyDensity Matrix is done");
        }

        constructBCDFromPSMatrix(EnergyDensity, comm_2D, desc, nrow, ncol, EDM);
        if(for_debug) 
        {
            //ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "EnergyDensity Matrix is converted to BCD format");            
            outlog("EnergyDensity Matrix is converted to BCD format");
            saveBCDMatrixToFile(comm_2D, desc, nrow, ncol, EDM, "EDM.dat");
        }
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
        if(for_debug) 
        {
            // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, 
            // "the BCD Matrix is converted to tripletList, non-zero elements are:", tripletList.GetSize());
            outlog("the BCD Matrix is converted to tripletList, non-zero elements are:", tripletList.GetSize());
            // std::string local_tripletList_filename="local_tripletList.txt";
            int myid;
            MPI_Comm_rank(comm_2D, &myid);
            std::string local_tripletList_filename="local_tripletList_"+std::to_string(myid)+".txt";
            saveTripletListToFile(tripletList, local_tripletList_filename);
        }
        // fill PSMatrix from tripletlist
        PSM.FillFromTripletList(tripletList);
        if(for_debug) // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, 
            // "the PSMatrix is filled from tripletList, size is", PSM.GetSize());
        {
            outlog("the PSMatrix is filled from tripletList, size is", PSM.GetSize());
        }
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
        if(for_debug) // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, 
            //"enter readTripletListFromBCD, initial tripletList size is", tripletList.GetSize());
        {
            outlog("enter readTripletListFromBCD, initial tripletList size is", tripletList.GetSize());
        }
        // read all non-zero values from BCD matrix into a tripletlist
        NTPoly::Triplet_r tmp_t;
        for(int i=0; i<ncol; ++i)
        {
            tmp_t.index_column=globalIndex(i, nblk, npcol, mypcol)+1;
            for(int j=0; j<nrow; ++j)
            {
                const int idx=i*nrow+j;
                const double val=M[idx];
                if(std::abs(val)<threshold) continue;
                tmp_t.index_row=globalIndex(j, nblk, nprow, myprow)+1;
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
        const int nFull=desc[2];  
        const int nblk=desc[4];
        int nprow, npcol, myprow, mypcol;
        Cblacs_gridinfo(blacs_context, &nprow, &npcol, &myprow, &mypcol);
        if(for_debug) //ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, 
            //"enter constructBCDFromPSMatrix, nblk is", nblk);
        {
            outlog("enter constructBCDFromPSMatrix, nblk is", nblk);
            outlog("nprow is "+std::to_string(nprow)+" npcol is " +std::to_string(npcol));
            outlog("myprow is "+std::to_string(myprow)+" mypcol is " +std::to_string(mypcol));
            // save PSM to local tripletlist file
            int myid;
            MPI_Comm_rank(comm_2D, &myid);
            NTPoly::TripletList_r local_DM_tripletList;
            PSM.GetTripletList(local_DM_tripletList);
            std::string local_tripletList_filename="local_DM_tripletList_"+std::to_string(myid)+".txt";
            saveTripletListToFile(local_DM_tripletList, local_tripletList_filename);
        }
        // count how many times the getMatrixBlock should be called
        int max_count=0;
        int my_count=((nrow-1)/nblk+1)*((ncol-1)/nblk+1);
        MPI_Allreduce(&my_count, &max_count, 1, MPI_INT, MPI_MAX, comm_2D);

        // gather matrix elements of current process to a tripletlist from the PSMatrix
        // and then fill the BCD matrix
        NTPoly::TripletList_r tripletList;
        for(int i=0; i<nrow; i+=nblk)
        {
            const int start_row=globalIndex(i, nblk, nprow, myprow); // start_row is the global index
            const int end_row=std::min(start_row+nblk, nFull);
            for(int j=0; j<ncol; j+=nblk)
            {
                const int start_col=globalIndex(j, nblk, npcol, mypcol); // start_col is the global index
                const int end_col=std::min(start_col+nblk, nFull);
                // The c++ interface of the NTPoly has already transformed the matrix index into Fortran format, which is one-based.
                // ref: PSMatrix.cc, line 132
                PSM.GetMatrixBlock(tripletList, start_row, end_row, start_col, end_col);                 
                if(for_debug)
                {
                    outlog("GetMatrixBlock, i is "+std::to_string(i)+" j is "+std::to_string(j));
                    outlog("GetMatrixBlock, start_row is "+std::to_string(start_row)+" end_row is "+std::to_string(end_row));
                    outlog("GetMatrixBlock, start_col is "+std::to_string(start_col)+" end_col is "+std::to_string(end_col));
                    outlog("GetMatrixBlock, number non-zero element is", tripletList.GetSize());
                }
                // fill the BCD matrix
                for(int k=0; k<tripletList.GetSize(); ++k)
                {
                    const NTPoly::Triplet_r tmp_t=tripletList.GetTripletAt(k);
                    const int gRow=tmp_t.index_row-1;
                    const int gCol=tmp_t.index_column-1;
                    const int lRow=localIndex(gRow, nblk, nprow, myprow);
                    const int lCol=localIndex(gCol, nblk, npcol, mypcol);
                    const int idx=lRow+lCol*nrow;
                    if(for_debug)
                    {
                        // outlog("tripletList tmp_t["+std::to_string(gRow)+", "+std::to_string(gCol)+"] = "+std::to_string(tmp_t.point_value));
                        // outlog("fill BCD Matrix, M("+std::to_string(lRow)+", "+std::to_string(lCol)+
                        //      ")("+std::to_string(idx)+")= "+std::to_string(tmp_t.point_value));
                        outlog("tmp_t["+std::to_string(gRow)+", "+std::to_string(gCol)+"] = "+std::to_string(tmp_t.point_value)+
                                " ==> M("+std::to_string(lRow)+", "+std::to_string(lCol)+")(idx: "+std::to_string(idx)+")= "+std::to_string(tmp_t.point_value));
                    }
                    M[idx]=tmp_t.point_value;
                }
            }
        }
        // some processes may still require some elements to be filled
        for(int i=0; i<max_count-my_count; ++i)
        {
            const int start_row=1;
            const int end_row=1;
            const int start_col=1;
            const int end_col=1;
            PSM.GetMatrixBlock(tripletList, start_row, end_row, start_col, end_col);
        }
        return 0;
    }
}
#endif // NTPOLY
