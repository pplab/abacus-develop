//#ifdef __NTPOLY
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
    namespace{ // Use an anonymous namespace for file-local toolkits
        
        MPI_Comm comm_2D_slice = MPI_COMM_NULL;
        std::vector<int> rank_slice;
        int nproc_slice = 0;
        bool require_init_NTPOLY=true;
        bool require_init_comm_slice = true;

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
        inline int globalIndex(const int localIndex, const int nblk, const int nprocs, const int myproc)
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
         * @param localProc The rank of the process row or column to which the local index belongs.
         * @return The local index corresponding to the given global index.
         */
        inline int localIndex(const int globalIndex, const int nblk, const int nprocs, int& localProc)
        {
            localProc = int((globalIndex % (nblk * nprocs)) / nblk);
            return int(globalIndex / (nblk * nprocs)) * nblk + globalIndex % nblk;
        }

        /**
         * @brief Calculates the global process number based on the BLACS grid layout and process coordinates.
         *
         * This function computes the global process number within a BLACS grid given the layout (row-major or column-major),
         * the BLACS context, and the process row and column coordinates. It also performs error checking on the layout
         * and the resulting process number.
         * 
         * *NOTE* The official BLACS function Cblacs_pnum seems not work correctly for the column-major layout.
         *        Therefore this function implements a workaround to calculate the process number.
         *
         * @param BLACS_LAYOUT The layout of the BLACS grid, 'R' or 'r' for row-major, 'C' or 'c' for column-major.
         * @param blacs_ctxt The BLACS context handle identifying the process grid.
         * @param prow The row coordinate of the process in the grid.
         * @param pcol The column coordinate of the process in the grid.
         * @return The global process number if successful, -1 if an error occurs.
         */
        inline int my_pnum(const char BLACS_LAYOUT, const int blacs_ctxt, const int prow, const int pcol)
        {
            int nprow, npcol, myprow, mypcol;
            Cblacs_gridinfo(blacs_ctxt, &nprow, &npcol, &myprow, &mypcol);
            int pnum=-1;
            if(BLACS_LAYOUT == 'R' || BLACS_LAYOUT == 'r')
                pnum= pcol + prow * npcol;
            else if(BLACS_LAYOUT == 'C' || BLACS_LAYOUT == 'c')
                pnum= prow + pcol * nprow;
            else
            {   
                std::cerr << "Error: Invalid BLACS layout specified. Use 'R' for row-major or 'C' for column-major." << std::endl;
                return -1; // Error code
            }
            if(pnum < 0 || pnum >= nprow * npcol)
            {
                std::cerr << "Error: Process number '" << pnum << "' is out of bounds." << std::endl;
                return -1; // Error code
            }
            return pnum;
        }

        /**
         * @brief Extracts information from a specific triplet in a TripletList.
         * 
         * This function retrieves a triplet from the given TripletList at the specified index,
         * calculates the global and local indices of the triplet, determines the receiving process rank,
         * and then obtains the rank within the slice corresponding to the receiving process.
         * 
         * @param blacs_context The BLACS context handle used to identify the process grid.
         * @param rank_slice An array mapping global process ranks to ranks within a slice.
         * @param nprow The number of process rows in the grid.
         * @param npcol The number of process columns in the grid.
         * @param nblk The number of elements in each block.
         * @param idx The index of the triplet in the TripletList to extract information from.
         * @param TL A constant reference to the TripletList object containing the triplets.
         * @param recv_rank_slice A reference to an integer where the rank within the slice of the receiving process will be stored.
         */
        inline void extractTripletInfo(const char BLACS_LAYOUT, const int blacs_context, const int rank_slice[], 
                const int nprow, const int npcol, const int nblk,
                const int idx, const NTPoly::TripletList_r & TL, 
                int& recv_rank_slice)
        {
            // Retrieve a triplet from the TripletList at the specified index
            NTPoly::Triplet_r tmp_t=TL.GetTripletAt(idx);
            // Calculate the global row index, adjusting for 0-based indexing
            int global_row_idx=tmp_t.index_row-1;
            // Calculate the global column index, adjusting for 0-based indexing
            int global_col_idx=tmp_t.index_column-1;
            // Variables to store the local process row and column numbers
            int local_prow, local_pcol;
            // Calculate the local row index and the local process row number
            int local_row_idx=localIndex(global_row_idx, nblk, nprow, local_prow);
            // Calculate the local column index and the local process column number
            int local_col_idx=localIndex(global_col_idx, nblk, npcol, local_pcol);
            // Determine the global rank of the receiving process using the BLACS context
            int recv_rank=my_pnum(BLACS_LAYOUT, blacs_context, local_prow, local_pcol);
            // Get the rank within the slice corresponding to the receiving process
            recv_rank_slice=rank_slice[recv_rank];
        }

        /**
         * @brief Extracts detailed information from a specific triplet in a TripletList.
         * 
         * This function retrieves a triplet from the given TripletList at the specified index.
         * It then calculates the global row and column indices of the triplet, 
         * determines the local process row and column numbers, finds the receiving process rank,
         * and extracts the value of the triplet.
         * 
         * @param blacs_context The BLACS context handle used to identify the process grid.
         * @param rank_slice An array mapping global process ranks to ranks within a slice.
         * @param nprow The number of process rows in the grid.
         * @param npcol The number of process columns in the grid.
         * @param nblk The number of elements in each block.
         * @param idx The index of the triplet in the TripletList to extract information from.
         * @param TL A constant reference to the TripletList object containing the triplets.
         * @param recv_rank_slice A reference to an integer where the rank within the slice of the receiving process will be stored.
         * @param global_row_idx A reference to an integer where the global row index of the triplet will be stored.
         * @param global_col_idx A reference to an integer where the global column index of the triplet will be stored.
         * @param val A reference to a double where the value of the triplet will be stored.
         */
        inline void extractTripletInfo(const char BLACS_LAYOUT, const int blacs_context, const int rank_slice[], 
                const int nprow, const int npcol, const int nblk,
                const int idx, const NTPoly::TripletList_r & TL, 
                int& recv_rank_slice, int& local_row_idx, int& local_col_idx, double& value)
        {
            // Retrieve a triplet from the TripletList at the specified index
            NTPoly::Triplet_r tmp_t=TL.GetTripletAt(idx);
            // Calculate the global row index of the triplet, adjusting for 0-based indexing
            int global_row_idx=tmp_t.index_row-1;
            // Calculate the global column index of the triplet, adjusting for 0-based indexing
            int global_col_idx=tmp_t.index_column-1;
            // Variables to store the local process row and column numbers
            int local_prow, local_pcol;
            // Calculate the local row index and the local process row number
            local_row_idx=localIndex(global_row_idx, nblk, nprow, local_prow);
            // Calculate the local column index and the local process column number
            local_col_idx=localIndex(global_col_idx, nblk, npcol, local_pcol);
            // Determine the global rank of the receiving process using the BLACS context
            int recv_rank=my_pnum(BLACS_LAYOUT, blacs_context, local_prow, local_pcol);
            
            // check process grid information
            static bool _check=true;
            if(_check)
            {
                for(int i=0; i<nprow; ++i)
                {
                    for(int j=0; j<npcol; ++j)
                    {
                        int rank=Cblacs_pnum(blacs_context, i, j);
                        outlog("Process at row " + std::to_string(i) + 
                               ", column " + std::to_string(j) + 
                               " has global rank " + std::to_string(rank));
                    }
                }
                _check=false;
            }
            // Get the rank within the slice corresponding to the receiving process
            recv_rank_slice=rank_slice[recv_rank];
            // Extract the value of the triplet
            value=tmp_t.point_value;
        }
        /**
         * @brief Creates a communicator for each process slice and establishes a mapping between the global communicator and local slice communicators.
         * 
         * This function either directly uses the input 2D communicator if the global number of slices is 1, 
         * or splits it into multiple slice communicators when the number of slices is greater than 1. 
         * It calculates the rank of each process within the slice and stores it in the `rank_slice` vector.
         * 
         * @param comm_2D The input 2D MPI communicator.
         * @param my_slice The slice number to which the current process belongs.
         * @param nproc The total number of processes in the input 2D communicator.
         * @param myid The rank of the current process in the input 2D communicator.
         */
        inline void createCommSlice(const MPI_Comm comm_2D, const int my_slice, const int nproc, const int myid)
        {  
            // Get the group of the input 2D communicator
            MPI_Group group_2D;
            MPI_Comm_group(comm_2D, &group_2D);

            // If the global number of slices is 1, directly use the input 2D communicator
            if(NTPoly::GetGlobalNumSlices()==1)
            {
                // Assign the 2D communicator to the slice communicator
                comm_2D_slice=comm_2D;
                // Get the number of processes in the slice communicator
                MPI_Comm_size(comm_2D_slice, &nproc_slice);
                // Resize the rank_slice vector to accommodate all processes
                rank_slice.resize(nproc);             
                // Initialize the rank_slice vector. Each process's slice rank equals its rank in the 2D communicator.
                for(int i=0; i<nproc; ++i)
                {
                    rank_slice[i]=i;
                }   
            }
            else{
                // When the global number of slices is greater than 1, create communicators and groups within each slice
                MPI_Group group_2D_slice = MPI_GROUP_NULL;
                // Split the input 2D communicator based on the slice number and process rank to create a slice communicator
                MPI_Comm_split(comm_2D, my_slice, myid, &comm_2D_slice);
                // Get the group of the slice communicator
                MPI_Comm_group(comm_2D_slice, &group_2D_slice);
                // Get the number of processes in the slice communicator
                MPI_Comm_size(comm_2D_slice, &nproc_slice);

                // Debug information, output relevant communicator and group information
                if(for_debug)
                {
                    outlog("Input comm_2D is", MPI_Comm_c2f(comm_2D));
                    outlog("Input my_slice is", my_slice);
                    outlog("Input myid is", myid);
                    outlog("comm_2D_slice is constructed, comm_2D_slice is", MPI_Comm_c2f(comm_2D_slice));
                    outlog("group_2D_slice is constructed, group_2D_slice is", MPI_Group_c2f(group_2D_slice));
                    outlog("group_2D is constructed, group_2D is", MPI_Group_c2f(group_2D));
                }

                // Initialize the rank_slice vector and set all elements to -1
                rank_slice.resize(nproc);
                for(int i=0; i<nproc; ++i)
                {
                    rank_slice[i]=-1;
                }

                // Store the ranks of processes within the slice
                int rank_list_slice[nproc_slice];
                // Store the ranks of processes within the slice in the global communicator
                int glocal_rank_list_slice[nproc_slice];
                // Initialize the rank list of processes within the slice
                for(int i=0; i<nproc_slice; ++i)
                {
                    rank_list_slice[i]=i;
                }            
                // Translate the ranks of processes within the slice to ranks in the global communicator
                MPI_Group_translate_ranks(group_2D_slice, nproc_slice, rank_list_slice, group_2D, glocal_rank_list_slice);

                // Debug information, indicating that the MPI_Group_translate_ranks function has been called
                if(for_debug) 
                {
                    outlog("MPI_Group_translate_ranks is called");
                }

                // Establish a mapping from global communicator ranks to ranks within the slice
                for(int i=0; i<nproc_slice; ++i)
                {
                    rank_slice[glocal_rank_list_slice[i]]=i;
                    outlog("rank_slice["+std::to_string(glocal_rank_list_slice[i])+ "] ="+ std::to_string(i));
                }

                // Debug information, output the number of processes in the current slice
                if(for_debug) 
                {
                    outlog("All processes are split to slices, nproc in current slice is", nproc_slice);
                }
            }
        }
    }

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
            outlog("nrow", nrow);
            outlog("ncol", ncol);
            int myid;
            MPI_Comm_rank(comm_2D, &myid);
            if(myid == 0)
            {
                saveParametersToFile("parameters.dat", nFull, nelec, nspin, converge_density, converge_overlap, threshold);
            }            
            // int saveBCDMatrixToFile(const MPI_Comm comm, const int* desc, const int nrow, const int ncol, double* matrix, const std::string& filename)
            saveBCDMatrixToFile(comm_2D, desc, nrow, ncol, H, "H.dat");
            saveBCDMatrixToFile(comm_2D, desc, nrow, ncol, S, "S.dat");
            // check desc array
            outlog("Descriptor array:");
            for(int i=0; i<9; ++i)
            {
                outlog("desc[" + std::to_string(i) + "] = ", desc[i]);
            }
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
            Hamiltonian.WriteToMatrixMarket("Hamiltonian.mtx");
            Overlap.WriteToMatrixMarket("Overlap.mtx");
        }

        // set permutation
        const int perm_dim = Hamiltonian.GetLogicalDimension();
        outlog("Permutation dimension: ", perm_dim);
        if (perm_dim <= 0) {
            outlog("Error: Permutation dimension is non-positive!", perm_dim);
            return -1;
        }
        NTPoly::Permutation permutation(perm_dim);
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
        
        if(for_debug) //ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Density Matrix is done");
        {
            outlog("Density Matrix is done");
            Density.WriteToMatrixMarket("DM.mtx");
        }
        // convert DM from the PSMatrix to a BCD matrix
        constructBCDFromPSMatrix(Density, comm_2D, 'C', desc, nrow, ncol, DM);

        if(for_debug)
        {
            //ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Density Matrix is converted to BCD format");
            outlog("Density Matrix is converted to BCD format");
            //saveBCDMatrixToFile(comm_2D, desc, nrow, ncol, DM, "DM.dat");
            //saveMatrixToFile("DM", DM, nrow, ncol);
            MPI_Barrier(comm_2D);
            outlog("DM is saved to file DM.dat");
        }

        // Solve the Energy Density Matrix
        if(for_debug) //ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "EnergyDensity Matrix is done");
        {
            outlog("threshold for EMD is", threshold);
            Hamiltonian.WriteToMatrixMarket("forEDM_Hamiltonian.mtx");
            Density.WriteToMatrixMarket("forEDM_Density.mtx");
        }
        NTPoly::DensityMatrixSolvers::EnergyDensityMatrix(Hamiltonian, Density, EnergyDensity, threshold);
        if(for_debug) //ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "EnergyDensity Matrix is done");
        {
            outlog("EnergyDensity Matrix is done");
            EnergyDensity.WriteToMatrixMarket("EDM.mtx");
        }

        constructBCDFromPSMatrix(EnergyDensity, comm_2D, 'C', desc, nrow, ncol, EDM);
        if(for_debug) 
        {
            //ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "EnergyDensity Matrix is converted to BCD format");            
            outlog("EnergyDensity Matrix is converted to BCD format");
            //saveBCDMatrixToFile(comm_2D, desc, nrow, ncol, EDM, "EDM.dat");
            //saveMatrixToFile("EDM", EDM, nrow, ncol);
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
        const MPI_Comm comm_2D, const char layout, const int desc[],
        const int nrow, const int ncol, double M[])
    {
        int blacs_context=desc[1];
        const int nFull=desc[2];
        const int nblk=desc[4];
        int nprow, npcol, myprow, mypcol;
        Cblacs_gridinfo(blacs_context, &nprow, &npcol, &myprow, &mypcol);

        // matrix are transformed within the process slice
        // setup slice parameters
        const int my_slice=NTPoly::GetGlobalMySlice();
        int myid, nproc;
        MPI_Comm_size(comm_2D, &nproc);
        MPI_Comm_rank(comm_2D, &myid);
        if(for_debug) // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, 
        {
            outlog("Enter constructBCDFromPSMatrix, my slice is", my_slice);
            outlog("M's size is", nrow*ncol);
        }
        if(require_init_comm_slice)
        {
            // init process slice communicator
            if(for_debug) 
            {
                outlog("Initializing process slice communicator");
            }
            createCommSlice(comm_2D, my_slice, nproc, myid);
            if(for_debug)
            {
                outlog("Process slice communicator initialized");
                outlog("comm_2D_slice is", MPI_Comm_c2f(comm_2D_slice));
            }
            require_init_comm_slice=false;
        }
        // transform PSMatrix to local tripletlist
        int n_send_element; // number of elements to be sent
        int n_recv_element; // number of elements to be received
        NTPoly::TripletList_r local_tripletList;
        PSM.GetTripletList(local_tripletList);
        n_send_element=local_tripletList.GetSize();
        if(for_debug) 
        {
            std::string local_tripletList_filename="local_input_tripletList_"+std::to_string(myid)+".txt";
            saveTripletListToFile(local_tripletList, local_tripletList_filename);
        }
        // count number of elements to be sent to each process
        std::vector<int> send_count(nproc_slice, 0);
        NTPoly::Triplet_r tmp_t;
        for(int i=0; i<n_send_element; ++i) 
        {
            int recv_rank_slice;
            extractTripletInfo(layout, blacs_context, rank_slice.data(),
                nprow, npcol, nblk, i, local_tripletList,
                recv_rank_slice);
            if(recv_rank_slice>=0) send_count[recv_rank_slice]++;
        }
        // build sender and receiver parameters for mpi_alltoallv
        std::vector<int> recv_count(nproc_slice, 0);
        // if(for_debug) 
        // {         
        //     outlog("nproc_slice is", nproc_slice); 
        //     outlog("comm_2D_slice is", MPI_Comm_c2f(comm_2D_slice));
        //     saveArrayToFile("send_count_before", send_count.data(), nproc_slice);
        //     saveArrayToFile("recv_count_before", recv_count.data(), nproc_slice);
        // }
        MPI_Alltoall(send_count.data(), 1, MPI_INT, recv_count.data(), 1, MPI_INT, comm_2D_slice);
        //MPI_Alltoall(&send_count[0], 1, MPI_INT, &recv_count[0], 1, MPI_INT, comm_2D_slice);
        if(for_debug) 
        {         
            outlog("nproc_slice is", nproc_slice); 
            saveArrayToFile("send_count", send_count.data(), nproc_slice);
            saveArrayToFile("recv_count", recv_count.data(), nproc_slice);
        }
        std::vector<int> send_displ(nproc_slice);
        std::vector<int> recv_displ(nproc_slice);
        send_displ[0]=0;
        recv_displ[0]=0;
        n_recv_element=recv_count[0];
        for(int i=1; i<nproc_slice; ++i)
        {
            send_displ[i]=send_displ[i-1]+send_count[i-1];
            recv_displ[i]=recv_displ[i-1]+recv_count[i-1];
            n_recv_element+=recv_count[i];
        }
        if(for_debug) 
        {
            outlog("n_send_element is", n_send_element);
            outlog("n_recv_element is", n_recv_element);            
            saveArrayToFile("send_displ", send_displ.data(), nproc_slice);
            saveArrayToFile("recv_displ", recv_displ.data(), nproc_slice);
        }

        // fill local elements and their index to send_data, send_row_index and send_col_index
        std::vector<double> send_data(n_send_element);
        std::vector<int> send_row_index(n_send_element);
        std::vector<int> send_col_index(n_send_element);
        std::vector<int> p_fill_to_send(nproc_slice, 0);
        
        // check the maximum index of each process grid
        std::vector<int> _max_row_idx(nproc_slice, -1);
        std::vector<int> _max_col_idx(nproc_slice, -1);
        for(int i=0; i<n_send_element; ++i)
        {
            int recv_rank_slice, local_row_idx, local_col_idx;
            double val;
            extractTripletInfo(layout,blacs_context, rank_slice.data(),
                nprow, npcol, nblk, i, local_tripletList,
                recv_rank_slice, local_row_idx, local_col_idx, val);
            if(recv_rank_slice>=0)
            {
                int t_idx=send_displ[recv_rank_slice]+p_fill_to_send[recv_rank_slice];
                p_fill_to_send[recv_rank_slice]++;
                send_data[t_idx]=val;
                send_row_index[t_idx]=local_row_idx;
                send_col_index[t_idx]=local_col_idx;

                if(local_row_idx>_max_row_idx[recv_rank_slice]) _max_row_idx[recv_rank_slice]=local_row_idx;
                if(local_col_idx>_max_col_idx[recv_rank_slice]) _max_col_idx[recv_rank_slice]=local_col_idx;
            }
        }
        outlog("max row index and col index in each process slice:");
        for(int i=0; i<nproc_slice; ++i)
        {
            outlog("Process slice " + std::to_string(i) 
                + ": max row index = " + std::to_string(_max_row_idx[i])
                + ": max col index = " + std::to_string(_max_col_idx[i]));
        }
        if(for_debug) 
        {
            outlog("send_data is filled");
            saveArrayToFile("send_data", send_data.data(), n_send_element);
            saveArrayToFile("send_row_index", send_row_index.data(), n_send_element);
            saveArrayToFile("send_col_index", send_col_index.data(), n_send_element);
        }
        // call MPI_Alltoallv to send elements to each process        
        std::vector<double> recv_data(n_recv_element);
        std::vector<int> recv_row_index(n_recv_element);
        std::vector<int> recv_col_index(n_recv_element);
        MPI_Alltoallv(send_data.data(), send_count.data(), send_displ.data(), MPI_DOUBLE,
            recv_data.data(), recv_count.data(), recv_displ.data(), MPI_DOUBLE, comm_2D_slice);
        MPI_Alltoallv(send_row_index.data(), send_count.data(), send_displ.data(), MPI_INT,
            recv_row_index.data(), recv_count.data(), recv_displ.data(), MPI_INT, comm_2D_slice);
        MPI_Alltoallv(send_col_index.data(), send_count.data(), send_displ.data(), MPI_INT,
            recv_col_index.data(), recv_count.data(), recv_displ.data(), MPI_INT, comm_2D_slice);
        if(for_debug) 
        {
            outlog("elements are exchanged between processes");
            saveArrayToFile("recv_data", recv_data.data(), n_recv_element);
            saveArrayToFile("recv_row_index", recv_row_index.data(), n_recv_element);
            saveArrayToFile("recv_col_index", recv_col_index.data(), n_recv_element);
        }
        // fill received elements index to BCD Matrix
        for(int i=0; i<n_recv_element; ++i)
        {
            int local_row_idx=recv_row_index[i];
            int local_col_idx=recv_col_index[i];
            int local_idx=local_row_idx+local_col_idx*nrow;
            M[local_idx]=recv_data[i];
        }
        if(for_debug) 
        {
            outlog("received elements are filled to BCD Matrix");
        }
        return 0;
    }

    /**
     * Constructs a Block Cyclic Distributed (BCD) matrix from a PSMatrix.
     * Use NTPoly transform function 
     * It has some issues, the matrix is not exactly the same as the original matrix,
     * and the performance is poor for big matrix with small nblk
     *
     * @param PSM The PSMatrix to be converted.
     * @param comm_2D The MPI communicator for the 2D grid.
     * @param desc The descriptor array for the BCD matrix.
     * @param nrow The number of rows in the local matrix.
     * @param ncol The number of columns in the local matrix.
     * @param M The BCD matrix to be filled.
     * @return Returns 0 if successful, or an error code if an error occurs.
     */
    int constructBCDFromPSMatrix_NTPoly(NTPoly::Matrix_ps& PSM, 
        const MPI_Comm comm_2D, const int desc[],
        const int nrow, const int ncol, double M[])
    {
        int blacs_context=desc[1];
        const int nFull=desc[2];  
        const int nblk=desc[4];
        int myid;
        MPI_Comm_rank(comm_2D, &myid);
        int nprow, npcol, myprow, mypcol;
        Cblacs_gridinfo(blacs_context, &nprow, &npcol, &myprow, &mypcol);
        if(for_debug) 
        {
            outlog("enter constructBCDFromPSMatrix, nblk is", nblk);
            outlog("nprow is "+std::to_string(nprow)+" npcol is " +std::to_string(npcol));
            outlog("myprow is "+std::to_string(myprow)+" mypcol is " +std::to_string(mypcol));
            // save PSM to local tripletlist file
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
    
    /**
     * @brief Cleans up MPI resources, specifically the 2D slice communicator.
     * 
     * This function checks if the 2D slice communicator is valid (not MPI_COMM_NULL).
     * If it is valid, it frees the communicator using MPI_Comm_free and sets it to MPI_COMM_NULL
     * to indicate that the communicator is no longer valid.
     */
    void cleanupMPIResources() 
    {
        // Check if the 2D slice communicator is valid
        if (MPI_COMM_NULL != comm_2D_slice) {
            // Free the 2D slice communicator
            MPI_Comm_free(&comm_2D_slice);
            // Set the 2D slice communicator to MPI_COMM_NULL to indicate it's no longer valid
            comm_2D_slice = MPI_COMM_NULL;
        }
    }
}
//#endif // NTPOLY
