// utility tools for debug
#include <string>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <TripletList.h>
#include <Triplet.h>
#include "Cblacs.h"
#include "utils.hpp"
extern "C"
{
    #include "scalapack.h" 
}

/**
 * @brief Initializes a BLACS grid and computes the number of local rows and columns for a matrix.
 * 
 * This function initializes a BLACS grid based on the given MPI communicator and parameters.
 * It calculates the number of process rows and columns in the grid, initializes the BLACS context,
 * and computes the number of local rows and columns for a matrix distributed in a block - cyclic fashion.
 * Finally, it initializes the matrix descriptor array.
 * 
 * @param comm The MPI communicator used for communication between processes.
 * @param nFull The total number of rows or columns in the global matrix.
 * @param nblk The block size used for block - cyclic distribution.
 * @param blacs_ctxt A reference to an integer where the BLACS context handle will be stored.
 * @param narows A reference to an integer where the number of local rows will be stored.
 * @param nacols A reference to an integer where the number of local columns will be stored.
 * @param desc A pointer to an array of integers used to store the matrix descriptor.
 */
void initBlacsGrid(MPI_Comm comm, int nFull, int nblk,
                   int& blacs_ctxt, int& narows, int& nacols, int* desc)
{
    char BLACS_LAYOUT='C';
    int ISRCPROC=0; 
    int nprows, npcols;
    int myprow, mypcol;
    int nprocs, myid;
    int info;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &myid);
    // set blacs parameters
    for(npcols=int(sqrt(double(nprocs))); npcols>=2; --npcols)
    {
        if(nprocs%npcols==0) break;
    }
    nprows=nprocs/npcols;
    outlog("nprows", nprows);
    outlog("npcols", npcols);

    //int comm_f = MPI_Comm_c2f(comm);
    blacs_ctxt=Csys2blacs_handle(comm);
    Cblacs_gridinit(&blacs_ctxt, &BLACS_LAYOUT, nprows, npcols);
    Cblacs_gridinfo(blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);

    narows=numroc_(&nFull, &nblk, &myprow, &ISRCPROC, &nprows);
    outlog("narows", narows);
    nacols=numroc_(&nFull, &nblk, &mypcol, &ISRCPROC, &npcols);
    outlog("nacols", nacols);
    descinit_(desc, &nFull, &nFull, &nblk, &nblk, &ISRCPROC, &ISRCPROC, &blacs_ctxt, &narows, &info);
    if(true)
    {
        outlog("BLACS context initialized", blacs_ctxt);
        outlog("Descriptor initialized with nFull " +std::to_string(nFull) + 
                                          ", nblk " + std::to_string(nblk) + 
                                          ", narows " + std::to_string(narows) + 
                                          ", blacs_ctxt " + std::to_string(blacs_ctxt));
        outlog("The desc array is: ");
        for(int i=0; i<9; ++i)
        {
            outlog(std::to_string(desc[i]));
        }
        outlog("Check process grid layout by pnum:");
        for(int i=0; i<nprows; ++i)
        {
            for(int j=0; j<npcols; ++j)
            {
                int rank=Cblacs_pnum(blacs_ctxt, i, j);
                outlog("process at row " + std::to_string(i) + 
                        ", column " + std::to_string(j) + 
                        " has global rank " + std::to_string(rank));
            }
        }
        outlog("Check process grid layout by pcoord:");
        for(int i=0; i<nprocs; ++i)
        {
            int prow, pcol;
            Cblacs_pcoord(blacs_ctxt, i, &prow, &pcol);
            outlog("process with global rank " + std::to_string(i) + 
                   " is at row " + std::to_string(prow) + 
                   ", column " + std::to_string(pcol));
        }
        outlog("Check current process rank and grid info:");
        outlog("My process rank is " + std::to_string(myid));
        outlog("My process row is " + std::to_string(myprow) + 
               ", my process column is " + std::to_string(mypcol));
    }
}

/**
 * @brief Saves the specified parameters to a file.
 * 
 * This function writes various parameters related to matrix calculations and simulations to a file.
 * Each parameter is written on a separate line with a descriptive key, making the file human - readable.
 * 
 * @param filename The name of the file to which the parameters will be saved.
 * @param nFull The total number of rows or columns in the matrix.
 * @param nelec The number of electrons.
 * @param nspin The number of spins.
 * @param converge_density The density convergence threshold.
 * @param converge_overlap The overlap convergence threshold.
 * @param threshold The threshold for other calculations.
 * @return Returns 0 if the file is successfully opened and written to; otherwise, returns 1.
 */
int saveParametersToFile(const std::string& filename, 
        const int nFull, const int nelec, const int nspin, 
        const double converge_density, const double converge_overlap, const double threshold)   
{
    std::ofstream outfile(filename);
    if (!outfile.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return 1;
    }
    outfile << "nFull: " << nFull << std::endl;
    outfile << "nelec: " << nelec << std::endl;
    outfile << "nspin: " << nspin << std::endl;
    outfile << "converge_density: " << converge_density << std::endl;
    outfile << "converge_overlap: " << converge_overlap << std::endl;
    outfile << "threshold: " << threshold << std::endl;
    outfile.close();
    return 0;
}

/**
 * @brief Loads parameters from a file into the provided variables.
 * 
 * This function reads a file containing parameters related to matrix calculations and simulations.
 * It parses each line of the file, extracts the parameter values based on their descriptive keys,
 * and stores these values in the provided variables.
 * 
 * @param filename The name of the file from which the parameters will be loaded.
 * @param nFull A reference to an integer where the total number of rows or columns in the matrix will be stored.
 * @param nelec A reference to an integer where the number of electrons will be stored.
 * @param nspin A reference to an integer where the number of spins will be stored.
 * @param converge_density A reference to a double where the density convergence threshold will be stored.
 * @param converge_overlap A reference to a double where the overlap convergence threshold will be stored.
 * @param threshold A reference to a double where the threshold for other calculations will be stored.
 * @return Returns 0 if the file is successfully opened and parsed; otherwise, returns 1.
 */
int loadParametersFromFile(const std::string& filename,
        int& nFull, int& nelec, int& nspin, 
        double& converge_density, double& converge_overlap, double& threshold)
{
    std::ifstream infile(filename);
    if (!infile.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return 1;
    }
    std::string line;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        std::string key;
        if (std::getline(iss, key, ':'))
        {
            if (key == "nFull")
            {
                iss >> nFull;
            }
            else if (key == "nelec")
            {
                iss >> nelec;
            }
            else if (key == "nspin")
            {
                iss >> nspin;
            }
            else if (key == "converge_density")
            {
                iss >> converge_density;
            }
            else if (key == "converge_overlap")
            {
                iss >> converge_overlap;
            }
            else if (key == "threshold")
            {
                iss >> threshold;
            }
        }
    }
    infile.close();
    return 0;
}

/**
 * @brief Saves a TripletList to a file.
 * 
 * This function writes each triplet in the provided TripletList to a specified file.
 * Each line in the file corresponds to a triplet, containing the row index, column index,
 * and the value of the triplet, separated by spaces.
 * 
 * @param tripletList A constant reference to the TripletList object to be saved.
 * @param filename The name of the file to which the TripletList will be saved.
 * @return Returns 0 if the file is successfully opened and written to; otherwise, returns 1.
 */
int saveTripletListToFile(const NTPoly::TripletList_r& tripletList, const std::string& filename)
{
    std::ofstream outfile(filename);
    if (!outfile.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return 1;
    }
    const int N=tripletList.GetSize();
    for (int i=0; i<N; ++i)
    {
        const NTPoly::Triplet_r t=tripletList.GetTripletAt(i);
        outfile << t.index_row << " " << t.index_column << " " << t.point_value << std::endl;
    }
    outfile.close();
    return 0;
}

/**
 * @brief Saves a local matrix (column major!) to a file.
 * 
 * This function writes a local column-major matrix to a specified file. Each row of the matrix
 * is written on a separate line, with elements separated by spaces.
 * 
 * @param narows The number of rows in the matrix.
 * @param nacols The number of columns in the matrix.
 * @param matrix A pointer to the local matrix data stored in column-major order.
 * @param filename The name of the file to which the matrix will be saved.
 * @return Returns 0 if the file is successfully opened and written to; otherwise, returns 1.
 */
int saveLocalMatrixToFile(const int narows, const int nacols, double* matrix, const std::string& filename)
{
    std::ofstream outfile(filename);
    if (!outfile.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return 1;
    }
    for (int i=0; i<narows; ++i) // i is the row index, j is the column index
    {
        for (int j=0; j<nacols; ++j)
        {
            outfile << matrix[i+j*narows] << " "; // attention: this is a column major matrix
        }
        outfile << std::endl;
    }
    outfile.close();
    return 0;
}

/**
 * @brief Saves a Block Cyclic Distributed (BCD) matrix to a file.
 * 
 * This function gathers a Block Cyclic Distributed (BCD) matrix from all processes and 
 * writes it to a specified file on the root process. The matrix is written in a row-major 
 * format, with elements separated by spaces and rows on separate lines.
 * 
 * @param comm The MPI communicator used for communication between processes.
 * @param desc The descriptor array for the BCD matrix, containing information about the matrix layout.
 * @param nrow The number of rows in the local matrix.
 * @param ncol The number of columns in the local matrix.
 * @param matrix A pointer to the local matrix data.
 * @param filename The name of the file to which the matrix will be saved.
 * @return Returns 0 if the file is successfully opened and written to; otherwise, returns 1.
 */
int saveBCDMatrixToFile(const MPI_Comm comm, const int* desc, const int nrow, const int ncol, const double* matrix, const std::string& filename)
{
    int return_val=0;
    // setup blacs environment
    int blacs_context=desc[1];
    const int nFull=desc[2];
    const int nblk=desc[4];
    int nprow, npcol, myprow, mypcol;
    Cblacs_gridinfo(blacs_context, &nprow, &npcol, &myprow, &mypcol);
    int myid;
    MPI_Comm_rank(comm, &myid);

    std::ofstream outfile;
    if(myid == 0)
    {
        outfile.open(filename);
        if (!outfile.is_open())
        {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return_val=1;
        }
    }
    MPI_Bcast(&return_val, 1, MPI_INT, 0, comm);
    if(return_val==1) return 1;

    double* a = const_cast<double*>(matrix);
    double* b; // buffer
    const int MAX_BUFFER_SIZE = 1e9; // max buffer size is 1GB

    int N = nFull;
    int M = std::max(1, std::min(nFull, (int)(MAX_BUFFER_SIZE / nFull / sizeof(double)))); // at lease 1 row, max size 1GB
    if (myid == 0)
        b = new double[M * N];
    else
        b = new double[1];

    int* desca=const_cast<int*>(desc);
    // set descb, which has all elements in the only block in the root process
    int descb[9] = {1, blacs_context, M, N, M, N, 0, 0, M};

    int ja = 1, ib = 1, jb = 1;
    for (int ia = 1; ia < nFull; ia += M)
    {
        int thisM = std::min(M, nFull - ia + 1); // nFull-ia+1 is the last few row to be saved
        // gather data rows by rows from all processes
        pdgemr2d_(&thisM, &N, a, &ia, &ja, desca, b, &ib, &jb, descb, &blacs_context);
        // write to the file
        if (myid == 0)
        {
            for (int i = 0; i < thisM; ++i)
            {
                for (int j = 0; j < N; ++j)
                {
                    outfile << b[i + j * M] << " ";
                }
                outfile << std::endl;
            }
        }
    }

    if (myid == 0)
        outfile.close();

    delete[] b;
    return return_val;
}

/**
 * @brief Loads a Block Cyclic Distributed (BCD) matrix from an ABACUS dumped file.
 * 
 * This function reads a matrix from a specified ABACUS - dumped file on the root process.
 * It then distributes the matrix in a Block Cyclic Distributed (BCD) format across all processes
 * using MPI and BLACS operations.
 * 
 * @param filename The name of the file from which the matrix will be loaded.
 * @param comm The MPI communicator used for communication between processes.
 * @param desc The descriptor array for the BCD matrix, containing information about the matrix layout.
 * @param matrix A pointer to the local matrix where the loaded data will be stored.
 * @return Returns 0 if the matrix is successfully loaded and distributed; otherwise, returns 1.
 */
int loadBCDMatrixFromABACUSFile(const std::string& filename, const MPI_Comm comm, const int* desc, double* matrix)
{
    int return_val=0;
    int blacs_context=desc[1];
    int N=desc[2];
    const int nblk=desc[4];
    int nprows, npcols, myprow, mypcol;
    Cblacs_gridinfo(blacs_context, &nprows, &npcols, &myprow, &mypcol);

    // read the matrix to b and send to all processes
    double* b=new double[N*N];

    int myid;
    MPI_Comm_rank(comm, &myid);

    std::ifstream matrixFile;
    if(myid == 0)
    {
        matrixFile.open(filename);
        if (!matrixFile.is_open())
        {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return_val=1;
        }
        // read the matrix to b
        int tmp;
        matrixFile>>tmp;
        std::cout<<"nFull="<<tmp<<"\n";
        double val;
        for(int i=0; i<N; ++i)
        {
            for(int j=i; j<N; ++j)
            {
                matrixFile>>val;
                b[i+j*N]=val;
                b[i*N+j]=val;
            }
        }
        matrixFile.close();
    }

    // set descb, which has all elements in the only block in the root process
    int descb[9] = {1, blacs_context, N, N, N, N, 0, 0, N};

    int isrc = 1;
    
    pdgemr2d_(&N, &N, b, &isrc, &isrc, descb, matrix, &isrc, &isrc, const_cast<int*>(desc), &blacs_context);

    delete[] b;
    return return_val;
}
