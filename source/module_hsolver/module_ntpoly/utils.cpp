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


void initBlacsGrid(MPI_Comm comm, int nFull, int nblk,
                   int& blacs_ctxt, int& narows, int& nacols, int* desc)
{
    char BLACS_LAYOUT='C';
    int ISRCPROC=0; // fortran array starts from 1
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
}

/**
 * Saves the parameters to a file.
 * 
 * @param filename The name of the file to which the parameters will be saved.
 * @param nFull The total number of rows or columns in the matrix.
 * @param nelec The number of electrons.
 * @param nspin The number of spins.
 * @param converge_density The density convergence threshold.
 * @param converge_overlap The overlap convergence threshold.
 * @param threshold The threshold for other calculations.
 * @return Returns 0 if successful, or an error code if an error occurs.
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
 * Saves a TripletList to a file.
 * 
 * @param tripletList The TripletList to be saved.
 * @param filename The name of the file to which the TripletList will be saved.
 * @return Returns 0 if successful, or an error code if an error occurs.
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
 * Saves a local matrix (column major) to a file.
 * 
 * @param N The dimension of the matrix.
 * @param matrix The local matrix to be saved.
 * @param filename The name of the file to which the matrix will be saved.
 * @return Returns 0 if successful, or an error code if an error occurs.
 */
int saveLocalMatrixToFile(const int N, double* matrix, const std::string& filename)
{
    std::ofstream outfile(filename);
    if (!outfile.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return 1;
    }
    for (int i=0; i<N; ++i) // i is the row index, j is the column index
    {
        for (int j=0; j<N; ++j)
        {
            outfile << matrix[i+j*N] << " "; // attention: this is a column major matrix
        }
        outfile << std::endl;
    }
    outfile.close();
    return 0;
}

/**
 * Saves a Block Cyclic Distributed (BCD) matrix to a file.
 * 
 * @param comm The MPI communicator.
 * @param desc The descriptor array for the BCD matrix.
 * @param nrow The number of rows in the local matrix.
 * @param ncol The number of columns in the local matrix.
 * @param matrix The local matrix to be saved.
 * @param filename The name of the file to which the matrix will be saved.
 * @return Returns 0 if successful, or an error code if an error occurs.
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
 * Loads a Block Cyclic Distributed (BCD) matrix from an ABACUS dumped file.
 * 
 * @param filename The name of the file from which the matrix will be loaded.
 * @param comm The MPI communicator.
 * @param desc The descriptor array for the BCD matrix.
 * @param matrix The local matrix to which the loaded data will be stored.
 * @return Returns 0 if successful, or an error code if an error occurs.
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
