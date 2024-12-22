// utility tools for debug
#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <TripletList.h>
#include <Triplet.h>
#include "Cblacs.h"
extern "C"
{
    #include "scalapack.h"
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
 * Saves a local matrix to a file.
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
    int retern_val=0;
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
            retern_val=1;
        }
    }
    MPI_Bcast(&retern_val, 1, MPI_INT, 0, comm);
    if(retern_val==1) return 1;

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
    return retern_val;
}

