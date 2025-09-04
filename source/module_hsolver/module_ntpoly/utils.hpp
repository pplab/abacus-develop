#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
#include <fstream>
#include <mpi.h>
#include <TripletList.h>
#include <Triplet.h>
extern "C"
{
    #include "scalapack.h" 
}

// Function to output log information
static inline void outlog(const std::string& str)
{
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    std::string logfilename="ntpoly_"+std::to_string(myid)+".log";
    std::ofstream outfile(logfilename, std::ios::app);
    if (!outfile.is_open())
    {
        std::cerr << "Error: Could not open file outlog.txt" << std::endl;
        return;
    }
    outfile << str << std::endl;
    outfile.close();
}

// Function to output log information with an number
template <typename T>
static inline void outlog(const std::string& str, const T num)
{
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    std::string logfilename="ntpoly_"+std::to_string(myid)+".log";
    std::ofstream outfile(logfilename, std::ios::app);
    if (!outfile.is_open())
    {
        std::cerr << "Error: Could not open file outlog.txt" << std::endl;
        return;
    }
    outfile << str <<" : "<<num<< std::endl;
    outfile.close();
}

static inline void saveArrayToFile(const std::string& file_prefix, const double* array, int size)
{
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    std::string filename = file_prefix + "_" + std::to_string(myid) + ".dat";
    std::ofstream outfile(filename);
    if (!outfile.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    for (int i = 0; i < size; ++i)
    {
        outfile << array[i] << std::endl;
    }
    outfile.close();    
}

static inline void saveArrayToFile(const std::string& file_prefix, const int* array, int size)
{
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    std::string filename = file_prefix + "_" + std::to_string(myid) + ".dat";
    std::ofstream outfile(filename);
    if (!outfile.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
    }
    for (int i = 0; i < size; ++i)
    {
        outfile << array[i] << std::endl;
    }
    outfile.close();
}

static inline void saveMatrixToFile(const std::string& file_prefix, const double* matrix, int rows, int cols)
{
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    std::string filename = file_prefix + "_" + std::to_string(myid) + ".dat";
    std::ofstream outfile(filename);
    if (!outfile.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            outfile << matrix[i * cols + j] << " ";
        }
    }
    outfile.close();
}

// Function to initialize the BLACS grid
void initBlacsGrid(MPI_Comm comm, int nFull, int nblk,
                   int& blacs_ctxt, int& narows, int& nacols, int* desc);

// Function to save parameters to a file
int saveParametersToFile(const std::string& filename, 
        const int nFull, const int nelec, const int nspin, 
        const double converge_density, const double converge_overlap, const double threshold);

// Function to load parameters from a file
int loadParametersFromFile(const std::string& filename,
        int& nFull, int& nelec, int& nspin, 
        double& converge_density, double& converge_overlap, double& threshold);

// Function to save a TripletList to a file
int saveTripletListToFile(const NTPoly::TripletList_r& tripletList, const std::string& filename);

// Function to save a local matrix to a file
int saveLocalMatrixToFile(const int narows, const int nacols, double* matrix, const std::string& filename);

// Function to save a Block Cyclic Distributed (BCD) matrix to a file
int saveBCDMatrixToFile(const MPI_Comm comm, const int* desc, const int nrow, const int ncol, const double* matrix, const std::string& filename);

// Loads a Block Cyclic Distributed (BCD) matrix from an ABACUS dumped file
int loadBCDMatrixFromABACUSFile(const std::string& filename, const MPI_Comm comm, const int* desc, double* matrix);
#endif // UTILS_HPP
