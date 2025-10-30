#pragma once

#include <mpi.h>
#include <iostream>
#include <string>

class MPITimer {
private:
    double start_time;
    bool is_running;

public:
    MPITimer() : is_running(false) {}

    void start() {
        start_time = MPI_Wtime();
        is_running = true;
    }

    double stop() {
        if (!is_running) {
            std::cerr << "Timer was not started!" << std::endl;
            return 0.0;
        }
        is_running = false;
        return MPI_Wtime() - start_time;
    }

    void print(const std::string& msg) {
        double elapsed = stop();
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        std::cout << "[MPI " << rank << "] " << msg << ": " << elapsed << "s" << std::endl;
    }
};
