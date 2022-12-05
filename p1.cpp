#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <vector>

using namespace std;

int main( int argc, char **argv ) {
    int rank, numprocs;

    // initiate MPI
    MPI_Init( &argc, &argv );

    // get size of the current communicator
    MPI_Comm_size( MPI_COMM_WORLD, &numprocs );

    // get current process rank
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    /*synchronize all processes*/
    MPI_Barrier( MPI_COMM_WORLD );
    double start_time = MPI_Wtime();

    // enter your code here
    freopen(argv[1], "r", stdin);
    
    int n, div ;

    if(rank == 0) {
        cin >> n ;
        int sq = (int) ceil(sqrt(n)) ; 

        div = (int) ceil(sqrt(n)/numprocs) ;
        div = max(div, 1);
    }
    MPI_Bcast(&div, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    int flag = 0 ; int others = 0 ; 
    int low = max(2, div * rank); int high = min(n-1, div * (rank + 1)); 
    
    // check if any number in this range [low, high] is a factor
    for(int i = low; i <= high; i++) {
        if(n % i == 0) {
            flag = 1;
            break ;
        }
    }
    // take max from all processes (ie, not prime even if one process finds a factor)
    MPI_Reduce(&flag, &others, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
 
    if(rank == 0)
    {
        ofstream out;
        out.open(argv[2]);
        flag = max(others, flag);
        if(flag == 0) out << "YES" << endl ;
        else out << "NO" << endl ; 
        out.close();
    }

    MPI_Barrier( MPI_COMM_WORLD );
    double end_time = MPI_Wtime() - start_time;
    double maxTime;
    // get max program run time for all processes
    MPI_Reduce( &end_time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    if ( rank == 0 ) {
        cout<<"Total time (s): "<<maxTime<<"\n";
    }

    // shut down MPI and close
    MPI_Finalize();
    return 0;
}