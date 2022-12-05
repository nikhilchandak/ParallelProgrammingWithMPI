#include <iostream>
#include <mpi.h>
#include <bits/stdc++.h>

using namespace std;

#define vi vector < int > 


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
    // freopen(argv[2], "w", stdout);
    int n ;
    double* mat ; 

    if(rank == 0) {
        cin >> n ;
        mat = new double[10101];
        for(int i = 0; i < n; i++)
            for(int j= 0; j <= n; j++)
                cin >> mat[i * (n + 1) + j] ;
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    int div = ceil( float(n) / numprocs ) ;
    int row_size = n + 1 ; // n for A and 1 for b (Ax = b)
    double sub[div * row_size] ; memset(sub, 0, sizeof(sub)); // sub-matrix of respective process 
    int perm[n] ; for(int i = 0; i < n; i++) perm[i] = i ;

    // MPI_Barrier(MPI_COMM_WORLD);
    for(int i = 0; i < div; i++) // demossable
        MPI_Scatter(mat + i * row_size * numprocs, row_size, MPI_DOUBLE, sub + i * row_size, row_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    vector < double > row(row_size + 1, 0); 


    for(int i = 0; i < n; i++) {
        int local_row = i/numprocs, row_rnk = i % numprocs ; 

        // my row 
        if(row_rnk == rank) {
            
            // find pivoting column
            int piv = i ; double mx = -1 ;
            for(int j = i; j < n; j++) {
                double val = fabs(sub[local_row * row_size + j]);
                if(val > mx) {
                    mx = val ;
                    piv = j ; 
                }
            }
            // update permutation of columns 
            swap(perm[i], perm[piv]);

            // swap the whole column 
            for(int row = 0; row < div; row++)
                swap(sub[row * row_size + i], sub[row * row_size + piv]);
            
            // normalize the given row
            double pivot = sub[local_row * row_size + i];
            sub[local_row * row_size + i] = 1;
            for (int k = i + 1; k < row_size; k++) 
                sub[local_row * row_size + k] /= pivot;

            // prepare row 
            for(int j = 0; j < row_size; j++) 
                row[j] = sub[local_row * row_size + j] ;
                
            row[row_size] = piv ; 

            // if current rank has the given row broadcast it other processes
            MPI_Bcast(row.data(), row_size + 1, MPI_DOUBLE, row_rnk, MPI_COMM_WORLD);

            // use this row to reduce other rows in your own sub
            for (int row = local_row + 1; row < div; row++) {
                double scale = sub[row * row_size + i];
                for (int k = i; k < row_size; k++) 
                    sub[row * row_size + k] -= scale * sub[local_row * row_size + k];
            }
        }
        else {
            MPI_Bcast(row.data(), row_size + 1, MPI_DOUBLE, row_rnk, MPI_COMM_WORLD);
            int piv = (int) row[row_size] ;

            // update permutation of columns 
            swap(perm[i], perm[piv]);
            
            // swap the whole column 
            for(int row = 0; row < div; row++)
                swap(sub[row * row_size + i], sub[row * row_size + piv]);

            for (int k = local_row; k < div; k++) {
                // this row is above the row received
                if( row_rnk >= rank and k == local_row)
                    continue;

                // subtract the row
                double scale = sub[k * row_size + i];
                for (int m = i; m < row_size; m++) 
                    sub[k * row_size + m] -= scale * row[m];
            }
        }
    }
    // save the original matrix
    double *mat2;
    if (rank == 0)
        mat2 = new double[div * numprocs * row_size];

    for(int i = 0; i < div; i++)
        MPI_Gather(sub + row_size * i, row_size, MPI_DOUBLE, mat2 + numprocs * row_size * i, row_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // back substitution
        vector<double> ans(n, 0), ans2(n, 0);

        for (int i = n - 1; i >= 0; i--) {
            ans[i] = mat2[i * row_size + n];
            for (int j = i + 1; j < n; j++) {
                ans[i] -= mat2[i * row_size + j] * ans[j];
            }
            ans[i] = ans[i] / mat2[i * row_size + i];
        }

        // map permutation 
        for(int p = 0; p < n; p++) ans2[perm[p]] = ans[p] ; 
        swap(ans, ans2);

        ofstream out;
        out.open(argv[2]);
        out << std::fixed << setprecision(7);
        for (double x : ans) 
            out << x << " ";

        // verify the answer is correct
        // bool flag = true;
        // for (int i = 0; i < n; i++) {
        // double sum = 0;
        // for (int j = 0; j < n; j++) {
        //     sum += mat[i * row_size + j] * ans[j];
        // }
        // // cout << sum << " " << mat[i * row_size + n] << endl;
        // if (fabs(sum - mat[i * row_size + n]) > 0.00001) {
        //     cout << "Wrong answer" << endl;
        //     flag = false;
        // }
        // }
        // if (flag) {
        // cout << "Correct answer" << endl;
        // }
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