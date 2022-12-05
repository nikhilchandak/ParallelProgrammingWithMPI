#include <iostream>
#include <mpi.h>
#include <bits/stdc++.h>

using namespace std;
#define vi vector < int > 

void process(vi& first, vi& second, int* mat, int* ans, int* head, int* tail, int& n, int& rnk) {
    
    for(int k = 0; k < (int) first.size(); k++) {
        int i = first[k], j = second[k] ; 
        if(i == j) continue ; 
        // fetch nodes from the two edges i and j
        int u = head[i], v = tail[i], x = head[j], y = tail[j] ; 

        // re-arrange if necessary 
        if(u > v) swap(u, v);
        if(x > y) swap(x, y);

        if(u == x and v == y) cout << "encountered same edge twice!" << endl ; // both edges are same!?
        int sum = -1 ; // NULL initially (stores the sum of the weights of the edges)

        // both endpoints different (ie., size 4 clique)
        if(u != x and v != y and u != y and x != v) {
            if(mat[u*n +v] >= 0 and mat[u*n +x] >= 0 and mat[u*n +y] >= 0 and mat[v*n +x] >= 0 and mat[v*n +y] >= 0 and mat[x*n +y] >= 0) {
                sum = mat[u*n +v] + mat[u*n +x] + mat[u*n +y] + mat[v*n +x] + mat[v*n +y] + mat[x*n +y] ; 
                ans[4 + sum]++ ; 
            }
        }
        else {
            // clique of size 3
            if(u == x or v == x) {
                if(mat[u*n +v] >= 0 and mat[u*n +y] >= 0 and mat[v*n +y] >= 0) {
                    sum = mat[u*n +v] + mat[u*n +y] + mat[v*n +y] ;
                    ans[sum]++ ;
                }
            }
            else if(v == y or u == y) {
                if(mat[v*n +u] >= 0 and mat[v*n +x] >= 0 and mat[u*n +x] >= 0) {
                    sum = mat[u*n +x] + mat[v*n +u] + mat[v*n +x] ;
                    ans[sum]++ ;
                }
            }
            else // both edges same?
                cout << "encountered repeated edge!" << endl ;
        }
    }
}

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

    int n, e ; int len, div, total ; 
    // stores the adjacency matrix but expanded in 1D array with row-major format (instead of 2D)
    int mat[10101] ; memset(mat, -1, sizeof(mat));
    
    int cnt[11] ; memset(cnt, 0, sizeof(cnt)); // store count of the different (size/weight) cliques
    int head[1001], tail[1001] ; // stores the head and tail of the edges 

    vi first, second, all1, all2 ;

    if(rank == 0) {
        cin >> n >> e; int u, v, w ; 
        for(int i = 0; i < e; i++) {
            cin >> u >> v >> w ; 
            u-- ; v-- ;
            mat[u*n + v] = w ; mat[v*n + u] = w ;
            head[i] = u ; tail[i] = v ; 
        }

        // create edge pairs 
        for(int i = 0; i < e; i++) {
            for(int j = i+1; j < e; j++) {
                all1.push_back(i);
                all2.push_back(j); 
            }
        }
        len = all1.size() ;
        div = ceil( float(len) / numprocs ) ;
        total = div * numprocs ;
        
        // add padding 
        for(int i = len, j = 0; i < total; i++, j++) {
            all1.push_back(0); all2.push_back(0) ;
        }
        // cout << len << " " << div << " " << total << " -- " << numprocs << endl ;
    }

    MPI_Bcast(&e, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&div, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    first.resize(div); second.resize(div); 

    MPI_Bcast(mat, 10101, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(head, 1001, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(tail, 1001, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(all1.data(), div, MPI_INT, first.data(), div, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(all2.data(), div, MPI_INT, second.data(), div, MPI_INT, 0, MPI_COMM_WORLD);

    process(first, second, mat, cnt, head, tail, n, rank);

    int ans[11] = {0};
    MPI_Reduce(cnt, ans, 11, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 

    // output 
    if (rank == 0){
        ofstream out;
        out.open(argv[2]);
        for(int i = 0; i < 11; i++){
            if (i < 4){
            out << 3 << " " << i << " " << ans[i] / 3 << endl;
            } else {
            out << 4 << " " << i - 4 << " " << ans[i] / 3 << endl;
            } 
        }
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