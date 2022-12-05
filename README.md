### Parallel Programming using MPI

## Q1 -- Prime Check

We find the square root as only numbers till there need to be checked. We divide sqrt by numprocs to mark the 
low/high limit for each process to check the numbers in its chunk whether they divide the given number.

Complexity : O(sqrt(N)/P)

## Q2 -- Clique Computation

For a given pair of edges <e1, e2> where e1 != e2 : either no node is common in
e1 and e2 which makes all the 4 nodes in total distinct thus we check whether they form a clique of size 4. In case a node
is common between e1 and e2 (can only be at most 1), then we have 3 distinct nodes overall thus acting as a candidate
for clique of size 3. 

Now for all such pairs of edges, we proportinately divide this to the processes (as in P1) and gather all their counts in the root.
Complexity : O(E^2/P)


## Q3 -- Gaussian Elimination

A root process takes input and distributes rows in a cyclic order among all the processes (first numproc rows where
each row is given to one process, like this further on with next rows).
Firstly we will make the A matrix upper triangular (with all diagonal elements = 1), and then 
back subsitute x values one by one to get all values.

To obtain upper trianlge formation :

- Each process upon receival of a row, sends it to further processes (if required) and updates all rows in its own quota by swapping columns along the received pivot index and making values below diagonal element of received row to 0.
- When a process gets all rows above its current row after their elimination, the current row is eliminated (pivot swapped and diagonal element=1) and sends the row, and pivot index to next process. All rows in current process' quota are also updated just like they are when a row is recieved from message.
- Each rows gets eliminated only when all other rows above it are eliminated.

Complexity : O(N^3/P)
