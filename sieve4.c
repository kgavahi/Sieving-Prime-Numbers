/*
 *   Sieve of Eratosthenes
 *
 *   Programmed by Michael J. Quinn
 *
 *   Modified by Kayhan Gavahi
 *
 *   Last modification: 24 March 2019
 */

#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN(a,b)  ((a)<(b)?(a):(b))

int main (int argc, char *argv[])
{
   unsigned long long   count;        /* Local prime count */
   double elapsed_time; /* Parallel execution time */
   unsigned long long    first;        /* Index of first multiple */
   unsigned long long    global_count; /* Global prime count */
   unsigned long long    high_value;   /* Highest value on this proc */
   unsigned long long    i;
   int    id;           /* Process ID number */
   unsigned long long    index;        /* Index of current prime */
   unsigned long long    low_value;    /* Lowest value on this proc */
   unsigned long long    low_value1;    /* Lowest value on this proc */
   char  *marked;       /* Portion of 2,...,'n' */
   unsigned long long    n;            /* Sieving from 2, ..., 'n' */
   int    p;            /* Number of processes */
   unsigned long long    proc0_size;   /* Size of proc 0's subarray */
   unsigned long long    prime;        /* Current prime */
   unsigned long long    size;         /* Elements in 'marked' */
   
   unsigned long long    size2;    /* size for local prime*/
   char                  *local_prime;          /* local_primes */
   unsigned long long    first2;     /* */
   
   unsigned long long    size3;     /* cache size */
   unsigned long long    end_block;     /*index of last   */
    

   MPI_Init (&argc, &argv);

   /* Start the timer */

      

   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   MPI_Comm_size (MPI_COMM_WORLD, &p);
   MPI_Barrier(MPI_COMM_WORLD);
   elapsed_time = -MPI_Wtime();



   n = atoll(argv[1]);


   /* Figure out this process's share of the array, as
      well as the integers represented by the first and
      last array elements */


   low_value = 3 +  2 * (id*((n/2)-1)/p);
   low_value1 = 3 +  2 * (0*((n/2)-1)/p);
   high_value = 1 + 2 * ((id+1)*((n/2)-1)/p);
   size = ((high_value - low_value)/2) + 1;

   size2 = (long long) sqrt((double) n);



   /* Allocate this process's share of the array. */

   marked = (char *) malloc (size);
	
   local_prime = (char *) malloc (size2);


   for (i = 0; i < size; i++) {
   marked[i] = 0;}

   for (i = 0; i < size2; i++) {
   local_prime[i] = 0; }  
   
   size3=10000;
   end_block=0;
   
do {
   index = 0;
   prime = 3;
   do {
      if (prime * prime > low_value)
         first = (prime * prime - low_value)/2;
      else {
         if (!(low_value % prime)) first = 0;
         else  {
            first = prime - (low_value % prime);
            if ((low_value+first)%2 == 0) {
               first = first+ prime;
            }
            first =first/2 ;
         }
      }

      first2 = (prime*prime -low_value1)/2; 
      

      for (i = first+end_block; i < MIN(end_block + size3,size); i += prime)
		  marked[i] = 1;

      for (i = first2; i<size2; i += prime) 
		  local_prime[i] = 1;  


      while (local_prime[++index]);
         prime = 2 * index + 3;

   } while (prime * prime <= n);
   	
	low_value += size3 * 2;
        end_block += size3; 
   } while (end_block<size);

   count = 0;
   global_count = 0;
   for (i = 0; i < size; i++)
      if (!marked[i]) count++;

   MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM,
                          0, MPI_COMM_WORLD);

   /* Stop the timer */

   elapsed_time += MPI_Wtime();


   /* Print the results */

   if (!id) {
      global_count; 
      printf ("There are %d primes less than or equal to %lld\n",
              global_count+1, n);
      printf ("SIEVE (%d) %10.6f\n", p, elapsed_time);
   }
   MPI_Finalize ();
   return 0;
}
