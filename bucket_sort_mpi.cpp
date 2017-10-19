#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

using namespace std;

const float xmin = 10.0;
const float xmax = 250000;
const float padding = 0.0;  
    
int numproc, myid;

/*****************************************************************************/

// Function to check whether the first n items of an array has been sorted.
int check(float *data,int nitems) {
  int sorted=1;
  int i;
  for(i=0;i<nitems;i++) {
     if(i && data[i]<data[i-1])
       {
         sorted=0;
         break;
       }
  }
  printf("totally %d items, sorted=%d\n",nitems, sorted);
}


// Create a data array to hold the given number of buckets for the
// given number of total data items. All buckets are held contiguously
float* create_buckets(int nbuckets, int nitems)
{
  int i;

  int ntotal = nbuckets * nitems;
  float* bucket = new float[(ntotal * sizeof(float))];
  for (i=0; i<ntotal; ++i)
      bucket[i] = 0;

  // return the address of the array of pointers to float arrays
  return bucket;
}



// Each process will have n/p numbers. Each value is examined and put into
// the appropriate small bucket.
void bucket_distribute(float *data, int ndata, float x1, float x2, int nbuckets,
                       float *bucket) 
{
  int i, count;

  // The range covered by one bucket
  float stepsize = (x2 - x1) / nbuckets;

  // The number of items thrown into each bucket. We would expect each
  // bucket to have a similar number of items, but they won't be
  // exactly the same. So we keep track of their numbers here.
  int* nitems = new int[nbuckets];
  for (i = 0; i < nbuckets; ++i)
    nitems[i] = 0;
  // Toss the data items into the correct bucket
  for (i = 0; i < ndata; ++i) {
    if(data[i] == padding)  // If it is a padding, all following data are paddings.
        break;
    // What bucket does this data value belong to?
    int bktno = (int)floor((data[i] - x1) / stepsize);
    int idx = bktno * ndata + nitems[bktno];
    // Put the data value into this bucket
    bucket[idx] = data[i];
    ++nitems[bktno];
  }
  delete nitems;
}



// The comparison function to use with the library qsort
// function. This will tell qsort to sort the numbers in ascending
// order.
int compare(const void* x1, const void* x2) {
  const float* f1 = (float *)x1;
  const float* f2 = (float *)x2;
  float diff = *f1 - *f2;

  return (diff < 0) ? -1 : 1;
}

// Main Entry
int main(int argc, char *argv[]) 
{
  /************  Phase 1 - Partitioning  **********/
  float *recvdata = NULL;
  float *data = NULL;
  double start, end;
  if(myid == 0)
    start = MPI_Wtime();
  //int numproc, myid;

  MPI_Init(&argc,&argv);
  
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  int items_count= 100000000;
  if(argc==2)
      items_count=atoi(argv[1]);
  int nitems_bkt = (int)(items_count / numproc); // The items number received by one processor via MPI_Scatter
  // It is very possible that there will be an remimnder after division, e.g. items_count is 100, divided by 7, the reminder is 2.
  // To make sure all reminders could be sent respectively, I add some paddings to data if needed.
  int mod = items_count % numproc;
  int paddingCount = 0;
  if(mod)
  {
      nitems_bkt = nitems_bkt + 1;
     // nitems_bkt++;
      paddingCount = numproc - mod;
  }
  int senddata_len = items_count + paddingCount;
  data = new float[senddata_len];
  recvdata = new float[nitems_bkt];
  if(myid == 0)
  { 
      for(int i=0;i<items_count;i++)
        data[i]=drand48()*(xmax-xmin-1)+xmin;
      for(int i = items_count;i < senddata_len; i ++)  // Add paddings
        data[i]=padding;
      check(data,items_count);
      printf("Begin bucket sorting...\n");
  }
  // Send data to all processors vai scatter, the last processor may receive paddings.
  MPI_Scatter(data,nitems_bkt,MPI_FLOAT,recvdata,nitems_bkt,MPI_FLOAT,0,MPI_COMM_WORLD);
  
  /************  Phase 2 - Small bucket distribution  **********/
  int small_bucket_count = numproc; // One processor is corresponding to one bucket.
  float *small_buckets = create_buckets(small_bucket_count,nitems_bkt);
  bucket_distribute(recvdata,nitems_bkt,xmin,xmax,small_bucket_count,small_buckets);
  /************  Phase 3 - Put in large buckets  **********/
  int large_bucket_size = numproc * nitems_bkt;
  float *large_bucket = new float[large_bucket_size];
  for (int i = 0; i < large_bucket_size; ++i)
    large_bucket[i] = 0;
  MPI_Alltoall(small_buckets, nitems_bkt, MPI_FLOAT, large_bucket, nitems_bkt, MPI_FLOAT, MPI_COMM_WORLD);
  /************  Phase 4 - Sort large buckets  **********/
  /*
   * Find how many of them are not 0 and move non-zero values to the front of the array.
   * e.g.: for array [0,0,0,1,2,4,3,0,6], after below operation, it will be [1,2,4,3,6,0,0,0,0], 
   * and new_index will be 5, this is same with the count of non-zero values.
   */
  int new_index = 0; // new_index will be the count of non-zero values. 
  for(int i = 0; i < large_bucket_size; i ++)
  {
      
      if(large_bucket[i] != 0)
      {
          if(new_index != i)
          {
              large_bucket[new_index] = large_bucket[i];  // Save the non-zero value to new_index and set current position be 0.
              large_bucket[i] = 0;
          }
	new_index++;
      }
  }
  qsort(large_bucket, new_index, sizeof(float), compare);
  /************  Phase 5 - Gather each large bucket  **********/
  int *displs = new int[numproc];
  int *large_bucket_nitems = new int[numproc]; // Items count in each bucket.
  for (int i = 0; i < numproc; i ++)
  {
    large_bucket_nitems[i] = 0;
    displs[i] = 0;
  }
  // Get the items count in each bucket.
  MPI_Gather(&new_index, 1, MPI_INT, large_bucket_nitems, 1, MPI_INT, 0, MPI_COMM_WORLD);
  for (int i = 0; i < numproc; i ++)
  {
    for(int j = 0; j < i; j ++)
      displs[i] += large_bucket_nitems[j];
  }
  // Get all sorted items and save to original array.
  MPI_Gatherv(large_bucket, new_index, MPI_FLOAT, data, large_bucket_nitems, displs, MPI_FLOAT, 0, MPI_COMM_WORLD);
  if(myid == 0)
  {
    end = MPI_Wtime();
    double time_used = end - start;
    check(data,items_count);
    printf("Bucket sort finished, time used: %f seconds\n", time_used);
  }
  delete displs;
  delete large_bucket_nitems;
  delete recvdata;
  delete data;
  delete large_bucket;
  delete small_buckets;
   
  MPI_Finalize();
}
