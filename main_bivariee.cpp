#include <cstdlib>
#include <cmath>
#include <iostream>
#include <cassert>
#include <omp.h>
#include <ctime>
#include <cstdlib>

#include "/usr/lib/x86_64-linux-gnu/openmpi/include/mpi.h"

#include "include/MpiBivariee.hpp"


int m = 1e+3;
const int n = 2; // n est fixé à 2
#define num_threads 8
//#define OMP


#ifdef OMP

double * Init(double * __restrict__ x,int size=1, int seed=0)
{
  //unsigned int Seed = 0;//time(NULL) ^ omp_get_thread_num() ^ getpt();
  // int offset = end - begin +1;
  srand(seed);
  
  #pragma parallel for
  for(int i = 0; i < size; i++) x[i] = (double) rand_r(&Seed) /((double) RAND_MAX);
  return x;
}

#endif


double * Init_X(double * __restrict__ X, int begin=0, int end=1, int n=1, int seed=0)
{
  srand(seed);
  int offset = end - begin +1;
  for(int i = 0; i < offset*n; i++)
  {
      X[i] =  (double) rand() /((double) RAND_MAX);
      //X[i*n+1] =  (double) rand() /((double) RAND_MAX);
  }
  return X;
}


double * Init(double * __restrict__ x,int begin=0, int end=1, int seed=0)
{
  srand(seed);
  int offset = end - begin +1;
  for(int i = 0; i < offset; i++) x[i] = (double) rand() /((double) RAND_MAX);
  return x;
}

void Print(double * X, int debut, int fin, int n,int rank)
{
  int offset = fin -debut +1;
  for(int i = 0; i < offset; i++) std::cout<<"rang ="<<rank<<"  i="<<i<<" xi="<<X[n*i]<<" xi+1="<<X[n*i+1]<<std::endl;
}

int main(int argc, char * argv[])
{
  MPI_Init(&argc , &argv);
  
  int rank = 0;
  int nproc = 0;

  omp_set_num_threads(num_threads);

  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int offset = m / nproc;

  // std::cout<<"nombre de proc="<<nproc<<std::endl;
   int begin = rank * offset;
  
  assert(m % nproc == 0); // Securité
  if (rank == nproc -1) offset += m % nproc;

   // int end = (rank +1) * offset - 1;
  int end = begin + offset;

  
  double * X = nullptr;
  size_t bytes_X = offset * n * sizeof(double);         //ceil((double)offset / ((double) 64))* 64* n* sizeof(double);
  int test1 = posix_memalign((void**) &X, 64, bytes_X);
  assert(!test1);

  double * Y=nullptr;
  size_t bytes_Y = ceil((double)offset / ((double) 64))* 64* sizeof(double);
  int test2 = posix_memalign((void**) &Y, 64, bytes_Y);
  assert(!test2);
  

  X = Init_X(X,begin,end,n,0);
  Y = Init(Y,begin,end,1);

//-------------- Definition de la classe-----------------------
  RegressionBivariee *  R = new RegressionBivariee(n); // Constructeur
  R->MaxIteration = 1000;
  R->m_loc = offset;
  R->begin = begin;
  R->end = end;
  R->offset = offset;
  R->m = m;
  R->n = n;
  R->rank = rank;
  R->nproc = nproc;

  // On Calcul la moyenne et l'erreur
  double moy = 0.0;
  double err = 0.0;
  double err_glob=0.;
  double moy_glob=0.;

  for(int i = 0; i < offset; i++)
    {
      moy += Y[i];
    }
  moy *= 1./((double) m);
  moy_glob =moy;

  if(nproc != 1)
    MPI_Allreduce(&moy, &moy_glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  R->moy = moy_glob;

  for(int i = 0; i < offset; i++)
    {
      double Resultat = Y[i]-moy;
      err += Resultat*Resultat;
    }
  err_glob = err;
  
  if(nproc != 1)
    MPI_Allreduce(&err, &err_glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   
  R->err2 =  err_glob;   
  R->eps = 1e-6;
  R->score = 0.5;
  R->alpha = 0.001;
  R->theta = Init(R->theta,0,n+1,1);
  R->theta0 = Init(R->theta0,0,n+1,1);  

  // Partie d'appel
  std::string test ="# iter J(theta) score theta0 theta1 theta2";
  std::cout<<"test1 ="<<test.length()<<std::endl;
  test = "+4.392934e-02 +3.896958e-01 +3.877523e-01 -5.805249e-02 +3.306638e-01";
  std::cout<<"test1 ="<<test.length()<<std::endl;
  
  R->fit(X,Y,R->alpha,R->m,R->n);
  

  free(Y);  free(X);  
  
  MPI_Finalize();
  return 0;
}
