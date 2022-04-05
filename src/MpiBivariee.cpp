#include <iostream>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <cassert>
#include "/usr/lib/x86_64-linux-gnu/openmpi/include/mpi.h"
#include <fstream>
#include <string>
#include <cstdlib> 

#include "../include/MpiBivariee.hpp"
#include <omp.h>
#include <immintrin.h>
#include <smmintrin.h>
#include <emmintrin.h>
#include <list>
#include <ctime>

#define OPTI_MPI_1
#define OPTI_MPI_2


// Script gnplot

void RegressionBivariee::Script_gnuplot3d(std::string & out_put_gnuplot_file, std::string & input_plot_file, std::string image_name)
{
  std::ofstream File_gnuplot;

  File_gnuplot.open(out_put_gnuplot_file);
  if(!File_gnuplot.is_open()) { std::cerr<<"Impossible d'ouvrir le fichier"<<std::endl; exit(EXIT_FAILURE);}

  File_gnuplot<<"set title \'Regression lineaire bivariee \' "<<std::endl;
  File_gnuplot<<"#set terminal eps"<<std::endl;
  File_gnuplot<<"set term postscript eps enhanced color"<<std::endl; // dirige les sorties vers un fichier postscript avec plusieurs options
  File_gnuplot<<"#set term x11"<<std::endl; // redirige les sorties vers une fenetre de l'ecran
  File_gnuplot<<"set output \""<<image_name<<"\""<<std::endl;
  File_gnuplot<<"set key right top "<<std::endl;
  File_gnuplot<<"set xlabel \'x1\' "<<std::endl;
  File_gnuplot<<"set ylabel \'x2\' "<<std::endl;
  File_gnuplot<<"set zlabel \'y\'  "<<std::endl;
  File_gnuplot<<"set view 10,30 "<<std::endl;
  File_gnuplot<<"#set view map "<<std::endl; // vue de dessus
  // surface
  File_gnuplot<<"set surface"<<std::endl; // trace la grille definissant la surface
  File_gnuplot<<"#set colorbox"<<std::endl; // Affiche la barre de couleur  

  // Plot
  std::string convert = "splot \" " + input_plot_file +"\" u 1:2:3";
   File_gnuplot<<convert<<std::endl;
  
}

// Split data train and data test 
void RegressionBivariee::Split_train_test(const double * __restrict__ X, const double * __restrict__ Y, double * __restrict__ X_train, double * __restrict__ X_test, double * __restrict__ Y_train, double * __restrict__ Y_test, double proba_train, int m_train, int m_test) const

{
  m_train = ceil(proba_train*m);
  m_test = m - m_train;

  // memcpy(X_train, X, 2*m_train);
  // memcpy(Y_train, Y, m_train);

  // memcpy(X_test, X+, 2*m_train);
  // memcpy(Y_train, Y, m_train);

  for(int i = 0; i <m_train; i++)
    {
      X_train[2*i] = X[2*i];
      X_train[2*i+1] = X[2*i+1];
      Y_train[i] = Y[i];
    }

  for(int i = 0; i < m_test; i++)
    {
      X_test[i*2] = X[2*i+m_train];
      X_test[i*2+1] = X[2*i+m_train+1];
      Y_test[i] = Y[i];
    }
  
}

// Init file

void RegressionBivariee::File_create(const std::string & file_name, int m,int seed)const
{
  std::ofstream File = std::ofstream();
  File.open(file_name);

  if(!File.is_open()) {std::cerr<<"Erreur d'ouverture de fichier"<<std::endl; exit(EXIT_FAILURE);}
  
}


// Output file for plot (after fit, the resultat will save with this member function)
void RegressionBivariee::DataFileBivariee(const double * __restrict__ X_train, const double * __restrict__ Y_train,int n_train, const double * __restrict__ X_test, const double * __restrict__ Y_test, int n_test, std::string & file_name) const
{
  // Ouverture du :fichier en écriture; Le format sera en colonne et il se peut que les dimmensions ne soient pas les mêmes

  assert(n_train >= n_test);

  std::ofstream File = std::ofstream();
  File.open(file_name);

  if(!File.is_open()) { std::cerr<<"Fichier non ouvert!!!"<<'\n'; exit(EXIT_FAILURE);}

  
  // if(rank == 0)
  //   {
  //     File<<"# Ce fichier contient un sous-ensemble test-train"<<'\n';
  //     File<<"# X_train, Y_train, X_test, Y_test "<<'\n';
  //     File <<'\n';
  //   }

  for(int i  = 0; i < n_train; i++)
    {
      File << X_train[i] <<'\t'<<Y_train[i];
      if(i < n_test) File<< X_test[i] <<'\t' << Y_test[i];

      File <<'\n';
    }
  File.close();
}



// Constructor
RegressionBivariee::RegressionBivariee(int n_ =2):n{n_}
{
  size_t bytes = 32*sizeof(double);
  
  int test = posix_memalign((void**)&theta,32,bytes);
  //if (test != 0) perror("Attention, problem de memoire");
  assert(!test);
  int test1 = posix_memalign((void**)&theta0,32,bytes);
  //if (test1 != 0) perror("Attention, problem de memoire");
  assert(!test1);
}

// Another Constructor
RegressionBivariee::RegressionBivariee(const std::string & file_name, double * __restrict__ X, double * __restrict__ Y,int m, int n, int size_of_one_line):n{2}
{
  // ouverture du fichier en mode lecture
  std::ifstream File = std::ifstream();
  File.open(file_name);

  if(!File.is_open()){ std::cerr<<"Impossible d'ouvrir le fichier "<<std::endl; exit(EXIT_FAILURE);}
  
  std::string x0,x1,y,read_line="";
  // On lit la taille, sur la première ligne du fichier
  File>>read_line;
  assert(RegressionBivariee::m == stoi(read_line));
  
  // while(++i < num_line) std::getline(File,read_line);

  // debut de lecture
  size_t beg_pos = begin * size_of_one_line + read_line.length();
  File.seekg(beg_pos,File.beg);

  for(int i = begin; i < end; i++)
    {
      File>>x0>>x1>>y;
      X[2*i] = std::stod(x0);
      X[2*i+1] = std::stod(x1);
      Y[i] = std::stod(y);
    }
  // assert(m == RegressionBivariee::m);
  // assert(n==RegressionBivariee::n);
}


// Calcule l'erreur de consistance : Somme_i (Ax - y)_i
double RegressionBivariee::ErreurDeConsistance(double *__restrict X, double *__restrict Y, double *__restrict theta, int m, int n)const
{
  double y = 0.;
  int offset = end - begin+1;
  int size = offset * n;

  int step = 4;
  int fin = size / step;
  int r = size % step;

  int warp_register = 4;
  int indice = 0;

  double ss = 0;
  assert(n == 2);
  assert(r==0);
  
  // Pour theta
  __m256d TH0 = _mm256_set1_pd(*theta);
  __m128d th = _mm_loadu_pd(theta+1);

  // On stock les deux valeurs dans un registre 256
  __m256d TH = _mm256_broadcast_pd(&th);
  
  //#pragma omp parallel for reduction(+:y)
  for(int i = 0; i < fin; i++)
    {
      indice = i*warp_register;
      __m256d X1 = _mm256_load_pd(X+indice);
      __m256d R1 = _mm256_mul_pd(X1,TH);
      __m256d R11 = _mm256_hadd_pd(R1,R1);
      __m256d R11_Theta0 = _mm256_add_pd(TH0,R11);
      double s[]={1,2,3,4};
      _mm256_storeu_pd(s, R11_Theta0);
      //std::cout<<"rank= "<< rank <<"somme, s="<<s[0]+s[2]<<std::endl;
       y += s[1] - Y[i*2]+ s[2]- Y[i*2+1];
    }

  return y;
}

double RegressionBivariee::JBivariee(double * __restrict X, double * __restrict Y, int m=1, int n=2)const
{
  double y = 0.;
  int offset = end - begin+1;
  //int size = offset * n;

  int step = 2;
  int fin = offset / step;
  int r = offset % step;

  const int warp_register = 4;
  int indice = 0;

  double ss = 0;
  assert(n == 2);
  assert(r==0);
  
  // Pour theta
  __m256d TH0 = _mm256_set1_pd(*RegressionBivariee::theta);
  __m128d th = _mm_loadu_pd(RegressionBivariee::theta+1);

  // On stock les deux valeurs dans un registre 256
  __m256d TH = _mm256_broadcast_pd(&th);
  
  //#pragma omp parallel for reduction(+:y) private(indice)
  //#pragma unrull(4)
  for(int i = 0; i < fin; i++)
    {
      indice = i*warp_register;
      __m256d X1 = _mm256_load_pd(X+indice);
      __m256d R1 = _mm256_mul_pd(X1,TH);
      __m256d R11 = _mm256_hadd_pd(R1,R1);
      __m256d R11_Theta0 = _mm256_add_pd(TH0,R11);
      double y1 = Y[i*step];double y2 = Y[i*step+1];
      __m256d Y1 = _mm256_setr_pd(y1,y1,y2,y2);
      __m256d R = _mm256_sub_pd(R11_Theta0,Y1);
      __m256d R_pow2 = _mm256_mul_pd(R,R);
      double s[]={1,2,3,4};
      _mm256_storeu_pd(s, R_pow2);
      //std::cout<<"rank= "<< rank <<"somme, s="<<s[0]+s[2]<<std::endl;
      y += s[1]+s[2];
    }

  return y* 0.5 * (1./((double)m));;
}


// Evalue le score (signé)

double RegressionBivariee::Score(double * __restrict__ X, double * __restrict__ Y, int m=1, int n=2)const
{
  double err1 = 0.;
  
  double err = 0.;
  int offset = end - begin + 1;

  int step = 2;
  int fin = offset / step;
  int r = offset % step;

  const int warp_register = 4;
  int indice = 0;

  double ss = 0;
  assert(n == 2);
  assert(r==0);
  
  // Pour theta
  __m256d TH0 = _mm256_set1_pd(*RegressionBivariee::theta);
  __m128d th = _mm_loadu_pd(RegressionBivariee::theta+1);

  // On stock les deux valeurs dans un registre 256
  __m256d TH = _mm256_broadcast_pd(&th);
  
  //#pragma omp parallel for reduction(+:y) private(indice)
  //#pragma unrull(4)
  for(int i = 0; i < fin; i++)
     {
       indice = i*warp_register;
       __m256d X1 = _mm256_load_pd(X+indice);
       __m256d R1 = _mm256_mul_pd(X1,TH);
       __m256d R11 = _mm256_hadd_pd(R1,R1);
       __m256d R11_Theta0 = _mm256_add_pd(TH0,R11);
       double y1 = Y[i*step];double y2 = Y[i*step+1];
       __m256d Y1 = _mm256_setr_pd(y1,y1,y2,y2);
      __m256d R = _mm256_sub_pd(R11_Theta0,Y1);
       __m256d R_pow2 = _mm256_mul_pd(R,R);
       double s[]={1,2,3,4};
       _mm256_storeu_pd(s, R_pow2);
       //std::cout<<"rank= "<< rank <<"somme, s="<<s[0]+s[2]<<std::endl;
       err1 += s[1]+s[2];
     }
  
  //  #pragma omp parallel for reduction(+:err1)
  // for(int i = 0; i < offset; i++)
  //   {
  //     err1 += pow(h(X+i*n,RegressionBivariee::n)-Y[i],2);
  //   }
  // ici 
  
  assert(RegressionBivariee::err2);
  err = err1/RegressionBivariee::err2;
  return err;
}

double RegressionBivariee::h(double * __restrict X_i, int n =2)const
{
  assert(n==2);
  return theta[0] + theta[1]*(*X_i) + theta[2]*(*(++X_i));
}

double * RegressionBivariee::H(const double * __restrict__ X, double * __restrict__ Y, int m, int n) const
{
  double y = 0.;
  int offset = end - begin+1;
  //int size = offset * n;

  int step = 2;
  int fin = offset / step;
  int r = offset % step;

  const int warp_register = 4;
  int indice = 0;

  double ss = 0;
  assert(n == 2);
  assert(r==0);
  
  // Pour theta
  __m256d TH0 = _mm256_set1_pd(*RegressionBivariee::theta);
  __m128d th = _mm_loadu_pd(RegressionBivariee::theta+1);

  // On stock les deux valeurs dans un registre 256
  __m256d TH = _mm256_broadcast_pd(&th);
  
  //#pragma omp parallel for reduction(+:y) private(indice)
  //#pragma unrull(4)
  for(int i = 0; i < fin; i++)
    {
      indice = i*warp_register;
      __m256d X1 = _mm256_load_pd(X+indice);
      __m256d R1 = _mm256_mul_pd(X1,TH);
      __m256d R11 = _mm256_hadd_pd(R1,R1);
      __m256d R11_Theta0 = _mm256_add_pd(TH0,R11);
      double s[]={1,2,3,4};
      _mm256_storeu_pd(s, R11_Theta0);
      //std::cout<<"rank= "<< rank <<"somme, s="<<s[0]+s[2]<<std::endl;
      Y[i] = s[1];
      Y[i*2+1]= s[2];
    }
  return Y;
}


void RegressionBivariee::fit(double * __restrict X, double* __restrict Y,double alpha =1.,int m=1, int n=2)
{
  assert(n == 2);

  int it = 0;
  double n_norme = 1.;
  double loc_n_norme = 1.;
  double last_n_norme = 0;
  int offset = end - begin +1;
  double glob_sc = 0.;
  std::ofstream file;


  loc_n_norme = JBivariee(X,Y,m,n);
  last_n_norme = loc_n_norme;
  
  if(nproc != 1)
    MPI_Allreduce(&loc_n_norme, &last_n_norme, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  double loc_sc = Score(X,Y,m,n);
  glob_sc = loc_sc;
   
  if(nproc != 1)
    MPI_Allreduce(&loc_sc, &glob_sc,1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  if(rank ==0){
     file.open("FichierDeConvergence.dat");

     if(!file.is_open()) { std::cerr<<"Fichier non ouvert!!!"<<'\n'; exit(EXIT_FAILURE);}
     
     file<<"# Fichierdeconvergence.dat"<<'\n';
    
     file<<"# iter"<<' '<<"J(theta)"<<' '<<"score"<<' '<<"theta0"<<' '<<"theta1"<<' '<<"theta2"<<' '<<'\n';
     file<<std::showpos<<std::scientific;

     file<<(double) it<<' '<<last_n_norme<<' '<<1-glob_sc<<' '<<theta0[0]<<' '<<theta0[1]<<' '<<theta0[2]<<'\n';
  }
  
   do
     {
       double s = 0.;
       double s_recv = 0.;

       s = ErreurDeConsistance(X,Y,theta0,m,n);
       s_recv = s;
       
       if(nproc != 1)
        MPI_Allreduce(&s,&s_recv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      
       // On copy la valeur dans theta0       
       for(int j = 0; j <= n; j++)
	 {
	   theta[j] = theta0[j] - alpha*(1./((double) m))*s_recv;
	   theta0[j] = theta[j];
	   //std::cout<<"i= "<<j<<" theta= "<<theta[j]<<std::endl;
	 }

       // On evalue la condition d'arrêt
       loc_n_norme = JBivariee(X,Y,m,n);
       n_norme = loc_n_norme;

       if(nproc != 1)
	 MPI_Allreduce(&loc_n_norme, &n_norme, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

       double tempo = n_norme;
       n_norme -= last_n_norme;
       last_n_norme = tempo;       
       if(nproc != 1)
        MPI_Allreduce(&s,&s_recv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      
       // On copy la valeur dans theta0       
       for(int j = 0; j <= n; j++)
	 {
	   theta[j] = theta0[j] - alpha*(1./((double) m))*s_recv;
	   theta0[j] = theta[j];
	   //std::cout<<"i= "<<j<<" theta= "<<theta[j]<<std::endl;
	 }

       // On evalue la condition d'arrêt

       loc_sc = Score(X,Y,m,n);
       glob_sc = loc_sc;
   
	if(nproc != 1)
	  MPI_Allreduce(&loc_sc, &glob_sc,1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      
       // std::cout<<" rank ="<< rank <<" norme= "<<n_norme<<std::endl;
       // std::cout<<" rank ="<<rank<<" score= "<<1-glob_sc<<std::endl;
	++it;
	if(rank == 0) { file<<(double) it<<' '<<last_n_norme<<' '<<1-glob_sc<<' '<<theta[0]<<' '<<theta[1]<<' '<<theta[2]<<'\n';}
       
       if(n_norme < 0.) n_norme = -1.* n_norme; 
       
     }while(( it < RegressionBivariee::MaxIteration) and (eps < n_norme) );

   file.close();
   
  RegressionBivariee::Iteration = it;
  if(rank==0)
    {
        std::cout<<" rank ="<<rank<<" iter= "<<Iteration<<std::endl;
	std::cout<<" rank ="<< rank <<" norme= "<<n_norme<<std::endl;
	std::cout<<" rank ="<<rank<<" score= "<<1-glob_sc<<std::endl;
    }
}
