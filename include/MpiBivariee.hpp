#ifndef MPIBIVARIEE_HPP
#include<string>
#include <list>
#include<vector>

struct Arguments
{
  int n_test;
  int n_train;
  double * __restrict__ X_train;
  double * __restrict__ Y_train;
  
  double * __restrict__ X_test;
  double * __restrict__ Y_test;
};


class RegressionBivariee
{
public:
  int Iteration=0;
  int MaxIteration = 0;
  int m_loc=1;
  int begin = 0;
  int end = 1;
  int offset = 0;
  int m=1;
  int n=2;
  int rank = 0;
  int nproc = 2;
  
  double moy = 0;
  double err2 = 0.;
  double eps = 1e-5;
  double score=0.;
  double alpha=1.;

  double* theta;
  double* theta0;
  

  // Fonctions membre
  
  double h(double * __restrict__ X_i, int n)const;

  double * H(const double * __restrict__ X, double * __restrict__ Y, int m, int n) const;

   double JBivariee(double * __restrict__ X, double * __restrict__ Y, int m, int n)const;

  double ErreurDeConsistance(double * __restrict__ X, double * __restrict__ Y, double * __restrict__ theta, int m, int n)const;


  void fit(double * __restrict__ X, double * __restrict__  Y, double alpha, int m, int n);
  
  double Score(double * __restrict__ X, double * __restrict__ y, int m, int n)const;
  
  void Split_train_test(const double * __restrict__ X, const double * __restrict__ Y, double * __restrict__ X_train, double * __restrict__ X_test, double * __restrict__ Y_train, double * __restrict__ Y_test, double proba_train = 0.70, int m_train = 1, int m_test=1) const;   // divise le dataset

  // Pour la sortie 
  void DataFileBivariee(const double * __restrict__ X_train, const double * __restrict__ Y_train,int n_train, const double * __restrict__ X_test, const double * __restrict__ Y_test, int n_test, std::string & file_name) const;

  void File_create(const std::string & file_name, int m,int seed)const;
  
  void Plot(const std::list<double *> & Arguments, std::string & gnplot_script ) const;
  
  void Print();

  // Scipt Gnuplot
  void Script_gnuplot3d(std::string & out_put_gnuplot_file, std::string & input_plot_file, std::string imp_name=std::string("Image_Reg_Bivariee.png"));

  
  // Constructeur
  RegressionBivariee(int n_);
  RegressionBivariee(const std::string & fil_name , double * __restrict__ X, double * __restrict__ Y,int m, int n=0,int size_of_one_line=69);
  
  // Destructeur
  ~RegressionBivariee() { delete [] theta; delete [] theta0;}
  
};












//double * init(int rank,double *x, int m);




#endif
