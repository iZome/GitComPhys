#ifndef ACTIVATOR
#define ACTIVATOR
#include<armadillo>


enum class LatticeType_type{
  SQU, TRI, HON
};

class Lattice;

class Activator{
public:
  Activator(int N, LatticeType_type);
  ~Activator();

protected:
  int n_of_sites;

  void convolution();
  void printConvolutionResult();
  void binomialComputer();

  std::vector<double> p_infValues;
  std::vector<double> w_clusterSizeValues;
  std::vector<double> chiValues;
  std::vector<double> P_inf_of_p;
  std::vector<double> P_infSquared_of_p;
  std::vector<double> s_of_p;


private:

}
