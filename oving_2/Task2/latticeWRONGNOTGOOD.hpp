#ifndef LATTICEWRONG
#define LATTICEWRONG
#include<armadillo>

class lattice{
public:
  lattice(unsigned int N_discret);

  void matrixBuilder();

  void save();

  ~lattice();

protected:
  unsigned int N{0};

  arma::vec values;
  arma::umat location;
  arma::sp_mat *matrix{NULL};


  void diagonal();

  void subDiagonal();

  void superDiagonal();

  void subsubDiagonal();

  void supersuperDiagonal();

  void initializeMatrix();

};

class DebugLattice: public lattice{
public:
  DebugLattice(unsigned int N):lattice(N){};
  void printMatrix();

};

#endif
