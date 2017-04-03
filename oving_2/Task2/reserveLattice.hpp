#ifndef LATTICE
#define LATTICE
#include<armadillo>


class Lattice{
public:
  Lattice(unsigned int numberOfBonds);

  ~Lattice();

  void clusterize();

  void save();

protected:
  unsigned int N;

  struct Bond{
    int startPos;
    int neighbour;
  };

  int M;
  int bondPerCite;
  int counter;
  int randomNumber;
  int biggestCluster;
  int numberOfBondsActivated;
  int writeCounter;
  double p_inf;
  double w_clusterSize;
  int average_s;
  double chi;
  double p_infAvg;
  double p_infSquaredAvg;
  //double B_mnp;
  double log_M{0};
  double log_n{0};
  double log_Mn{0};
  int p_count{0};

  void generateSquareBond(int position);
  void generateTriangularBond(int position);
  void generateNeighbours();
  int randomNumberGenerator(int lowerPos);
  void shuffleList();
  int findRootNode(int nodePos);
  void activator();
  void printClusters();
  void findBiggestCluster(int potentiallyBigger);
  void biggestClusterToFile();
  void giantComponent();
  void weightedAvgClusterSize();
  void susceptibility(int i);
  void convolution();
  void printConvolutionResult();
  void binomialComputer();

  arma::vec binomial;

  std::vector<Bond> bonds;
  arma::vec sites;
  std::vector<int> clusters;
  std::vector<int> currentBiggestCluster;
  std::vector<double> p_infValues;
  std::vector<double> w_clusterSizeValues;
  std::vector<double> chiValues;
  std::vector<double> P_inf_of_p;
  std::vector<double> P_infSquared_of_p;
  std::vector<double> s_of_p;

};

class DebugLattice: public Lattice{
public:
  DebugLattice(unsigned int N):Lattice(N){};
  void printBonds();
  void printStructure(Bond bond);
  void printAllStructures();
  void writeGiantComponentsToFile();
  void checkBmnp();

};

#endif
