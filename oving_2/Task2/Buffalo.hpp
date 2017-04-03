#ifndef BUFFALO
#define BUFFALO

class Buffalo{
public:
  Buffalo(unsigned int numberOfSites, double occupationProb);  // set number of sites in X and Y direction and occupation probability
  ~Buffalo();
  void initialize();

protected:
  unsigned int N;                       // number of sites in X and Y direction
  double p;                             // occupation probability
  bool** occupied;                      // NxN array to store state of lattice
                                        // <type>** means pointer to a <type> pointer
  int clusters;                         // number of clusters
  int** label;                          // NxN array of cluster labels for each of the sites
  int currentLabel;
  int spanningLabel;

  void newSample();
  void addNewNeighbour(int i, int j);

};

class DebugBuffalo: public Buffalo{
public:
  DebugBuffalo(unsigned int N, double p):Buffalo(N,p){};

  void debug();

};

#endif
