#include<iostream>
#include<armadillo>
#include"Lattice.hpp"
#include<iterator>
#include<algorithm>
#include<fstream>
extern "C" double gsl_sf_lnfact(unsigned int n);

using namespace std;

Lattice::Lattice(unsigned int numberOfBonds):N(numberOfBonds),clusters(N*N){
  bondPerCite = 2;
  M = bondPerCite*N*N;
  bonds.reserve(M);
  fill(clusters.begin(), clusters.end(),-1);
  //random seed
  //srand(time(NULL));

  //set seed
  srand(7);

  p_count = 10000;
  P_inf_of_p.resize(p_count);
  P_infSquared_of_p.resize(p_count);
  s_of_p.resize(p_count);

  binomial.set_size(M);
  writeCounter = 0;
  biggestCluster = 0;
  numberOfBondsActivated = 0;
  average_s = N*N;
  incorrect_s = 0;

}

Lattice::~Lattice(){
}

void Lattice::generateNeighbours(){
  for(int i = 0; i < N*N; i++){
    generateSquareBond(i);
  }
}

void Lattice::generateSquareBond(int position){
  Bond bond;

  if((position+1)%N == 0){
    bond.startPos = position; bond.neighbour = position-N+1;
    bonds.push_back(bond);
  }
  else{
    bond.startPos = position; bond.neighbour = position+1;
    bonds.push_back(bond);
  }

  if(position > N*N-N-1){
    bond.startPos = position; bond.neighbour = position-N*N+N;
    bonds.push_back(bond);
  }
  else{
    bond.neighbour = position + N;
    bonds.push_back(bond);
  }
}

void Lattice::generateTriangularBond(int position){
  Bond bond;

  bond.startPos = position;


  if(!((position+1)%N == 0 || position > N*N-N-1)){
    bond.neighbour = position+1; bonds.push_back(bond);         //Bond to site east
    bond.neighbour = position + N; bonds.push_back(bond);       //bond to site south
    bond.neighbour = position + N + 1; bonds.push_back(bond);   //bond to site to south-east
  }
  else if((position+1)%N == 0 && position != N*N-1){
    bond.neighbour = position-N+1; bonds.push_back(bond);       //Bond to site east
    bond.neighbour = position + N; bonds.push_back(bond);       //bond to site south
    bond.neighbour = position + 1; bonds.push_back(bond);       //bond to site to south-east
  }

  else if(position > N*N-N-1 && position != N*N-1){
    bond.neighbour = position+1; bonds.push_back(bond);         //Bond to site east
    bond.neighbour = position + N - N*N; bonds.push_back(bond);       //bond to site south
    bond.neighbour = position + N + 1 - N*N; bonds.push_back(bond);   //bond to site to south-east
  }
  else{
    bond.neighbour = (N-1)*N; bonds.push_back(bond);         //Bond to site east
    bond.neighbour = N-1; bonds.push_back(bond);       //bond to site south
    bond.neighbour = 0; bonds.push_back(bond);   //bond to site to south-east
  }

}


//Not in use
int Lattice::randomNumberGenerator(int lowerPos){
  randomNumber = rand() % (N*N-lowerPos) + lowerPos; //generate a random number from lowerPos to N-1
  return randomNumber;
}

void Lattice::shuffleList(){
  random_shuffle(bonds.begin(), bonds.end());
}

int Lattice::findRootNode(int nodePos){
  int rootPos = nodePos;
  if (clusters[rootPos] >= 0){
    rootPos =findRootNode(clusters[nodePos]);
    clusters[nodePos] = rootPos;
  }
  return rootPos;
}

void Lattice::activator(){
  double percentageTransformer = 100.0/(double)bonds.size();
  for (int i = 0; i < bonds.size(); i++){
    int siteFrom = bonds[i].startPos;
    int siteTo = bonds[i].neighbour;
    int rootFrom = findRootNode(siteFrom);
    int rootTo = findRootNode(siteTo);

    printf("\r[%.2f%%]",(double)i*percentageTransformer);

    if(rootFrom == rootTo){
    }
    else if(clusters[rootFrom] > clusters[rootTo]){
      average_s -= (pow(clusters[rootTo],2) + pow(clusters[rootFrom],2));
      clusters[rootTo] += clusters[rootFrom];
      average_s += pow(clusters[rootTo],2);
      clusters[rootFrom] = rootTo;
      findBiggestCluster(rootTo);
      /*if(average_s < 0){
        cout << "rootTo largest" << endl;
        cout << rootFrom << "   " << rootTo << endl;
        cout << clusters[rootTo] << "   " << clusters[rootFrom] << endl;
        perror("Incorrect average_s.");
        incorrect_s += 1;
        cout << average_s << endl;
      }*/

    }
    else{
      average_s -= (pow(clusters[rootTo],2) + pow(clusters[rootFrom],2));
      clusters[rootFrom] += clusters[rootTo];
      average_s += pow(clusters[rootFrom],2);
      clusters[rootTo] = rootFrom;
      findBiggestCluster(rootFrom);

    }
    numberOfBondsActivated++;
    if(i == 1000000){                   //numberOfBondsActivated%(2*N*N/(5))
      biggestClusterToFile();
    }
    giantComponent();
    weightedAvgClusterSize();
    susceptibility(i);
    //printClusters();
  }
  convolution();
  printConvolutionResult();
}

void Lattice::findBiggestCluster(int potentiallyBigger){
  if(clusters[biggestCluster] > clusters[potentiallyBigger]){
    biggestCluster = potentiallyBigger;
  }
}

void Lattice::biggestClusterToFile(){
  currentBiggestCluster.clear();
  writeCounter++;

  currentBiggestCluster.push_back(N);
  for(int i = 0; i < N*N; i++){
    if (clusters[i] >= 0){
      if (findRootNode(clusters[i]) == biggestCluster){
        currentBiggestCluster.push_back(i);
      }
    }
  }

  stringstream filename;
  filename << "data/test" <<static_cast<int>(writeCounter)<< ".csv";        //5*numberOfBondsActivated/(2*N*N)
  ofstream output_file(filename.str().c_str());

  ostream_iterator<int> output_iterator(output_file, "\n");
  copy(currentBiggestCluster.begin(), currentBiggestCluster.end(), output_iterator);

}

void Lattice::giantComponent(){
  p_inf = -((double)clusters[biggestCluster])/((double)(N*N));
  p_infValues.push_back(p_inf);
}

void Lattice::weightedAvgClusterSize(){
  double N_d_squared = (double)(N*N);
  if(p_inf >= 1.0){
    w_clusterSize = 0;
  }
  else{
    w_clusterSize = ((double)average_s-pow(N_d_squared*p_inf,2))/(N_d_squared*(1.0-p_inf));
  }

  w_clusterSizeValues.push_back(w_clusterSize);
}

void Lattice::susceptibility(int i){
  p_infAvg = (p_infAvg*i+p_inf)/(i+1);
  p_infSquaredAvg = (p_infSquaredAvg*i+pow(p_inf,2))/(i+1);
  chi = sqrt(p_infSquaredAvg-pow(p_infAvg,2));
  chiValues.push_back(chi);
}

void Lattice::convolution(){
  //log_M = gsl_sf_lnfact(M);
  cout << incorrect_s << endl;
  binomialComputer();

  cout << log_M << ", " << log_n << ", " << log_Mn << p_infValues.size() << endl;
  #pragma omp parallel for
  for (int j = 1; j < p_count; j++){

    printf("\r[%.2f%%]",(double)j/(p_count/100));
    double p = (double)j/p_count;
    double temp_p_inf_of_p = 0;
    double temp_p_infSquared_of_p = 0;
    double temp_s_of_p = 0;

    for (int i = 1; i < M; i++){
      double B_mnp = exp(binomial(i) + i*log(p) + (M-i)*log(1.0-p));

      temp_p_inf_of_p += B_mnp * p_infValues[i];
      temp_p_infSquared_of_p += B_mnp * pow(p_infValues[i],2);
      temp_s_of_p += B_mnp * w_clusterSizeValues[i];

    }
    P_inf_of_p[j] = temp_p_inf_of_p;
    P_infSquared_of_p[j] = temp_p_infSquared_of_p;
    s_of_p[j] =temp_s_of_p;

    //cout << P_inf_of_p[j] << ", " << P_infSquared_of_p[j] << ", " << s_of_p[j] << endl;
  }
  /*for (int i = 0; i < P_inf_of_p.size(); i++){
    cout << "P_inf: "<< P_inf_of_p[i] << " , P_inf² " << P_infSquared_of_p[i] << " , s_of_p: " << s_of_p[i] << endl;
  }*/
}


void Lattice::binomialComputer(){
  log_M = gsl_sf_lnfact(M);
  log_n = 0;
  log_Mn = log_M;
  for (int i = 1; i < M; ++i){
    log_n += log(i);
    log_Mn -= log(M-i+1);
    binomial(i) = log_M - log_n - log_Mn;
  }
}

void Lattice::printClusters(){
  counter = 0;
  for (vector<int>::iterator it=clusters.begin(); it!=clusters.end(); ++it){
    if((counter)%N == 0 && counter > 0){
      std::cout << '\n';
    }
    std::cout << ' ' << *it;
    counter++;
  }
  std::cout << '\n' << "Biggest cluster: " << biggestCluster << '\n';
  std::cout << '\n' << "P: " << (float)numberOfBondsActivated / (float)(bonds.size()) << '\n';
  std::cout << '\n';
}

void Lattice::printConvolutionResult(){
  cout << "hi" << endl;
  cout << P_inf_of_p.size() << endl;
  cout << p_infValues.size() << endl;
  ofstream output_file("data/P_inf_of_p.csv");

  ostream_iterator<double> output_iterator(output_file, "\n");
  copy(P_inf_of_p.begin(), P_inf_of_p.end(), output_iterator);

  output_file.close();
  output_file.clear();

  output_file.open("data/P_infSquared_of_p.csv");
  ostream_iterator<double> output_iterator2(output_file , "\n");
  copy(P_infSquared_of_p.begin(), P_infSquared_of_p.end(), output_iterator2);

  output_file.close();
  output_file.clear();

  output_file.open("data/s_of_p.csv");
  ostream_iterator<double> output_iterator3(output_file , "\n");
  copy(s_of_p.begin(), s_of_p.end(), output_iterator3);
}


/////////////DEBUG////////////////

void DebugLattice::printStructure(Bond bond){
  cout <<"Start posimaition: "<< bond.startPos << "; ";
  cout <<"Bond to : " << bond.neighbour << endl;
}

void DebugLattice::printAllStructures(){
  for(int i = 0; i < bonds.size(); i++){
    printStructure(bonds[i]);
  }
}

void DebugLattice::writeGiantComponentsToFile(){
  ofstream output_file("data/giantComponents.csv");

  ostream_iterator<double> output_iterator(output_file, "\n");
  copy(p_infValues.begin(), p_infValues.end(), output_iterator);

  output_file.close();
  output_file.clear();

  output_file.open("data/weightedAvgClusterSize.csv");
  ostream_iterator<double> output_iterator2(output_file , "\n");
  copy(w_clusterSizeValues.begin(), w_clusterSizeValues.end(), output_iterator2);

  output_file.close();
  output_file.clear();

  output_file.open("data/chiValues.csv");
  ostream_iterator<double> output_iterator3(output_file , "\n");
  copy(chiValues.begin(), chiValues.end(), output_iterator3);
}

void DebugLattice::checkBmnp(){
  arma::vec bmnp_gsl(M);
  arma::vec bmnp_me(M);

  double m = gsl_sf_lnfact(M);
  double n = 0;
  double mn = gsl_sf_lnfact(M-1);

  for (int i = 1; i < M; i++){

    bmnp_gsl(i) = m - gsl_sf_lnfact(i) - gsl_sf_lnfact(M-i);
    bmnp_me(i) = m - n - mn;

    if(bmnp_gsl(i) - bmnp_me(i) > 1E-5){
      cout <<"gsl: "<< bmnp_gsl(i) << ". me: " << bmnp_me(i) << endl;
    }
    n += log(i+1);
    mn -= log(M-i);
  }
}

void DebugLattice::printBonds(){
  generateNeighbours();
  //printAllStructures();
  shuffleList();
  //printAllStructures();
  activator();
  //biggestClusterToFile();
  writeGiantComponentsToFile();
  //checkBmnp();
  /*for (int i = 0; i < 100; i++) {
    cout << randomNumberGenerator(2) << endl;
  }*/
}
