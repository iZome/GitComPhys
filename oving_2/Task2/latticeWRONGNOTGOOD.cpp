#include<iostream>
#include<armadillo>
#include"lattice.hpp"

using namespace std;

lattice::lattice(unsigned int N_discret):N(N_discret){
  unsigned int length = 3*N*N;
  location.set_size(2,length);
  values.set_size(length);
}

lattice::~lattice(){
  if(matrix != NULL){
    delete matrix;
  }
}

void lattice::diagonal(){
  for(int i = 0; i < N*N; i++){
    location(0,i) = i;
    location(1,i) = i;
    values(i) = 1;
  }
}

void lattice::subDiagonal(){
  int startPos = N*N;
  int counter = 0;
  for(int i = startPos; i < 2*startPos-N+1; i++){
    int index_location = startPos+counter;
    if((i)%N != 0){
      continue;
    }



    location(0,index_location) = i+N-1-startPos;
    location(1,index_location) = i-startPos;
    values(index_location) = 1;
    counter++;
  }
}

void lattice::superDiagonal(){
  int startPos = N*N+N;
  int counter = 0;
  for(int i = startPos; i < startPos + N*N -N+1; i++){
    int index_location = startPos+counter;
    if((i+1)%N == 0){
      continue;
    }

    location(0,index_location) = i-startPos;
    location(1,index_location) = i+1-startPos;
    values(index_location) = 1;
    counter++;
  }
}

void lattice::subsubDiagonal(){
  int startPos = 2*N*N;
  for(int i = startPos; i < startPos + N; i++){
    location(0,i) = i - startPos + N*N - N;
    location(1,i) = i - startPos;
    values(i) = 1;
  }
}

void lattice::supersuperDiagonal(){
  int startPos = 2*N*N+N;
  for(int i = startPos; i < startPos + N*N -N; i++){
    location(0,i) = i - startPos ;
    location(1,i) = i - startPos + N;
    values(i) = 1;
  }
}

void lattice::initializeMatrix(){
  diagonal();
  subDiagonal();
  superDiagonal();
  subsubDiagonal();
  supersuperDiagonal();
  cout << "hi" << endl;
}

void lattice::matrixBuilder(){
  initializeMatrix();
  if(matrix != NULL){
    delete matrix;
  }
  matrix = new arma::sp_mat(location,values,true);
  cout << "Thomas liker deg!" << endl;
}


/////////////DEBUG////////////////

void DebugLattice::printMatrix(){
  matrixBuilder();
  arma::mat denseMatrix(*matrix);
  cout << denseMatrix << endl;
}
