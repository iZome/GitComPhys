#include<iostream>
#include"Buffalo.hpp"

using namespace std;

Buffalo::Buffalo(unsigned int numberOfSites, double occupationProb):N(numberOfSites){
  p = occupationProb;
}

Buffalo::~Buffalo(){

}

void Buffalo::addNewNeighbour(int i, int j){
  if(occupied[i][j] && label[i][j] == 0){
    label[i][j] = currentLabel;
    if(i < N-1){
      addNewNeighbour(i+1,j);
    }
    if(i > 0){
      addNewNeighbour(i-1,j);
    }
    if(j > 0){
      addNewNeighbour(i,j-1);
    }
    if(j < N-1){
      addNewNeighbour(i,j+1);
    }
  }
}

void Buffalo::newSample(){
  //visit all sites and occupy it with probability p
  for(int i = 0; i < N; i++){
    for(int j = 0; j < N; j++){
      occupied[i][j] = (rand() / (RAND_MAX + 1.0)) < p;
      label[i][j] = 0;
    }
  }

  //find and label all cluters of occupied sites
  clusters = 0;
  for(int i = 0; i < N; i++){
    for(int j = 0; j < N; j++){
      if(occupied[i][j] && label[i][j] == 0){
        currentLabel = ++clusters;
        addNewNeighbour(i,j);
      }
    }
  }

  //checking each cluster for percolation
  spanningLabel = 0;
  for(int cluster = 1; cluster <= clusters; ++cluster){
    bool west = false;
    for(int j = 0; j < N; j++){
      if(label[0][j] == cluster){west = true; break;}
    }

    bool east = false;
    for(int j = 0; j < N; j++){
      if(label[N-1][j] == cluster){east = true; break;}
    }

    bool north = false;
    for(int i = 0; i < N; i++){
      if(label[i][0] == cluster){north = true; break;}
    }

    bool south = false;
    for(int i = 0; i < N; i++){
      if(label[i][N-1] == cluster){south = true; break;}
    }
    // checks if cluster goes to all boundaries
    if(west && east && north && south){
        spanningLabel = cluster;
        break;
    }
  }
}

void Buffalo::initialize(){
  //allocate memory for arrays occupied and label
  occupied = new bool* [N];
  label = new int* [N];
  for (int i = 0; i < N; i++){
    occupied[i] = new bool [N];
    label [i] = new int [N];
  }
}

///////DEBUG///////

void DebugBuffalo::debug(){
  initialize();
  newSample();
  /*for(int i = 0; i < N; i++){
    for(int j = 0; j < N; j++){
      cout << "Occupied: " << occupied[i][j] << endl;
      cout << "Label: " << label[i][j] << endl;
    }
  }*/
  cout << spanningLabel << endl;
}
