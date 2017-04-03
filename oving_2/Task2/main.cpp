#include<iostream>
#include"Lattice.hpp"
#define DEBUG

using namespace std;

int main(){
Lattice test(5);

  #ifdef DEBUG
    DebugLattice debug(1000);
    debug.printBonds();
  #endif

  //test.generateBonds(1);
  return 0;
}
