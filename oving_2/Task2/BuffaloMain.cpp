#include <iostream>
#include"Buffalo.hpp"
#define DEBUG

using namespace std;

int main(){
  Buffalo percolation(3, 0.65);

  #ifdef DEBUG
    DebugBuffalo debug(500,0.80);
    debug.debug();
  #endif

  percolation.initialize();

  return 0;
}
