#include "Activator.hpp"
#include "Lattice.hpp"
#include <aramdillo>

extern "C" double gsl_sf_lnfact(unsigned int n);

using namespace std;

Activator::Activator(unsigned int numberOfBonds, LatticeType_type lattice):N(numberOfBonds),lattice(lattice){
  n_of_sites = N*N;
  int bonds;
  switch(lattice){
    case LatticeType_type::SQU:
      bonds = 2*n_of_sites;
      break;
    case LatticeType_type::SQU:
      bonds = 3*n_of_sites;
      break;
    case LatticeType_type::SQU:
      //bonds = 2*n_of_sites; TO DO
      break;
  }
  p_infValues.resize(bonds);
  w_clusterSizeValues.resize(bonds);
  chiValues.resize(bonds);

  P_inf_of_p.resize(bonds);
  P_infSquared_of_p.resize(bonds);
  s_of_p.resize(bonds);

}
