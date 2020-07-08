//
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Created 2018
//authors: rudraa (2018)
//

#ifndef DEFORMATIONMAP_H_
#define DEFORMATIONMAP_H_
#include <deal.II/fe/fe_values.h>
#include "functionEvaluations.h"
using namespace dealii;

//Forward declarations
template <class T, int dim>
  void evaluateVectorFunctionGradient(FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, Table<3, T>& gradU);
template <class T, int dim>
  void getInverse(Table<2, T>& matrix, Table<2, T>& invMatrix, T& det);

//deformationMap structure
template <class T, int dim>
struct deformationMap{
deformationMap(unsigned int n_q_points): F(n_q_points, dim, dim), invF(n_q_points, dim, dim), detF(n_q_points){}
  Table<3, T> F, invF;
  Table<1, T> detF;
};

//Compute deformation map
template <class T, int dim>
void getDeformationMap(FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, deformationMap<T, dim>& defMap, const unsigned int iteration){
  unsigned int n_q_points= fe_values.n_quadrature_points;
  //evaluate dx/dX
  Table<3, T> gradU(n_q_points, dim, dim);
  evaluateVectorFunctionGradient<T, dim>(fe_values, DOF, ULocal, gradU);

  //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){
    Table<2, T > Fq(dim, dim), invFq(dim, dim); T detFq;
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
	defMap.F[q][i][j] = Fq[i][j] = (i==j) + gradU[q][i][j]; //F (as double value)
      }
    }
    getInverse<T, dim>(Fq, invFq, detFq); //get inverse(F)
    defMap.detF[q] = detFq;
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
	defMap.invF[q][i][j] = invFq[i][j];
      }
    }
    //detF
    if (defMap.detF[q].val()<=1.0e-15 && iteration==0){
      printf("**************Non positive jacobian detected**************. Value: %12.4e\n", defMap.detF[q].val());
      for (unsigned int i=0; i<dim; ++i){
	for (unsigned int j=0; j<dim; ++j) printf("%12.6e  ", defMap.F[q][i][j].val());
	printf("\n"); exit(-1);
      }
      //throw "Non positive jacobian detected";
    }
  }
}

#endif /*DEFORMATIONMAP_H_ */

