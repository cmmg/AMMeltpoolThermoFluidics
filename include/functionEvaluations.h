//
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Created 2011
//authors: rudraa (2011, 2018)
//

#ifndef FUNCTIONEVALUATIONS_H_
#define FUNCTIONEVALUATIONS_H_
#include <deal.II/fe/fe_values.h>
#include "deformationMap.h"
using namespace dealii;

template <class T, int dim>
void evaluateScalarFunction(FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, Table<1, T>& U){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  
  //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){
    U[q]=0.0; //U
    for (unsigned int k=0; k<dofs_per_cell; ++k){
      if (fe_values.get_fe().system_to_component_index(k).first==DOF){
	U[q]+=ULocal[k]*fe_values.shape_value(k, q); //U
      }
    }
  }
}

template <class T, int dim>
void evaluateScalarFunctionGradient(FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, Table<2, T>& gradU, deformationMap<T, dim>& defMap, bool gradientInCurrentConfiguration){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  Table<1, T> refGradU(dim);
  //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){
    for (unsigned int i=0; i<dim; ++i){refGradU[i]=0.0;}
    for (unsigned int k=0; k<dofs_per_cell; ++k){
      if (fe_values.get_fe().system_to_component_index(k).first==DOF){
	for (unsigned int i=0; i<dim; ++i){
	  refGradU[i]+=ULocal[k]*fe_values.shape_grad(k, q)[i]; //gradU
	}
      }
    }
    //Transform gradient to current configuration. gradW=(F^-T)*GradW
    for (unsigned int i=0; i<dim; ++i){
      if (gradientInCurrentConfiguration==false) gradU[q][i]=refGradU[i];
      else{
	gradU[q][i]=0.0;
	for (unsigned int j=0; j<dim; ++j){
	  gradU[q][i]+=defMap.invF[q][j][i]*refGradU[j];
	}
      }
    }
  }
}

template <class T, int dim>
void evaluateVectorFunction(FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, Table<2, T>& U){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){
    for (unsigned int i=0; i<dim; ++i){
      U[q][i]=0.0;
    }
    for (unsigned int k=0; k<dofs_per_cell; ++k){
      unsigned int ck = fe_values.get_fe().system_to_component_index(k).first - DOF;
      if (ck>=0 && ck<dim){
	U[q][ck]+=ULocal[k]*fe_values.shape_value(k, q); //U
      }
    }
  }
}

template <class T, int dim>
void evaluateVectorFunctionGradient(FEValues<dim>& fe_values, unsigned int DOF, Table<1, T>& ULocal, Table<3, T>& gradU){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;

  //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
	gradU[q][i][j]=0.0;
      }
    }
    for (unsigned int k=0; k<dofs_per_cell; ++k){
      unsigned int ck = fe_values.get_fe().system_to_component_index(k).first - DOF;
      if (ck>=0 && ck<dim){
	for (unsigned int i=0; i<dim; ++i){
	  gradU[q][ck][i]+=ULocal[k]*fe_values.shape_grad(k, q)[i]; //gradU
	}
      }
    }
  }
}

template <class T, int dim>
T determinantOfMinor(unsigned int theRowHeightY, unsigned int theColumnWidthX, Table<2, T>& matrix){
  unsigned int x1 = theColumnWidthX == 0 ? 1 : 0;  /* always either 0 or 1 */
  unsigned int x2 = theColumnWidthX == 2 ? 1 : 2;  /* always either 1 or 2 */
  unsigned int y1 = theRowHeightY   == 0 ? 1 : 0;  /* always either 0 or 1 */
  unsigned int y2 = theRowHeightY   == 2 ? 1 : 2;  /* always either 1 or 2 */
  return matrix[y1][x1]*matrix[y2][x2] - matrix[y1][x2]*matrix[y2][x1];
}


template <class T, int dim>
void getInverse(Table<2, T>& matrix, Table<2, T>& invMatrix, T& det){
  if (dim==1){
    det=matrix[0][0];
    invMatrix[0][0]=1.0/det;
  }
  else if(dim==2){
    det=matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0];
    invMatrix[0][0]=matrix[1][1]/det;
    invMatrix[1][0]=-matrix[1][0]/det;
    invMatrix[0][1]=-matrix[0][1]/det;
    invMatrix[1][1]=matrix[0][0]/det;
  }
  else if(dim==3){
    det=  matrix[0][0]*determinantOfMinor<T, dim>(0, 0, matrix) - matrix[0][1]*determinantOfMinor<T, dim>(0, 1, matrix) +  matrix[0][2]*determinantOfMinor<T, dim>(0, 2, matrix);
    for (int y=0;  y< dim;  y++){
      for (int x=0; x< dim;  x++){
	invMatrix[y][x] = determinantOfMinor<T, dim>(x, y, matrix)/det;
	if( ((x + y) % 2)==1){invMatrix[y][x]*=-1;}
      }
    }
  }
  else throw "dim>3";
  if (std::abs(det)< 1.0e-15){
    printf("**************Near zero determinant in Matrix inversion***********************\n"); throw "Near zero determinant in Matrix inversion";
  }
}

#endif /* FUNCTIONEVALUATIONS_H_ */

