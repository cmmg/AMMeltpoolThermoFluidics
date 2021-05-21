//
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Created 2012
//authors: rudraa (2012, 2018)
//

#ifndef PROJECTION_H_
#define PROJECTION_H_
#include "functionEvaluations.h"
#include "supplementaryFunctions.h"
#include <math.h>     

//Chemistry residual implementation
template <int dim>
void residualForProjection(FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1, double >& ULocal, dealii::Table<1, Sacado::Fad::DFad<double> >& Pr_ULocal, dealii::Table<1, Sacado::Fad::DFad<double> >& Rp, double currentTime, double totalTime) {
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  
  unsigned int n_q_points_face= fe_face_values.n_quadrature_points;
  const unsigned int faces_per_cell = GeometryInfo<dim>::faces_per_cell;
  
  //Velocity and Pressure  
  dealii::Table<2,Sacado::Fad::DFad<double> > phi_j(n_q_points, dim);
  dealii::Table<3,double > vel_j(n_q_points, dim,dim);


  
  //Interpolate on all cells 
  for (unsigned int q=0; q<n_q_points; ++q) {    
    for (unsigned int j=0; j<dim; j++) {
      phi_j[q][j]=0;
      for (unsigned int k=0; k<dim; k++) {
	vel_j[q][j][k]=0;
      }      
    }
    
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;

      if (ck>=3 && ck < 6) {
	for (unsigned int j=0; j<dim; j++) {
	  vel_j[q][ck-dim][j]+=fe_values.shape_grad_component(i, q, ck)[j]*ULocal[i];
	}
      }

      else if (ck==6) {
	for (unsigned int j=0; j<dim; j++) {
	  phi_j[q][j]+=fe_values.shape_grad_component(i, q, ck)[j]*Pr_ULocal[i];
	  
	}	
      }
      
    }
        
  }


  //evaluate Residual on cell
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;        
    for (unsigned int q=0; q<n_q_points; ++q) {              
      if (ck==6) {
	for (unsigned int j=0; j<dim; j++) {
	  Rp[i]+=-(1.0)*fe_values.shape_grad_component(i, q, ck)[j]*(phi_j[q][j])*fe_values.JxW(q);
	}
	//if kinetic pressure is used
	Rp[i]+=-(1.5/dt)*fe_values.shape_value_component(i, q, ck)*(vel_j[q][0][0])*fe_values.JxW(q);
        Rp[i]+=-(1.5/dt)*fe_values.shape_value_component(i, q, ck)*(vel_j[q][1][1])*fe_values.JxW(q);
	Rp[i]+=-(1.5/dt)*fe_values.shape_value_component(i, q, ck)*(vel_j[q][2][2])*fe_values.JxW(q);
	//if dynamic pressure is used
	//Rp[i]+=-(RHO)*(1.5/dt)*fe_values.shape_value_component(i, q, ck)*(vel_j[q][0][0])*fe_values.JxW(q);
        //Rp[i]+=-(RHO)*(1.5/dt)*fe_values.shape_value_component(i, q, ck)*(vel_j[q][1][1])*fe_values.JxW(q);

			
      }               
    }
  }

  
}

#endif /* PROJECTION_H_ */
