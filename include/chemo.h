//
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Created 2012
//authors: rudraa (2012, 2018)
//

#ifndef CHEMO_H_
#define CHEMO_H_
#include "functionEvaluations.h"
#include "supplementaryFunctions.h"
#include <math.h>     

//Chemistry residual implementation
template <int dim>
void residualForChemo(FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1, Sacado::Fad::DFad<double> >& ULocal, dealii::Table<1, double>& ULocalConv, dealii::Table<1, double>& ULocalConvConv, dealii::Table<1, double>& Pr_ULocalConv,dealii::Table<1, double>& Pr_ULocalConvConv, dealii::Table<1, Sacado::Fad::DFad<double> >& R, double currentTime, double totalTime) {
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  
  unsigned int n_q_points_face= fe_face_values.n_quadrature_points;
  const unsigned int faces_per_cell = GeometryInfo<dim>::faces_per_cell;

  dealii::Table<2,Sacado::Fad::DFad<double> > vel(n_q_points,dim);
  dealii::Table<3,Sacado::Fad::DFad<double> > vel_j(n_q_points,dim,dim);
  dealii::Table<1,double> press_conv(n_q_points),phi_conv(n_q_points),phi_conv_conv(n_q_points);      
  dealii::Table<2,double> vel_conv(n_q_points,dim),vel_conv_conv(n_q_points,dim),vel_star(n_q_points,dim);
  dealii::Table<3,double> vel_conv_j(n_q_points,dim,dim),vel_conv_conv_j(n_q_points,dim,dim),vel_star_j(n_q_points,dim,dim);
  
  //Interpolate on all cells 
  for (unsigned int q=0; q<n_q_points; ++q) {  
    press_conv[q]=0; phi_conv[q]=0; phi_conv_conv[q]=0;

    //initialization : filling it with zero and later updated
    for (unsigned int j=0; j<dim; j++) {
      vel[q][j]=0; vel_conv[q][j]=0; vel_conv_conv[q][j]=0; vel_star[q][j]=0;
      for (unsigned int k=0; k<dim; k++) {
	vel_j[q][j][k]=0;
	vel_star_j[q][j][k]=0;
	vel_conv_j[q][j][k]=0;
	vel_conv_conv_j[q][j][k]=0;
      }      
    }
    
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;      
      if (ck>=0 && ck<2) {
	vel[q][ck]+=fe_values.shape_value_component(i, q, ck)*ULocal[i];
	vel_conv[q][ck]+=fe_values.shape_value_component(i, q, ck)*ULocalConv[i];
	vel_conv_conv[q][ck]+=fe_values.shape_value_component(i, q, ck)*ULocalConvConv[i];     	
	for (unsigned int j=0; j<dim; j++) {
	  vel_j[q][ck][j]+=fe_values.shape_grad_component(i, q, ck)[j]*ULocal[i];
	  vel_conv_j[q][ck][j]+=fe_values.shape_grad_component(i, q, ck)[j]*ULocalConv[i];
	  vel_conv_conv_j[q][ck][j]+=fe_values.shape_grad_component(i, q, ck)[j]*ULocalConv[i];
	}	
      }
      else if (ck==2) {
	press_conv[q]+=fe_values.shape_value_component(i, q, ck)*ULocalConv[i];
	phi_conv[q]+=fe_values.shape_value_component(i, q, ck)*Pr_ULocalConv[i];
	phi_conv_conv[q]+=fe_values.shape_value_component(i, q, ck)*Pr_ULocalConvConv[i];	
      }
           
    }

     for (unsigned int j=0; j<dim; j++) {
       vel_star[q][j]=2.0*vel_conv[q][j]-vel_conv_conv[q][j];
       for (unsigned int k=0; k<dim; k++) {
	vel_star_j[q][j][k]=2.0*vel_conv_j[q][j][k]-vel_conv_conv_j[q][j][k];
      }      
    }
          
  }
 
  //evaluate Residual on cell
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;         
    for (unsigned int q=0; q<n_q_points; ++q) {        
      if(ck>=0 && ck < 2) {
	//Massterm
       R[i]+=(0.5/dt)*fe_values.shape_value_component(i, q, ck)*(3.0*vel[q][ck]-4.0*vel_conv[q][ck]+vel_conv_conv[q][ck])*fe_values.JxW(q);

       //Laplace term and pressure
	for (unsigned int j = 0; j < dim; j++){
	  R[i]+=(nu)*fe_values.shape_grad_component(i, q, ck)[j]*(vel_j[q][ck][j])*fe_values.JxW(q);		 	  
	}
	R[i]+=-(1.0)*fe_values.shape_grad_component(i, q, ck)[ck]*(press_conv[q])*fe_values.JxW(q);
	R[i]+=-(4.0/3.0)*fe_values.shape_grad_component(i, q, ck)[ck]*(phi_conv[q])*fe_values.JxW(q);
	R[i]+=-(-1.0/3.0)*fe_values.shape_grad_component(i, q, ck)[ck]*(phi_conv_conv[q])*fe_values.JxW(q);
	
	
	//advection term
	//first part	
	for (unsigned int j = 0; j < dim; j++){	  
	  R[i]+=fe_values.shape_value_component(i, q, ck)*(vel_star[q][j]*vel_j[q][ck][j])*fe_values.JxW(q);
	}
	//second part
	R[i]+=0.5*fe_values.shape_value_component(i, q, ck)*(vel_star_j[q][0][0]*vel[q][ck])*fe_values.JxW(q);
	R[i]+=0.5*fe_values.shape_value_component(i, q, ck)*(vel_star_j[q][1][1]*vel[q][ck])*fe_values.JxW(q);  	
      }
          
                                 
    }
  }



  
}

#endif /* CHEMO_H_ */
