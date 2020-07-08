//
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Created 2012
//authors: rudraa (2012, 2018)
//

#ifndef CHEMO_H_
#define CHEMO_H_
#include "functionEvaluations.h"
#include "supplementaryFunctions.h"

//Chemistry residual implementation
template <int dim>
void residualForChemo(FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1, Sacado::Fad::DFad<double> >& ULocal, dealii::Table<1, double>& ULocalConv, dealii::Table<1, Sacado::Fad::DFad<double> >& R, double currentTime, double totalTime, dealii::Table<1,double>& c_conv,  dealii::Table<1,double>& phi_conv ) {
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
    
  //evaluate gradients 
  dealii::Table<1,Sacado::Fad::DFad<double> > c(n_q_points), phi(n_q_points),  mu(n_q_points);
  dealii::Table<2,Sacado::Fad::DFad<double> > c_j(n_q_points, dim), phi_j(n_q_points,dim),mu_j(n_q_points, dim);
  dealii::Table<2,double> c_conv_j(n_q_points, dim) , phi_conv_j(n_q_points,dim) ;
  for (unsigned int q=0; q<n_q_points; ++q) {
    c[q]=0.0; c_conv[q]=0; phi_conv[q]=0.0 ; mu[q]=0.0; 
    for (unsigned int j=0; j<dim; j++) {c_j[q][j]=0.0; c_conv_j[q][j]=0.0;  phi_j[q][j]=0.0 ;phi_conv_j[q][j]=0.0 ;mu_j[q][j]=0.0;}
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
      if (ck==0) { c[q]+=fe_values.shape_value(i, q)*ULocal[i]; c_conv[q]+=fe_values.shape_value(i, q)*ULocalConv[i];}
      else if (ck==1){ phi[q]+=fe_values.shape_value(i, q)*ULocal[i]; phi_conv[q]+=fe_values.shape_value(i, q)*ULocalConv[i];  }
      else if (ck==2){ mu[q]+=fe_values.shape_value(i, q)*ULocal[i];  }
   
      
      for (unsigned int j=0; j<dim; j++) {
	if (ck==0) {
	  c_j[q][j]+=fe_values.shape_grad(i, q)[j]*ULocal[i];
	  // c_conv_j[q][j]+=fe_values.shape_grad(i, q)[j]*ULocalConv[i];
	}
	else if (ck==1)
	  { phi_j[q][j]+=fe_values.shape_grad(i, q)[j]*ULocal[i];
	    //  phi_conv_j[q][j]+=fe_values.shape_grad(i, q)[j]*ULocalConv[i];
	}
	
      }
      
      // theta[q]+= atan2(mu_j[q][1],mu_j[q][0]) ; 
    }

    
    
  }
  
  //evaluate Residual
  double Kappa[] =InterfaceEnergyParameter;
  Sacado::Fad::DFad<double> M= Mobility;
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    for (unsigned int q=0; q<n_q_points; ++q) {
      Sacado::Fad::DFad<double> dfdc  = 1*(phi[q]- lam*c[q] + lam*c[q]*phi[q]*phi[q])*(1-phi[q]*phi[q]);
      Sacado::Fad::DFad<double>  theta = atan2(phi_j[q][1],phi_j[q][0]) ;
      //std::cout <<theta.val()<<std::endl;     
      Sacado::Fad::DFad<double> W = W0*(1+ epm*cos(mm*(theta-theta0))) ;
      Sacado::Fad::DFad<double> tau = tau0*(1+ epm*cos(mm*(theta-theta0)))  ;
      dealii::Table<1,Sacado::Fad::DFad<double> > BIGTERM(dim) ;
      BIGTERM[0]=W*W*phi_j[q][0] + W0*epm*mm*W*sin(mm*(theta-theta0))*phi_j[q][1] ;
      BIGTERM[1]= W*W*phi_j[q][1] - W0*epm*mm*W*sin(mm*(theta-theta0))*phi_j[q][0] ;
      // std::cout <<BIGTERM[0].val()<<std::endl;
      
	
      if (ck==0) {
	R[i] +=  (1/dt)*fe_values.shape_value(i, q)*(c[q]-c_conv[q])*fe_values.JxW(q);  
	R[i] -= (0.5/tau)*(fe_values.shape_value(i, q)*(mu[q])*fe_values.JxW(q));   
	for (unsigned int j = 0; j < dim; j++) {	  
	  R[i] += fe_values.shape_grad(i, q)[j]*D*c_j[q][j]*fe_values.JxW(q);  
	}
      }
      
      else if(ck==1){
	R[i] +=  fe_values.shape_value(i, q)*(phi[q] - phi_conv[q])*fe_values.JxW(q); 
	R[i] -= (1.0/tau)*fe_values.shape_value(i, q)*(mu[q])*fe_values.JxW(q);  
      }

      else if (ck==2) {
	R[i] +=  fe_values.shape_value(i, q)*(mu[q])*fe_values.JxW(q);
	R[i] -=  fe_values.shape_value(i, q)*(dfdc)*fe_values.JxW(q);

	for (unsigned int j = 0; j < dim; j++) {	 
	   R[i]+= fe_values.shape_grad(i, q)[j]*BIGTERM[j]*fe_values.JxW(q) ;
	   //  R[i]+= (1.0/tau)*(fe_values.shape_grad(i, q)[j]*phi_j[q][j]*fe_values.JxW(q));
	}

	

	
      }

      
    }
  } 

  
  
}

#endif /* CHEMO_H_ */
