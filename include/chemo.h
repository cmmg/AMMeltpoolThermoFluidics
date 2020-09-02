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
void residualForChemo(FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1, Sacado::Fad::DFad<double> >& ULocal, dealii::Table<1, double>& ULocalConv, dealii::Table<1, double>& ULocalConvConv, dealii::Table<1, double>& Pr_ULocalConv,dealii::Table<1, double>& Pr_ULocalConvConv, dealii::Table<1, double>& T_ULocalConv,dealii::Table<1, double>& T_ULocalConvConv, dealii::Table<1, Sacado::Fad::DFad<double> >& R, double currentTime, double totalTime) {
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  
  unsigned int n_q_points_face= fe_face_values.n_quadrature_points;
  const unsigned int faces_per_cell = GeometryInfo<dim>::faces_per_cell;
  
  dealii::Table<1,Sacado::Fad::DFad<double> > T(n_q_points), liquid(n_q_points);
  dealii::Table<1,double> T_conv(n_q_points), liquid_conv(n_q_points);
  dealii::Table<1,double> LiqfaceConv(n_q_points);
  dealii::Table<2,Sacado::Fad::DFad<double> > Tfaceconv_j(n_q_points,dim);

  dealii::Table<2,Sacado::Fad::DFad<double> > vel(n_q_points,dim);
  dealii::Table<3,Sacado::Fad::DFad<double> > vel_j(n_q_points,dim,dim),gamma_tense(n_q_points,dim,dim);
  dealii::Table<1,double> press_conv(n_q_points),phi_conv(n_q_points),phi_conv_conv(n_q_points);      
  dealii::Table<2,double> vel_conv(n_q_points,dim),vel_conv_conv(n_q_points,dim),vel_star(n_q_points,dim);
  dealii::Table<3,double> vel_conv_j(n_q_points,dim,dim),vel_conv_conv_j(n_q_points,dim,dim),vel_star_j(n_q_points,dim,dim);

  Sacado::Fad::DFad<double>  RHO_T;
  //Interpolate on all cells 
  for (unsigned int q=0; q<n_q_points; ++q) {  
    press_conv[q]=0; phi_conv[q]=0; phi_conv_conv[q]=0;
    T[q]=0; T_conv[q]=0;liquid[q]=0; liquid_conv[q]=0;
    //initialization : filling it with zero and later updated
    for (unsigned int j=0; j<dim; j++) {
      vel[q][j]=0; vel_conv[q][j]=0; vel_conv_conv[q][j]=0; vel_star[q][j]=0;
      for (unsigned int k=0; k<dim; k++) {
	vel_j[q][j][k]=0;
	vel_star_j[q][j][k]=0;
	vel_conv_j[q][j][k]=0;
	vel_conv_conv_j[q][j][k]=0;
	gamma_tense[q][j][k]=0;
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
      
      else if (ck==3) {
	T_conv[q]+=fe_values.shape_value_component(i, q, ck)*T_ULocalConv[i];
     }
    
      else if (ck==4) {
	liquid_conv[q]+=fe_values.shape_value_component(i, q, ck)*T_ULocalConv[i];
      }
                
    }

     for (unsigned int j=0; j<dim; j++) {
       vel_star[q][j]=2.0*vel_conv[q][j]-vel_conv_conv[q][j];
       for (unsigned int k=0; k<dim; k++) {
	vel_star_j[q][j][k]=2.0*vel_conv_j[q][j][k]-vel_conv_conv_j[q][j][k];
	gamma_tense[q][j][k]=0.5*vel_j[q][j][k]+0.5*vel_j[q][k][j];
      }      
    }
     
     
  }
  

   //Interpolate over faces
  for (unsigned int f=0; f < faces_per_cell; f++) { 
    fe_face_values.reinit (cell, f);
    // double CHECKL=cell->face(f)->center()[0];
      double CHECKH=cell->face(f)->center()[1];
      //double CHECKW=cell->face(f)->center()[2];      
    //    if(CHECKL ==0 ||CHECKL== problem_Length||CHECKW ==0 ||CHECKW== problem_Width||CHECKH== problem_Height) {
    if(cell->face(f)->center()[1] == problemHeight) {
      for (unsigned int q=0; q<n_q_points_face; ++q) {
	LiqfaceConv[q]=0;
	for (unsigned int j=0; j<dim ; ++j) { Tfaceconv_j[q][j]=0; }
	for (unsigned int i=0; i<dofs_per_cell; ++i) {
	  const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
	  if (ck==3) {
	    for (unsigned int j=0; j<dim ; ++j) {
	      Tfaceconv_j[q][j]+=fe_face_values.shape_grad_component(i, q,ck)[j]*T_ULocalConv[i];
	    }	    
	  }
	  else if (ck==4) { LiqfaceConv[q]+=fe_face_values.shape_value_component(i, q,ck)*T_ULocalConv[i];}
	  
	}	
      }
    }    
  }

 
  //evaluate Residual on cell
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;         
    for (unsigned int q=0; q<n_q_points; ++q) {        
      if(ck>=0 && ck<2) {
	
	//Massterm
	R[i]+=(0.5/dt)*fe_values.shape_value_component(i, q, ck)*(3.0*vel[q][ck]-4.0*vel_conv[q][ck]+vel_conv_conv[q][ck])*fe_values.JxW(q);
	       
	//viscous term
	for (unsigned int j = 0; j < dim; j++){
	  //R[i]+= (mu/RHO)*fe_values.shape_grad_component(i, q, ck)[j]*(vel_j[q][ck][j])*fe_values.JxW(q);
	  R[i]+= (mu/RHO)*fe_values.shape_grad_component(i, q, ck)[j]*(gamma_tense[q][ck][j])*fe_values.JxW(q);
	  R[i]+= (mu/RHO)*fe_values.shape_grad_component(i, q, j)[ck]*(gamma_tense[q][ck][j])*fe_values.JxW(q);
	}

	 // pressure (kineitc)
	R[i]+=-(1.0)*fe_values.shape_grad_component(i, q, ck)[ck]*(press_conv[q])*fe_values.JxW(q);
	R[i]+=-(4.0/3.0)*fe_values.shape_grad_component(i, q, ck)[ck]*(phi_conv[q])*fe_values.JxW(q);
	R[i]+=-(-1.0/3.0)*fe_values.shape_grad_component(i, q, ck)[ck]*(phi_conv_conv[q])*fe_values.JxW(q);
	
	//free convection
	if (ck==1) R[i]+=-fe_values.shape_value_component(i, q, ck)*(gravity*BETA*(T_conv[q]-TSS))*fe_values.JxW(q);
	
	double AA;
	//presure drop due to mush zone		
	AA=(1.0e+08)*((1.0-liquid_conv[q])*(1.0-liquid_conv[q]))/(liquid_conv[q]*liquid_conv[q]*liquid_conv[q]+1.0e-03); 	
	R[i]+=fe_values.shape_value_component(i, q, ck)*((AA)*vel[q][ck])*fe_values.JxW(q);
      

	
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

   //surface integral
  for (unsigned int f=0; f < faces_per_cell; f++) { 
    fe_face_values.reinit (cell, f);
    //double CHECKL=cell->face(f)->center()[0];
    double CHECKH=cell->face(f)->center()[1];
    //double CHECKW=cell->face(f)->center()[2];      
    //    if(CHECKL ==0 ||CHECKL== problem_Length||CHECKW ==0 ||CHECKW== problem_Width||CHECKH== problem_Height) {
      if(cell->face(f)->center()[1] == problemHeight /*&& cell->face(f)->center()[2]==0.5*problem_Width*/) {      
      //evaluate Residual on face
      for (unsigned int i=0; i<dofs_per_cell; ++i) {
	const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
	
	for (unsigned int q=0; q<n_q_points_face; ++q) {
	  if (ck==0 && LiqfaceConv[q]>0.05) {  
	    //R[i] +=-fe_face_values.shape_value_component(i, q, ck)*(dGammadT*Tface_j[q][ck])*fe_face_values.JxW(q);
	  }
	}
      }

    }

  }
  
}

#endif /* CHEMO_H_ */
