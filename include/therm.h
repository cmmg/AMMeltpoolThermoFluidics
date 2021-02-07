//
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Created 2012
//authors: rudraa (2012, 2018)
//

#ifndef THERM_H_
#define THERM_H_
#include "functionEvaluations.h"
#include "supplementaryFunctions.h"
#include <math.h>     

//Chemistry residual implementation
template <int dim>
void residualForTherm(FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1, Sacado::Fad::DFad<double> >& T_ULocal, dealii::Table<1, double>& T_ULocalConv,dealii::Table<1, double>& T_ULocalConvConv, dealii::Table<1, Sacado::Fad::DFad<double> >& ULocal, dealii::Table<1, double>& ULocalConv,dealii::Table<1, double>& ULocalConvConv, dealii::Table<1, Sacado::Fad::DFad<double> >& R, double currentTime, double totalTime) {
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  
  unsigned int n_q_points_face= fe_face_values.n_quadrature_points;
  const unsigned int faces_per_cell = GeometryInfo<dim>::faces_per_cell;

  dealii::Table<1,double> T_conv(n_q_points),Tface_conv(n_q_points);
  dealii::Table<1,double> T_convconv(n_q_points);

  dealii::Table<1,double> liquid_conv(n_q_points),liquidface_conv(n_q_points);
  dealii::Table<1,double> liquid_convconv(n_q_points);
  
  //evaluate gradients 
  dealii::Table<1,Sacado::Fad::DFad<double> > T(n_q_points), Tface(n_q_points_face), liquid(n_q_points); 
  dealii::Table<2,Sacado::Fad::DFad<double> > T_j(n_q_points, dim), liquid_j(n_q_points,dim);
  dealii::Table<2,Sacado::Fad::DFad<double> > vel(n_q_points,dim);

  //Interpolate on all cells 
  for (unsigned int q=0; q<n_q_points; ++q) {
    T[q]=0.0; T_conv[q]=0.0; T_convconv[q]=0.0; liquid[q]=0.0; liquid_conv[q]=0.0; liquid_convconv[q]=0.0;    
    for (unsigned int j=0; j<dim; j++) {T_j[q][j]=0.0; liquid_j[q][j]=0.0; vel[q][j]=0.0; }
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;     
       if (ck>=0 && ck<3) { 
	 vel[q][ck]+=fe_values.shape_value_component(i, q, ck)*ULocalConv[i];
      }
           
       else if (ck==4) { 
	T[q]+=fe_values.shape_value_component(i, q, ck)*T_ULocal[i]; T_conv[q]+=fe_values.shape_value_component(i, q, ck)*T_ULocalConv[i];
	T_convconv[q]+=fe_values.shape_value_component(i, q, ck)*T_ULocalConvConv[i];
	for (unsigned int j=0; j<dim; j++) {
	  T_j[q][j]+=fe_values.shape_grad_component(i, q, ck)[j]*T_ULocal[i];	  
	}
      }
      else if (ck==5) {
	liquid[q]+=fe_values.shape_value_component(i, q,ck)*T_ULocal[i];
	liquid_conv[q]+=fe_values.shape_value_component(i, q,ck)*T_ULocalConv[i];
	liquid_convconv[q]+=fe_values.shape_value_component(i, q,ck)*T_ULocalConvConv[i];
	for (unsigned int j=0; j<dim; j++) {
	  liquid_j[q][j]+=fe_values.shape_grad_component(i, q,ck)[j]*T_ULocal[i];
	}	
      }
    }
    
  }

  //Interpolate over faces
  for (unsigned int f=0; f < faces_per_cell; f++) { 
    fe_face_values.reinit (cell, f);
    //double CHECKL=cell->face(f)->center()[0];
    double CHECKH=cell->face(f)->center()[1];
    //double CHECKW=cell->face(f)->center()[2];      
    //    if(CHECKL ==0 ||CHECKL== problem_Length||CHECKW ==0 ||CHECKW== problem_Width||CHECKH== problem_Height) {
    if(cell->face(f)->center()[1] == problemHeight) {      
      for (unsigned int q=0; q<n_q_points_face; ++q) {
	Tface[q]=0.0; Tface_conv[q]=0.0;   liquid_conv[q]=0;
	for (unsigned int i=0; i<dofs_per_cell; ++i) {
	  const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
	  if (ck==4) {
	    Tface[q]+=fe_face_values.shape_value_component(i, q, ck)*T_ULocal[i]; 
	    Tface_conv[q]+=fe_face_values.shape_value_component(i, q, ck)*T_ULocalConv[i];	    
	  }
	  else if (ck==5) {
	    liquidface_conv[q]+=fe_face_values.shape_value_component(i, q, ck)*T_ULocalConv[i];	    
	  }	  
	}       
      }
    }
    
  }

  //evaluate Residual on cell
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    Sacado::Fad::DFad<double>  KK_T,CC_T;      
    for (unsigned int q=0; q<n_q_points; ++q) {     
      if (ck==4) {		
	KK_T=(1-liquid_conv[q])*KKS+(liquid_conv[q])*KKL;
	CC_T=(1-liquid_conv[q])*CCS+(liquid_conv[q])*CCL;

	//Mass term
	R[i]+=(0.5/dt)*fe_values.shape_value_component(i, q, ck)*(3.0*T[q]-4.0*T_conv[q]+T_convconv[q])*fe_values.JxW(q);

	//Advection term : temperature and phi
	for (unsigned int j = 0; j < dim; j++){
	  R[i] +=(LATENT/CC_T)*fe_values.shape_value_component(i, q,ck)*(vel[q][j])*(liquid_j[q][j])*fe_values.JxW(q);
	  R[i] +=fe_values.shape_value_component(i, q,ck)*(vel[q][j])*(T_j[q][j])*fe_values.JxW(q);
	}

	//sink term latent
	R[i] +=(LATENT/CC_T)*(0.5/dt)*fe_values.shape_value_component(i, q, ck)*(3.0*liquid[q]-4.0*liquid_conv[q]+liquid_convconv[q])*fe_values.JxW(q);
	
	//diffusion terms
	for (unsigned int j = 0; j < dim; j++){	
	  R[i] +=(KK_T/RHO/CC_T)*fe_values.shape_grad_component(i, q,ck)[j]*T_j[q][j]*fe_values.JxW(q);
	}

	//Laser sinl Term
	//Sink term laser
	Point<dim> qPoint=fe_values.quadrature_point(q);  
	Sacado::Fad::DFad<double>  LASER =(ABSORB)*DD*(PP/3.1416/spotRadius/spotRadius/LAYER);
	LASER*=std::exp(-(BB/spotRadius/spotRadius)*((qPoint[0]-VV*currentTime)*(qPoint[0]-VV*currentTime)))  ;
	LASER*=std::exp(-(BB/spotRadius/spotRadius)*((2.0*qPoint[2]-problemWidth)*(2.0*qPoint[2]-problemWidth) ))  ;
	LASER*=std::exp(-(BB/LAYER/LAYER)*((qPoint[1]-problemHeight)*(qPoint[1]-problemHeight) ))  ;
	R[i] +=-(1.0/RHO/CC_T)*fe_values.shape_value_component(i, q, ck)*(LASER)*fe_values.JxW(q);
	
      }

      else if (ck==5) {
	Sacado::Fad::DFad<double> FRACTION;
	FRACTION=std::tanh((5.0/deltaT)*(T[q]-0.5*(TSS+TLL)));
	FRACTION=(1+FRACTION)*0.5;
	R[i] += fe_values.shape_value_component(i, q,ck)*(liquid[q]-FRACTION )*fe_values.JxW(q);	
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
	  if (ck==4) {  
	    //Point<dim> qPoint=fe_face_values.quadrature_point(q);	  
	    Sacado::Fad::DFad<double>  dTRAD= Tface[q]*Tface[q]*Tface[q]*Tface[q] - Tamb*Tamb*Tamb*Tamb;	 	  
	    Sacado::Fad::DFad<double>  KK_T=0,CC_T=0;
	    
	   
	    KK_T=(1-liquidface_conv[q])*KKS+(liquidface_conv[q])*KKL;
	    CC_T=(1-liquidface_conv[q])*CCS+(liquidface_conv[q])*CCL;

	    
	    R[i] += (1.0/RHO/CC_T)*(HH)*fe_face_values.shape_value_component(i, q,ck)*(Tface[q]-Tamb)*fe_face_values.JxW(q);
	    R[i] += (1.0/RHO/CC_T)*(SIG*em)*fe_face_values.shape_value_component(i, q,ck)*(dTRAD)*fe_face_values.JxW(q);
	  
	  }
	}
      }

    }

  }

  
}

#endif /* THERM_H_ */
