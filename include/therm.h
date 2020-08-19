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
void residualForTherm(FEValues<dim>& fe_values, unsigned int DOF, FEFaceValues<dim>& fe_face_values, const typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1, Sacado::Fad::DFad<double> >& ULocal, dealii::Table<1, double>& ULocalConv, dealii::Table<1, Sacado::Fad::DFad<double> >& R, double currentTime, double totalTime) {
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  
  unsigned int n_q_points_face= fe_face_values.n_quadrature_points;
  const unsigned int faces_per_cell = GeometryInfo<dim>::faces_per_cell;

  dealii::Table<1,double> c_conv(n_q_points);
  //evaluate gradients 
  dealii::Table<1,Sacado::Fad::DFad<double> > c(n_q_points), cface(n_q_points_face), liquid(n_q_points),liquid_conv(n_q_points); 
  dealii::Table<2,Sacado::Fad::DFad<double> > c_j(n_q_points, dim), liquid_j(n_q_points,dim);
  dealii::Table<2,Sacado::Fad::DFad<double> > vel(n_q_points,dim);

  //Interpolate on all cells 
  for (unsigned int q=0; q<n_q_points; ++q) {
    c[q]=0.0; c_conv[q]=0.0; liquid[q]=0.0; liquid_conv[q]=0.0;    
    for (unsigned int j=0; j<dim; j++) {c_j[q][j]=0.0; vel[q][j]=0; }
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
      const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;     
       if (ck>=0 && ck< 2) { 
	 vel[q][ck]+=fe_values.shape_value_component(i, q, ck)*ULocal[i];
      }
           
      if (ck==3) { 
	c[q]+=fe_values.shape_value(i, q)*ULocal[i]; c_conv[q]+=fe_values.shape_value(i, q)*ULocalConv[i];
	for (unsigned int j=0; j<dim; j++) {
	  c_j[q][j]+=fe_values.shape_grad(i, q)[j]*ULocal[i];
	}
      }
      else if (ck==4) {
	liquid[q]+=fe_values.shape_value(i, q)*ULocal[i]; liquid_conv[q]+=fe_values.shape_value(i, q)*ULocalConv[i];
	for (unsigned int j=0; j<dim; j++) {
	  liquid_j[q][j]+=fe_values.shape_grad(i, q)[j]*ULocal[i];
	}	
      }
    }
    
  }

  //Interpolate over faces
  for (unsigned int f=0; f < faces_per_cell; f++) { 
    fe_face_values.reinit (cell, f);
    double CHECKL=cell->face(f)->center()[0];
    double CHECKH=cell->face(f)->center()[1];
    double CHECKW=cell->face(f)->center()[2];      
    //    if(CHECKL ==0 ||CHECKL== problem_Length||CHECKW ==0 ||CHECKW== problem_Width||CHECKH== problem_Height) {
    if(cell->face(f)->center()[1] == problemHeight) {
      
      for (unsigned int q=0; q<n_q_points_face; ++q) {
	cface[q]=0.0;    
	for (unsigned int i=0; i<dofs_per_cell; ++i) {
	  const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
	  if (ck==3) {cface[q]+=fe_face_values.shape_value(i, q)*ULocal[i]; }
	}
        
      }


    }
    
  }

  //evaluate Residual on cell
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    Sacado::Fad::DFad<double>  KK_T,CC_T;      
    for (unsigned int q=0; q<n_q_points; ++q) {     
      if (ck==3) {		
	R[i] += (1.0/dt)*fe_values.shape_value(i, q)*(c[q]-c_conv[q])*fe_values.JxW(q);  
	//R[i] += (VV)*fe_values.shape_value(i, q)*c_j[q][0]*fe_values.JxW(q);
	if (c_conv[q]<=TLL && c_conv[q]>=0.0)	  {
	    //KK_T=1.57 + (1.6*pow(10.0,-2.0)*c_conv[q])- (pow(10.0,-6)*c_conv[q]*c_conv[q]);
	    //CC_T= 492.4 + (0.025*c_conv[q]) -(4.18*pow(10.0,-6)*c_conv[q]*c_conv[q]) ;	    	    
	    //SS316
	    KK_T=11.82+(1.06*pow(10.0,-2))*c_conv[q];
	    CC_T=330.9+0.563*c_conv[q]-(4.015*pow(10.0,-4))*c_conv[q]*c_conv[q]+ (9.465*pow(10.0,-8))*c_conv[q]*c_conv[q]*c_conv[q];
	}
	else if (c_conv[q]>=TLL) {KK_T=KK ; CC_T=CC;}
	//sink term latent
	R[i] +=(LATENT/CC_T)*(1.0/dt)*fe_values.shape_value(i, q)*(liquid[q]-liquid_conv[q])*fe_values.JxW(q);
	
	for (unsigned int j = 0; j < dim; j++){
	  R[i] +=(LATENT/CC_T)*fe_values.shape_value(i, q)*(vel[q][j])*(liquid_j[q][j])*fe_values.JxW(q);
	}

	//Sink term laser
	Point<dim> qPoint=fe_values.quadrature_point(q);	  
	Sacado::Fad::DFad<double>  LASER =(ABSORB)*DD*(PP/3.1416/spotRadius/spotRadius/LAYER);
	LASER*=std::exp(-(BB/spotRadius/spotRadius)*((qPoint[0]-VV*currentTime)*(qPoint[0]-VV*currentTime)))  ;
	//LASER*=std::exp(-(BB/spotRadius/spotRadius)*((2.0*qPoint[2]-problem_Width)*(2.0*qPoint[2]-problem_Width) ))  ;
	LASER*=std::exp(-(BB/LAYER/LAYER)*((qPoint[1]-problemHeight)*(qPoint[1]-problemHeight) ))  ;	
	R[i] +=-(1.0/RHO/CC_T)*fe_values.shape_value(i, q)*(LASER)*fe_values.JxW(q);		
	for (unsigned int j = 0; j < dim; j++){
	  R[i] += (KK_T/RHO/CC_T)*fe_values.shape_grad(i, q)[j]*c_j[q][j]*fe_values.JxW(q);
	}
	
      }

      else if (ck==4) {
	Sacado::Fad::DFad<double> FRACTION;
	FRACTION=std::tanh(2.5*(c[q]-0.5*(TLL+TSS))/(TLL-TSS)/2.0) ;
	FRACTION=(1+FRACTION)*0.5;
	R[i] += fe_values.shape_value(i, q)*(liquid[q]-FRACTION )*fe_values.JxW(q);	
      }
     
    }
  }


  //surface integral
  for (unsigned int f=0; f < faces_per_cell; f++) { 
    fe_face_values.reinit (cell, f);
    double CHECKL=cell->face(f)->center()[0];
    double CHECKH=cell->face(f)->center()[1];
    double CHECKW=cell->face(f)->center()[2];      
    if(CHECKL ==0 ||CHECKL== problem_Length||CHECKW ==0 ||CHECKW== problem_Width||CHECKH== problem_Height) {
    //  if(cell->face(f)->center()[1] == problem_Height /*&& cell->face(f)->center()[2]==0.5*problem_Width*/) {      
      //evaluate Residual on face
      for (unsigned int i=0; i<dofs_per_cell; ++i) {
	const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
	
	for (unsigned int q=0; q<n_q_points_face; ++q) {
	  if (ck==3) {  
	  Point<dim> qPoint=fe_face_values.quadrature_point(q);	  
	  Sacado::Fad::DFad<double>  Tambient =Tamb;
	  Sacado::Fad::DFad<double>  dTRAD= cface[q]*cface[q]*cface[q]*cface[q] - Tambient*Tambient*Tambient*Tambient;
	  //Sacado::Fad::DFad<double>  LASER = DD*(PP/3.1416/spot/spot)*std::exp(-BB*(qPoint[0]-VV*currentTime)*(qPoint[0]-VV*currentTime)/spot/spot) ;
	  // Sacado::Fad::DFad<double>  LASER = DD*(PP/3.1416/spot/spot)*std::exp(-(BB/spot/spot)*((qPoint[0]-VV*currentTime)*(qPoint[0]-VV*currentTime) + (2.0*qPoint[2]-problem_Width)*(2.0*qPoint[2]-problem_Width) ))  ;
	  
	  
	  Sacado::Fad::DFad<double>  KK_T,CC_T;
	  if (c_conv[q]<=TLL && c_conv[q]>=0.0)
	  {
	    //   KK_T=1.57 + (1.6*pow(10.0,-2.0)*c_conv[q])- (pow(10.0,-6)*c_conv[q]*c_conv[q]);
	    //CC_T= 492.4 + (0.025*c_conv[q]) -(4.18*pow(10.0,-6)*c_conv[q]*c_conv[q]) ;

	    //SS316
	    KK_T=11.82+(1.06*pow(10.0,-2))*c_conv[q];
	    CC_T=330.9+0.563*c_conv[q]-(4.015*pow(10.0,-4))*c_conv[q]*c_conv[q]+ (9.465*pow(10.0,-8))*c_conv[q]*c_conv[q]*c_conv[q];
	   
	    

	  }
	else if (c_conv[q]>=TLL) {KK_T=KK ; CC_T=CC;}
	  R[i] += (1.0/RHO/CC_T)*(HH)*fe_face_values.shape_value(i, q)*(cface[q]-Tambient)*fe_face_values.JxW(q);
	  R[i] += (1.0/RHO/CC_T)*(SIG*em)*fe_face_values.shape_value(i, q)*(dTRAD)*fe_face_values.JxW(q);

	  //  R[i] += -(1.0/RHO/CC_T)*fe_face_values.shape_value(i, q)*(LASER)*fe_face_values.JxW(q);
	  
	  }
	}
      }

    }

  }

  
}

#endif /* THERM_H_ */
