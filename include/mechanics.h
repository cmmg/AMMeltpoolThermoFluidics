//
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Created 2011
//authors: rudraa (2011, 2018)
//
#ifndef MECHANICS_H_
#define MECHANICS_H_
#include "functionEvaluations.h"
#include "supplementaryFunctions.h"
#include "deformationMap.h"


//Mechanics implementation
template <class T, int dim>
  void evaluateStress(const FEValues<dim>& fe_values, const unsigned int DOF, const Table<1, T>& M_ULocal, const Table<1, double >& T_ULocalConv, Table<3, T>& P, Table <1,T>& W, const typename DoFHandler<dim>::active_cell_iterator &cell){

  //number of quadrature poits
  unsigned int n_q_points= fe_values.n_quadrature_points;
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;

  // unsigned int n_q_points_face= fe_face_values.n_quadrature_points;
  //const unsigned int faces_per_cell = GeometryInfo<dim>::faces_per_cell;
  
  Table<3, Sacado::Fad::DFad<double>> gradU(n_q_points, dim, dim);
  Table<1,Sacado::Fad::DFad<double>> Temp(n_q_points);
  
   //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q) {
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
	gradU[q][i][j]=0.0;
      }
    }
    for (unsigned int k=0; k<dofs_per_cell; ++k){
      unsigned int ck = fe_values.get_fe().system_to_component_index(k).first - DOF;
      if (ck>=0 && ck<3 /*3*/){
	for (unsigned int i=0; i<dim; ++i){
	  gradU[q][ck][i]+=M_ULocal[k]*fe_values.shape_grad_component(k, q, ck)[i]; //gradU
	}
      }
    }
    
  }              

 
  //loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q) { 
    //Fe
    Table<2, Sacado::Fad::DFad<double> > Fe (dim, dim);
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
	//Fe[i][j]=defMap.F[q][i][j];
      }
    }
    //E
    Table<2, Sacado::Fad::DFad<double> > E (dim, dim);
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){	
	//small strain: E is epsilon
	E[i][j] = 0.5*(gradU[q][i][j]+gradU[q][j][i]);
      }
    }

    
    //S
    //Material moduli
    double Y=elasticModulus, nu=PoissonsRatio;
    //Lame parameters
    double lambda=(nu*Y)/((1+nu)*(1-2*nu)), mu=Y/(2*(1+nu));
    //double lambda=LAM, mu=MU;
    Table<2, Sacado::Fad::DFad<double> > S (dim, dim);
  
    double C11=lambda+2*mu, C12=lambda, C44=mu;    
    //Constitutive relations for 3D
    S[0][0]=C11*E[0][0]+C12*E[1][1]+C12*E[2][2];
    S[1][1]=C12*E[0][0]+C11*E[1][1]+C12*E[2][2];                                                                                          
    S[2][2]=C12*E[0][0]+C12*E[1][1]+C11*E[2][2]; 
    S[0][1]=S[1][0]=C44*E[0][1];
    S[0][2]=S[2][0]=C44*E[0][2];
    S[1][2]=S[2][1]=C44*E[1][2];
    
    //P
    //Calculate W= 0.5*e:C:e '
    
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
	//small strain
	P[q][i][j]=S[i][j];
	W[q]+=0.5*E[i][j]*S[i][j];
      }
    }
    
  }

}

//mechanics residual implementation
template <int dim>
void residualForMechanics(FEValues<dim>& fe_values, unsigned int DOF, typename DoFHandler<dim>::active_cell_iterator &cell, double dt, dealii::Table<1, Sacado::Fad::DFad<double> >& M_ULocal, dealii::Table<1, double>& T_ULocalConv,  dealii::Table<1, Sacado::Fad::DFad<double> >& R, double currentTime, double totalTime) {

  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;

  //unsigned int n_q_points_face= fe_face_values.n_quadrature_points;
  //const unsigned int faces_per_cell = GeometryInfo<dim>::faces_per_cell;
  
  dealii::Table<2,Sacado::Fad::DFad<double> > ux(n_q_points,dim) ;
  dealii::Table<2,Sacado::Fad::DFad<double> > uy(n_q_points,dim) ;
  dealii::Table<2,Sacado::Fad::DFad<double> > uz(n_q_points,dim) ;
  dealii::Table<1,double> TempConv(n_q_points) ;
  
  for (unsigned int q=0; q<n_q_points; ++q) {
    TempConv[q]=0.0;
    for (unsigned int i=0; i<dim; ++i){
	ux[q][i]=0.0;
	uy[q][i]=0.0;
	uz[q][i]=0.0;
    }


     for (unsigned int i=0; i<dofs_per_cell; ++i) {
       const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;       
       
       if (ck==0) {
	 //ux[q]+=fe_values.shape_grad(i, q)[ck]*ULocal[i];	 
	 for (unsigned int j=0; j<dim; ++j) {
	   ux[q][j]+=fe_values.shape_grad_component(i, q, ck)[j]*M_ULocal[i];	 
	 }
       }
       
       if (ck==1) {
	 //uy[q]+=fe_values.shape_grad(i, q)[ck]*ULocal[i];
	 for (unsigned int j=0; j<dim; ++j) {
	   uy[q][j]+=fe_values.shape_grad_component(i, q, ck)[j]*M_ULocal[i];	 
	 }
	 
       }
       if (ck==2/*ck==2*/) {	   
	 for (unsigned int j=0; j<dim; ++j) {
	   uz[q][j]+=fe_values.shape_grad_component(i, q, ck)[j]*M_ULocal[i];	 
	 }	 
       }
       
       if (ck==7) {
	 for (unsigned int j=0; j<dim; ++j) {
	   //TempConv_j[q][j]+=fe_values.shape_grad_component(i, q, ck)[j]*T_ULocalConv[i];	 
	 }
	 TempConv[q]+=fe_values.shape_value_component(i, q, ck)*T_ULocalConv[i];	 
       }
       
     }
  }
  
  
  // double PP=0;
  //if (CURRENT > 0) {PP=PressureMin+(PressureMax-PressureMin)*(currentTime/TotalTime);}
  //temporary arrays
  Table<3,Sacado::Fad::DFad<double> > P (n_q_points, dim, dim);
  //Table<3,Sacado::Fad::DFad<double> > PFace (n_q_points, dim, dim);
  Table<1, Sacado::Fad::DFad<double> > W (n_q_points);
    
  //evaluate stress
  evaluateStress<Sacado::Fad::DFad<double>, dim>(fe_values, DOF, M_ULocal, T_ULocalConv, P, W,/*defMap*/ cell);
  //Material moduli
  double Y=elasticModulus, nu=PoissonsRatio,ALPHA=expCoeff; 
  //Lame parameters
  double lambda=(nu*Y)/((1+nu)*(1-2*nu)), mu=Y/(2*(1+nu));
  double thermalConst=(3.0*lambda+2.0*mu)*ALPHA; 

  //evaluate Residual
  for (unsigned int i=0; i<dofs_per_cell; ++i) {
    const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - DOF;
    
    for (unsigned int q=0; q<n_q_points; ++q){
      if (ck>=0 && ck<3 /*3*/) {  
	for (unsigned int d = 0; d < dim; d++){
	  R[i] +=-fe_values.shape_grad_component(i, q, ck)[d]*P[q][ck][d]*fe_values.JxW(q);	  
	}
	
	if (ck==1) R[i]+=-fe_values.shape_value_component(i, q, ck)*(gravity*RHO)*fe_values.JxW(q);
	

	R[i] +=(thermalConst)*fe_values.shape_grad_component(i, q, ck)[ck]*(TempConv[q]-Tamb)*fe_values.JxW(q);	  
	//if (R[i].val()!=0) std::cout <<"value of R is non zero " << R[i].val()<<std::endl;
      }               
    }
       
  }
  


}

#endif /* MECHANICS_H_ */
