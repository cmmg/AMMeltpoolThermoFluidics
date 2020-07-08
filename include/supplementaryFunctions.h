//
//Computational Mechanics and Multiphysics Group @ UW-Madison
//Created 2011
//authors: rudraa (2011, 2018)
//

#ifndef SUPPLEMENTARYFUNCTIONS_H_
#define SUPPLEMENTARYFUNCTIONS_H_
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/numerics/vector_tools.h>

void solveSystem(SparseMatrix<double>& system_matrix, Vector<double>& system_rhs, Vector<double>& dU, double tolerance){
    //Iterative solve
  dU=0; SolverControl solver_control1 (2000, tolerance, true, true); solver_control1.set_failure_criterion(1.0e8); solver_control1.log_frequency(100); solver_control1.log_result(true);
  try{
    //BICG, Jacobi
    printf("BICG, Jacobi:  ");
    SolverGMRES<> solver (solver_control1); PreconditionJacobi<> preconditioner; preconditioner.initialize(system_matrix, 1.0); solver.solve (system_matrix, dU, system_rhs, preconditioner);
    printf("iterative solve complete (Steps:%4u, Tol:%11.4e, InitialRes:%11.4e, Res:%11.4e).\n", solver_control1.last_step(), tolerance,  solver_control1.initial_value(), solver_control1.last_value());
  }
  catch(...){
    dU=0; SolverControl solver_control2 (2000, tolerance, true, true); solver_control2.set_failure_criterion(1.0e8); solver_control2.log_frequency(100); solver_control2.log_result(true);
    try{
      //GMRES, Jacobi
      printf("failed (Steps:%4u, Res:%11.4e) \nGMRES, Jacobi: ", solver_control1.last_step(), solver_control1.last_value());
      SolverBicgstab<> solver (solver_control2); PreconditionJacobi<> preconditioner; preconditioner.initialize(system_matrix, 1.0); solver.solve (system_matrix, dU, system_rhs, preconditioner);
      printf("iterative solve complete (Steps:%4u, Tol:%11.4e, InitialRes:%11.4e, Res:%11.4e).\n", solver_control2.last_step(), tolerance, solver_control2.initial_value(), solver_control2.last_value());
    }
    catch(...){
      dU=0; SolverControl solver_control3 (2000, tolerance, true, true); solver_control3.set_failure_criterion(1.0e8); solver_control3.log_frequency(100); solver_control3.log_result(true);
      try{
	//BiCG, SOR
	printf("failed (Steps:%4u, Res:%11.4e) \nBiCG, SOR:     ", solver_control2.last_step(), solver_control2.last_value());
	SolverBicgstab<> solver (solver_control3); PreconditionSOR<> preconditioner; preconditioner.initialize(system_matrix, 1.0); solver.solve (system_matrix, dU, system_rhs, preconditioner);
	printf("iterative solve complete (Steps:%4u, Tol:%11.4e, InitialRes:%11.4e, Res:%11.4e).\n", solver_control3.last_step(), tolerance, solver_control3.initial_value(), solver_control3.last_value());
      }
      catch(...){
	dU=0; SolverControl solver_control4 (2000, tolerance, true, true); solver_control4.set_failure_criterion(1.0e8); solver_control4.log_frequency(100); solver_control4.log_result(true);
	try{
	  //GMRES, SOR
	  printf("failed (Steps:%4u, Res:%11.4e) \nGMRES, SOR:    ", solver_control3.last_step(), solver_control3.last_value());
	  SolverGMRES<> solver (solver_control4); PreconditionSOR<> preconditioner; preconditioner.initialize(system_matrix, 1.0); solver.solve (system_matrix, dU, system_rhs, preconditioner);
	  printf("iterative solve complete (Steps:%4u, Tol:%11.4e, InitialRes:%11.4e, Res:%11.4e).\n", solver_control4.last_step(), tolerance, solver_control4.initial_value(), solver_control4.last_value());
	}
	catch(...){
	  printf("failed (Steps:%4u, Res:%11.4e) \nDirect Solve: ", solver_control4.last_step(), solver_control4.last_value());
	  dU=0; SparseDirectUMFPACK  A_direct; A_direct.initialize(system_matrix);  A_direct.vmult (dU, system_rhs);
	  printf("direct solve complete.\n");
	}
      }
    }
  }
}

#endif /* SUPPLEMENTARYFUNCTIONS_H_ */

