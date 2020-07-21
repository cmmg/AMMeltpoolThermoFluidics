//Computational Mechanics and Multiphysics Group @ UW-Madison
//Basic framework for Phase Field (Cahn-Hilliard mixed formulation)
//Created May 2018
//authors: rudraa (2018)
//

//deal.II headers
#include "include/headers.h"
//input parameter headers
#include "parameters.h"
//physics headers
#include "include/chemo.h"
#include "include/projection.h"
//Namespace
namespace phaseField1
{
  using namespace dealii;

  //Initial conditions
  template <int dim>
  class InitalConditions: public Function<dim>{
  public:
    InitalConditions (): Function<dim>(DIMS){ }
    void vector_value (const Point<dim>   &p, Vector<double>   &values) const{
        Assert (values.size() == DIMS, ExcDimensionMismatch (values.size(), DIMS));     	
	values(0)=0.0;
	values(1)=0.0;
	values(2)=25-p[0];
	values(3)=0;
    }
    
  };

                                    
  template <int dim>
  class phaseField{
  public:
    phaseField ();
    ~phaseField ();
    void run ();

  private:
    void applyBoundaryConditions(const unsigned int increment);
    void setup_system ();
    void assemble_system ();
    void solveIteration (const bool lever);
    void solve (const bool lever);
    void refine_grid ();
    void output_results (const unsigned int increment);
    MPI_Comm                                  mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;
    FESystem<dim>                             fe;
    DoFHandler<dim>                           dof_handler;
    IndexSet                                  locally_owned_dofs;
    IndexSet                                  locally_relevant_dofs;
    ConstraintMatrix                          constraints, constraintsZero;
    LA::MPI::SparseMatrix                     system_matrix;
    LA::MPI::Vector                           locally_relevant_solution, U, Un, UGhost, UnGhost, dU;
    LA::MPI::Vector                           Unn, UnnGhost; 
    LA::MPI::Vector                           system_rhs;
    ConditionalOStream                        pcout;
    TimerOutput                               computing_timer;

    void update_pressure ();
    void applyBoundaryConditions_projection(const unsigned int increment);
    void setup_system_projection ();
    void assemble_system_projection();

    ConstraintMatrix                          Pr_constraints, Pr_constraintsZero;
    LA::MPI::SparseMatrix                     Pr_system_matrix;
    LA::MPI::Vector                           Pr_U, Pr_Un, Pr_UGhost, Pr_UnGhost, Pr_dU;
    LA::MPI::Vector                           Pr_system_rhs;
    
    //solution variables
    unsigned int currentIncrement, currentIteration;
    double totalTime, currentTime, dt;
    std::vector<std::string> nodal_solution_names; std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;
  };

  template <int dim>
  phaseField<dim>::phaseField ():
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator,
                   typename Triangulation<dim>::MeshSmoothing
                   (Triangulation<dim>::smoothing_on_refinement |
                    Triangulation<dim>::smoothing_on_coarsening)),
    fe(FE_Q<dim>(2),2,FE_Q<dim>(1),2),
    dof_handler(triangulation),
    pcout (std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator)== 0)),
    computing_timer (mpi_communicator, pcout, TimerOutput::summary, TimerOutput::wall_times){
    //solution variables
    dt=TimeStep; totalTime=TotalTime;
    currentIncrement=0; currentTime=0;

    //nodal Solution names
    //nodal_solution_names.push_back("ux"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
    //nodal_solution_names.push_back("uy"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);

   for (unsigned int i=0; i<dim; ++i){
      nodal_solution_names.push_back("velocity"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
    }

   
   nodal_solution_names.push_back("pressure"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
   nodal_solution_names.push_back("PHI"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  }
  
  template <int dim>
  phaseField<dim>::~phaseField (){
    dof_handler.clear ();   
  }

  //Apply boundary conditions for diffusion step
  template <int dim>
  void phaseField<dim>::applyBoundaryConditions(const unsigned int increment){
    constraints.clear (); 
    constraints.reinit (locally_relevant_dofs);
    //    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    constraintsZero.clear (); 
    constraintsZero.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, constraintsZero);
    
    //Setup boundary conditions
    //No Dirchlet BC are necessary for the parabolic problem
    std::vector<bool> uBCI (DIMS, false); uBCI[0]=true;uBCI[1]=true;
    std::vector<bool> uBCW (DIMS, false); uBCW[0]=true;uBCW[1]=true;
    std::vector<bool> uBCO (DIMS, false); uBCO[0]=true; uBCO[1]=true; 
    
    // 1 : walls top and bowttom , 2 : inlet 3: outlet 4: cavity walls
    
    //top wall and bottom wall
    VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(DIMS), constraints,uBCW);
    VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(DIMS), constraintsZero,uBCW);

    //cavity walls
    VectorTools::interpolate_boundary_values (dof_handler, 4, ZeroFunction<dim>(DIMS) , constraints,uBCW);
    VectorTools::interpolate_boundary_values (dof_handler, 4, ZeroFunction<dim>(DIMS) , constraintsZero,uBCW);

    //outlet
    //VectorTools::interpolate_boundary_values (dof_handler, 3, ZeroFunction<dim>(DIMS) , constraints,uBCO);
    //VectorTools::interpolate_boundary_values (dof_handler, 3, ZeroFunction<dim>(DIMS) , constraintsZero,uBCO);
    

    //inlet
    std::vector<double> inletValue (DIMS);
    inletValue[0]=1.0;
    inletValue[1]=0.0;
    inletValue[2]=0.0;
    inletValue[3]=0.0;
    VectorTools::interpolate_boundary_values (dof_handler, 2, ConstantFunction<dim>(inletValue), constraints,uBCI);
    VectorTools::interpolate_boundary_values (dof_handler, 2, ZeroFunction<dim>(DIMS) , constraintsZero,uBCI);
    
    constraints.close ();
    constraintsZero.close ();
  }

  template <int dim>
  void phaseField<dim>::applyBoundaryConditions_projection(const unsigned int increment){
    Pr_constraints.clear (); 
    Pr_constraints.reinit (locally_relevant_dofs);
    //    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    Pr_constraintsZero.clear (); 
    Pr_constraintsZero.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dof_handler, Pr_constraintsZero);
    
    //Setup boundary conditions
    //No Dirchlet BC are necessary for the parabolic problem
    std::vector<bool> uBC (DIMS, false); uBC[3]=true; //uBC[2]=true;
    
    // 1 : walls top and bowttom , 2 : inlet 3: outlet 4: cavity walls
   
    //outlet
    VectorTools::interpolate_boundary_values (dof_handler, 3, ZeroFunction<dim>(DIMS) , Pr_constraints,uBC);
    VectorTools::interpolate_boundary_values (dof_handler, 3, ZeroFunction<dim>(DIMS) , Pr_constraintsZero,uBC);       

    //    VectorTools::interpolate_boundary_values (dof_handler, 4, ZeroFunction<dim>(DIMS) , Pr_constraints,uBC);
    //VectorTools::interpolate_boundary_values (dof_handler, 4, ZeroFunction<dim>(DIMS) , Pr_constraintsZero,uBC);       

    //VectorTools::interpolate_boundary_values (dof_handler, 2, ZeroFunction<dim>(DIMS) , Pr_constraints,uBC);
    //VectorTools::interpolate_boundary_values (dof_handler, 2, ZeroFunction<dim>(DIMS) , Pr_constraintsZero,uBC);       

    //VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(DIMS) , Pr_constraints,uBC);
    //VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(DIMS) , Pr_constraintsZero,uBC);       

    Pr_constraints.close ();
    Pr_constraintsZero.close ();
  }

  
  
  //Setup
  template <int dim>
  void phaseField<dim>::setup_system (){
    TimerOutput::Scope t(computing_timer, "setup");
    dof_handler.distribute_dofs (fe);
    
    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler,
                                             locally_relevant_dofs);
    
    locally_relevant_solution.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

    //Non-ghost vectors
    system_rhs.reinit (locally_owned_dofs, mpi_communicator);
    U.reinit (locally_owned_dofs, mpi_communicator);
    Un.reinit (locally_owned_dofs, mpi_communicator);
    dU.reinit (locally_owned_dofs, mpi_communicator);
    
    Unn.reinit (locally_owned_dofs, mpi_communicator);
    
    //Ghost vectors
    UGhost.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    UnGhost.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

    UnnGhost.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

    UnGhost.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);


    //call applyBoundaryConditions to setup constraints matrix needed for generating the sparsity pattern
    applyBoundaryConditions(0);
    
    DynamicSparsityPattern dsp (locally_relevant_dofs);
    DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern (dsp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
    system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);

    //log info to screen
    pcout << "   Number of active cells:       "
	  << triangulation.n_global_active_cells()
	  << std::endl
	  << "   Number of degrees of freedom: "
	  << dof_handler.n_dofs()
	  << std::endl;
  }


  //Setup
  template <int dim>
  void phaseField<dim>::setup_system_projection (){
    TimerOutput::Scope t(computing_timer, "setup");
    
    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler,
                                             locally_relevant_dofs);
    
    locally_relevant_solution.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

    //Non-ghost vectors
    Pr_system_rhs.reinit (locally_owned_dofs, mpi_communicator);
    Pr_U.reinit (locally_owned_dofs, mpi_communicator);
    Pr_Un.reinit (locally_owned_dofs, mpi_communicator);
    Pr_dU.reinit (locally_owned_dofs, mpi_communicator);
   
    //Ghost vectors
    Pr_UGhost.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    Pr_UnGhost.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
     
    //call applyBoundaryConditions to setup constraints matrix needed for generating the sparsity pattern
    applyBoundaryConditions_projection(0);
    
    DynamicSparsityPattern Pr_dsp (locally_relevant_dofs);
    DoFTools::make_sparsity_pattern (dof_handler, Pr_dsp, Pr_constraints, false);
    SparsityTools::distribute_sparsity_pattern (Pr_dsp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
    Pr_system_matrix.reinit (locally_owned_dofs, locally_owned_dofs, Pr_dsp, mpi_communicator);
   
    
  }

  

  

  //Assembly
  template <int dim>
  void phaseField<dim>::assemble_system (){
    TimerOutput::Scope t(computing_timer, "assembly");
    system_rhs=0.0; system_matrix=0.0;
    const QGauss<dim>  quadrature_formula(FEOrder+2);
    const QGauss<dim-1>	face_quadrature_formula (FEOrder+1);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points |
                             update_JxW_values);
    FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors);
    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       local_rhs (dofs_per_cell);
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    unsigned int n_q_points= fe_values.n_quadrature_points;

    unsigned int n_q_points_face= fe_face_values.n_quadrature_points;
    const unsigned int faces_per_cell = GeometryInfo<dim>::faces_per_cell;
  
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned()){
	fe_values.reinit (cell);
	local_matrix = 0; local_rhs = 0; 
	cell->get_dof_indices (local_dof_indices);
	 //AD variables
	Table<1, Sacado::Fad::DFad<double> > ULocal(dofs_per_cell); Table<1, double > ULocalConv(dofs_per_cell);

	Table<1, double > ULocalConvConv(dofs_per_cell);
	
	for (unsigned int i=0; i<dofs_per_cell; ++i){
	  const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - 0;
	  if (std::abs(UGhost(local_dof_indices[i]))<1.0e-16){ULocal[i]=0.0;}
	  else{ULocal[i]=UGhost(local_dof_indices[i]);}
	  ULocal[i].diff (i, dofs_per_cell);
	  ULocalConv[i]= UnGhost(local_dof_indices[i]);	 
	  ULocalConvConv[i]= UnnGhost(local_dof_indices[i]);
	}

	//setup residual vector
	Table<1, Sacado::Fad::DFad<double> > R(dofs_per_cell); 
	for (unsigned int i=0; i<dofs_per_cell; ++i) {R[i]=0.0;}
	
	 
	//populate residual vctor 
	//residualForChemo(fe_values, 0, fe_face_values, cell, dt, ULocal, ULocalConv, R, currentTime, totalTime);
	residualForChemo(fe_values, 0, fe_face_values, cell, dt, ULocal, ULocalConv, ULocalConvConv, R, currentTime, totalTime);
	       		
	//evaluate Residual(R) and Jacobian(R')
	for (unsigned int i=0; i<dofs_per_cell; ++i) {
	  for (unsigned int j=0; j<dofs_per_cell; ++j){
	    // R' by AD
	    //local_matrix(i,j)= R[i].fastAccessDx(j);
	    local_matrix(i,j)= R[i].dx(j);
	  }
	  //R
	  local_rhs(i) = -R[i].val();
	}
	if ((currentIteration==0)&&(currentIncrement==1)){
	  constraints.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, system_matrix, system_rhs);	
	}
	else{
	  constraintsZero.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, system_matrix, system_rhs);
	}
      }
    system_matrix.compress (VectorOperation::add);
    system_rhs.compress (VectorOperation::add);
    if(isnan(system_matrix.frobenius_norm())) pcout << "\nK norm diffusion has a nan: " << system_matrix.frobenius_norm() << std::endl; 
  }


   //Assembly
  template <int dim>
  void phaseField<dim>::assemble_system_projection (){
    TimerOutput::Scope t(computing_timer, "assembly");
    Pr_system_rhs=0.0; Pr_system_matrix=0.0;
    const QGauss<dim>  quadrature_formula(FEOrder+1);
    const QGauss<dim-1>	face_quadrature_formula (FEOrder);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points |
                             update_JxW_values);
    FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, update_values | update_quadrature_points | update_JxW_values | update_normal_vectors);
    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       local_rhs (dofs_per_cell);
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    unsigned int n_q_points= fe_values.n_quadrature_points;

    unsigned int n_q_points_face= fe_face_values.n_quadrature_points;
    const unsigned int faces_per_cell = GeometryInfo<dim>::faces_per_cell;
  
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned()){
	fe_values.reinit (cell);
	local_matrix = 0; local_rhs = 0; 
	cell->get_dof_indices (local_dof_indices);
	 //AD variables
	Table<1, Sacado::Fad::DFad<double> > ULocal(dofs_per_cell); Table<1, double > ULocalConv(dofs_per_cell);
	Table<1, Sacado::Fad::DFad<double> > Pr_ULocal(dofs_per_cell); Table<1, double > Pr_ULocalConv(dofs_per_cell);
	
	for (unsigned int i=0; i<dofs_per_cell; ++i){
	  const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - 0;
	  if (std::abs(Pr_UGhost(local_dof_indices[i]))<1.0e-16){Pr_ULocal[i]=0.0;}	  
	  else{Pr_ULocal[i]=Pr_UGhost(local_dof_indices[i]);}
	  Pr_ULocal[i].diff (i, dofs_per_cell);
	  Pr_ULocalConv[i]= Pr_UnGhost(local_dof_indices[i]);
	  
	  //if (std::abs(UGhost(local_dof_indices[i]))<1.0e-16){ULocal[i]=0.0;}
	  //else{ULocal[i]=UGhost(local_dof_indices[i]);}
	  //ULocal[i].diff (i, dofs_per_cell);
	  ULocalConv[i]= UnGhost(local_dof_indices[i]);	  
	}


	//setup residual vector
	Table<1, Sacado::Fad::DFad<double> > Rp(dofs_per_cell); 
	for (unsigned int i=0; i<dofs_per_cell; ++i) {Rp[i]=0.0;}
		 
	//populate residual vctor 
	//residualForChemo(fe_values, 0, fe_face_values, cell, dt, ULocal, ULocalConv, R, currentTime, totalTime);
	residualForProjection(fe_values, 0, fe_face_values, cell, dt, /*ULocal,*/  Pr_ULocal, ULocalConv, Pr_ULocalConv, Rp, currentTime, totalTime);
	       		
	//evaluate Residual(R) and Jacobian(R')
	for (unsigned int i=0; i<dofs_per_cell; ++i) {
	  for (unsigned int j=0; j<dofs_per_cell; ++j){
	    // R' by AD
	    //local_matrix(i,j)= R[i].fastAccessDx(j);
	    local_matrix(i,j)= Rp[i].dx(j);
	  }
	  //R
	  local_rhs(i) = -Rp[i].val();
	}
	if ((currentIteration==0)&&(currentIncrement==1)){
	  Pr_constraints.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, Pr_system_matrix, Pr_system_rhs);	
	}
	else{
	  Pr_constraintsZero.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, Pr_system_matrix, Pr_system_rhs);
	}
      }
    Pr_system_matrix.compress (VectorOperation::add);
    Pr_system_rhs.compress (VectorOperation::add);
    if(isnan(Pr_system_matrix.frobenius_norm())) pcout << "\nK norm projection has a nan: " << Pr_system_matrix.frobenius_norm() << std::endl;
  }

 
  
  //Solve
  template <int dim>
  void phaseField<dim>::solveIteration(const bool lever){
    TimerOutput::Scope t(computing_timer, "solve");
    LA::MPI::Vector completely_distributed_solution (locally_owned_dofs, mpi_communicator);    
      if(lever) {				
	//Iterative solvers from Petsc and Trilinos
    SolverControl solver_control (dof_handler.n_dofs(), 1e-12);
#ifdef USE_PETSC_LA
    LA::SolverGMRES solver(solver_control, mpi_communicator);
#else
    LA::SolverGMRES solver(solver_control);
#endif
    LA::MPI::PreconditionAMG preconditioner;
    LA::MPI::PreconditionAMG::AdditionalData data;
#ifdef USE_PETSC_LA //dof_handler.get_fe().n_components()
    //data.symmetric_operator = true;
#else
    // Trilinos defaults are good 
#endif
    //if lever is true run for diffusion part 
    preconditioner.initialize(system_matrix, data);
    solver.solve (system_matrix, completely_distributed_solution, system_rhs, preconditioner);
    pcout << "   Solved in " << solver_control.last_step()
          << " iterations." << std::endl;
	
    
	//Direct solver MUMPS
    //	SolverControl cn;
    //	PETScWrappers::SparseDirectMUMPS solver(cn, mpi_communicator);
    //	solver.set_symmetric_mode(false);
    //	solver.solve(system_matrix, completely_distributed_solution, system_rhs);
    
  	if ((currentIteration==0)&&(currentIncrement==1)){
	  constraints.distribute (completely_distributed_solution);
	}
	else{
	  constraintsZero.distribute (completely_distributed_solution);
	}
	locally_relevant_solution = completely_distributed_solution;
        dU = completely_distributed_solution; 	
      }
      
  
    //if lever is false run for projection parr   
      if (!lever) {
	//Iterative solvers from Petsc and Trilinos
	SolverControl solver_control (dof_handler.n_dofs(), 1e-12);
#ifdef USE_PETSC_LA
	LA::SolverGMRES solver(solver_control, mpi_communicator);
#else
    LA::SolverGMRES solver(solver_control);
#endif
    LA::MPI::PreconditionAMG preconditioner;
    LA::MPI::PreconditionAMG::AdditionalData data;
#ifdef USE_PETSC_LA //dof_handler.get_fe().n_components()
    //data.symmetric_operator = true;
#else
    // Trilinos defaults are good 
#endif
    //if lever is true run for diffusion part 
    preconditioner.initialize(system_matrix, data);
    solver.solve (Pr_system_matrix, completely_distributed_solution, Pr_system_rhs, preconditioner);
    pcout << "   Solved in " << solver_control.last_step()
          << " iterations." << std::endl;     
	//Direct solver MUMPS	
	//SolverControl cn;
	//PETScWrappers::SparseDirectMUMPS solver(cn, mpi_communicator);
	//solver.set_symmetric_mode(false);
	//solver.solve(Pr_system_matrix, completely_distributed_solution, Pr_system_rhs);   
	if ((currentIteration==0)&&(currentIncrement==1)){
	  Pr_constraints.distribute (completely_distributed_solution);
	}
	else{
	    Pr_constraintsZero.distribute (completely_distributed_solution);
	}
	locally_relevant_solution = completely_distributed_solution;
	Pr_dU=completely_distributed_solution;
      }
              
  }   
 
  
  /*    
 //Solve iteration
  template <int dim>
  void phaseField<dim>::solveIteration (const bool lever)
  {
    TimerOutput::Scope t(computing_timer, "solve");
    PETScWrappers::MPI::Vector
      completely_distributed_solution (locally_owned_dofs,mpi_communicator);

    PETScWrappers::PreconditionSSOR preconditioner(Pr_system_matrix), Pr_preconditioner(system_matrix);
     
    if (lever) {
    //distributed_incremental_displacement = incremental_displacement;
    SolverControl           solver_control (dof_handler.n_dofs(),
                                            1e-16*system_rhs.l2_norm());
    PETScWrappers::SolverGMRES cg (solver_control,
                                mpi_communicator);
    // PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix);
    //PETScWrappers::PreconditionILU preconditioner(Pr_system_matrix);
    cg.solve (system_matrix,  completely_distributed_solution, system_rhs,
              preconditioner);

    if ((currentIteration==0)&&(currentIncrement==1)){
      constraintsZero.distribute (completely_distributed_solution);
    }
    else{
      constraints.distribute (completely_distributed_solution);
    }    
    locally_relevant_solution = completely_distributed_solution;
    dU = completely_distributed_solution;
    pcout << "   Solved in " << solver_control.last_step()
	  << " iterations." << std::endl;
    }
    
    if (!lever) {
          //distributed_incremental_displacement = incremental_displacement;
    SolverControl           solver_control (dof_handler.n_dofs(),
                                            1e-16*Pr_system_rhs.l2_norm());
    PETScWrappers::SolverGMRES cg (solver_control,
                                mpi_communicator);
    //PETScWrappers::PreconditionBlockJacobi preconditioner(Pr_system_matrix);
    //PETScWrappers::PreconditionILU preconditioner(Pr_system_matrix);
    cg.solve (Pr_system_matrix,  completely_distributed_solution, Pr_system_rhs,
              Pr_preconditioner);

    if ((currentIteration==0)&&(currentIncrement==1)){
      constraintsZero.distribute (completely_distributed_solution);
    }
    else{
      constraints.distribute (completely_distributed_solution);
    }    
    locally_relevant_solution = completely_distributed_solution;
    Pr_dU = completely_distributed_solution;
    pcout << "   Solved in " << solver_control.last_step()
	  << " iterations." << std::endl;
    
    }
     
  }
  */
  

  
  //Solve
  template <int dim>
  void phaseField<dim>::solve(const bool lever) {
    double res=1, tol=1.0e-12, abs_tol=1.0e-14, initial_norm=0, current_norm=0;
    double machineEPS=1.0e-15;
    currentIteration=0;
    char buffer[200];    
    
    if(lever) {
      
      pcout << "Solving for diffusion "<<std::endl;
    while (true){
      if (currentIteration>=10){sprintf(buffer, "maximum number of iterations reached without convergence. \n"); pcout<<buffer; break;exit (1);}
      if (current_norm>1/std::pow(tol,2)){sprintf(buffer, "\n norm is too high. \n\n"); pcout<<buffer; break; exit (1);}
      assemble_system();       
      current_norm=system_rhs.l2_norm();     
      initial_norm=std::max(initial_norm, current_norm);
      res=current_norm/initial_norm;
      sprintf(buffer,"inc:%3u (time:%10.3e, dt:%10.3e), iter:%2u, abs-norm: %10.2e, rel-norm: %10.2e\n", currentIncrement, currentTime, dt,  currentIteration, current_norm, res); pcout<<buffer; 
      if ((currentIteration>1) && ((res<tol) || (current_norm<abs_tol))){sprintf(buffer,"residual converged in %u iterations.\n\n", currentIteration); pcout<<buffer; break;}
      solveIteration(true);
      U+=dU; UGhost=U; 
      ++currentIteration;
    }
    Unn=Un;UnnGhost=Unn;  //previous to previous time step converged solution  
    Un=U;  //this step converged solution copied and this will get updated in update pressure
    UnGhost=Un;
    }

    if(!lever) {
      pcout << "Solving for projection "<<std::endl;	    
      while (true){
      if (currentIteration>=10){sprintf(buffer, "maximum number of iterations reached without convergence. \n"); pcout<<buffer; break;exit (1);}
      if (current_norm>1/std::pow(tol,2)){sprintf(buffer, "\n norm is too high. \n\n"); pcout<<buffer; break; exit (1);}
      assemble_system_projection();       
      current_norm=Pr_system_rhs.l2_norm();
      //     if(isnan(current_norm)) pcout << "rhs has a nan: " << current_norm << std::endl;      
      initial_norm=std::max(initial_norm, current_norm);
      res=current_norm/initial_norm;
      sprintf(buffer,"inc:%3u (time:%10.3e, dt:%10.3e), iter:%2u, abs-norm: %10.2e, rel-norm: %10.2e\n", currentIncrement, currentTime, dt,  currentIteration, current_norm, res); pcout<<buffer; 
      if ((currentIteration>1) && ((res<tol) || (current_norm<abs_tol))){sprintf(buffer,"residual converged in %u iterations.\n\n", currentIteration); pcout<<buffer; break;}
      solveIteration(false);
      Pr_U+=Pr_dU; Pr_UGhost=Pr_U; 
      ++currentIteration;
    }     
     Pr_Un=Pr_U; Pr_UnGhost=Pr_Un;

    }
  
  }


  
  //Adaptive grid refinement
  template <int dim>
  void phaseField<dim>::update_pressure ()  {
    //Un+=Pr_Un;
    //UnGhost.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    //UnGhost=Un;
    //applyBoundaryConditions(0);

    
    // TimerOutput::Scope t(computing_timer, "adaptiveRefinement");
    const QGauss<dim>  quadrature_formula(FEOrder+2);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points);
    unsigned int dofs_per_cell= fe_values.dofs_per_cell;
    unsigned int n_q_points= fe_values.n_quadrature_points;

    std::vector<Vector<double> > quadSolDiffusion,quadSolProjection;
    for (unsigned int q=0; q<n_q_points; ++q){
      quadSolDiffusion.push_back(dealii::Vector<double>(DIMS)); //2 since there are two degree of freedom per cell
      quadSolProjection.push_back(dealii::Vector<double>(DIMS)); //2 since there are two degree of freedom per cell
    }
   
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    typename parallel::distributed::Triangulation<dim>::active_cell_iterator t_cell = triangulation.begin_active();
   
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    for (;cell!=endc; ++cell) {
      
      if (cell->is_locally_owned() ) {

	fe_values.reinit (cell);
        
        cell->get_dof_indices (local_dof_indices);
         //AD variables	
	std::vector<double> temp (dofs_per_cell);
	//double temp2=0;
	for (unsigned int i=0; i<dofs_per_cell; ++i) {	   
	  const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - 0;	   
	  temp[i]=0;
	  if (ck==3) {
	      //note value of phi
	      temp[i]=Pr_UnGhost(local_dof_indices[i]);			    
	    }
	    
	}
	
	for (unsigned int i=0; i<dofs_per_cell; ++i) {	   
	  const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - 0;
	    if (ck==2) {
	      //update value of pressure
	      Un(local_dof_indices[i])+=temp[i];
	      
	    }	    
	    
	}
				
      }
      ++t_cell;

      } 

    Un.compress(VectorOperation::add);  	    	              	  	 	    
	    
  }

  
  //Output
  template <int dim>
  void phaseField<dim>::output_results (const unsigned int cycle) {
    TimerOutput::Scope t(computing_timer, "output");
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (UnGhost, nodal_solution_names, DataOut<dim>::type_dof_data, nodal_data_component_interpretation);    

    Vector<float> subdomain (triangulation.n_active_cells());
    for (unsigned int i=0; i<subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector (subdomain, "subdomain");
    
    data_out.build_patches (FEOrder);
    const std::string filename = ("solution-" +
                                  Utilities::int_to_string (cycle, 2) +
                                  "." +
                                  Utilities::int_to_string
                                  (triangulation.locally_owned_subdomain(), 4));
    std::ofstream output ((filename + ".vtu").c_str());
    data_out.write_vtu (output);
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
      std::vector<std::string> filenames;
      for (unsigned int i=0;
	   i<Utilities::MPI::n_mpi_processes(mpi_communicator);
	   ++i)
	filenames.push_back ("solution-" +
			     Utilities::int_to_string (cycle, 2) +
			     "." +
			     Utilities::int_to_string (i, 4) +
			     ".vtu");
      
      std::ofstream master_output (("solution-" +
				    Utilities::int_to_string (cycle, 2) +
				    ".pvtu").c_str());
      data_out.write_pvtu_record (master_output, filenames);
    }
  }

   

    //Solve problem
  template <int dim>
  void phaseField<dim>::run (){   
    GridIn<dim> grid_in;
    
    grid_in.attach_triangulation(triangulation);
      {
	std::string   filename = "nsbench2.inp";
	std::ifstream file(filename);
	Assert(file, ExcFileNotOpen(filename.c_str()));
	grid_in.read_ucd(file);
      }  
    
      //boundary_ids = triangulation.get_boundary_ids(); 
      

      /*  
    std::vector<unsigned int> numRepetitions;
    numRepetitions.push_back(XSubRf); // x refinement
    numRepetitions.push_back(YSubRf); // y refinement
    // numRepetitions.push_back(ZSubRf); // Z refinement
 
    Point<2> p1 (0,0);
    Point<2> p2 (problemWidth,problemHeight);
   
    GridGenerator::subdivided_hyper_rectangle (triangulation, numRepetitions, p1, p2, true);
    //GridGenerator::channel_with_cylinder(triangulation,0.03,2,2.0,true)
    //GridGenerator::hyper_cube_with_cylindrical_hole(triangulation,.25,.5,5.0,1,true ); 	
    */
    
    
    triangulation.refine_global (globalRefinementFactor); //global refinement
    setup_system(); //initial setup
    setup_system_projection();
    
    //setup initial conditions
    VectorTools::interpolate(dof_handler, InitalConditions<dim>(), U); Un=U;
    Unn=Un;
    //  VectorTools::interpolate(dof_handler, InitalConditions<dim>(), Pr_U); Pr_Un=Pr_U;    
   
    //sync ghost vectors to non-ghost vectors
    UGhost=U;
    UnGhost=Un;
    UnnGhost=Unn; 
    Pr_UGhost=Pr_U; Pr_UnGhost=Pr_Un;
    output_results (0);

    //Time stepping
    currentIncrement=0;
    for (currentTime=0; currentTime<totalTime; currentTime+=dt){
      currentIncrement++;
      //pcout << std::endl;      
      solve(true); //for diffuse solve
      // pcout << std::endl;
      //solve(false);// for projection solve
      //update_pressure();
      //UnGhost=Un;
      int NSTEP=(currentTime/dt);
      if (NSTEP%250==0) output_results(currentIncrement);      
      pcout << std::endl;
     
    }
    //computing_timer.print_summary ();
  }
}


int main(int argc, char *argv[]){
  try
    {
      using namespace dealii;
      using namespace phaseField1;
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      phaseField<2> problem;
      problem.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
