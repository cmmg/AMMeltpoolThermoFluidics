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
#include "include/therm.h"
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
	values(0)=0.0;  //ux
	values(1)=0.0; //uy
	values(2)=0.0; //pressure	
	//Initial conditon for T
	values(3)=0.0; 	
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
    void solveIteration ();
    void solve ();
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
    LA::MPI::Vector                           UItm, UItmold, UItmGhost, DIFFU;     
    LA::MPI::Vector                           system_rhs;
    ConditionalOStream                        pcout;
    TimerOutput                               computing_timer;

    void L2_projection();
    void update_pressure (); 
    void applyBoundaryConditions_projection(const unsigned int increment);
    void setup_system_projection ();
    void assemble_system_projection();
    void solveIteration_Pr ();
    void solve_Pr ();
    void L2_solveIteration ();
    
    IndexSet                                  Pr_locally_owned_dofs;
    IndexSet                                  Pr_locally_relevant_dofs;
    DoFHandler<dim>                           Pr_dof_handler;
    ConstraintMatrix                          Pr_constraints, Pr_constraintsZero;
    AffineConstraints<double>                 Pr_mean_value_constraints;

    LA::MPI::SparseMatrix                     Pr_system_matrix;
    LA::MPI::Vector                           Pr_locally_relevant_solution,Pr_U, Pr_Un, Pr_UGhost, Pr_UnGhost, Pr_dU;
    LA::MPI::Vector                           Pr_system_rhs;
    LA::MPI::Vector                           Pr_UnnGhost, Pr_Unn;
    LA::MPI::Vector                           Pr_UItmGhost, Pr_UItm, Pr_UItmold;

    ConstraintMatrix                          L2_constraints;
    LA::MPI::SparseMatrix                     L2_Mass_matrix;
    LA::MPI::Vector                           L2_locally_relevant_solution, L2_U,L2_UGhost;    
    LA::MPI::Vector                           L2_system_rhs;


    DoFHandler<dim>                           T_dof_handler;
    IndexSet                                  T_locally_owned_dofs;
    IndexSet                                  T_locally_relevant_dofs;
    ConstraintMatrix                          T_constraints, T_constraintsZero;
    LA::MPI::SparseMatrix                     T_system_matrix;
    LA::MPI::Vector                           T_locally_relevant_solution, T_U, T_Un, T_UGhost, T_UnGhost, T_dU;
    LA::MPI::Vector                           T_Unn, T_UnnGhost;
    LA::MPI::Vector                           T_UItm,T_UItmold, T_UItmGhost, DIFFT;
    LA::MPI::Vector                           T_system_rhs;
   
    
    void applyBoundaryConditions_temp(const unsigned int increment);
    void setup_system_temp ();
    void assemble_system_temp();
    void solveIteration_temp ();
    void solve_temp ();
    void update_TintoU (); 
    
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
    fe(FE_Q<dim>(2),2,FE_Q<dim>(1),1,FE_Q<dim>(2),1),
    dof_handler(triangulation),
    Pr_dof_handler(triangulation),
    T_dof_handler(triangulation),
    pcout (std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator)== 0)),
    computing_timer (mpi_communicator, pcout, TimerOutput::summary, TimerOutput::wall_times){
    //solution variables
    dt=TimeStep; totalTime=TotalTime;
    currentIncrement=0; currentTime=0;

    //nodal Solution names

   for (unsigned int i=0; i<dim; ++i){
      nodal_solution_names.push_back("velocity"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
    }
   
   nodal_solution_names.push_back("pressure"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
   nodal_solution_names.push_back("Temp"); nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
 
  }
  
  template <int dim>
  phaseField<dim>::~phaseField (){
    dof_handler.clear ();
    Pr_dof_handler.clear();
    T_dof_handler.clear();
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
    std::vector<bool> uB (DIMS, false); uB[0]=true; uB[1]=true; 
    std::vector<bool> uBT (DIMS, false); uBT[0]=true; uBT[1]=true; 
 
    // 1 : walls top and bowttom , 2 : inlet 3: outlet 4: cavity walls
    
    //left
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(DIMS), constraints,uB);
    VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(DIMS), constraintsZero,uB);

    //right
    VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(DIMS) , constraints,uB);
    VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<dim>(DIMS) , constraintsZero,uB);

    //bottom
    VectorTools::interpolate_boundary_values (dof_handler, 2, ZeroFunction<dim>(DIMS), constraints,uB);
    VectorTools::interpolate_boundary_values (dof_handler, 2, ZeroFunction<dim>(DIMS), constraintsZero,uB);

    //top
    VectorTools::interpolate_boundary_values (dof_handler, 3, ZeroFunction<dim>(DIMS) , constraints,uBT);
    VectorTools::interpolate_boundary_values (dof_handler, 3, ZeroFunction<dim>(DIMS) , constraintsZero,uBT);    
    
    constraints.close ();
    constraintsZero.close ();
  }

  
  template <int dim>
  void phaseField<dim>::applyBoundaryConditions_projection(const unsigned int increment){
    Pr_constraints.clear (); L2_constraints.clear(); 
    Pr_constraints.reinit (Pr_locally_relevant_dofs);  L2_constraints.reinit (Pr_locally_relevant_dofs); 
    //DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    Pr_constraintsZero.clear (); 
    Pr_constraintsZero.reinit (Pr_locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (Pr_dof_handler, Pr_constraintsZero);
    
    //Setup boundary conditions
    std::vector<bool> uBC (DIMS, false);  uBC[2]=true;

    FEValuesExtractors::Scalar pressure(2); //dim is 2

    typename DoFHandler<dim>::active_cell_iterator
      cell = Pr_dof_handler.begin_active(), endc = Pr_dof_handler.end();

    //std::vector<bool> boundary_dofs(Pr_dof_handler.n_dofs(), false);
    IndexSet boundary_dofs;
    
    //constraints first dof of pressure to zero
    DoFTools::extract_boundary_dofs(Pr_dof_handler, fe.component_mask(pressure),
    boundary_dofs);

   const unsigned int first_boundary_dof = std::distance(boundary_dofs.begin(),
							 std::find (boundary_dofs.begin(), boundary_dofs.end(), true));
   
   // 1 : walls top and bowttom , 2 : inlet 3: outlet 4: cavity walls       
   //VectorTools::interpolate_boundary_values (Pr_dof_handler, 0, ZeroFunction<dim>(DIMS) , Pr_constraints,uBC);
   //VectorTools::interpolate_boundary_values (Pr_dof_handler, 0, ZeroFunction<dim>(DIMS) , Pr_constraintsZero,uBC);       

   //VectorTools::interpolate_boundary_values (Pr_dof_handler, 1, ZeroFunction<dim>(DIMS) , Pr_constraints,uBC);
   //VectorTools::interpolate_boundary_values (Pr_dof_handler, 1, ZeroFunction<dim>(DIMS) , Pr_constraintsZero,uBC);       

   //VectorTools::interpolate_boundary_values (Pr_dof_handler, 2, ZeroFunction<dim>(DIMS) , Pr_constraints,uBC);
   //VectorTools::interpolate_boundary_values (Pr_dof_handler, 2, ZeroFunction<dim>(DIMS) , Pr_constraintsZero,uBC);       

   //VectorTools::interpolate_boundary_values (Pr_dof_handler, 3, ZeroFunction<dim>(DIMS) , Pr_constraints,uBC);
   //VectorTools::interpolate_boundary_values (Pr_dof_handler, 3, ZeroFunction<dim>(DIMS) , Pr_constraintsZero,uBC);       
    
   Pr_mean_value_constraints.clear();
   Pr_mean_value_constraints.add_line(first_boundary_dof);

   for (unsigned int i = first_boundary_dof + 1; i < Pr_dof_handler.n_dofs(); ++i)
     //     if (boundary_dofs[i] == true)
     if (boundary_dofs.is_element(i) == true)
       Pr_mean_value_constraints.add_entry(first_boundary_dof, i, -1);
   Pr_mean_value_constraints.close();

   //Pr_constraints.add_line(first_boundary_dof);
   Pr_constraints.close ();
   Pr_constraintsZero.close ();
   L2_constraints.close ();
  }


  //Apply boundary conditions for diffusion step
  template <int dim>
  void phaseField<dim>::applyBoundaryConditions_temp(const unsigned int increment){
    T_constraints.clear (); 
    T_constraints.reinit (T_locally_relevant_dofs);
    //    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    T_constraintsZero.clear (); 
    T_constraintsZero.reinit (T_locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (T_dof_handler, T_constraintsZero);
    
   
    //Setup boundary conditions
    //std::vector<bool> uBL (DIMS, false); uBL[3]=true;    
    std::vector<bool> uBB (DIMS, false); uBB[3]=true;   
  
    // 1 : walls top and bowttom , 2 : inlet 3: outlet 4: cavity walls
    //right cold  wall
    VectorTools::interpolate_boundary_values (T_dof_handler, 1, ZeroFunction<dim>(DIMS), T_constraints,uBB);
    // VectorTools::interpolate_boundary_values (T_dof_handler, 1, ConstantFunction<dim>(0,DIMS), T_constraints,uBB);
    VectorTools::interpolate_boundary_values (T_dof_handler, 1, ZeroFunction<dim>(DIMS), T_constraintsZero,uBB);
    
    //left hot wall
    VectorTools::interpolate_boundary_values (T_dof_handler, 0, ConstantFunction<dim>(1.0,DIMS), T_constraints,uBB);
    VectorTools::interpolate_boundary_values (T_dof_handler, 0, ZeroFunction<dim>(DIMS) , T_constraintsZero,uBB);

    T_constraints.close ();
    T_constraintsZero.close ();
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
    
    UItm.reinit (locally_owned_dofs, mpi_communicator);
    UItmold.reinit (locally_owned_dofs, mpi_communicator);
    DIFFU.reinit(T_locally_owned_dofs, mpi_communicator);
    //Ghost vectors
    UGhost.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    UnGhost.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    UnnGhost.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    
    UItmGhost.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
    
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
    Pr_dof_handler.distribute_dofs (fe);
    Pr_locally_owned_dofs = Pr_dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (Pr_dof_handler,
                                             Pr_locally_relevant_dofs);
    
    Pr_locally_relevant_solution.reinit (Pr_locally_owned_dofs, Pr_locally_relevant_dofs, mpi_communicator);

    L2_locally_relevant_solution.reinit (Pr_locally_owned_dofs, Pr_locally_relevant_dofs, mpi_communicator);
    
    //Non-ghost vectors
    Pr_system_rhs.reinit (Pr_locally_owned_dofs, mpi_communicator);
    Pr_U.reinit (Pr_locally_owned_dofs, mpi_communicator);
    Pr_Un.reinit (Pr_locally_owned_dofs, mpi_communicator);
    Pr_dU.reinit (Pr_locally_owned_dofs, mpi_communicator);
    Pr_Unn.reinit (Pr_locally_owned_dofs, mpi_communicator);


    Pr_UItm.reinit (Pr_locally_owned_dofs, mpi_communicator);
    Pr_UItmold.reinit (Pr_locally_owned_dofs, mpi_communicator);

    L2_U.reinit (Pr_locally_owned_dofs, mpi_communicator);
    L2_system_rhs.reinit (Pr_locally_owned_dofs, mpi_communicator);
    
    //Ghost vectors
    Pr_UGhost.reinit (Pr_locally_owned_dofs, Pr_locally_relevant_dofs, mpi_communicator);
    Pr_UnGhost.reinit (Pr_locally_owned_dofs, Pr_locally_relevant_dofs, mpi_communicator);
    Pr_UnnGhost.reinit (Pr_locally_owned_dofs, Pr_locally_relevant_dofs, mpi_communicator);

    Pr_UItmGhost.reinit (Pr_locally_owned_dofs, Pr_locally_relevant_dofs, mpi_communicator);

    L2_UGhost.reinit (Pr_locally_owned_dofs, Pr_locally_relevant_dofs, mpi_communicator);
    //call applyBoundaryConditions to setup constraints matrix needed for generating the sparsity pattern
    applyBoundaryConditions_projection(0);
    
    DynamicSparsityPattern Pr_dsp (Pr_locally_relevant_dofs);
    DoFTools::make_sparsity_pattern (Pr_dof_handler, Pr_dsp, Pr_constraints, false);
    
    Pr_mean_value_constraints.condense(Pr_dsp);

SparsityTools::distribute_sparsity_pattern (Pr_dsp, Pr_dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, Pr_locally_relevant_dofs);
    Pr_system_matrix.reinit (Pr_locally_owned_dofs, Pr_locally_owned_dofs, Pr_dsp, mpi_communicator);


    L2_Mass_matrix.reinit (Pr_locally_owned_dofs, Pr_locally_owned_dofs, Pr_dsp, mpi_communicator);
    DynamicSparsityPattern L2_dsp (Pr_locally_relevant_dofs);
    DoFTools::make_sparsity_pattern (Pr_dof_handler, L2_dsp, L2_constraints, false);
    SparsityTools::distribute_sparsity_pattern (L2_dsp, Pr_dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, Pr_locally_relevant_dofs);
    L2_Mass_matrix.reinit (Pr_locally_owned_dofs, Pr_locally_owned_dofs, L2_dsp, mpi_communicator);    
  }

    //Setup
  template <int dim>
  void phaseField<dim>::setup_system_temp (){
    TimerOutput::Scope t(computing_timer, "setup");
    T_dof_handler.distribute_dofs (fe);
    
    T_locally_owned_dofs = T_dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (T_dof_handler,
                                             T_locally_relevant_dofs);  
    T_locally_relevant_solution.reinit (T_locally_owned_dofs, T_locally_relevant_dofs, mpi_communicator);

    //Non-ghost vectors
    T_system_rhs.reinit (T_locally_owned_dofs, mpi_communicator);
    T_U.reinit (T_locally_owned_dofs, mpi_communicator);
    T_Un.reinit (T_locally_owned_dofs, mpi_communicator);
    T_dU.reinit (T_locally_owned_dofs, mpi_communicator);    
    T_Unn.reinit (T_locally_owned_dofs, mpi_communicator);

    T_UItm.reinit (T_locally_owned_dofs, mpi_communicator);
    T_UItmold.reinit (T_locally_owned_dofs, mpi_communicator);
    DIFFT.reinit(T_locally_owned_dofs, mpi_communicator);

    //Ghost vectors
    T_UGhost.reinit (T_locally_owned_dofs, T_locally_relevant_dofs, mpi_communicator);
    T_UnGhost.reinit (T_locally_owned_dofs, T_locally_relevant_dofs, mpi_communicator);
    T_UnnGhost.reinit (T_locally_owned_dofs, T_locally_relevant_dofs, mpi_communicator);
    
    T_UItmGhost.reinit (T_locally_owned_dofs, T_locally_relevant_dofs, mpi_communicator);

    //call applyBoundaryConditions to setup constraints matrix needed for generating the sparsity pattern
    applyBoundaryConditions_temp(0);
    
    DynamicSparsityPattern T_dsp (T_locally_relevant_dofs);
    DoFTools::make_sparsity_pattern (T_dof_handler, T_dsp, T_constraints, false);
    SparsityTools::distribute_sparsity_pattern (T_dsp, T_dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, T_locally_relevant_dofs);
    T_system_matrix.reinit (T_locally_owned_dofs, T_locally_owned_dofs, T_dsp, mpi_communicator);
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
    //    FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, update_values | update_gradients | update_quadrature_points | update_JxW_values | update_normal_vectors);

    FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula,
				      update_values    |  update_gradients |
				      update_quadrature_points |
				      update_JxW_values  |
				      update_normal_vectors);

    
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

	Table<1, double > Pr_ULocalConv(dofs_per_cell);	Table<1, double > Pr_ULocalConvConv(dofs_per_cell);

	Table<1, double > T_ULocalConv(dofs_per_cell); Table<1, double > T_ULocalConvConv(dofs_per_cell);
	
	for (unsigned int i=0; i<dofs_per_cell; ++i){
	  const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - 0;

	  //if (std::abs(UGhost(local_dof_indices[i]))<1.0e-16){ULocal[i]=0.0;}
	  //else{ULocal[i]=UGhost(local_dof_indices[i]);}
	  //ULocal[i].diff (i, dofs_per_cell);
	  
	  if (std::abs(UItmGhost(local_dof_indices[i]))<1.0e-16){ULocal[i]=0.0;}
	  else{ULocal[i]=UItmGhost(local_dof_indices[i]);}
	  ULocal[i].diff (i, dofs_per_cell);
	  
	  ULocalConv[i]= UnGhost(local_dof_indices[i]);	 
	  ULocalConvConv[i]= UnnGhost(local_dof_indices[i]);

	  Pr_ULocalConv[i]= Pr_UnGhost(local_dof_indices[i]);	 
	  Pr_ULocalConvConv[i]= Pr_UnnGhost(local_dof_indices[i]);

	  //T_ULocalConv[i]= T_UnGhost(local_dof_indices[i]);	 
	  //T_ULocalConvConv[i]= T_UnnGhost(local_dof_indices[i]);
	  
	  //use T_conv from T intermediate
	  T_ULocalConv[i]= T_UItmGhost(local_dof_indices[i]);	 
	  T_ULocalConvConv[i]= T_UnnGhost(local_dof_indices[i]);

	}

	//setup residual vector
	Table<1, Sacado::Fad::DFad<double> > R(dofs_per_cell); 
	for (unsigned int i=0; i<dofs_per_cell; ++i) {R[i]=0.0;}
	
	residualForChemo(fe_values, 0, fe_face_values, cell, dt, ULocal, ULocalConv, ULocalConvConv,Pr_ULocalConv,Pr_ULocalConvConv,T_ULocalConv,T_ULocalConvConv, R,currentTime,totalTime);
	
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
  
    typename DoFHandler<dim>::active_cell_iterator cell = Pr_dof_handler.begin_active(), endc = Pr_dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned()){
	fe_values.reinit (cell);
	local_matrix = 0; local_rhs = 0; 
	cell->get_dof_indices (local_dof_indices);
	 //AD variables
	Table<1, double>  ULocal(dofs_per_cell); Table<1, double > ULocalConv(dofs_per_cell);
	Table<1, Sacado::Fad::DFad<double> > Pr_ULocal(dofs_per_cell); Table<1, double > Pr_ULocalConv(dofs_per_cell);
	
	for (unsigned int i=0; i<dofs_per_cell; ++i){
	  const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - 0;
	  
	  if (std::abs(Pr_UGhost(local_dof_indices[i]))<1.0e-16){Pr_ULocal[i]=0.0;}	  
	  else{Pr_ULocal[i]=Pr_UGhost(local_dof_indices[i]);}
	  
	  //if (std::abs(Pr_UItmGhost(local_dof_indices[i]))<1.0e-16){Pr_ULocal[i]=0.0;}	  
	  //else{Pr_ULocal[i]=Pr_UItmGhost(local_dof_indices[i]);}

	  Pr_ULocal[i].diff (i, dofs_per_cell);
	  Pr_ULocalConv[i]= Pr_UnGhost(local_dof_indices[i]);
	  //here UGhost has converged values from diffusion step
	  ULocal[i]=UGhost(local_dof_indices[i]);
	}


	//setup residual vector
	Table<1, Sacado::Fad::DFad<double> > Rp(dofs_per_cell); 
	for (unsigned int i=0; i<dofs_per_cell; ++i) {Rp[i]=0.0;}
		 
	//populate residual vctor 
	residualForProjection(fe_values, 0, fe_face_values, cell, dt, ULocal, Pr_ULocal, Rp, currentTime, totalTime);
	       		
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

 

    //Assembly
  template <int dim>
  void phaseField<dim>::assemble_system_temp (){
    TimerOutput::Scope t(computing_timer, "assembly");
    T_system_rhs=0.0; T_system_matrix=0.0;
    const QGauss<dim>  quadrature_formula(FEOrder+2);
    const QGauss<dim-1>	face_quadrature_formula (FEOrder+1);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points |
                             update_JxW_values);
    //    FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, update_values | update_gradients | update_quadrature_points | update_JxW_values | update_normal_vectors);

    FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula,
				      update_values    |  update_gradients |
				      update_quadrature_points |
				      update_JxW_values  |
				      update_normal_vectors);

    
    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       local_rhs (dofs_per_cell);
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    unsigned int n_q_points= fe_values.n_quadrature_points;

    unsigned int n_q_points_face= fe_face_values.n_quadrature_points;
    const unsigned int faces_per_cell = GeometryInfo<dim>::faces_per_cell;
  
    typename DoFHandler<dim>::active_cell_iterator cell = T_dof_handler.begin_active(), endc = T_dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned()){
	fe_values.reinit (cell);
	local_matrix = 0; local_rhs = 0; 
	cell->get_dof_indices (local_dof_indices);
	 //AD variables
	Table<1, Sacado::Fad::DFad<double> > T_ULocal(dofs_per_cell); Table<1, double > T_ULocalConv(dofs_per_cell);
	Table<1, double > T_ULocalConvConv(dofs_per_cell);

	Table<1, Sacado::Fad::DFad<double> > ULocal(dofs_per_cell); Table<1, double > ULocalConv(dofs_per_cell);
	Table<1, double > ULocalConvConv(dofs_per_cell);

	for (unsigned int i=0; i<dofs_per_cell; ++i){
	  const unsigned int ck = fe_values.get_fe().system_to_component_index(i).first - 0;

	  //  if (std::abs(T_UGhost(local_dof_indices[i]))<1.0e-16){T_ULocal[i]=0.0;}
	  //else{T_ULocal[i]=T_UGhost(local_dof_indices[i]);}

	  if (std::abs(T_UItmGhost(local_dof_indices[i]))<1.0e-16){T_ULocal[i]=0.0;}
	  else{T_ULocal[i]=T_UItmGhost(local_dof_indices[i]);}
		  

	  T_ULocal[i].diff (i, dofs_per_cell);
	  T_ULocalConv[i]= T_UnGhost(local_dof_indices[i]);	 
	  T_ULocalConvConv[i]= T_UnnGhost(local_dof_indices[i]);

	  //momentum equation terms
	  ULocal[i]= UItmGhost(local_dof_indices[i]);	 
	  //ULocal[i]= UGhost(local_dof_indices[i]);	 
	  ULocalConv[i]= UnGhost(local_dof_indices[i]);	 
	  ULocalConvConv[i]= UnnGhost(local_dof_indices[i]);	  
	}

	//setup residual vector
	Table<1, Sacado::Fad::DFad<double> > R(dofs_per_cell); 
	for (unsigned int i=0; i<dofs_per_cell; ++i) {R[i]=0.0;}
			
	residualForTherm(fe_values, 0, fe_face_values, cell, dt, T_ULocal, T_ULocalConv, T_ULocalConvConv,ULocal, ULocalConv, ULocalConvConv, R, currentTime, totalTime) ;
	
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
	  T_constraints.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, T_system_matrix, T_system_rhs);	
	}
	else{
	  T_constraintsZero.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, T_system_matrix, T_system_rhs);
	}
      }
    T_system_matrix.compress (VectorOperation::add);
    T_system_rhs.compress (VectorOperation::add);
    if(isnan(T_system_matrix.frobenius_norm())) pcout << "\nK norm diffusion has a nan: " << T_system_matrix.frobenius_norm() << std::endl; 
  }



  
  //Solve
  template <int dim>
  void phaseField<dim>::solveIteration(){

    TimerOutput::Scope t(computing_timer, "solve");
    LA::MPI::Vector completely_distributed_solution (locally_owned_dofs, mpi_communicator);          
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
	
    //std::cout <<"sys matx  is " <<system_matrix.frobenius_norm()<<std::endl;

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
 

 //Solve
  template <int dim>
  void phaseField<dim>::solveIteration_Pr(){
    //Iterative solvers from Petsc and Trilinos
    TimerOutput::Scope t(computing_timer, "solve");
    LA::MPI::Vector completely_distributed_solution (Pr_locally_owned_dofs, mpi_communicator);    	
    SolverControl solver_control (Pr_dof_handler.n_dofs(), 1e-12);
#ifdef USE_PETSC_LA
    LA::SolverGMRES solver(solver_control, mpi_communicator);
#else
    LA::SolverGMRES solver(solver_control);
#endif
    LA::MPI::PreconditionAMG preconditioner;
    LA::MPI::PreconditionAMG::AdditionalData data;
#ifdef USE_PETSC_LA //dof_handler.get_fe().n_components()
    // data.symmetric_operator = true;
#else
    // Trilinos defaults are good 
#endif
    //if lever is true run for diffusion part 
    preconditioner.initialize(Pr_system_matrix, data);
    
    solver.solve (Pr_system_matrix, completely_distributed_solution, Pr_system_rhs, preconditioner);
    pcout << "   Solved in " << solver_control.last_step()
          << " iterations." << std::endl;     
    		

    	//Direct solver MUMPS
    //SolverControl cn;
    //PETScWrappers::SparseDirectMUMPS solver(cn, mpi_communicator);
    //solver.set_symmetric_mode(true);
    //solver.solve(system_matrix, completely_distributed_solution, system_rhs);


    if ((currentIteration==0)&&(currentIncrement==1)){
	  Pr_constraints.distribute (completely_distributed_solution);
    }
	else{
	  Pr_constraintsZero.distribute (completely_distributed_solution);
	}
    
	Pr_locally_relevant_solution = completely_distributed_solution;
	Pr_dU=completely_distributed_solution;
                    
  }   

  
  //Solve
  template <int dim>
 void phaseField<dim>::L2_solveIteration(){
    //Iterative solvers from Petsc and Trilinos
    TimerOutput::Scope t(computing_timer, "solve");
    LA::MPI::Vector completely_distributed_solution (Pr_locally_owned_dofs, mpi_communicator);    	
    SolverControl solver_control (Pr_dof_handler.n_dofs(), 1e-12);
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
    preconditioner.initialize(L2_Mass_matrix, data);
    
    solver.solve (L2_Mass_matrix, completely_distributed_solution, L2_system_rhs, preconditioner);  
    //    pcout << "   Solved in " << solver_control.last_step()
    //    << " iterations." << std::endl;         	       

    L2_constraints.distribute (completely_distributed_solution);
    L2_locally_relevant_solution = completely_distributed_solution;
    L2_U=completely_distributed_solution;
    if (isnan(L2_U.l2_norm())) std:: cout << "Nan in L2projection" << std::endl;
    if (isinf(L2_U.l2_norm())) std:: cout << "inf in L2projection" << std::endl;
    
    L2_UGhost=L2_U;
	
  }
    

    //Solve
  template <int dim>
  void phaseField<dim>::solveIteration_temp(){

    TimerOutput::Scope t(computing_timer, "solve");
    LA::MPI::Vector completely_distributed_solution (T_locally_owned_dofs, mpi_communicator);          
      //Iterative solvers from Petsc and Trilinos
    SolverControl solver_control (T_dof_handler.n_dofs(), 1e-12);
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
    preconditioner.initialize(T_system_matrix, data);
    solver.solve (T_system_matrix, completely_distributed_solution, T_system_rhs, preconditioner);
    pcout << "   Solved in " << solver_control.last_step()
          << " iterations." << std::endl;
	
    //std::cout <<"sys matx  is " <<system_matrix.frobenius_norm()<<std::endl;

	//Direct solver MUMPS
    //	SolverControl cn;
    //	PETScWrappers::SparseDirectMUMPS solver(cn, mpi_communicator);
    //	solver.set_symmetric_mode(false);
    //	solver.solve(system_matrix, completely_distributed_solution, system_rhs);
    
  	if ((currentIteration==0)&&(currentIncrement==1)){
	  T_constraints.distribute (completely_distributed_solution);
	}
	else{
	  T_constraintsZero.distribute (completely_distributed_solution);
	}
	T_locally_relevant_solution = completely_distributed_solution;
        T_dU = completely_distributed_solution; 	
                         
  }   


  
  //Solve
  template <int dim>
  void phaseField<dim>::solve() {     
    pcout << "Solving for intermediate toggle  "<<std::endl<<std::endl;      
    double error_U=10; double error_T=10; double MAXX=10; unsigned int toggleSteps=0;
    while (true) {    
      //for (unsigned int tog=0 ; tog<20 ; tog++) {
      double res=1, tol=1.0e-12, abs_tol=1.0e-11, initial_norm=0, current_norm=0;
      double machineEPS=1.0e-15;
      currentIteration=0;
      char buffer[200];          
      pcout << "Solving for momentum "<<std::endl;      
      while (true){
	if (currentIteration>=10){sprintf(buffer, "maximum number of iterations reached without convergence. \n \n"); pcout<<buffer; break;exit (1);}
	if (current_norm>1/std::pow(tol,2)){sprintf(buffer, "\n norm is too high. \n\n"); pcout<<buffer; break; exit (1);}
	assemble_system();       
	current_norm=system_rhs.l2_norm();     
	initial_norm=std::max(initial_norm, current_norm);
	res=current_norm/initial_norm;
	sprintf(buffer,"inc:%3u (time:%10.3e, dt:%10.3e), iter:%2u, abs-norm: %10.2e, rel-norm: %10.2e\n", currentIncrement, currentTime, dt,  currentIteration, current_norm, res); pcout<<buffer; 
	if ((currentIteration>1) && ((res<tol) || (current_norm<abs_tol))){sprintf(buffer,"residual converged in %u iterations.\n\n", currentIteration); pcout<<buffer; break;}
	solveIteration();
	//U+=dU; UGhost=U;
	UItm+=dU; UItmGhost=UItm; 
	++currentIteration;
      }
      //Solve for T using UItm
      solve_temp();     
      
      //calculate tolerance 
      ++toggleSteps;
      DIFFU=UItm; DIFFU-=UItmold; DIFFT=T_UItm; DIFFT-=T_UItmold;
      error_U=DIFFU.l2_norm();error_T=DIFFT.l2_norm(); 
      MAXX=std::max(error_U,error_T);
      //  sprintf(buffer,"time:%10.3e, dt:%10.3e, togg:%u, error: %10.2e \n \n", currentTime, dt,  toggleSteps, MAXX); pcout<<buffer;       
      if (MAXX<TOLERANCE ) {	
	pcout << "Value of MAXX is:  " << MAXX << std::endl; 
	sprintf(buffer,"intermediate steps converged in %u steps, error is: %10.2e \n \n", toggleSteps, MAXX); pcout<<buffer; break;
      }
      
      
      if (1/*MAXX<TOLERANCE*/ ) {
	//	error_U=10; error_T=10;
	//pcout << "Value of MAXX is:  " << MAXX << std::endl; 
	sprintf(buffer,"intermediate step no. is  %u steps, error is: %10.2e \n \n", toggleSteps, MAXX); pcout<<buffer; //break;
	}
      
      
      if (toggleSteps>20){sprintf(buffer, "maximum number of intermediate iterations reached without convergence . \n \n"); pcout<<buffer; break;exit(1);}           
      error_U=10; error_T=10; MAXX=10;
      UItmold=UItm; T_UItmold=T_UItm;
    }
    
    //Do update after satisfying tolerance
    
    U=UItm; UGhost=U;     
    T_U=T_UItm; T_UGhost=T_U; 
    pcout << std::endl;   
    solve_Pr();// for projection solve : solve based on converged UItm 
    // L2_projection();
    update_pressure(); //update values of U and UGhost    
    Unn=Un; UnnGhost=Unn;        
    T_Unn=T_Un; T_UnnGhost=T_Unn;
    T_Un=T_U; T_UnGhost=T_Un;
    update_TintoU();
    Un=U; UnGhost=Un; // copy updated values in Un and UnGhost;
  }



   //Solve
  template <int dim>
  void phaseField<dim>::solve_Pr() {     
      double res=1, tol=1.0e-12, abs_tol=1.0e-11, initial_norm=0, current_norm=0;
      double machineEPS=1.0e-15;
      currentIteration=0;
      char buffer[200];          
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
      solveIteration_Pr();
      Pr_U+=Pr_dU; Pr_UGhost=Pr_U;
      ++currentIteration;
    }
      //Pr_U.operator*=(1.5/dt);
      Pr_Unn=Pr_Un; Pr_UnnGhost=Pr_Unn;
      Pr_Un=Pr_U; Pr_UnGhost=Pr_Un;
     //update values   
  }


  
   //Solve
  template <int dim>
  void phaseField<dim>::solve_temp() {     
      double res=1, tol=1.0e-12, abs_tol=1.0e-11, initial_norm=0, current_norm=0;
      double machineEPS=1.0e-15;
      currentIteration=0;
      char buffer[200];          
      pcout << "Solving for temperature "<<std::endl;	    
      while (true){
      if (currentIteration>=10){sprintf(buffer, "maximum number of iterations reached without convergence. \n"); pcout<<buffer; break;exit (1);}
      if (current_norm>1/std::pow(tol,2)){sprintf(buffer, "\n norm is too high. \n\n"); pcout<<buffer; break; exit (1);}
      assemble_system_temp();       
      current_norm=T_system_rhs.l2_norm();
      //     if(isnan(current_norm)) pcout << "rhs has a nan: " << current_norm << std::endl;      
      initial_norm=std::max(initial_norm, current_norm);
      res=current_norm/initial_norm;
      sprintf(buffer,"inc:%3u (time:%10.3e, dt:%10.3e), iter:%2u, abs-norm: %10.2e, rel-norm: %10.2e\n", currentIncrement, currentTime, dt,  currentIteration, current_norm, res); pcout<<buffer; 
      if ((currentIteration>1) && ((res<tol) || (current_norm<abs_tol))){sprintf(buffer,"residual converged in %u iterations.\n\n", currentIteration); pcout<<buffer; break;}
      solveIteration_temp();
      //T_U+=T_dU; T_UGhost=T_U; 
      T_UItm+=T_dU; T_UItmGhost=T_UItm; 
      ++currentIteration;
    }
    
  }

  
  
  //Adaptive grid refinement
  template <int dim>
  void phaseField<dim>::L2_projection ()  {    
    // TimerOutput::Scope t(computing_timer, "adaptiveRefinement");
    L2_system_rhs=0.0; L2_Mass_matrix=0.0;
    const QGauss<dim>  quadrature_formula(FEOrder+2);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points |update_JxW_values);
    unsigned int dofs_per_cell= fe_values.dofs_per_cell;
    unsigned int n_q_points= fe_values.n_quadrature_points;    
    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       local_rhs (dofs_per_cell);
    std::vector<Vector<double> > quadSolutions,Pr_quadSolutions;
    std::vector< std::vector< Tensor< 1, dim, double >>> quadGradU;  
    for (unsigned int q=0; q<n_q_points; ++q){
      Pr_quadSolutions.push_back(dealii::Vector<double>(DIMS)); //2 since there are two degree of freedom per cell
      quadSolutions.push_back(dealii::Vector<double>(DIMS));
      quadGradU.push_back(std::vector< Tensor< 1, dim, double >>(DIMS));
    }
        
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    typename parallel::distributed::Triangulation<dim>::active_cell_iterator t_cell = triangulation.begin_active();   
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    std::vector<double> press_update(dofs_per_cell);
    //get values of p,phi and nu time div v 
    for (;cell!=endc; ++cell)      
      if (cell->is_locally_owned() ) {
	fe_values.reinit (cell);        
        cell->get_dof_indices (local_dof_indices);

	fe_values.get_function_values(UnnGhost, quadSolutions);
	fe_values.get_function_values(Pr_UnGhost, Pr_quadSolutions);
	fe_values.get_function_gradients(UGhost, quadGradU);	

	//rotational pressure update
	for (unsigned int q=0; q<n_q_points; ++q) {
	     for (unsigned int i=0; i<dofs_per_cell ; ++i) {
	       //RHS
	       unsigned int ci = fe_values.get_fe().system_to_component_index(i).first - 0;
	       if(ci==2) {
		 for(unsigned int d=0; d<dim ; ++d) {
		   local_rhs[i] +=fe_values.shape_value_component(i, q, ci)*(quadGradU[q][d][d])*fe_values.JxW(q);		   
		 }
		 //local_rhs[i] +=fe_values.shape_value_component(i, q, ci)*(Pr_quadSolutions[q][2])*fe_values.JxW(q);
		 //local_rhs[i] +=fe_values.shape_value_component(i, q, ci)*(quadSolutions[q][2])*fe_values.JxW(q);
	       }
	       
	       for (unsigned int j=0; j<dofs_per_cell ; ++j) {
		 unsigned int cj = fe_values.get_fe().system_to_component_index(j).first - 0;
		 if(ci==cj && ci==2) { 
		   local_matrix[i][j]+=fe_values.shape_value_component(i, q,dim)*fe_values.shape_value_component(j, q,dim)*fe_values.JxW(q);	    }
	       }
	     }
	}
	
	L2_constraints.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, L2_Mass_matrix, L2_system_rhs);        	     	
      }
       
    L2_Mass_matrix.compress (VectorOperation::add);
    L2_system_rhs.compress (VectorOperation::add);
    L2_solveIteration();
  }



    //Adaptive grid refinement
  template <int dim>
  void phaseField<dim>::update_pressure ()  {    
    // TimerOutput::Scope t(computing_timer, "adaptiveRefinement");
    const QGauss<dim>  quadrature_formula(FEOrder+2);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points |update_JxW_values);
    unsigned int dofs_per_cell= fe_values.dofs_per_cell;
    unsigned int n_q_points= fe_values.n_quadrature_points;
    
       
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    typename parallel::distributed::Triangulation<dim>::active_cell_iterator t_cell = triangulation.begin_active();
   
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    std::vector<double> press_update(dofs_per_cell);

    //get values of p,phi and nu time div v 
    for (;cell!=endc; ++cell) {     
      if (cell->is_locally_owned() ) {
	fe_values.reinit (cell);        
        cell->get_dof_indices (local_dof_indices);
	
	//pressure update
	for (unsigned int i=0; i<dofs_per_cell; ++i) {
	  press_update[i]=0;
	  const unsigned int ci = fe_values.get_fe().system_to_component_index(i).first - 0;	   
	  //note value of pressue L2_ UGhost	  
	  if (ci==2) {
	    press_update[i]=Pr_UnGhost(local_dof_indices[i]); //phi_k+1 is noted
	    press_update[i]+=UnGhost(local_dof_indices[i]); //press_k is noted and added to phi_k+1
	    //    press_update[i]+=-(nu)*L2_UGhost(local_dof_indices[i]); //nu. div.v
	    U(local_dof_indices[i])=press_update[i];	    

	  }	  
	}
	     	
      }
    }
    
    
    U.compress(VectorOperation::insert);	
  }


  
  //Adaptive grid refinement
  template <int dim>
  void phaseField<dim>::update_TintoU ()  {    
    // TimerOutput::Scope t(computing_timer, "adaptiveRefinement");
    const QGauss<dim>  quadrature_formula(FEOrder+2);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points |update_JxW_values);
    unsigned int dofs_per_cell= fe_values.dofs_per_cell;
    unsigned int n_q_points= fe_values.n_quadrature_points;
    
       
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    typename parallel::distributed::Triangulation<dim>::active_cell_iterator t_cell = triangulation.begin_active();
   
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    std::vector<double> T_update(dofs_per_cell), Liq_update(dofs_per_cell);

    //get values of p,phi and nu time div v 
    for (;cell!=endc; ++cell) {     
      if (cell->is_locally_owned() ) {
	fe_values.reinit (cell);        
        cell->get_dof_indices (local_dof_indices);
	
	//pressure update
	for (unsigned int i=0; i<dofs_per_cell; ++i) {
	  T_update[i]=0;Liq_update[i]=0;
	  const unsigned int ci = fe_values.get_fe().system_to_component_index(i).first - 0;	   
	  //note value of pressue L2_ UGhost	  
	  if (ci==3) {
	    T_update[i]=T_UnGhost(local_dof_indices[i]); //phi_k+1 is noted
	    U(local_dof_indices[i])=T_update[i];	    
	  }
	  if (ci==5) {
	    //Liq_update[i]=T_UnGhost(local_dof_indices[i]); //phi_k+1 is noted
	    //U(local_dof_indices[i])=Liq_update[i];	    
	  }
	}
	     	
      }
    }
    
    
    U.compress(VectorOperation::insert);	
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
    
    
    /*GridIn<dim> grid_in; 
    grid_in.attach_triangulation(triangulation);
      {
	std::string   filename = "nsbench2.inp";
	std::ifstream file(filename);
	Assert(file, ExcFileNotOpen(filename.c_str()));
	grid_in.read_ucd(file);
	}*/  
      
    std::vector<unsigned int> numRepetitions;
    numRepetitions.push_back(XSubRf); // x refinement
    numRepetitions.push_back(YSubRf); // y refinement
    //numRepetitions.push_back(ZSubRf); // z refinement
    Point<dim> p1 (0,0);
    Point<dim> p2 (problemLength,problemHeight);
    GridGenerator::subdivided_hyper_rectangle (triangulation, numRepetitions, p1, p2, true);
   


    triangulation.refine_global (globalRefinementFactor); //global refinement
    setup_system(); //initial setup
    setup_system_projection();   
    setup_system_temp();

    
    //setup initial conditions
    VectorTools::interpolate(dof_handler, InitalConditions<dim>(), U); Un=U; Unn=Un;
    Pr_U=U;Pr_Un=Pr_U;Pr_Unn=Pr_Un;
    T_U=U; T_Un=T_U; T_Unn=T_Un ;
    //sync ghost vectors to non-ghost vectors
    UGhost=U; UnGhost=Un;UnnGhost=Unn;
    Pr_UGhost=Pr_U;Pr_UnGhost=Pr_Un;Pr_UnnGhost=Pr_Unn;

    T_UGhost=U; T_UnGhost=Un;T_UnnGhost=Unn;
   
    output_results (0);
    
    //Time stepping
    currentIncrement=0;
    for (currentTime=2*dt; currentTime<totalTime; currentTime+=dt){
      currentIncrement++;
      pcout << std::endl;         
            
      if (currentTime>4*dt) {
	solve();
      } //for diffuse solve                  
      else {
      solve_temp();
      update_TintoU();
      Unn=Un; UnnGhost=Unn; // copy updated values in Un and UnGhost;
      Un=U; UnGhost=Un; // copy updated values in Un and UnGhost;
      } 

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
