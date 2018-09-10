# NeuronTransportIGA
NeuronTransportIGA performs material transport simulation in complex geometry of neurons using isogeometric analysis (IGA)

## Dependencies
* *[TREES toolbox](http://www.treestoolbox.org/)*
* *[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)*
* *[METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)*
* *[PETSc 3.6.1](https://www.mcs.anl.gov/petsc/)*

## User guide
The package contains four subfolders for different stages of the simulation workflow. 
User should follow the steps to finish the whole simulation workflow:


### 1. meshgeneration (Matlab)

* **Description:** this code is used to smooth the neuron skeleton and generate hexahedral control mesh for IGA simulation.
* **Input:**
    * \*.swc (neuron skeleton information)   
    * mesh_parameter.txt (parameter settings for generating the mesh)    
* **Output:**
    * controlmesh.vtk (control mesh file)
* **To run:**
1. Download *[TREES toolbox](http://www.treestoolbox.org/)* and put it inside the meshgeneration directory.
2. Use *TreeSmooth.m* as the code for smoothing neuron skeleton and generate *\*_smooth.swc* file. User need to set smooth parameters in *mesh_parameter.txt* and input/output path in *TreeSmooth.m*.
3. Use *Hexmesh_main.m* as the code to generate the mesh for smoothed skeleton. User need to set bifurcation refinement parameters in *mesh_parameter.txt* and input/output path in *Hexmesh_main.m*.
* **Notes:**
1. The neuron skeleton dataset can be downloaded from *[NeuroMorpho.org](http://neuromorpho.org/)*. Before using the downloaded dataset, user needs to visualize the skeleton using TREES Toolbox and check the skeleton to make sure the geometry only has bifurcation structure and has no overlap. The bad geometry needs manually modifying before smoothing.
2. In *mesh_parameter.txt*, user need to set five parameters:
    * n_noisesmooth:      set iteration steps for noise smoothing, default: 100
    * ratio_bifur_node:   set bifurcation nodes smooth ratio, range: 0 to 1
    * ratio_noisesmooth:  set noise smooth ratio, range: 0 to 1
    * seg_length:         set the Bezier smooth segments length, the value mainly depends on the neuron scale
    * ratio_refine:       set the parameter used to calculate the refinement around bifurcation region range: 0 to 1 
3. To obtain a good quality mesh, user needs to try several parameter settings for mesh generation. User is recommended to check the scaled Jacobian of the mesh in Paraview before using it in simulation. It's better to make sure the minimum scaled Jacobian is bigger than 0.1.
          
### 2. spline_src (C++)

* **Description:** this code is used to construct truncated hierarchical tricubic spline and extract Bezier information for IGA based on the input controlmesh.
* **Input:**
    * controlmesh.vtk (control mesh file)
* **Output:**
    * bzpt.txt (Bezier point information)
    * cmat.txt (The matrix from Bezier extraction)
    * bzmeshinfo.txt (Bezier element connectivity, used for mesh partition)         
* **To compile:** (requires *[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)*)

    ` >> make`

* **To run:**

   ` >> ./spline <meshfilepath>` 

   `meshfilepath` is the path that contains *controlmesh.vtk*

   Example: 

   `>> ./spline ../example/bifurcation/`

### 3. METIS (C++)

* **Description:**
    The open source library *[METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)* is used to partition the mesh for parallel computing.

* **Input:**
    * bzmeshinfo.txt  (Bezier mesh information)
* **Output:**
    * bzmeshinfo.txt.epart.nparts (mesh partition file)
* **To run：**
    User can use the stand-alone program mpmetis in METIS library to partition a mesh into a specified number of parts.
    
   ` >> mpmetis <meshfile> process_num` 
     
   `process_num` is the number of parts that the mesh will be partitioned info
    
    Example:
    
    ` >> mpmetis bzmeshinfo.txt 28`
        
### 4. nsvms_src (C++)

* **Description:** this code is used to derive velocity field by solving incompressible steady-state Navier-Stokes equation. The code is paralleled using MPI to accelerate the computation.
                
* **Input:**
    * controlmesh.vtk  (control mesh file)
    * bzmeshinfo.txt.epart.nparts (mesh partition file)
    * bzpt.txt & cmat.txt (Bezier extraction information)
    * initial_velocity.txt (initial velocity information)
    * simulation_parameter.txt (parameter settings for the model)
* **Output:**
    * controlmesh_VelocityPressure_0.vtk (Velocity field on control mesh for visualization)
    * final_physics_VelocityPressureBezier.vtk (Velocity field in physical domain for visualization)
    * velocityfield.txt (Velocity information for transport simulation)

* **To compile:** (requires *[PETSc 3.6.1](https://www.mcs.anl.gov/petsc/)*)

   ` >> make` 

* **To run:**

   ` >> mpiexec -np process_num ./nsvms <meshfilepath> process_num`
   
   `process_num` is the number of processors used in simulation

   Example:

   ` >> mpiexec -np 28 ./nsvms ../example/bifurcation/ 28`
        
### 5. transport_src (C++)
* **Description:** this code is used to perform transport simulation and obtain concentration result.
        The code is paralleled using MPI to accelerate the computation. 

	User can define the initial condition in `UserSetting::SetInitialCondition` function in *UserSetting.cpp* and compile again to apply settings.

                
* **Input:**
    * controlmesh.vtk  (control mesh file)
    * bzmeshinfo.txt.epart.nparts (mesh partition file)
    * bzpt.txt & cmat.txt (Bezier extraction information)
    * velocity.txt (velocity field information)
    * simulation_parameter.txt (parameter settings for the model)
* **Output:**
    * controlmesh_allparticle_\*.vtk (Concentration for visualization)
    
    the output path name is set using all parameters value
            
* **To compile:** (requires *[PETSc 3.6.1](https://www.mcs.anl.gov/petsc/)*)

   ` >> make`  

* **To run:**

   ` >> mpiexec -np process_num ./transport <meshfilepath> process_num`
   
   `process_num` is the number of processors used in simulation

   Example: 

   ` >> mpiexec -np 28 ./transport ../example/bifurcation/ 28`

The solver code uses parallel computation (MPI) to accelerate computation. 
All result .vtk files are tested and can be opened in *[Paraview 5.4.0](https://www.paraview.org/)*

In `./example` folder, we provide all the input files for generating the result in our paper. For the single bifurcation neurite model, we provide all the intermediate output files in `./example/bifurcation/`. User can use this model to learn how to run the whole package.

