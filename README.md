# NeuronTransportIGA
NeuronTransportIGA performs material transport simulation in complex geometry of neurons using isogeometric analysis (IGA)

# User guide
The package contains four subfolders for different stages of the simulation pipeline. 
User should follow the steps to finish the whole simulation workflow:

### 1. meshgeneration (Matlab)
* **Description:** this code is used to smooth the neuron skeleton and generate hexahedral mesh for IGA simulation.
					 
* **Input:**
    *.swc (neuron skeleton information)	
    mesh_parameter.txt (parameter settings for generating the mesh)
  **Output:**
    controlmesh.vtk (control mesh file)
* **To run:**
1. Download [TREES toolbox](http://www.treestoolbox.org/) and put it inside the meshgeneration directory.
2. Use *TreeSmooth.m* as the code for smoothing neuron skeleton and generate **_smooth.swc* file. User need to set smooth parameters in *mesh_parameter.txt* and input/output path in *TreeSmooth.m*.
3. Use *Hexmesh_main.m* as the code to generate the mesh for smoothed skeleton. User need to set bifurcation refinement parameters in *mesh_parameter.txt* and input/output path in *Hexmesh_main.m*.
* **Notes:**
1. The neuron skeleton dataset can be downloaded from [NeuroMorpho.org](http://neuromorpho.org/). Before using the downloaded dataset, user needs to visualize the skeleton using TREES Toolbox and check the skeleton to make sure the geometry only has bifurcation structure and has no overlap. The bad geometry needs manual modifying right away.
2. In *mesh_parameter.txt*, user need to set five parameters:
    * n_noisesmooth:      set iteration steps for noise smooth
    * ratio_bifur_node:   set bifurcation nodes smooth ratio, range:0 to 1
	* ratio_noisesmooth:  set noise smooth ratio, range:0 to 1
	* seg_length:         set the bezier smooth segments length, the value mainly depends on the neuron scale
	* ratio_refine:       set the parameter used to calculate the refinement around bifurcation region range:0 to 1 
3. To obtain a good quality mesh, user needs to try several parameter settings for mesh generation. User is recommended to check the scaled Jacobian of the mesh in Paraview before using it in simulation. It's better to make sure the minimum scaled Jacobian is bigger than 0.1.
		  
### 2.spline_src (C++)
* **Description:** this code is used to construct B-spline and extract Bezier information for IGA based on the input controlmesh.
					 
* **Input:** 
    controlmesh.vtk (control mesh file)
  **Output:**
    bzpt.txt (Bezier point information)
    cmat.txt (The matrix from Bezier extraction)
    bzmeshinfo.txt (Bezier element connectivity, used for mesh partition)
		   
* **To compile:** (requires [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page))
    `>> make`


* **To run:**
`>> ./spline <meshfilepath>` 
(`meshfilepath` is the path that contains *controlmesh.vtk*)
Example: 
`>> ./spline ../example/bifurcation/`



