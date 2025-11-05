# RN\.py

A tool for running resistor network simulations based on simple voxel structures

## Table of Contents 

- [System requirements](#system-requirements)
- [Installation guide](#installation-guide)
- [Demos and user instructions](#demos-and-user-instructions)
    - [Example 1: series connected materials](#example-1---series-connected-materials-runtime--30-s)
    - [Example 2: parallel connected materials](#example-2---parallel-connected-materials-runtime--10-s)
    - [Example 3: mixed materials with clustering](#example-3---mixed-materials-with-clustering-runtime--10-min)
    - [Addition: Parallel script execution](#addition-parallel-script-execution)
- [References](#references)


## System requirements

RN\.py only requires packages from the standard python library as well as the numpy and the matplotlib package. The versions currently used in our lab are listed below:

- **Python**: 3.12.3
- **NumPy**: 1.26.4
- **Matplotlib**: 3.8.4

Small resistor networks like in the examples below can well be calculated on a standard PC. Computing large resistor networks however can require the usage of high performance computing infrastructure as long runtimes and large RAM are necessary.

## Installation guide

Installation time: < 5 min

- Download RN\.py as well as the libraries composite_lib.py and hotplate_lib.py from https://doi.org/10.17879/16948580876.
- Make sure you have installed the Python, NumPy and Matplotlib versions listed above.
- Specify the location of composite_lib.py and hotplate_lib.py in the RN\.py script. 

## Demos and user instructions

The first 30 lines of the code import the necessary libraries and define the important parameters for running the simulation. The following lines then generate the resistor network, approximate the steady state and calculate the final effective conductivity. .pdf files showing the composite microstructure as well as the progression of the residual value are generated, besides .txt files containing the effective conductivity and the residual value as a function of iteration respectively. Example inputs (lines 9 to 30 in RN\.py) and expected simulation results for effective conductivities are given below:

### Example 1 - series connected materials (runtime < 30 s)

For series connected materials 

$$ 
\kappa_{\text{eff}} = \left(\frac{\varphi_{1}}{\kappa_{1}} + \frac{\varphi_{2}}{\kappa_{2}}\right)^{-1} 
$$

can be expected [1]. Connecting equal volume fractions ($\varphi$) of materials with thermal conductivities of 0.71 $\text{W} \cdot (\text{m} \cdot \text{K})^{-1}$ and 0.32 $\text{W} \cdot (\text{m} \cdot \text{K})^{-1}$ should hence result in an effective thermal conductivitiy of 0.4412 $\text{W} \cdot (\text{m} \cdot \text{K})^{-1}$. To simulate a series-connected composite with materials showing thermal conductivities of 0.71 $\text{W} \cdot (\text{m} \cdot \text{K})^{-1}$ and 0.32 $\text{W} \cdot (\text{m} \cdot \text{K})^{-1}$, the following inputs can be specified.

```python
'''input parameters to run the calculation'''
name = 'series_connection'
Lvox = 20 # length of total structure in voxels
disp_volfrac = 50 # %, volumefraction of the material with conductivity kap_disp

kap_cont = 0.71 # W/(m*K), thermal cond (continuous Phase)  
kap_disp = 0.32 # W/(m*K), thermal cond (dispersed Phase)  
R_int = 0 # K*m^2/W, interfacial thermal resistance

x_len = 600*10**-6 # m, total length of the structure
dx = x_len/Lvox # m, distance between two voxels in m

T_left = 50 # °C, left hotplate Temperature
T_right = 0 # °C, right hotplate Temperature
T_comp = 25 # °C, initial composite Temperature

###############################################################################

'''generating and plotting the inital voxel structure''' 
voxelstructure = comp.seriesconnection3D(Lvox, disp_volfrac) # generate the voxelstructure 
comp.compositefigure(voxelstructure, show=True, save=True, name=name) # plot and save the figure of the voxelstructure

```

The expected simulation output is 0.4412 $\text{W} \cdot (\text{m} \cdot \text{K})^{-1}$, hence well agreeing with the analytical result.

### Example 2 - parallel connected materials (runtime < 10 s)

For parallel connected materials 

$$ 
\kappa_{\text{eff}} = \varphi_{1} \cdot \kappa_{1} + \varphi_{2} \cdot \kappa_{2}
$$

can be expected [1]. Connecting equal volume fractions ($\varphi$) of materials with thermal conductivities of 0.71 $\text{W} \cdot (\text{m} \cdot \text{K})^{-1}$ and 0.32 $\text{W} \cdot (\text{m} \cdot \text{K})^{-1}$ should hence result in an effective thermal conductivitiy of 0.515 $\text{W} \cdot (\text{m} \cdot \text{K})^{-1}$. Taking the input lines provided in Example 1, but changing the voxelstructure to 

```python
voxelstructure = comp.parallelconnection3D(Lvox, disp_volfrac) # generate the voxelstructure
```

leads to 0.5150 $\text{W} \cdot (\text{m} \cdot \text{K})^{-1}$ as expected simulation output, hence well agreeing with the analytical solution.


### Example 3 - mixed materials with clustering (runtime < 10 min)

To simulate the effective electronic conductivity of a composite with 60x60x60 total voxels (as shown section S21 of this work) containing 80 % of voxels representing material with an electronic conductivity of 5.22 $\text{mS} \cdot \text{cm}^{-1}$ and 20 % of insulating inclusions with a clustersize of 4 voxels the following input could be specified

```python
name = 'mixed_materials_with_clustering'
Lvox = 60 # length of total structure in voxels
disp_volfrac = 20 # %, volume fraction of the dispersed Phase
sclust = 4 # number of voxels per cluster of the dispersed Phase

kap_cont = 5.22 # mS/cm, electron cond (continuous Phase)  
kap_disp = 10**-100 # mS/cm, electron cond (dispersed Phase)  
R_int = 0 # interfacial resistance

x_len = 600*10**-6 # m, total length of the structure
dx = x_len/Lvox # m, distance between two voxels in m

T_left = 2.02 # V, Potential of the left side
T_right = 2.00 # V, Potential of the right side
T_comp = 2.01 # V, initial Potential

###############################################################################

'''generating and plotting the inital voxel structure''' 
voxelstructure = comp.blobs3D(Lvox, disp_volfrac, sclust) # generate the voxelstructure 
comp.compositefigure(voxelstructure, show=True, save=True, name=name) # plot and save the figure of the voxelstructure

```
The expected output is 3.02 $\text{mS} \cdot \text{cm}^{-1}$, however slight deviations ($\pm$ 0.01 $\text{mS} \cdot \text{cm}^{-1}$) can occur due to the randomness of the generated microstructures.

### Addition: Parallel script execution

To start multiple, similar RN\.py calculations at once the executing_parallel.py script can be used. To run it, it should be placed in the same folder as the RN\.py script. For such a use case, RN\.py could be prepared like this:

```python
name = input('samplename: ')
Lvox = 60 # length of total structure in voxels
disp_volfrac = int(input('volumefraction of the dispersed phase / %: '))
sclust = 4 # number of voxels per cluster of the dispersed Phase

kap_cont = 0.71 # W/(m*K), thermal cond (continuous Phase)  
kap_disp = 0.32 # W/(m*K), thermal cond (dispersed Phase)  
R_int = 0 # K*m^2/W, interfacial thermal resistance

x_len = 600*10**-6 # m, total length of the structure
dx = x_len/Lvox # m, distance between two voxels in m

T_left = 50 # °C, left hotplate Temperature
T_right = 0 # °C, right hotplate Temperature
T_comp = 25 # °C, initial composite Temperature

###############################################################################

'''generating and plotting the inital voxel structure''' 
voxelstructure = comp.blobs3D(Lvox, disp_volfrac, sclust) # generate the voxelstructure 
comp.compositefigure(voxelstructure, show=True, save=True, name=name) # plot and save the figure of the voxelstructure
```

while, the inputs in the executing_parallel.py script could be defined as

```python
'''Inputs needed to run the script'''
path = r'RN.py' 
disp_volfracs = [0,10,20,30,40,50,60,70,80,90,100] # composites with volumefractions to be calculated
name_appendage = 'Lvox60_sclust4' # choose a general file_name appendage. Each file name automatically starts with the volumefraction of the dispersed phase
```

11 processes are then started in parallel calculating effective transport of composites with different volume fractions.

## References

1. Carson, J.K. Modelling thermal diffusivity of heterogeneous materials based on thermal diffusivities of components with implications for thermal diffusivity and thermal conductivity measurement. *Int. J. Thermophys.* **43**, 108 (2022). https://doi.org/10.1007/s10765-022-03037-6   