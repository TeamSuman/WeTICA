# WeTICA

WeTICA is a binless weighted ensemble based enhanced sampling method to study rare event kinetics on a predefined linear CV space. It uses WEPY weighted ensemble software & OpenMM MD engine to run and store simulation data.

WEPY (https://adicksonlab.github.io/wepy/index.html)

OpenMM (https://openmm.org)


## Input Data Preparation

To start WeTICA simulation, at least 4 major files are required.

1) An equilibrated, neutrilize and solvated initial structure file in GROMACS .gro format.
2) System topology file in GROMACS .top format.
3) Target structure file in GROMACS .gro format.
4) An eigenvector file in .txt format where each column represents one eigenvector.
5) A file containing the pairwise selected atom ids in .txt format (Optional).

Example files are provided in the "Systems" folder of this repository.

#### Note 1 
If AmberTools is used to prepare system and topology files, use the following command to convert the amber structure and topology to GROMACS format.

      amb2gro_top_gro.py -p <name>.prmtop -c <name>.inpcrd -t <name>.top -g <name>.gro

#### Note 2
If CHARMM-GUI is used to prepare system and topology, files in GROMACS format can be generated from there itself.

#### Note 3
Currently WeTICA framework supports two featurization schemes: 1) CA-CA atom pairwise distances ('allCA'), 2) User selected atom pairwise distances ('sel_pair'). To maintain consistency in calculations, both the eigenvector generation and the simulation should be performed using the same feature.

#### Note 4
0 based selected atom ids (MDTraj convention) are required for the 'sel_pair' featurization scheme. This file is not required for the 'allCA' feature.





## Usage

#### step-1
First install OpenMM and other dependencies.

      conda create -n openmm python==3.9 numpy==1.25.2 pandas==1.3.5
      conda activate openmm
      conda install -c conda-forge openmm=8.0.0 mdtraj=1.9.9
      pip install h5py==3.12.1 networkx==3.2.1 tabulate==0.9.0 Jinja2==3.1.4 pint==0.24.3 eliot==1.14.0

#### step-2
Clone this repository in the local machine or HPC. 

      git clone https://github.com/TeamSuman/WeTICA.git


#### step-3
Go to the "Scripts" folder of this repository and run the WeTICA.py script with the necessary command line arguments as showcased below:

      python WeTICA.py -h
  
      -dire DIRE          Full PATH of the directory containing system related input files
  
      -num NUM            Number of walkers [Optional, default=16]
  
      -rid RID            Run index
  
      -steps STEPS        Number of MD steps between two consecutive resampling processes [optional, default=10000 -> 20 ps]
  
      -cyc CYC            Total number of cycles in each run [optional, default=20000]
  
      -init INIT          Starting structure file (GROMACS .gro format)
  
      -top TOP            System topology file (GROMACS .top format)
  
      -ref REF            Target structure file (GROMACS .gro format)
  
      -vec VEC            Eigenvector file
  
      -vid VID [VID ...]  Eigenvectors to be used e.g 0-> 1st vector, 1-> 2nd vector, ...
  
      -feat FEAT          Select feature e.g 'sel_pair' or 'allCA'.'sel_pair'-> selected pairwise atom distances & 'allCA'-> pairwise CA-CA distances
  
      -pair PAIR          File containing indices of selected atoms [Optional, but required for 'sel_pair' feature]
  
      -dist1 DIST1        Merge distance cut-off (unitless) [Optional, default=1.0]
  
      -dist2 DIST2        Warped distance cut-off (unitless) [Optional, default=0.25]
  
      -tem TEM            Temperature in Kelvin [Optional, default=300]
    
      -gid GID [GID ...]  GPU ids

##### Example

This example showcases the use of WeTICA.py script to study the unfolding of Protein G using the first two eigenvectors (specified as: -vid 0 1) as the collective variables (CVs).

    python WeTICA.py -dire Systems/protein_G -rid 1 -init folded.gro -top topol.top -ref unfolded.gro -vec eigenvectors.txt -vid 0 1
    -feat sel_pair  -pair atom_pairs.txt -gid 0 1 2 3 -tem 350

## On-the-fly monitoring of the simulation progress

While the WE simulation is running, a new file named "Info_*.txt" will be generated alongside the main HDF5 file inside the same folder that contains the system related input files. Informations from this file can be used to monitor the progress of the simulation on-the-fly without opening the HDF5 file. Gradual decrease of the 'Clst walk. dist' value indicates the progress of the simulation. Is “Clst walk. dist” ≈ d_warp? 
If yes, there is posibility that some walkers reach the target. Now its time for convergence check as described below.



## Data Analysis

To compute the Mean first passage time (MFPT), a jupyter notebook is provided in the "Analysis notebooks" folder of this repository. Three major analysis tutorials are provided in this notebook as follows:

    A) Convergence check
    B) MFPT calculation
    C) Generation of WE productive trajectory

For further information, visit wepy documentation (https://adicksonlab.github.io/wepy/_source/tutorials/data_analysis/index.html).

