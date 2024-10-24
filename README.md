# WeTICA

WeTICA is a binless weighted ensemble based enhanced sampling method to study rare event kinetics on a predefined linear CV space. It uses WEPY weighted ensemble software & OpenMM MD engine to run and store simulation data.

WEPY (https://adicksonlab.github.io/wepy/index.html)

OpenMM (https://openmm.org)


## Input Data Preparation

To start WeTICA simulation, at least 4 major files are required.

1) An equilibrated, neutrilize and solvated $\bf{initial\ structure\ file}$ (GROMACS .gro format).
2) $\bf{System\ topology\ file}$ (GROMACS .top format).
3) $\bf{Target\ structure\ file}$ (GROMACS .gro format).
4) An $\bf{eigenvector\ file}$ (.txt format). Each column represents one eigenvector.
5) Pairwise selected atom id file (.txt format, Optional).

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
Follow the instructions step-by-step.

#### step-1 

      git clone https://github.com/TeamSuman/WeTICA.git

#### step-2

      cd WeTICA
      conda create --name wepy
      conda activate wepy
      conda install pip==21.2.2 wheel==0.37.1 python==3.7.13
      pip install -r requirements.txt
      conda install -c conda-forge openmm==7.5.1


#### step-3

      cd Scripts
      
      python WeTICA.py -h
  
             -dire DIRE          Full PATH of the directory containing system related input files
  
             -num NUM            Number of walkers [Optional, default=24]
  
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
  
             -dist1 DIST1        Merge distance cut-off (unitless) [Optional, default=0.5]
  
             -dist2 DIST2        Warped distance cut-off (unitless) [Optional, default=0.75]
  
             -tem TEM            Temperature in Kelvin [Optional, default=300]
    
             -gid GID [GID ...]  GPU ids

##### Example

This example showcases the use of WeTICA.py script to study the unfolding of Protein G using the first two eigenvectors (specified as: -vid 0 1) as the collective variables (CVs).

    python WeTICA.py -dire Systems/protein_G -rid 1 -init folded.gro -top topol.top -ref unfolded.gro -vec eigenvectors.txt -vid 0 1
    -feat sel_pair  -pair atom_pairs.txt -gid 0 1 2 3 -tem 350

## On-the-fly monitoring of the simulation progress

Open the $\bf{"Info.txt"}$ file generated inside the same folder that contains the system related input files. Gradual decrease of the $\bf{"Clst\ walk. dist"}$ value indicates the progress of the simulation. Is $\bf{Clst\ walk. dist ≈ d_{warp}?}$ If yes, there is posibility that some walkers reach the target. Now its time for convergence check as described below.


## Data Analysis

To compute the Mean first passage time (MFPT), a jupyter notebook is provided in the "Analysis notebooks" folder of this repository. Three major analysis tutorials are provided in this notebook as follows:

    A) Convergence check
    B) MFPT calculation
    C) Generation of WE productive trajectory

For further information, visit wepy documentation (https://adicksonlab.github.io/wepy/_source/tutorials/data_analysis/index.html).

