# WeTICA


## Data Preparation

To start WeTICA weighted ensemble simulation, at least 4 major files are required. Example files are provided in the "Systems" directory.

1) An equilibrated neutrilize solvated initial structure file in GROMACS .gro format.
2) System topology file in GROMACS .top format.
3) Target structure file in GROMACS .gro format.
4) An eigenvector file in .txt format where each column represents one eigenvector.
5) A file containing the pairwise selected 0-based atom ids in .txt format (Optional).

Note 1a. If AmberTools is used to prepare system and topology files, "amb2gro_top_gro.py" python script available in AmberTools can be used to convert the amber structure and topology to GROMACS format.

Note 1b. If CHARMM-GUI is used to prepare system and topology, files in GROMACS formats can be generated from there itself.

Note 1c. ParmEd package is generally used to convert files between different format. For more information visit:
https://github.com/ParmEd/ParmEd




## Usage

First download the "wepy" weighted ensemble software. Open the link provided below and follow the installation guide.

      https://adicksonlab.github.io/wepy/_source/installation.html

Inside the "Scripts" folder of this repository, several python files are there. The main script is the WeTICA.py. Download all the python files from this folder and put them together inside a folder in the local machine or HPC. Activate the wepy environment and run the WeTICA.py script with the necessary command line arguments as described below. Note that all the available arguments are not needed to specify explicitly. Many optional arguments are there with default values. 

      python WeTICA.py -h
  
      -dire DIRE          Directory containing system related input files
  
      -num NUM            Number of walkers [Optional, default=24]
  
      -rid RID            Run index
  
      -steps STEPS        Number of MD steps between two consecutive resampling processes [optional, default=10000 -> 20 ps]
  
      -cyc CYC            Total number of cycles in each run [optional, default=20000]
  
      -init INIT          Starting structure file (GROMACS .gro format)
  
      -top TOP            System topology file (GROMACS .top format)
  
      -ref REF            Target structure file (GROMACS .gro format)
  
      -vec VEC            Eigenvector file
  
      -vid VID [VID ...]  Eigenvectors to be used e.g 0-> 1st vector, 1-> 2nd vector, ...
  
      -feat FEAT          Select feature e.g 'sel_pair' or 'allCA'.'sel_pair'-> selected pairwise atom distance & 'allCA'-> pairwise all CA-CA distance
  
      -pair PAIR          File containing indices of selected atoms [Optional, but required for 'sel_pair' feature]
  
      -dist1 DIST1        Merge distance cut-off (unitless) [Optional, default=1.0]
  
      -dist2 DIST2        Warped distance cut-off (unitless) [Optional, default=0.25]
  
      -tem TEM            Temperature in Kelvin [Optional, default=300]
  
      -ngpu NGPU          Number of available GPU cards
  
      -gid GID [GID ...]  GPU ids

#### Example

This example showcases the use of WeTICA.py script to study the unfolding of Protein G using the first two eigenvectors (specified as: -vid 0 1) as the collective variables (CVs).

    python WeTICA.py -dire Systems/protein_G -rid 1 -init folded.gro -top topol.top -ref unfolded.gro -vec eigenvectors.txt -vid 0 1
    -feat sel_pair  -pair atom_pairs.txt -ngpu 4 -gid 0 1 2 3 -tem 350

#### On-the-fly monitoring of the simulation progress

While the WE simulation is running, a new file named "Info_*.txt" will be generated alongside the main HDF5 file inside the same folder that contains the system related input files. Informations from this file can be used to monitor the progress of the simulation on-the-fly without opening the HDF5 file. Gradual decrease of the 'Clst walk. dist' value indicates the progress of the simulation. Is “Clst walk. dist” ≈ d_warp? 
If yes, there is posibility that some walkers reach the target. Now its time for convergence check as described below.



## Data Analysis

To compute the kinetic observables, a jupyter notebook is provided in the "Analysis notebooks" folder of this repository. Three major analysis tutorials are provided in this notebook as follows:

    A) Convergence check
    B) Mean first passage time (MFPT) calculation
    C) Generation of WE productive trajectory

Open the notebook and follow the instructions to calculate the above mentioned quantities. No major a priori knowledge about the HDF5 file is required to run the notebook.


## Workflow Diagram



![WE_TICA_protocol](https://github.com/user-attachments/assets/d60c2ea9-e2ce-445f-bda6-c1c66eaa5fef)

