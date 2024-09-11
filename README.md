# WeTICA


## Usage

Inside the "Scripts" folder, several python files are there. The main script is the WeTICA.py. Download all the python files from this folder and put them together inside a folder in the local machine or HPC. Activete the wepy environment and run the WeTICA.py script with the necessary command line arguments as described below. Note that all the available arguments are not needed to specify explicitly. Many optional arguments are there with default values.

      python WeTICA.py -h
  
      -dire DIRE          Directory containing input files
  
      -num NUM            Number of walkers [Optional, default=24]
  
      -rid RID            Run index
  
      -steps STEPS        Number of MD steps between two consecutive resampling processes [optional, default=10000]
  
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

### Example

    python WeTICA.py -dire Systems/protein_G -rid 1 -init folded.gro -top topol.top -ref unfolded.gro -vec eigenvectors.txt -vid 0 1
    -feat sel_pair  -pair atom_pairs.txt -ngpu 4 -gid 0 1 2 3 -tem 350
