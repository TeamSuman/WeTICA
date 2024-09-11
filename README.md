# WeTICA

python WeTICA.py -h

arguments:

      -h, --help          show this help message and exit
  
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
