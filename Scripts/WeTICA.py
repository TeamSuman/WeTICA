
import os
import os.path as osp
import pickle as pkl

import simtk.openmm.app as omma
import simtk.openmm as omm
import simtk.unit as unit

import mdtraj as mdj

from wepy.sim_manager import Manager
from wepy.runners.openmm import OpenMMGPUWalkerTaskProcess, OpenMMRunner, OpenMMWalker, OpenMMState, gen_sim_state
from wepy.reporter.hdf5 import WepyHDF5Reporter
from wepy.work_mapper.task_mapper import TaskMapper
from wepy.util.mdtraj import mdtraj_to_json_topology

from walker_pkl_reporter import WalkersPickleReporter

from wepy.reporter.dashboard import DashboardReporter
from wepy.reporter.openmm import OpenMMRunnerDashboardSection


import pandas as pd
import numpy as np
import argparse

from Features import SelectedAtomsPairwise_distance, All_CACAdistance
from DistanceMetric import InvDist_FromRef 
from Resampler_new import REVOResampler
from warped import TargetBC





##------------------------------------------------ User defined inputs ----------------------------------------------##

##-------------------------------------------------------------------------------------------------------------------##


parser = argparse.ArgumentParser(description="Usage: python WeTICA.py [Options]\n")
parser.add_argument("-dire", type=str, help="Directory containing system related input files",required=True)
parser.add_argument("-num", type=int, help="Number of walkers [Optional, default=24]",required=False, default=24)
parser.add_argument("-rid", type=int, help="Run index",required=True)
parser.add_argument("-steps",  type=int, help="Number of MD steps between two consecutive resampling processes [optional, default=10000 -> 20 ps]",required=False, default=10000)
parser.add_argument("-cyc", type=int, help="Total number of cycles in each run [optional, default=20000]",required=False, default=20000)
parser.add_argument("-init", type=str, help="Starting structure file (GROMACS .gro format)",required=True)
parser.add_argument("-top", type=str, help="System topology file (GROMACS .top format)",required=True)
parser.add_argument("-ref", type=str, help="Target structure file (GROMACS .gro format)",required=True)
parser.add_argument("-vec", type=str, help="Eigenvector file",required=True)
parser.add_argument("-vid", type=int, help="Eigenvectors to be used e.g 0-> 1st vector, 1-> 2nd vector, ...",required=True, nargs='+')
parser.add_argument("-feat", type=str, help="Select feature e.g 'sel_pair' or 'allCA'.'sel_pair'-> selected pairwise atom distances & 'allCA'-> pairwise CA-CA distances",required=True)
parser.add_argument("-pair", type=str, help="File containing indices of selected atoms [Optional, but required for 'sel_pair' feature]", required=False)
parser.add_argument("-dist1", type=float, help="Merge distance cut-off (unitless) [Optional, default=1.0]",required=False, default=1.0)
parser.add_argument("-dist2", type=float, help="Warped distance cut-off (unitless) [Optional, default=0.25]",required=False, default=0.25)
parser.add_argument("-tem", type=float, help="Temperature in Kelvin [Optional, default=300]",required=False, default=300)
parser.add_argument("-ngpu", type=int, help="Number of available GPU cards", required=True)
parser.add_argument("-gid", type=int, help="GPU ids", required=True, nargs='+')
args = parser.parse_args()

DIR = args.dire
num_walkers = args.num                             
run =  args.rid                                      
n_steps =  args.steps                              
n_cycles =  args.cyc                             
start =  args.init                         
topol =  args.top                          
target =  args.ref                    
eigen =  args.vec                   
vec_list = args.vid 
sel_feat =  args.feat                      
atompairs = args.pair                  
d_merge =  args.dist1                                
d_warped = args.dist2 
temp = args.tem                              
n_gpu    = args.ngpu                                
gpu_ids  =  args.gid 


##-------------------------------------------------------- End ------------------------------------------------------##

##-------------------------------------------------------------------------------------------------------------------##







##-------------------------------------------- Main program using WEPY ----------------------------------------------##

##-------------------------------------------------------------------------------------------------------------------##


#### Paths: Set it to your preference
inp_path = DIR
gro_path = f'{inp_path}/'+start 
top_path = f'{inp_path}/'+topol
tar_path = f'{inp_path}/'+target
outputs_dir = f'{inp_path}/simdata_run{run}_steps{n_steps}_cycs{n_cycles}'
os.makedirs(outputs_dir, exist_ok=True)
os.mknod(f'{inp_path}/Info_{run}.txt')
####


if __name__ == "__main__":

    # Load atoms and topology objects 
    gro = omma.gromacsgrofile.GromacsGroFile(gro_path)
    top = omma.gromacstopfile.GromacsTopFile(top_path, periodicBoxVectors=gro.getPeriodicBoxVectors())
        
    # Get positions from gro file
    pos = gro.getPositions()

    # Make system
    system = top.createSystem(nonbondedMethod=omma.PME, nonbondedCutoff=1.0*unit.nanometer, constraints=omma.HBonds)
    system.addForce(omm.openmm.MonteCarloBarostat(1*unit.bar, temp*unit.kelvin)) # NPT ensemble
    integrator = omm.openmm.LangevinMiddleIntegrator(temp*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
    
    # Generate a new simtk "state"
    new_simtk_state = gen_sim_state(pos, system, integrator)
    
    # set up the OpenMMRunner with your system
    runner = OpenMMRunner(system, top.topology, integrator, platform='CUDA')
    
    # Load the reference state
    ref = mdj.load(tar_path)

    # Eigenvectors
    eigenvector_file = pd.read_csv(f'{inp_path}/'+eigen, sep='\t', header=None)
    eigenvectors = []
    for i in vec_list:
        eigenvectors.append(np.array(eigenvector_file.iloc[:,i]))

    if sel_feat == "sel_pair":

        # Read atom pairs file
        pair_file = pd.read_csv(f'{inp_path}/'+atompairs, sep='\t', header=None)

        # Atom pair indices
        atom_id1 = np.array(pair_file.iloc[:,0])
        atom_id2 = np.array(pair_file.iloc[:,1])
        atom_indexs = [atom_id1, atom_id2]

        # Projection of the reference state
        ref_pos = SelectedAtomsPairwise_distance(ref.xyz[0], atom_id1, atom_id2, eigenvectors)

    elif sel_feat == "allCA":

        # Select CA atoms
        atom_id = ref.topology.select('name CA')
        atom_indexs = atom_id

        # Projection of the reference state
        ref_pos = All_CACAdistance(ref.xyz[0], atom_id, eigenvectors)

    else:
        raise ValueError("Unrecognized feature selection")



    # Get the walker topology in a json format
    json_top = mdtraj_to_json_topology(mdj.load(gro_path).top)
    
    # Set up parameters for running the simulation
    init_weight = 1.0 / num_walkers
    
    
    
    
    
    #--------------------------------
    # Building wepy objects
    #--------------------------------
    print('Creating the wepy objects...')
    # Generate the walker state in wepy format
    walker_state = OpenMMState(new_simtk_state)
        
    # Make a list of the initial walkers
    init_walkers = [OpenMMWalker(walker_state, init_weight) for i in range(num_walkers)]
    
    
    # Distance metric to be used in resampling
    proj_distance = InvDist_FromRef(atom_indexs, ref_pos, eigenvectors, sel_feat)
    
    # Set up the Resampler with the parameters
    resampler = REVOResampler(distance=proj_distance,
                              init_state=walker_state,
                              merge_dist=d_merge,
                              run_id=run,
                              path=inp_path)
    
    # Set up the boundary conditions for a non-eq ensemble
    tbc = TargetBC(cutoff_distance=d_warped,
                   initial_state=walker_state,
                   idxs=atom_indexs,
                   target_pos=ref_pos,
                   eigenvecs=eigenvectors,
                   feat=sel_feat)
    
    
    # Set up the HDF5 reporter
    hdf5_reporter = WepyHDF5Reporter(save_fields=('positions','box_vectors'),
                                file_path=osp.join(outputs_dir,f'wepy.results.h5') ,
                                resampler=resampler,
                                boundary_conditions=tbc,
                                topology=json_top)

    # Set up the pickle reporter (Essential for restarts)
    out_folder_pkl = osp.join(outputs_dir,f'pkls')
    pkl_reporter = WalkersPickleReporter(save_dir = out_folder_pkl,
                                     freq = 1,
                                     num_backups = 2)

    # Set up the dashboard reporter
    dashboard_path = osp.join(outputs_dir,f'wepy.dash.org')
    openmm_dashboard_sec = OpenMMRunnerDashboardSection(runner)
    dashboard_reporter = DashboardReporter(file_path = dashboard_path,
                                        runner_dash = openmm_dashboard_sec)


    # Create a work mapper for NVIDIA GPUs for a GPU cluster
    mapper = TaskMapper(walker_task_type=OpenMMGPUWalkerTaskProcess,
                        num_workers=n_gpu,
                        platform='CUDA',
                        device_ids=gpu_ids)

    # Build the simulation manager
    sim_manager = Manager(init_walkers,
                          runner=runner,
                          resampler=resampler,
                          boundary_conditions=tbc,
                          work_mapper=mapper,
                          reporters=[hdf5_reporter, pkl_reporter, dashboard_reporter])

    #------------------------------
    # Run the simulation
    #------------------------------
    print('Running the simulation...')
    # run a simulation with the manager for 'n_cycles' with 'n_steps' of integrator steps in each
    steps_list = [n_steps for i in range(n_cycles)]

    # and..... go!
    sim_manager.run_simulation(n_cycles,
                               steps_list)


##-------------------------------------------------------- End ------------------------------------------------------##

##-------------------------------------------------------------------------------------------------------------------##
