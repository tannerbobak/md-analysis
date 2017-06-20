"""
This module defines some methods for trimming down large trajectory files.
Each method accepts GROMACS .xtc files for trajectories and GROMACS .gro
files for topologies.
"""
import mdtraj as md
import time


def clip(trajectory_file, topology_file, start, end):
    """
    Truncates a trajectory to be between the frames specified by start and
    end. It then saves the truncated trajectory to file with the suffix '_clipped'.
    :param trajectory_file: File path to load the trajectory from.
    :param topology_file: File path to load the topology from.
    :param start: Start of frames to truncate to (inclusive).
    :param end: End of frames to truncate to (exclusive).
    :return: The truncated trajectory.
    """
    trajectory = md.load(trajectory_file, top=topology_file)
    trajectory = trajectory[start:end]
    trajectory.save_xtc(trajectory_file[:-4] + '_clipped.xtc')
    return trajectory


def remove_solvent(trajectory_file, topology_file):
    """
    Removes solvent from a trajectory and saves the modified trajectory and topology
    to disk as an XTC and GRO file, respectively. Each will have the suffix
    '_no_solvent'
    :param trajectory_file: File path to load the trajectory from.
    :param topology_file: File path to load the topology from.
    :return: The updated trajectory.
    """
    my_time = time.time()
    # Load only frame 1 to select the protein atoms
    print('Loading protein atoms...')
    atoms = md.load_xtc(trajectory_file, top=topology_file, frame=1).topology.select('protein')
    print('Loading main trajectory...')
    trajectory = md.load_xtc(trajectory_file, top=topology_file, atom_indices=atoms)
    # No longer needed, by just loading protein there will be no solvent. Not calling this will save some time.
    # print('Removing solvent...')
    # trajectory.remove_solvent(inplace=True)
    print('Writing new trajectory to disk...')
    trajectory.save_xtc(trajectory_file[:-4] + '_no_solvent.xtc')
    print('Writing new topology to disk...')
    trajectory[0].save_gro(trajectory_file[:-4] + '_no_solvent.gro') # Save frame 1 to a gro file
    print('Done. Time elapsed: ' + str(time.time()-my_time) + ' seconds.')
    return trajectory

"""
# Stuff for running from file. Commented out to use this file as a module.

def main():
    trajectory_file = '../output/no-peptide/sh2b1-trjconv.xtc'
    topology_file = '../output/no-peptide/sh2b1.gro'
    # start = 0
    # end = 10
    # clip(trajectory_file, topology_file, start, end)
    remove_solvent(trajectory_file, topology_file)

if __name__ == "__main__":
    main()
"""