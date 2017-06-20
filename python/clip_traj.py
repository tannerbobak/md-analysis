import mdtraj as md
import time


def clip(trajectory_file, topology_file, start, end):
    trajectory = md.load(trajectory_file, top=topology_file)
    trajectory[start:end].save_xtc(trajectory_file[:-4] + '_clipped.xtc')


def remove_solvent(trajectory_file, topology_file):
    my_time = time.time()
    # Load only frame 1 to select the protein atoms
    print('Loading protein atoms...')
    atoms = md.load_xtc(trajectory_file, top=topology_file, frame=1).topology.select('protein')
    print('Loading main trajectory...')
    trajectory = md.load_xtc(trajectory_file, top=topology_file, atom_indices=atoms)
    print('Removing solvent...')
    trajectory.remove_solvent(inplace=True)
    print('Writing trajectory to disk...')
    trajectory.save_xtc(trajectory_file[:-4] + '_no_solvent.xtc')
    print('Writing topology to disk...')
    trajectory[0].save_gro(trajectory_file[:-4] + '_no_solvent.gro') # Save frame 1 to a gro file
    print('Done. Time elapsed: ' + str(time.time()-my_time) + ' seconds.')


def main():
    trajectory_file = '../output/no-peptide/sh2b1-trjconv.xtc'
    topology_file = '../output/no-peptide/sh2b1.gro'
    # start = 0
    # end = 10
    # clip(trajectory_file, topology_file, start, end)
    remove_solvent(trajectory_file, topology_file)

if __name__ == "__main__":
    """
        summary: Clip a trajectory down to fewer frames and atoms.
    """
    main()