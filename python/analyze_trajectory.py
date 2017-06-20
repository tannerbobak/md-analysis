import mdtraj as md
import itertools
from interaction import interaction_residues_pool

def main():
    """ Analyze a trajectory for relaxation times. """

    # Residues of the peptide in the crystal structure that we are interested in.
    peptide_residues = range(102,115)

    # Standard capture rate is every 0.2ps
    ps_per_frame = 100
    stride_len = int(ps_per_frame/0.2)
    cutoff = 0.25  # nm

    """
    Where to output the distance matrix.
    """
    output_path = '../output/peptide/distanceData_' + str(cutoff*1000) + 'pm.txt'

    """
    The output files from GROMACS
    """
    trajectory_file = '../output/peptide/sh2b1-trjconv_no_solvent.xtc'
    topology_file = '../output/peptide/sh2b1-trjconv_no_solvent.gro'

    """
    These files should have their first frame contain the bound crystal structure so 
    relevant residues near the peptide can be identified.
    """
    bound_crystal_trajectory_file = '../output/peptide/sh2b1-trjconv_no_solvent.xtc'
    bound_crystal_topology_file = '../output/peptide/sh2b1-trjconv_no_solvent.gro'
    # bound_crystal_trajectory_file = trajectory_file
    # bound_crystal_topology_file = topology_file

    # Just need the first frame
    print('Reading trajectories...')
    bound_crystal_trajectory = md.load_frame(bound_crystal_trajectory_file, index=0, top=bound_crystal_topology_file)
    trajectory = md.load(trajectory_file, top=topology_file, stride=stride_len) #stride*2fs is the time between frames here

    bound_crystal_trajectory.remove_solvent(inplace=True)
    trajectory.remove_solvent(inplace=True)
    trajectory.superpose(trajectory, frame=0)

    print('Determining interacting residues...')
    inter_atoms, inter_res = interaction_residues_pool(bound_crystal_trajectory[0], bound_crystal_trajectory[0], peptide_residues, cutoff)

    """
    Few options:
        1. Cut out just the interacting residues and measure how the RMSD between that 
           residue and the others evolves with time. ***
           - use every 5th frame
        2. Measure the rate of change of the RMSD of the total protein w.r.t time
    """

    print('Computing contacts...')
    contacts_array = list(itertools.product(inter_res, inter_res))  # There are n^2 possible contacts
    distances, residue_pairs = md.compute_contacts(trajectory, contacts=contacts_array,
                                                   scheme='closest-heavy', ignore_nonprotein=True)

    contacts_array = residue_pairs
    contact_map = md.geometry.squareform(distances, contacts_array)

    # Slice the array so only every nth element is used
    # n = 50
    # contact_map = contact_map[::n]

    # Plot a table
    n_residues = len(inter_res)
    n_frames = contact_map.shape[0]

    print('Writing data file...')
    outfile = open(output_path, 'w')

    # Write the residues that are in this file.
    header = 'Interacting Residues (Count: ' + str(n_residues) + '): '
    for i in range(0, n_residues):
        header = header + str(inter_res[i]) + ','
    header = header[:-1] # Truncate trailing comma.
    outfile.write(header + '\n')

    for frame_number in range(0, n_frames):
        # data = np.zeros((n_residues, n_residues))
        for i in range(0, n_residues):
            row = ''
            for j in range(0, n_residues):
                # data[i, j] = contact_map[frame_number, inter_res[i], inter_res[j]]
                row = row + str(contact_map[frame_number, inter_res[i], inter_res[j]]) + ','
            row = row[:-1] # Truncate trailing comma.
            outfile.write(row + '\n')

if __name__ == "__main__":
    main()
