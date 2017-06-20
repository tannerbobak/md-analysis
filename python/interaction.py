import numpy as np
import mdtraj as md


def interaction_matrix(traj, residues, cutoff):
    interaction_numbers = np.empty((traj.n_frames, traj.n_frames))
    inter_numb_pair = np.empty((traj.n_frames, traj.n_frames))
    for i in range(traj.n_frames):
        for j in range(traj.n_frames):
            atom_temp, res_temp = interaction_residues_pool(traj[i], traj[j], residues, cutoff)
            interaction_numbers[i, j] = len(res_temp)
            atom_sum, res_sum = interaction_residues_pair(traj[i], traj[j], residues, cutoff)
            inter_numb_pair[i, j] = res_sum

    return interaction_numbers, inter_numb_pair


def interaction_residues_pool(frame1, frame2, residues, cutoff):
    """
    summary:
        compares frame1 and frame2 to determine how the chemical interactions present in each frame are different.
    input:
        frame1 - first frame to be compared
        frame2 - second frame to be compared
        residues - index of residues to be used in analysis
    output:
        inter_atom - number of atoms interacting
        inter_res - number of residues interacting
    """

    # determining atoms to find neighbors for
    atoms = None
    for res in residues:
        string = "resid " + str(res)
        a_temp = frame1.topology.select(string)
        if atoms is None:
            atoms = np.array(a_temp)
        else:
            atoms = np.append(atoms, a_temp)

    # determining neighbors and excluding atoms in listed residues
    nei1 = md.compute_neighbors(frame1, cutoff, atoms, haystack_indices=None)
    nei2 = md.compute_neighbors(frame2, cutoff, atoms, haystack_indices=None)
    union = (set(nei1[0]) & set(nei2[0]))
    exclude = list(union - set(atoms))
    inter_atom = exclude

    # determining residues of interacting atoms
    inter_res = None
    for a in exclude:
        inter_res_temp = frame1.topology.atom(a).residue.index
        if inter_res is None:
            inter_res = np.array(inter_res_temp)
        else:
            inter_res = np.append(inter_res, inter_res_temp)
    try:
        inter_res = list(set(inter_res))
    except:
        inter_res = []

    if inter_atom is None:
        inter_atom = []

    return inter_atom, inter_res


def interaction_residues_pair(frame1, frame2, residues, cutoff):
    """
    summary:
        compares frame1 and frame2 to determine how the chemical interactions present in each frame are different.
    input:
        frame1 - first frame to be compared
        frame2 - second frame to be compared
        residues - index of residues to be used in analysis
    output:

    """

    # determining atoms for exclude
    atoms_exlcude = None
    for res in residues:
        string = "resid " + str(res)
        a_temp = frame1.topology.select(string)
        if atoms_exlcude is None:
            atoms_exlcude = np.array(a_temp)
        else:
            atoms_exlcude = np.append(atoms_exlcude, a_temp)

    # determining atoms to find neighbors for
    sum_atoms = 0
    sum_residues = 0
    sum_atom_inter = 0
    sum_residue_inter = 0
    for res in residues:
        atoms = None
        string = "resid " + str(res)
        a_temp = frame1.topology.select(string)
        if atoms is None:
            atoms = np.array(a_temp)

        # determining neighbors and excluding atoms in listed residues
        nei1 = md.compute_neighbors(frame1, cutoff, atoms, haystack_indices=None)
        nei2 = md.compute_neighbors(frame2, cutoff, atoms, haystack_indices=None)
        intersection = list(set(nei1[0]) & set(nei2[0]))
        exclude = list(set(intersection) - set(atoms_exlcude))
        inter_atom = exclude

        # determining potential residue interactions
        union_temp = list(set(nei1[0]) | set(nei2[0]))
        union = list(set(union_temp) - set(atoms_exlcude))
        inter_res_pot = None
        for a in union:
            inter_res_temp = frame1.topology.atom(a).residue.index
            if inter_res_pot is None:
                inter_res_pot = np.array(inter_res_temp)
            else:
                inter_res_pot = np.append(inter_res_pot, inter_res_temp)
        if inter_res_pot is None:
            inter_res_pot = 0
        else:
            try:
                inter_res_pot = len(list(set(inter_res_pot)))
            except:
                inter_res_pot = 1
        sum_res_inter = inter_res_pot

        # determining residues of interacting atoms
        inter_res = None
        for a in exclude:
            inter_res_temp = frame1.topology.atom(a).residue.index
            if inter_res is None:
                inter_res = np.array(inter_res_temp)
            else:
                inter_res = np.append(inter_res, inter_res_temp)

        # determining the number of interactions
        if inter_res is None:
            inter_res = 0
        else:
            try:
                inter_res = len(list(set(inter_res)))
            except:
                inter_res = 1

        if inter_atom is None:
            inter_atom = 0
        else:
            inter_atom = len(list(inter_atom))

        # summing the intersection
        sum_atoms += inter_atom
        sum_residues += inter_res

        # summing the number of potential residue interactions
        sum_residue_inter += sum_res_inter

    try:
        percent_res = 1 - float(sum_residues) / sum_residue_inter
    except ZeroDivisionError:
        percent_res = 1

    return sum_atoms, percent_res


def main():
    """used for testing the algorithm"""

    # files and input parameters
    traj_file = 'PLCpep7.xtc'
    gro_file = 'PLCpep7.gro'
    peptide_residues = [102, 103, 104, 105, 106, 107, 108, 109]
    cutoff = 0.25  # nm

    trj = md.load(traj_file, top=gro_file)
    traj = trj[0:10]
    traj.remove_solvent(inplace=True)
    traj.superpose(traj, frame=0)
    frame1 = traj[1]
    frame2 = traj[-1]

    # running analysis
    interaction_residues_pool(frame1, frame2, peptide_residues, cutoff)
    interaction_numbers, inter_numb_pair = interaction_matrix(traj, peptide_residues, cutoff)
    print(inter_numb_pair)
    distances = np.empty((traj.n_frames, traj.n_frames))
    for i in range(traj.n_frames):
        distances[i] = md.rmsd(traj, traj, i)
    print('Max pairwise rmsd: %f nm' % np.max(distances))
    print(distances)