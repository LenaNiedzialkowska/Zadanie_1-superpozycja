import Bio.PDB
import numpy as np
from Bio.SVDSuperimposer import SVDSuperimposer


def calc_rmsd(coords1, coords2):
    diff = coords1 - coords2
    return np.sqrt(np.mean(np.sum(diff**2, axis=1)))


pdb_parser = Bio.PDB.PDBParser(QUIET=True)

sample_structure = pdb_parser.get_structure("sample", "R1107TS081.pdb")
ref_structure = pdb_parser.get_structure("reference", "R1107_reference.pdb")

lref_atoms = ref_structure.get_atoms()
lsample_atoms = sample_structure.get_atoms()

ref_atoms = []
sample_atoms = []

for a in lref_atoms:
    ref_atoms.append(a)
for a in lsample_atoms:
    sample_atoms.append(a)

ref_models = [ref_atoms[0:1461], ref_atoms[0:1464], ref_atoms[0:1461], ref_atoms[0:1461]]
sample_models = [sample_atoms[3:1464], sample_atoms[1467:2928], sample_atoms[2931:4392], sample_atoms[4395:5856]]

ref_model = ref_structure[0]
sample_model = sample_structure[0]

i = 0
while i < 4:
    # print(ref_models[i])
    # print(sample_models[i])

    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_models[i], sample_models[i])
    super_imposer.apply(sample_model.get_atoms())

    aligned_coords_ref = np.array(([atom.get_coord() for atom in ref_models[i]]))
    aligned_coords_sample = np.array(([atom.get_coord() for atom in sample_models[i]]))

    sup = SVDSuperimposer()
    sup.set(aligned_coords_ref, aligned_coords_sample)
    sup.run()
    ref_coords = sup.reference_coords
    transformed = sup.get_transformed()
    rmsd = calc_rmsd(transformed, ref_coords)
    print("RMSD dla modelu {0}: ".format(i+1), rmsd)

    i += 1

