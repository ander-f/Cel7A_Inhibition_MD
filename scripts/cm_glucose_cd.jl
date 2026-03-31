#LOADING SOFTWARES
using PDBTools, ComplexMixtures, DelimitedFiles

#READING PDB AND SELECTING ATOMS
atoms = readPDB("./solvated.pdb")
protein = select(atoms,"index >= 1 and index <= 6615")
glu = select(atoms,"resname BGL")

solute = ComplexMixtures.AtomSelection(protein,nmols=1)
glucose = ComplexMixtures.AtomSelection(glu,natomspermol=24)

#READING TRAJECTORY
trajectory = ComplexMixtures.Trajectory("./final_traj/prod_final.xtc",solute,glucose)

#GETTING MDDFs
options = ComplexMixtures.Options(bulk_range=(10.0, 15.0), firstframe=1, lastframe=5000)
results = ComplexMixtures.mddf(trajectory,options)

#SAVING RESULTS
ComplexMixtures.save(results,"./final_traj/results-glucose_protein_glycan.json")
