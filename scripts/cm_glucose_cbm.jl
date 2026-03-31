
#LOADING SOFTWARES
using PDBTools, ComplexMixtures, DelimitedFiles

#READING PDB AND SELECTING ATOMS
atoms = readPDB("./solvated.pdb")
protein = select(atoms,"index >= 1 and index <= 537")
glu = select(atoms,"resname BGL")

solute = ComplexMixtures.AtomSelection(protein,nmols=1)
glucose = ComplexMixtures.AtomSelection(glu,natomspermol=24)

#READING TRAJECTORY
trajectory = ComplexMixtures.Trajectory("./prod.xtc",solute,glucose)

#GETTING MDDFs
options = ComplexMixtures.Options(bulk_range=(10.0, 15.0), firstframe=2001, lastframe=10000)
results = ComplexMixtures.mddf(trajectory,options)

#SAVING RESULTS
ComplexMixtures.save(results,"./results-glucose_protein_glycan_last80ns.json")
