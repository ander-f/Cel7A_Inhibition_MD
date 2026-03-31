#LOADING SOFTWARES
using PDBTools, ComplexMixtures, DelimitedFiles

#READING PDB AND SELECTING ATOMS
atoms = readPDB("./solvated.pdb")
protein = select(atoms,"index >= 1 and index <= 537")
water = select(atoms,"resname TIP3")

solute = ComplexMixtures.AtomSelection(protein,nmols=1)
solvent = ComplexMixtures.AtomSelection(water,natomspermol=3)

#READING TRAJECTORY
trajectory = ComplexMixtures.Trajectory("./prod.xtc",solute,solvent)

#GETTING MDDFs
options = ComplexMixtures.Options(bulk_range=(10.0, 15.0), firstframe=2001, lastframe=10000)
results = ComplexMixtures.mddf(trajectory,options)

#SAVING RESULTS
ComplexMixtures.save(results,"./results-water_protein_glycan_last80ns.json")
