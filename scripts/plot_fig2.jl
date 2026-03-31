import Pkg
Pkg.activate(".")
Pkg.instantiate()

using DelimitedFiles, LaTeXStrings, Plots, EasyFit, ColorSchemes, ComplexMixtures, PDBTools, Formatting, JSON3

conc = "1.5"
chain = "cbm"
dir="./.."

# Will use moving averages for more pretty graphs
ma(data) = movingaverage(data,10).x

# Plot defaults
plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=1.5,
    framestyle=:box,
    label=nothing,
    grid=false,
    dpi=300,
    legend=:topleft,
    top_margin=1Plots.mm,
    bottom_margin=3.4Plots.mm,
    left_margin=3.6Plots.mm, 
    #right_margin=margin,
    xlim=[0,8.0],
    ylim=occursin("cbm", chain) ? [-0.08,3.5] : [-0.08,3.8],
    minorticks=false)

plot(layout = (1,3))
plot(
    layout = @layout([
        [ a b c ] ; [ d{0.6h} ; e{0.6h} ]    ])
)

divide=8
c = [ get(ColorSchemes.brg,(8)/divide),
      get(ColorSchemes.brg,(7)/divide),
      get(ColorSchemes.brg,(6)/divide),
      get(ColorSchemes.brg,(5)/divide),
      get(ColorSchemes.brg,(4)/divide),
      get(ColorSchemes.brg,(3)/divide),
      get(ColorSchemes.brg,(2)/divide),
      get(ColorSchemes.brg,(1)/divide),
      ]

######## Chain in the presence of the sugars ########
# Calculate the protein group contributions
sp=1
cosolvent="glucose"
atoms = readPDB("$dir/pdbs/solvated_$(cosolvent)0.5M_$(chain).pdb")
cm = occursin("cbm", chain) ? ComplexMixtures.load("$dir/analyses/results-$(cosolvent)_$(conc)M_cbm_last80ns.json") :
                              ComplexMixtures.load("$dir/analyses/results-$(cosolvent)_$(conc)M_cd.json")
contrib = occursin("cbm", chain) ? contributions(cm, SoluteGroup(PDBTools.select(atoms,"index >= 1 and index <= 537"))) :
                                   contributions(cm, SoluteGroup(PDBTools.select(atoms,"index >= 1 and index <= 6615")))
plot!(cm.d,(ma(contrib)),
      xaxis=("r/Å",0:2:8.0),ylabel=L"\mathsf{g_{ps} \ (r)}",
      color=c[7],
      label="Total",
      legend=:topright,
      titlefontsize=10,
      subplot=sp
      )

polar = contributions(cm, SoluteGroup(PDBTools.select(atoms,"protein and polar")))
cn_polar = contributions(cm, SoluteGroup(PDBTools.select(atoms,"protein and polar")); type=:coordination_number)
plot!(cm.d,ma(polar),
label="Polar",
color=:green,
subplot=sp)

nonpolar = contributions(cm, SoluteGroup(PDBTools.select(atoms,"protein and not polar")))
cn_nonpolar = contributions(cm, SoluteGroup(PDBTools.select(atoms,"protein and not polar")); type=:coordination_number)
plot!(cm.d,ma(nonpolar),
label="Nonpolar",
color=:goldenrod,
subplot=sp)
hline!([1.0],linestyle=:dash,color=:black,alpha=0.6,subplot=sp)

###
glycans = occursin("cbm", chain) ? contributions(cm, SoluteGroup(PDBTools.select(atoms,"index >= 494 and index <= 537"))) :
                                   contributions(cm, SoluteGroup(PDBTools.select(atoms,"index >= 6136 and index <= 6615")))
plot!(cm.d,ma(glycans),
label="Glycans",
color=:black,
subplot=sp)

sp=2
cosolvent="mannose"
atoms = readPDB("$dir/pdbs/solvated_$(cosolvent)0.5M_$(chain).pdb")
cm = occursin("cbm", chain) ? ComplexMixtures.load("$dir/analyses/results-$(cosolvent)_$(conc)M_cbm_last80ns.json") :
                              ComplexMixtures.load("$dir/analyses/results-$(cosolvent)_$(conc)M_cd.json")
contrib = occursin("cbm", chain) ? contributions(cm, SoluteGroup(PDBTools.select(atoms,"index >= 1 and index <= 537"))) :
                                   contributions(cm, SoluteGroup(PDBTools.select(atoms,"index >= 1 and index <= 6615")))

plot!(cm.d,(ma(contrib)),
      xaxis=("r/Å",0:2:8.0),ylabel=L"\mathsf{g_{ps} \ (r)}",
      color=c[7],
      label="Total",
      legend=:topright,
      titlefontsize=10,
      subplot=sp
      )

polar = contributions(cm, SoluteGroup(PDBTools.select(atoms,"protein and polar")))
cn_polar = contributions(cm, SoluteGroup(PDBTools.select(atoms,"protein and polar")); type=:coordination_number)
plot!(cm.d,ma(polar),
label="Polar",
color=:green,
subplot=sp)

nonpolar = contributions(cm, SoluteGroup(PDBTools.select(atoms,"protein and not polar")))
cn_nonpolar = contributions(cm, SoluteGroup(PDBTools.select(atoms,"protein and not polar")); type=:coordination_number)
plot!(cm.d,ma(nonpolar),
label="Nonpolar",
color=:goldenrod,
subplot=sp)
hline!([1.0],linestyle=:dash,color=:black,alpha=0.6,subplot=sp)

###
glycans = occursin("cbm", chain) ? contributions(cm, SoluteGroup(PDBTools.select(atoms,"index >= 494 and index <= 537"))) :
                                   contributions(cm, SoluteGroup(PDBTools.select(atoms,"index >= 6136 and index <= 6615")))
plot!(cm.d,ma(glycans),
label="Glycans",
color=:black,
subplot=sp)

sp=3
cosolvent="xylose"
atoms = readPDB("$dir/pdbs/solvated_$(cosolvent)0.5M_$(chain).pdb")
cm = occursin("cbm", chain) ? ComplexMixtures.load("$dir/analyses/results-$(cosolvent)_$(conc)M_cbm_last80ns.json") :
                              ComplexMixtures.load("$dir/analyses/results-$(cosolvent)_$(conc)M_cd.json")
contrib = occursin("cbm", chain) ? contributions(cm, SoluteGroup(PDBTools.select(atoms,"index >= 1 and index <= 537"))) :
                                   contributions(cm, SoluteGroup(PDBTools.select(atoms,"index >= 1 and index <= 6615")))

plot!(cm.d,(ma(contrib)),
      xaxis=("r/Å",0:2:8.0),ylabel=L"\mathsf{g_{ps} \ (r)}",
      color=c[7],
      label="Total",
      legend=:topright,
      titlefontsize=10,
      subplot=sp
      )

polar = contributions(cm, SoluteGroup(PDBTools.select(atoms,"protein and polar")))
cn_polar = contributions(cm, SoluteGroup(PDBTools.select(atoms,"protein and polar")); type=:coordination_number)
plot!(cm.d,ma(polar),
label="Polar",
color=:green,
subplot=sp)

nonpolar = contributions(cm, SoluteGroup(PDBTools.select(atoms,"protein and not polar")))
cn_nonpolar = contributions(cm, SoluteGroup(PDBTools.select(atoms,"protein and not polar")); type=:coordination_number)
plot!(cm.d,ma(nonpolar),
label="Nonpolar",
color=:goldenrod,
subplot=sp)
hline!([1.0],linestyle=:dash,color=:black,alpha=0.6,subplot=sp)

###
glycans = occursin("cbm", chain) ? contributions(cm, SoluteGroup(PDBTools.select(atoms,"index >= 494 and index <= 537"))) :
                                   contributions(cm, SoluteGroup(PDBTools.select(atoms,"index >= 6136 and index <= 6615")))
plot!(cm.d,ma(glycans),
label="Glycans",
color=:black,
subplot=sp)

sp=4
# PDB file of the system simulated
pdb = readPDB("./../pdbs/solvated_glucose0.5M_cbm.pdb")

# Inform which is the solute
protein = PDBTools.select(pdb, "protein")
solute = AtomSelection(protein, nmols=1)

# Obtain pretty labels for the residues in the x-axis
residues = collect(eachresidue(protein))
labels = PDBTools.oneletter.(resname.(residues)).*format.(resnum.(residues))

# We will plot only the range of ARGS[1] of residues
irange_arg = "1:36"
start, stop = parse.(Int, split(irange_arg, ":"))
irange = start:stop
# Use irange as UnitRange{Int64}
println(irange)
# Load the results
    results = load("./../analyses/results-glucose_1.5M_cbm_last80ns.json")
    # columns equal to the number of residues
    rescontrib = zeros(length(results.mddf), length(residues))

    # Each column is then filled up with the contributions of each residue
    for (ires, residue) in pairs(residues)
        rescontrib[:, ires] .= contributions(results, SoluteGroup(residue))
    end

    rescontrib_ensemble1 = copy(rescontrib)

    # Load the results of the ensemble
    results = load("./../analyses/results-mannose_1.5M_cbm_last80ns.json")
    # columns equal to the number of residues
    rescontrib = zeros(length(results.mddf), length(residues))

    # Each column is then filled up with the contributions of each residue
    for (ires, residue) in pairs(residues)
        rescontrib[:, ires] .= contributions(results, SoluteGroup(residue))
    end

    rescontrib_ensemble2 = copy(rescontrib)

    rescontrib = rescontrib_ensemble1 .- rescontrib_ensemble2

    clims = [  minimum(rescontrib),
               maximum(rescontrib) ]

    println(clims)

    # Plot only for distances within 1.5 and 3.5:
    idmin = findfirst(d -> d > 1.5, results.d)
    idmax = findfirst(d -> d > 3.5, results.d)

    contourf!(irange, results.d[idmin:idmax], rescontrib[idmin:idmax, irange],
              clims=(-0.0405, 0.0405),
              color=cgrad(:RdBu,rev=true),
              linewidth=0.1, linecolor=:black,
              colorbar=true,#:none, 
              levels=6,
              xlabel="", ylabel=L"\mathrm{r/\AA}",
              xticks=nothing,
              xlim=(1,36),
              ylim=(1.5,3.5),
              subplot=sp
            )

sp=5
# Load the results
    results = load("./../analyses/results-glucose_1.5M_cbm_last80ns.json")
    # columns equal to the number of residues
    rescontrib = zeros(length(results.mddf), length(residues))

    # Each column is then filled up with the contributions of each residue
    for (ires, residue) in pairs(residues)
        rescontrib[:, ires] .= contributions(results, SoluteGroup(residue))
    end

    rescontrib_ensemble1 = copy(rescontrib)

    # Load the results of the ensemble
    results = load("./../analyses/results-xylose_1.5M_cbm_last80ns.json")
    # columns equal to the number of residues
    rescontrib = zeros(length(results.mddf), length(residues))

    # Each column is then filled up with the contributions of each residue
    for (ires, residue) in pairs(residues)
        rescontrib[:, ires] .= contributions(results, SoluteGroup(residue))
    end

    rescontrib_ensemble2 = copy(rescontrib)

    rescontrib = rescontrib_ensemble1 .- rescontrib_ensemble2

    clims = [  minimum(rescontrib),
               maximum(rescontrib) ]

    println(clims)

    # Plot only for distances within 1.5 and 3.5:
    idmin = findfirst(d -> d > 1.5, results.d)
    idmax = findfirst(d -> d > 3.5, results.d)

    contourf!(irange, results.d[idmin:idmax], rescontrib[idmin:idmax, irange],
              clims=(-0.0405, 0.0405),
              color=cgrad(:RdBu,rev=true),
              linewidth=0.1, linecolor=:black,
              colorbar=true,#:none, 
              levels=6,
              xlabel="Residue", ylabel=L"\mathrm{r/\AA}",
              xticks=(irange, labels[irange]), xrotation=90,
              xtickfont=font(8, plot_font),
              xlim=(1,36),
              ylim=(1.5,3.5),
              bottommargin=8Plots.mm,
              subplot=sp
            )

for sp in 1:3
  annotate!(4.3, 1.7, text(latexstring("\$\\mathrm{$(uppercase(chain)) \\ in }\$"), :left, 9), subplot=sp)
  annotate!(3.6, 1.35, text(latexstring("\$\\mathrm{\\ $(conc) \\ M \\ $(cosolvent)}\$"), :left, 9), subplot=sp)
end

plot!(size=(675,500))

# Save figure
savefig("./../figures/plot_fig2.svg")
