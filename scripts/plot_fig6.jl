import Pkg
Pkg.activate(".")
Pkg.instantiate()

using DelimitedFiles, LaTeXStrings, Plots, EasyFit, ColorSchemes, DataFrames, CSV, JSON3, PDBTools, ComplexMixtures, Formatting

function parse_array_column_to_float_vectors!(df::DataFrame, column_name::Symbol)
  df[!, column_name] = [
      typeof(x) <: AbstractString ? 
          Float64.(collect(JSON3.read(strip(x, ['"', ' '])))) : 
          Float64.(collect(x))
      for x in df[!, column_name]
  ]
end

# Load data
protein="cd"
cosolvent="glucose"
cm = occursin("cbm", protein) ? 
     CSV.read("./../analyses/average_MDDF_KBI_par_$(cosolvent)_$(protein)_glycan_rep1_20_last80ns_end.csv", DataFrame)    :
     CSV.read("./../analyses/average_MDDF_KBI_par_$(cosolvent)_$(protein)_glycan_rep1_5_end.csv", DataFrame)

parse_array_column_to_float_vectors!(cm, :Distance)
parse_array_column_to_float_vectors!(cm, :Avg_MDDF_wat)
parse_array_column_to_float_vectors!(cm, :Err_MDDF_wat)
parse_array_column_to_float_vectors!(cm, :Avg_MDDF_cos)
parse_array_column_to_float_vectors!(cm, :Err_MDDF_cos)

# Will use moving averages for more pretty graphs
ma(data) = movingaverage(data,10).x

# Plot defaults
margin = 1.5Plots.mm
plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=1.5,
    framestyle=:box,
    label=nothing,
    grid=false,
    dpi=300,
    legend=:topleft,
    bottom_margin=margin,
    top_margin=4margin, 
    left_margin=1.8*margin, 
    right_margin=4margin,
    alpha=0.8,
    minorticks=false)

scalefontsizes(); scalefontsizes(1.3)

plot(layout=@layout [ [ a ; b ] [ c ; d ] [ e ; f ] ; [ g{0.3h} ; h{0.3h} ; i{0.3h} ] ])

divide=8
c = occursin("cbm", protein) ? 
    [ get(ColorSchemes.brg,(8)/divide),
      get(ColorSchemes.brg,(7)/divide),
      get(ColorSchemes.brg,(6)/divide),
      get(ColorSchemes.brg,(5)/divide),
      get(ColorSchemes.brg,(4)/divide),
      get(ColorSchemes.brg,(3)/divide),
      get(ColorSchemes.brg,(2)/divide),
      get(ColorSchemes.brg,(1)/divide),
      ] :
    [ get(ColorSchemes.brg,(8)/divide),
      get(ColorSchemes.brg,(4.5)/divide),
      get(ColorSchemes.brg,(1)/divide),
      ]

conc = occursin("cbm", protein) ? 
       [ "0.1",
         "0.2",
         "0.3",
         "0.4",
         "0.5",
         "1.0",
         "1.5",
         "2.0"
       ] :
       [ #"0.3",
         #"0.4",
         "0.5",
         "1.0",
         #"1.5",
         "2.0"
       ]
       
sp=1

  for i in 1:length(conc)
    
    # Cosolvent  
    plot!(xlabel="",ylabel=L"\mathrm{g^{md}_{pc} \ (r)}",
          bottom_margin = -4Plots.mm,
          subplot=sp)
    
    plot!(cm.Distance[1],(ma(cm.Avg_MDDF_cos[i])),
          yticks=range(2, 4, length=3),
          color=c[i],
          xlim=[1.5,3.],
          ylim=[2.5,3.6],
          subplot=sp,
          )

  end
annotate!(2.5, 3.4, text(latexstring("\$\\mathrm{Glucose}\$"), :center, 16), subplot=sp)

sp=2

  for i in 1:length(conc)

    # Water
    plot!(xlabel=L"\mathrm{r/Å}",ylabel=L"\mathrm{g^{md}_{pw} \ (r)}",subplot=sp)
    plot!(cm.Distance[1],(ma(cm.Avg_MDDF_wat[i])),
          yticks=range(0.8, 2., length=3),
          color=c[i],
          xlim=[1.5,3.],
          ylim=occursin("cbm", protein) ? [0.8,2.2] : [0.8,2.5],
          #bottom_margin = -20Plots.mm,
          subplot=sp
          )

  end

annotate!(2.5, 2.3, text(latexstring("\$\\mathrm{Water}\$"), :center, 16), subplot=sp)
annotate!(2.5, 2.0, text(latexstring("\$\\mathrm{in \\ Glucose}\$"), :center, 16), subplot=sp)

#################
cosolvent="mannose"
cm = occursin("cbm", protein) ? 
     CSV.read("./../analyses/average_MDDF_KBI_par_$(cosolvent)_$(protein)_glycan_rep1_20_last80ns_end.csv", DataFrame)    :
     CSV.read("./../analyses/average_MDDF_KBI_par_$(cosolvent)_$(protein)_glycan_rep1_5_end.csv", DataFrame)

parse_array_column_to_float_vectors!(cm, :Distance)
parse_array_column_to_float_vectors!(cm, :Avg_MDDF_wat)
parse_array_column_to_float_vectors!(cm, :Err_MDDF_wat)
parse_array_column_to_float_vectors!(cm, :Avg_MDDF_cos)
parse_array_column_to_float_vectors!(cm, :Err_MDDF_cos)

sp=3

  for i in 1:length(conc)
    
    # Cosolvent  
    plot!(xlabel="",ylabel=L"\mathrm{g^{md}_{pc} \ (r)}",
          bottom_margin = -4Plots.mm,
          subplot=sp)
    
    plot!(cm.Distance[1],(ma(cm.Avg_MDDF_cos[i])),
          yticks=range(2, 4, length=3),
          color=c[i],
          xlim=[1.5,3.],
          ylim=[2.5,3.6],
          subplot=sp,
          )

  end

annotate!(2.5, 3.4, text(latexstring("\$\\mathrm{Mannose}\$"), :center, 16), subplot=sp)

sp=4

  for i in 1:length(conc)

    # Water
    plot!(xlabel=L"\mathrm{r/Å}",ylabel=L"\mathrm{g^{md}_{pw} \ (r)}",subplot=sp)
    plot!(cm.Distance[1],(ma(cm.Avg_MDDF_wat[i])),
          yticks=range(0.8, 2., length=3),
          color=c[i],
          xlim=[1.5,3.],
          ylim=occursin("cbm", protein) ? [0.8,2.2] : [0.8,2.5],
          #bottom_margin = -20Plots.mm,
          subplot=sp
          )

  end

annotate!(2.5, 2.3, text(latexstring("\$\\mathrm{Water}\$"), :center, 16), subplot=sp)
annotate!(2.5, 2.0, text(latexstring("\$\\mathrm{in \\ Glucose}\$"), :center, 16), subplot=sp)
  
#################
cosolvent="xylose"
cm = occursin("cbm", protein) ? 
     CSV.read("./../analyses/average_MDDF_KBI_par_$(cosolvent)_$(protein)_glycan_rep1_20_last80ns_end.csv", DataFrame)    :
     CSV.read("./../analyses/average_MDDF_KBI_par_$(cosolvent)_$(protein)_glycan_rep1_5_end.csv", DataFrame)

parse_array_column_to_float_vectors!(cm, :Distance)
parse_array_column_to_float_vectors!(cm, :Avg_MDDF_wat)
parse_array_column_to_float_vectors!(cm, :Err_MDDF_wat)
parse_array_column_to_float_vectors!(cm, :Avg_MDDF_cos)
parse_array_column_to_float_vectors!(cm, :Err_MDDF_cos)

sp=5

  for i in 1:length(conc)
    
    # Cosolvent  
    plot!(xlabel="",ylabel=L"\mathrm{g^{md}_{pc} \ (r)}",
          bottom_margin = -4Plots.mm,
          subplot=sp)
    
    plot!(cm.Distance[1],(ma(cm.Avg_MDDF_cos[i])),
          yticks=range(2, 4, length=3),
          color=c[i],
          xlim=[1.5,3.],
          ylim=[2.5,3.6],
          subplot=sp,
          )

  end
annotate!(2.5, 3.4, text(latexstring("\$\\mathrm{Xylose}\$"), :center, 16), subplot=sp)

sp=6

  for i in 1:length(conc)

    # Water
    plot!(xlabel=L"\mathrm{r/Å}",ylabel=L"\mathrm{g^{md}_{pw} \ (r)}",subplot=sp)
    plot!(cm.Distance[1],(ma(cm.Avg_MDDF_wat[i])),
          yticks=range(0.8, 2., length=3),
          color=c[i],
          xlim=[1.5,3.],
          ylim=occursin("cbm", protein) ? [0.8,2.2] : [0.8,2.5],
          #bottom_margin = -20Plots.mm,
          subplot=sp
          )

  end

annotate!(2.5, 2.3, text(latexstring("\$\\mathrm{Water}\$"), :center, 16), subplot=sp)
annotate!(2.5, 2.0, text(latexstring("\$\\mathrm{in \\ Xylose}\$"), :center, 16), subplot=sp)

sp=7
# PDB file of the system simulated
pdb = readPDB("./../pdbs/solvated_glucose2.0M_cd.pdb")

# Inform which is the solute
protein = PDBTools.select(pdb, "protein")
solute = AtomSelection(protein, nmols=1)
residue_within_5A_substract = readdlm("./../pdbs/residues_at_5A_from_GLC.dat")
residue_within_5A_substract = Int.(residue_within_5A_substract[:,2])

# We will plot only the range of ARGS[1] of residues
irange_arg = "1:54"
start, stop = parse.(Int, split(irange_arg, ":"))
irange = start:stop
# Use irange as UnitRange{Int64}
println(irange)

# Obtain pretty labels for the residues in the x-axis
residues = collect(eachresidue(protein))[residue_within_5A_substract .- 1] # minus 1 because of the sequence starting at S2
labels = PDBTools.oneletter.(resname.(residues)).*format.(resnum.(residues))

# Load the results
results = load("./../analyses/results-glucose_2.0M_cd.json")
# columns equal to the number of residues
rescontrib_glucose = zeros(length(results.mddf), length(residues))
# Each column is then filled up with the contributions of each residue
for (ires, residue) in pairs(residues)
    rescontrib_glucose[:, ires] .= contributions(results, SoluteGroup(residue))
end

# Plot only for distances within 1.5 and 3.5:
idmin = findfirst(d -> d > 1.5, results.d)
idmax = findfirst(d -> d > 3.5, results.d)
heatmap!(irange, results.d[idmin:idmax], rescontrib_glucose[idmin:idmax, irange],
         #clims=(clims[1], clims[2]),
          clims=(0.0, 0.025),
          color=cgrad(:RdBu,rev=true)[0.5:0.01:1.0],
          linewidth=0.1, linecolor=:black,
          colorbar=true,#:none, 
          levels=10,
          xlabel="", ylabel=L"\mathrm{r/\AA}",
          #xticks=(irange[1:1:end], labels[irange[1:1:end]]), xrotation=90,
          #xtickfont=font(8, plot_font),
          xticks=nothing,
          subplot=sp
        )
annotate!(1, 3.8, text(latexstring("\$\\mathrm{Glucose \\ (2.0 \\ M)}\$"), :left, 12), subplot=sp)

sp=8
# Load the results
results = load("./../analyses/results-mannose_2.0M_cd.json")
# columns equal to the number of residues
rescontrib_mannose = zeros(length(results.mddf), length(residues))
# Each column is then filled up with the contributions of each residue
for (ires, residue) in pairs(residues)
    rescontrib_mannose[:, ires] .= contributions(results, SoluteGroup(residue))
end

rescontrib_diff_glu_less_man = rescontrib_glucose - rescontrib_mannose
# Plot only for distances within 1.5 and 3.5:
idmin = findfirst(d -> d > 1.5, results.d)
idmax = findfirst(d -> d > 3.5, results.d)
heatmap!(irange, results.d[idmin:idmax], rescontrib_diff_glu_less_man[idmin:idmax, irange],
          clims=(-0.006, 0.006),
          color=cgrad(:RdBu,rev=true),
          linewidth=0.1, linecolor=:black,
          colorbar=true,#:none, 
          levels=10,
          xlabel="", ylabel=L"\mathrm{r/\AA}",
          #xticks=(irange[1:1:end], labels[irange[1:1:end]]), xrotation=90,
          #xtickfont=font(8, plot_font),
          xticks=nothing,
          subplot=sp
        )

annotate!(1, 3.8, text(latexstring("\$\\mathrm{Glucose - mannose \\ (2.0 \\ M)}\$"), :left, 12), subplot=sp)

sp=9
# Load the results
results = load("./../analyses/results-xylose_2.0M_cd.json")
# columns equal to the number of residues
rescontrib_xylose = zeros(length(results.mddf), length(residues))
# Each column is then filled up with the contributions of each residue
for (ires, residue) in pairs(residues)
    rescontrib_xylose[:, ires] .= contributions(results, SoluteGroup(residue))
end

rescontrib_diff_glu_less_xyl = rescontrib_glucose - rescontrib_xylose

# Plot only for distances within 1.5 and 3.5:
idmin = findfirst(d -> d > 1.5, results.d)
idmax = findfirst(d -> d > 3.5, results.d)
heatmap!(irange, results.d[idmin:idmax], rescontrib_diff_glu_less_xyl[idmin:idmax, irange],
          clims=(-0.006, 0.006),
          color=cgrad(:RdBu,rev=true),
          linewidth=0.1, linecolor=:black,
          colorbar=true,#:none, 
          levels=10,
          xlabel="Residue", ylabel=L"\mathrm{r/\AA}",
          xticks=(irange[1:1:end], labels[irange[1:1:end]]), xrotation=90,
          xtickfont=font(8, plot_font),
          subplot=sp
        )

annotate!(1, 3.8, text(latexstring("\$\\mathrm{Glucose - xylose \\ (2.0 \\ M)}\$"), :left, 12), subplot=sp)

plot!(size=(700*1.3,780*1.3))

# Save figure
savefig("./../figures/plot_fig6.svg")
