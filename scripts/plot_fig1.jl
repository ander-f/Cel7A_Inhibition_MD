import Pkg
Pkg.activate(".")
Pkg.instantiate()

using DelimitedFiles, LaTeXStrings, Plots, EasyFit, ColorSchemes, DataFrames, CSV, JSON3, Formatting, PDBTools, ComplexMixtures

cosolvent = "glucose"
protein = "cbm"

function parse_array_column_to_float_vectors!(df::DataFrame, column_name::Symbol)
  df[!, column_name] = [
      typeof(x) <: AbstractString ? 
          Float64.(collect(JSON3.read(strip(x, ['"', ' '])))) : 
          Float64.(collect(x))
      for x in df[!, column_name]
  ]
end

# Load data
cm = occursin("cbm", protein) ? 
     CSV.read("./../analyses/average_MDDF_KBI_par_$(cosolvent)_$(protein)_glycan_rep1_20_last80ns_end.csv", DataFrame)    :
     CSV.read("./../analyses/average_MDDF_KBI_par_$(cosolvent)_$(protein)_glycan_rep1_5_end.csv", DataFrame)

parse_array_column_to_float_vectors!(cm, :Distance)
parse_array_column_to_float_vectors!(cm, :Avg_MDDF_wat)
parse_array_column_to_float_vectors!(cm, :Err_MDDF_wat)
parse_array_column_to_float_vectors!(cm, :Avg_KBI_wat)
parse_array_column_to_float_vectors!(cm, :Err_KBI_wat)
parse_array_column_to_float_vectors!(cm, :Avg_MDDF_cos)
parse_array_column_to_float_vectors!(cm, :Err_MDDF_cos)
parse_array_column_to_float_vectors!(cm, :Avg_KBI_cos)
parse_array_column_to_float_vectors!(cm, :Err_KBI_cos)

# Print the \Gamma (cm.Avg_par) +/- error (cm.Err_par) for each concentration
for i in 1:length(cm.Avg_par)
    println("Concentration: $(round(cm.Avg_bulk_conc[i], digits=2)) M, Γ: $(round(cm.Avg_par[i], digits=3)) ± $(round(cm.Err_par[i], digits=3))")
end

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
    top_margin=margin, 
    left_margin=2*margin, 
    right_margin=margin,
    alpha=0.8,
    minorticks=false)

scalefontsizes(); scalefontsizes(1.3)

plot(
    layout = @layout([
        [ a [ b ; c ] ; [ c d ] e ] a{0.1w}
        s{0.2h}
    ])
)


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
    [ #get(ColorSchemes.brg,(8)/divide),
      #get(ColorSchemes.brg,(7)/divide),
      #get(ColorSchemes.brg,(6)/divide),
      #get(ColorSchemes.brg,(5)/divide),
      get(ColorSchemes.brg,(4)/divide),
      get(ColorSchemes.brg,(3)/divide),
      #get(ColorSchemes.brg,(2)/divide),
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
    plot!(xaxis = ("r/Å",0:2:8),ylabel=L"\mathrm{g^{md}_{pw \ or \ pc}}",
          subplot=sp)
    plot!(cm.Distance[1],(ma(cm.Avg_MDDF_cos[i])),
          color=c[i],
          xlim=[0-0.2,8],
          ylim=occursin("cbm", protein) ? [-0.1,3.6] : [-0.1,3.6],
          subplot=sp)

    plot!(cm.Distance[1],(ma(cm.Avg_MDDF_wat[i])),
          color=c[i],
          xlim=[0-0.2,8],
          subplot=sp)

  end

sp=2

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

sp=3

  for i in 1:length(conc)

    # Water
    plot!(xlabel=L"\mathrm{r/Å}",ylabel=L"\mathrm{g^{md}_{pw} \ (r)}",subplot=sp)
    plot!(cm.Distance[1],(ma(cm.Avg_MDDF_wat[i])),
          yticks=range(0.8, 2., length=3),
          color=c[i],
          xlim=[1.5,3.],
          ylim=occursin("cbm", protein) ? [0.8,2.2] : [0.8,2.5],
          subplot=sp
          )

  end

sp=4

  for i in 1:length(conc)
    
    # Water
    plot!(xaxis = ("r/Å",0:2:12),
                  ylabel=L"\mathrm{G_{pw \ or \ pc}} / \mathrm{L\ mol^{-1}}",
                  subplot=sp)
    plot!(cm.Distance[1],(ma(cm.Avg_KBI_wat[i] ./ 1000)),
          color=c[i], xlim=(0-0.2,12),
          bottom_margin= -4Plots.mm,
          right_margin=-6*margin,
          ylim=occursin("cbm", protein) ? [-7.5,0.5] : [-70,1.0],
          subplot=sp)

  end

sp=5
  for i in 1:length(conc)
    
    # Cosolvent
    plot!(xaxis = ("r/Å",0:2:12),
                   ylabel="",
                   subplot=sp)
    plot!(cm.Distance[1],(ma(cm.Avg_KBI_cos[i] ./ 1000)),
          color=c[i],xlim=[0-0.2,12],
          ylim=occursin("cbm", protein) ? [-7.5,0.5] : [-70,1.0],
          subplot=sp)

  end

sp=6

ticks = string.(round.(cm.Avg_bulk_conc, digits=1))

divide=8
df = occursin("cbm", protein) ? 
     DataFrame(data=cm.Avg_par, err=cm.Err_par,
              cs=[get(ColorSchemes.brg,(8)/divide),
                  get(ColorSchemes.brg,(7)/divide),
                  get(ColorSchemes.brg,(6)/divide),
                  get(ColorSchemes.brg,(5)/divide),
                  get(ColorSchemes.brg,(4)/divide),
                  get(ColorSchemes.brg,(3)/divide),
                  get(ColorSchemes.brg,(2)/divide),
                  get(ColorSchemes.brg,(1)/divide)
      ] ) :
     DataFrame(data=cm.Avg_par, err=cm.Err_par,
              cs=[get(ColorSchemes.brg,(8)/divide),
                  get(ColorSchemes.brg,(4.5)/divide),
                  get(ColorSchemes.brg,(1)/divide)
      ] )

bars=occursin("cbm", protein) ?  collect(1:8) : collect(1:3)

  bar!((bars)', df.data', yerr=df.err',
              xlabel=latexstring("\$\\mathrm{$(cosolvent) \\ conc. \\ (mol \\ L^{-1})}\$"),
              ylabel=latexstring("\$\\mathrm{\\Gamma_{$(cosolvent)}}\$"),
              ylim = (occursin("cbm", protein) ? (-0.35,1.15) : (-29,0.5)),
              label="",
              bottom_margin=-8Plots.mm,
	            minorticks=false,
              xticks=(1:1:8,ticks),
              xtickfontsize=10,
              subplot=sp,
              linecolor=:black,#:match,
              lw=0.6,
              color=[ [i] for i in df.cs]')

sp=7

  for i in 1:length(ticks)
      plot!([1],[1],
  	        xlabel=L"",
            ylabel="",
  	        border=:none,
            lw=4,
            label="$(ticks[i])",
            subplot=sp,color=c[i])
  end

sp=8
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
results = load("./../analyses/results-$(cosolvent)_1.5M_cbm_last80ns.json")
# columns equal to the number of residues
rescontrib = zeros(length(results.mddf), length(residues))

# Each column is then filled up with the contributions of each residue
for (ires, residue) in pairs(residues)
    rescontrib[:, ires] .= contributions(results, SoluteGroup(residue))
end

    clims = [  minimum(rescontrib),
               maximum(rescontrib) ]

    println(clims)

    # Plot only for distances within 1.5 and 3.5:
    idmin = findfirst(d -> d > 1.5, results.d)
    idmax = findfirst(d -> d > 3.5, results.d)

    contourf!(irange, results.d[idmin:idmax], rescontrib[idmin:idmax, irange],
              clims=(0.0, 0.3), # limit of urea, just to test
              color=cgrad(:RdBu,rev=true)[0.5:0.01:1.0],
              linewidth=0.1, linecolor=:black,
              colorbar=true,#:none, 
              levels=10,
              xlabel="Residue", ylabel=L"\mathrm{r/\AA}",
              xticks=(irange, labels[irange]), xrotation=90,
              xtickfont=font(8, plot_font),
              xlim=(1,36),
              ylim=(1.5,3.5),
              #bottommargin=30Plots.mm,
              subplot=sp
            )

plot!(size=(910,754))

savefig("./../figures/plot_fig1.svg")
