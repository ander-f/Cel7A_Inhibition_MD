import Pkg
Pkg.activate(".")
Pkg.instantiate()

using DelimitedFiles, Plots, StatsPlots, LaTeXStrings, ColorSchemes, PDBTools, Formatting

asa = "wasa"

#W38: atoms from 502 to 525
#W40: atoms from 550 to 573
#W367: atoms from 5177 to 5200
#W376: atoms from 5327 to 5350

base_dir = "./.."
cosolvents = ["glucose", "mannose", "xylose"]
concentrations = ["0.5", "1.0", "2.0"]

wat_sasa = readdlm("$base_dir/analyses/average_$(asa)_by_atom_water.dat")
println("SASA in water:")
println("W38: ", sum(wat_sasa[502:525, 2]))
println("W40: ", sum(wat_sasa[550:573, 2]))
println("W367: ", sum(wat_sasa[5177:5200, 2]))
println("W376: ", sum(wat_sasa[5327:5350, 2]))

println("")

for conc in concentrations
    for cosol in cosolvents
        cos_sasa = readdlm("$base_dir/analyses/average_$(asa)_by_atom_$(cosol)_$(conc).dat")
        println("SASA in $(cosol) at $(conc) M:")
        println("W38: ", sum(cos_sasa[502:525, 2]))
        println("W40: ", sum(cos_sasa[550:573, 2]))
        println("W367: ", sum(cos_sasa[5177:5200, 2]))
        println("W376: ", sum(cos_sasa[5327:5350, 2]))
        println("")
    end
end

# Plot for the residues W38, W40, W367, W376 in water and cosolvents
res_labels = ["W38", "W40", "W367", "W376"]
indices = [502:525, 550:573, 5177:5200, 5327:5350]

# values for water
wat_vals = [ sum(wat_sasa[r, 2]) for r in indices ]

#conc = ARGS[1]
# values for cosolvents
sasa_cos_05M = Dict{String, Vector{Float64}}()
for cosol in cosolvents
    f = "$base_dir/analyses/average_$(asa)_by_atom_$(cosol)_0.5.dat"
    if !isfile(f)
        @warn "SASA file not found, skipping" file=f
        sasa_cos_05M[cosol] = fill(NaN, length(indices))
        continue
    end
    data = readdlm(f)
    sasa_cos_05M[cosol] = [ sum(data[intersect(1:size(data,1), collect(r)), 2]) for r in indices ]
end

sasa_cos_1M = Dict{String, Vector{Float64}}()
for cosol in cosolvents
    f = "$base_dir/analyses/average_$(asa)_by_atom_$(cosol)_1.0.dat"
    if !isfile(f)
        @warn "SASA file not found, skipping" file=f
        sasa_cos_1M[cosol] = fill(NaN, length(indices))
        continue
    end
    data = readdlm(f)
    sasa_cos_1M[cosol] = [ sum(data[intersect(1:size(data,1), collect(r)), 2]) for r in indices ]
end

sasa_cos_2M = Dict{String, Vector{Float64}}()
for cosol in cosolvents
    f = "$base_dir/analyses/average_$(asa)_by_atom_$(cosol)_2.0.dat"
    if !isfile(f)
        @warn "SASA file not found, skipping" file=f
        sasa_cos_2M[cosol] = fill(NaN, length(indices))
        continue
    end
    data = readdlm(f)
    sasa_cos_2M[cosol] = [ sum(data[intersect(1:size(data,1), collect(r)), 2]) for r in indices ]
end

# assemble matrix: columns = water, glucose, mannose, xylose (in that order if available)
cols = ["water"]
append!(cols, cosolvents)
mat_05M = zeros(length(indices), length(cols))
for (j, col) in enumerate(cols)
    if col == "water"
        mat_05M[:, j] = wat_vals
    else
        mat_05M[:, j] = sasa_cos_05M[col]
    end
end

mat_1M = zeros(length(indices), length(cols))
for (j, col) in enumerate(cols)
    if col == "water"
        mat_1M[:, j] = wat_vals
    else
        mat_1M[:, j] = sasa_cos_1M[col]
    end
end

mat_2M = zeros(length(indices), length(cols))
for (j, col) in enumerate(cols)
    if col == "water"
        mat_2M[:, j] = wat_vals
    else
        mat_2M[:, j] = sasa_cos_2M[col]
    end
end

# Configuration for plots
plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    framestyle=:box,
    label=nothing,
    grid=false,
    lw=1.5,
    dpi=300,
    legend=:topleft,
    xlabel="Residue",
    ylabel=latexstring("\\mathrm{" * uppercase(asa) * " \\ (\\AA^2)}"),
    #L"\mathrm{SASA \ (\AA^2)}",
    bar_width=0.7,
    leftmargin=1Plots.Measures.cm,
    topmargin=0.4Plots.Measures.cm,
    bottommargin=1Plots.Measures.cm,
    rightmargin=0.2Plots.Measures.cm,
    ylim=(600, 9800),
    markersize=5,
    markerstrokewidth=0,
    minorticks=false,
)

plot(
    layout = @layout([
        [ a{0.7h} b{0.7h} c{0.7h} ] ; d{0.4h} ]), size=(1200*0.7, 900*0.7)
)

colors = [ColorSchemes.tab10[1], ColorSchemes.tab10[2], ColorSchemes.tab10[3], ColorSchemes.tab10[4], ColorSchemes.tab10[5]]

# use grouped bar plot so legend maps correctly to columns. mat is residues x conditions (water and cosolvents)
sp=1
scatter!(res_labels, mat_05M[:,1], subplot=sp, label="Water", title="0.5 M Cosolvent", color=:black)
scatter!(res_labels, mat_05M[:,2], subplot=sp, label="Glucose", color=colors[1])
scatter!(res_labels, mat_05M[:,3], subplot=sp, label="Mannose", color=colors[2])
scatter!(res_labels, mat_05M[:,4], subplot=sp, label="Xylose", color=colors[3])

sp = 2
scatter!(res_labels, mat_1M[:,1], subplot=sp, label="Water", title="1.0 M Cosolvent", color=:black)
scatter!(res_labels, mat_1M[:,2], subplot=sp, label="Glucose", color=colors[1])
scatter!(res_labels, mat_1M[:,3], subplot=sp, label="Mannose", color=colors[2])
scatter!(res_labels, mat_1M[:,4], subplot=sp, label="Xylose", color=colors[3])

sp = 3
scatter!(res_labels, mat_2M[:,1], subplot=sp, label="Water", title="2.0 M Cosolvent", color=:black)
scatter!(res_labels, mat_2M[:,2], subplot=sp, label="Glucose", color=colors[1])
scatter!(res_labels, mat_2M[:,3], subplot=sp, label="Mannose", color=colors[2])
scatter!(res_labels, mat_2M[:,4], subplot=sp, label="Xylose", color=colors[3])

sp=4
pdb = readPDB("./../pdbs/solvated_glucose2.0M_cd.pdb")
protein = PDBTools.select(pdb, "protein")
residues = collect(eachresidue(protein))
labels = PDBTools.oneletter.(resname.(residues)).*format.(resnum.(residues))

### avg RMSF in water ###
res_rmsf = readdlm("./../analyses/average_rmsf_water.dat")
plot!(res_rmsf[:,1], res_rmsf[:,2],
      xlabel=L"\mathrm{Residue \ Number}",
      ylabel=L"\mathrm{Average \ RMSF \ (\AA)}",
      label=L"\mathrm{Water}",
      color=:black,
      legend=:topleft,
      xticks = (2:10:434, string.(labels[1:10:end])),
      xrotation=90,
      xlim=(1, 433),
      ylim=(0, 8.8),
      subplot=sp
     )

for cosol in [ "glucose", "mannose", "xylose" ]
    res_rmsf = readdlm("./../analyses/average_rmsf_$(cosol)_2.0.dat")
    plot!(res_rmsf[:,1], res_rmsf[:,2],
          xlabel=L"\mathrm{Residue \ Number}",
          ylabel=L"\mathrm{Average \ RMSF \ (\AA)}",
          label=latexstring("\$\\mathrm{$(cosol) \\ 2.0 \\ M}\$"),
          color=colors[findfirst(==(cosol), cosolvents)],
          legend=:topleft,
          subplot=sp
         )
    # RMSFs of the W38, W40, W367, W376
    w38 = res_rmsf[37,2]
    w40 = res_rmsf[39,2]
    w367 = res_rmsf[366,2]
    w376 = res_rmsf[375,2]
    println("RMSF for W38: $w38")
    println("RMSF for W40: $w40")
    println("RMSF for W367: $w367")
    println("RMSF for W376: $w376")
    println("")
end

sp=4
# B1 loop
vline!([51],linestyle=:dash,color=:black,alpha=0.6, subplot=sp)
vline!([56],linestyle=:dash,color=:black,alpha=0.6, subplot=sp)
# A1 loop
vline!([98],linestyle=:dash,color=:black,alpha=0.6, subplot=sp)
vline!([102],linestyle=:dash,color=:black,alpha=0.6, subplot=sp)
# B2 loop
vline!([194],linestyle=:dash,color=:black,alpha=0.6, subplot=sp)
vline!([201],linestyle=:dash,color=:black,alpha=0.6, subplot=sp)
# B3 loop
vline!([244],linestyle=:dash,color=:black,alpha=0.6, subplot=sp)
vline!([253],linestyle=:dash,color=:black,alpha=0.6, subplot=sp)
# B4 loop
vline!([339],linestyle=:dash,color=:black,alpha=0.6, subplot=sp)
vline!([342],linestyle=:dash,color=:black,alpha=0.6, subplot=sp)
# A3 loop
vline!([369],linestyle=:dash,color=:black,alpha=0.6, subplot=sp)
vline!([373],linestyle=:dash,color=:black,alpha=0.6, subplot=sp)
# A4 loop
vline!([383],linestyle=:dash,color=:black,alpha=0.6, subplot=sp)
vline!([392],linestyle=:dash,color=:black,alpha=0.6, subplot=sp)
# A2 loop
vline!([399],linestyle=:dash,color=:black,alpha=0.6, subplot=sp)
vline!([411],linestyle=:dash,color=:black,alpha=0.6, subplot=sp)

plot!()

savefig("./../figures/plot_fig4.svg")
