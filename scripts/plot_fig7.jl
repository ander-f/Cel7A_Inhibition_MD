import Pkg
Pkg.activate(".")
Pkg.instantiate()

using DelimitedFiles, LaTeXStrings, Plots, EasyFit, ColorSchemes, DataFrames, CSV, JSON3

function parse_array_column_to_float_vectors!(df::DataFrame, column_name::Symbol)
  df[!, column_name] = [
      typeof(x) <: AbstractString ? 
          Float64.(collect(JSON3.read(strip(x, ['"', ' '])))) : 
          Float64.(collect(x))
      for x in df[!, column_name]
  ]
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
    legend=:bottomleft,
    bottom_margin=3*margin,
    top_margin=margin, 
    left_margin=5*margin, 
    right_margin=margin,
    alpha=0.8,
    markerstrokewidth=2.5,
    minorticks=false)

scalefontsizes(); scalefontsizes(1.3)

#plot(layout=(2,2))
plot(layout=@layout [ [ a b ] [ d e ] ; c f ] )

cosolvents = ["glucose", "mannose", "xylose" ]
c = [ColorSchemes.tab10[1], ColorSchemes.tab10[2], ColorSchemes.tab10[3]]

# CBM
# Catalytic domain (CD)
for i in 1:length(cosolvents)

  # Load data
  cm_ADD_ff = CSV.read("./../analyses/average_MDDF_KBI_par_$(cosolvents[i])_cbm_glycan_rep1_20_last80ns_end.csv", DataFrame)

  parse_array_column_to_float_vectors!(cm_ADD_ff, :Distance)
  parse_array_column_to_float_vectors!(cm_ADD_ff, :Avg_MDDF_wat)
  parse_array_column_to_float_vectors!(cm_ADD_ff, :Err_MDDF_wat)
  parse_array_column_to_float_vectors!(cm_ADD_ff, :Avg_KBI_wat)
  parse_array_column_to_float_vectors!(cm_ADD_ff, :Err_KBI_wat)
  parse_array_column_to_float_vectors!(cm_ADD_ff, :Avg_MDDF_cos)
  parse_array_column_to_float_vectors!(cm_ADD_ff, :Err_MDDF_cos)
  parse_array_column_to_float_vectors!(cm_ADD_ff, :Avg_KBI_cos)
  parse_array_column_to_float_vectors!(cm_ADD_ff, :Err_KBI_cos)

  # Extrair o último valor (na maior distância) de cada vetor
  last_KBI_wat = [ (typeof(v) <: AbstractVector && !isempty(v)) ? v[end] : NaN for v in cm_ADD_ff.Avg_KBI_wat ] / 1000
  last_Err_KBI_wat = [ (typeof(v) <: AbstractVector && !isempty(v)) ? v[end] : NaN for v in cm_ADD_ff.Err_KBI_wat ] / 1000

  last_KBI_cos = [ (typeof(v) <: AbstractVector && !isempty(v)) ? v[end] : NaN for v in cm_ADD_ff.Avg_KBI_cos ] / 1000
  last_Err_KBI_cos = [ (typeof(v) <: AbstractVector && !isempty(v)) ? v[end] : NaN for v in cm_ADD_ff.Err_KBI_cos ] / 1000

  # plot usando apenas o último valor (final distance)
  scatter!(cm_ADD_ff.Avg_bulk_conc, last_KBI_wat, yerr=last_Err_KBI_wat,
             label="",
             xlabel=L"\mathrm{[Cosolvent] \ (M)}",
             ylabel=L"\mathrm{G_{pw}} / \mathrm{L\ mol^{-1} - CBM}",
             xlim=(0,2.2),
             ylim=(-3.4, -1.15),
             markerstrokecolor=c[i],
             marker=:diamond,
             color=:white,
             subplot=1
             )
  scatter!(cm_ADD_ff.Avg_bulk_conc, last_KBI_cos, yerr=last_Err_KBI_cos,
             label="",
             xlabel=L"\mathrm{[Cosolvent] \ (M)}",
             ylabel=L"\mathrm{G_{pc}} / \mathrm{L\ mol^{-1} - CBM}",
             xlim=(0,2.2),
             ylim=(-3.4, -1.15),
             markerstrokecolor=c[i],
             marker=:diamond,
             color=:white,
             subplot=2
             )
    scatter!(cm_ADD_ff.Avg_bulk_conc, cm_ADD_ff.Avg_par, yerr=cm_ADD_ff.Err_par,
             label="",
             xlabel=L"\mathrm{[Cosolvent] \ (M)}",
             ylabel=L"\mathrm{Γ_{pc} - CBM}",
             xlim=(0,2.2),
             ylim=(-0.4,1.2),
             markerstrokecolor=c[i],
             marker=:diamond,
             color=:white,
             subplot=5
             )

end
hline!([0.0], color=:grey, linestyle=:dash, label="", subplot=5)

# Catalytic domain (CD)
for i in 1:length(cosolvents)

  # Load data
  cm_ADD_ff = CSV.read("./../analyses/average_MDDF_KBI_par_$(cosolvents[i])_cd_glycan_rep1_5_end.csv", DataFrame)

  parse_array_column_to_float_vectors!(cm_ADD_ff, :Distance)
  parse_array_column_to_float_vectors!(cm_ADD_ff, :Avg_MDDF_wat)
  parse_array_column_to_float_vectors!(cm_ADD_ff, :Err_MDDF_wat)
  parse_array_column_to_float_vectors!(cm_ADD_ff, :Avg_KBI_wat)
  parse_array_column_to_float_vectors!(cm_ADD_ff, :Err_KBI_wat)
  parse_array_column_to_float_vectors!(cm_ADD_ff, :Avg_MDDF_cos)
  parse_array_column_to_float_vectors!(cm_ADD_ff, :Err_MDDF_cos)
  parse_array_column_to_float_vectors!(cm_ADD_ff, :Avg_KBI_cos)
  parse_array_column_to_float_vectors!(cm_ADD_ff, :Err_KBI_cos)

  # Extrair o último valor (na maior distância) de cada vetor
  last_KBI_wat = [ (typeof(v) <: AbstractVector && !isempty(v)) ? v[end] : NaN for v in cm_ADD_ff.Avg_KBI_wat ] / 1000
  last_Err_KBI_wat = [ (typeof(v) <: AbstractVector && !isempty(v)) ? v[end] : NaN for v in cm_ADD_ff.Err_KBI_wat ] / 1000

  last_KBI_cos = [ (typeof(v) <: AbstractVector && !isempty(v)) ? v[end] : NaN for v in cm_ADD_ff.Avg_KBI_cos ] / 1000
  last_Err_KBI_cos = [ (typeof(v) <: AbstractVector && !isempty(v)) ? v[end] : NaN for v in cm_ADD_ff.Err_KBI_cos ] / 1000

  # plot usando apenas o último valor (final distance)
  scatter!(cm_ADD_ff.Avg_bulk_conc, last_KBI_wat, yerr=last_Err_KBI_wat,
             label="",
             xlabel=L"\mathrm{[Cosolvent] \ (M)}",
             ylabel=L"\mathrm{G_{pw}} / \mathrm{L\ mol^{-1} - CD}",
             xlim=(0,2.2),
             ylim=(-47,-31),
             markerstrokecolor=c[i],
             marker=:diamond,
             color=:white,
             subplot=3
             )
  scatter!(cm_ADD_ff.Avg_bulk_conc, last_KBI_cos, yerr=last_Err_KBI_cos,
             label="",
             xlabel=L"\mathrm{[Cosolvent] \ (M)}",
             ylabel=L"\mathrm{G_{pc}} / \mathrm{L\ mol^{-1} - CD}",
             xlim=(0,2.2),
             ylim=(-47,-31),
             markerstrokecolor=c[i],
             marker=:diamond,
             color=:white,
             subplot=4
             )
    scatter!(cm_ADD_ff.Avg_bulk_conc, cm_ADD_ff.Avg_par, yerr=cm_ADD_ff.Err_par,
             label="$(cosolvents[i])",
             xlabel=L"\mathrm{[Cosolvent] \ (M)}",
             ylabel=L"\mathrm{Γ_{pc} - CD}",
             xlim=(0,2.2),
             ylim=(-30.0,0),
             markerstrokecolor=c[i],
             marker=:diamond,
             color=:white,
             subplot=6
             )

end

plot!(size=(1360*0.9, 720*0.9))

savefig("./../figures/plot_fig7.svg")
