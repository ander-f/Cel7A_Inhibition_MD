import Pkg
Pkg.activate(".")
Pkg.instantiate()

using StatsPlots, Plots, LaTeXStrings, DelimitedFiles, ColorSchemes

# 1. Configuração Global (Estilo Publication-Ready)
plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    framestyle=:box,
    grid=false,
    lw=0.5,
    dpi=600,
    alpha=0.4,
    leftmargin=0.4Plots.Measures.cm,
    xlim=(0.0, 12),
    ylim=(-0.02, 1),
    size=(550*0.9, 750*0.9), # Proporção ideal para 3 subplots verticais
    minorticks=false
)

# 2. Definição de Bins (IMPORTANTE: Garante larguras iguais)
# Ajuste o intervalo (0.0 a 6.0) e o passo (0.05) conforme seus dados
step_size = 0.1
bin_edges = 0.0:step_size:13.2 # 13.2 is the maximum RMSD value for water.
norm_type = :pdf # Normaliza a área total para 1 (Soma das alturas * 0.1 = 1)
#norm_type = :probability # Normaliza a soma das alturas para 1 (Soma das alturas = 1, mas área total pode ser diferente de 1)

base_dir = "/home/ander/Desktop/dev_projects/proj_cel7A_inhibition/ADD_ff"
cosolvents = ["glucose", "mannose", "xylose"]
concentrations = ["0.5", "1.0", "2.0"]
colors = [ColorSchemes.tab10[1], ColorSchemes.tab10[2], ColorSchemes.tab10[3]]

# Inicializa o layout
p = plot(layout=(3,1), link=:x) # link=:x garante que o eixo X seja o mesmo para todos

for (i, conc) in enumerate(concentrations)
    
    # --- Água (Referência em preto no fundo) ---
    res_water = readdlm("./../analyses/all_rmsd_loopB3_water.dat")[1:2:end, 1]
    println(maximum(res_water))
    println("")

    density!(p[i], res_water,
           color=:black,
           lw=2,
           alpha=1.0,
           #fill=(0, 0.12, :black),
           label=L"\mathrm{Water}",
           #bandwidth = 0.1, trim = false
    )

    # --- Cossolventes ---
    for (j, cosol) in enumerate(cosolvents)
        file_path = "$base_dir/analyses/all_rmsd_loopB3_$(cosol)_$(conc).dat"
        if isfile(file_path)
            data = readdlm(file_path)[1:1:end, 1]
            println(maximum(data))
            println("")

            density!(p[i], data,
                   color=colors[j],
                   lw=2,
                   alpha=1.0,
                   label=latexstring("\\mathrm{$(uppercasefirst(cosol)) \\ $(conc) \\ M}"),
            )
         else
            @warn "File not found, skipping" file=file_path
        end
    end

    # 3. Ajustes de Estética por Subplot
    ylabel!(p[i], L"\mathrm{Probability \ Density}")
    
    # Remove x-label e ticks dos subplots superiores para economizar espaço
    if i < 3
        xlabel!(p[i], "")
    else
        xlabel!(p[i], L"\mathrm{RMSD \ of \ B3 \ Loop \ (\AA)}")
    end
    
    # Legenda compacta e organizada
    plot!(p[i], legend=:topright, legend_foreground_color=:white, legend_font_pointsize=8)
end

indexes = [ "A)", 
            "B)",
            "C)"
          ]
for sp in 1:3
    annotate!(-1.8, 1, text(latexstring("\$\\mathrm{$(indexes[sp])}\$"), :left, 12), subplot=sp)
end

plot!()
# 4. Salvar
savefig("./../figures/plot_fig5.svg")
