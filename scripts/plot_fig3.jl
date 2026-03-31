import Pkg
Pkg.activate(".")
Pkg.instantiate()

using DelimitedFiles, LaTeXStrings, Plots, EasyFit
using ColorSchemes, ComplexMixtures, PDBTools

#cosolvent = ARGS[1]

# =======================
# Moving average
# =======================
ma(data) = movingaverage(data, 10).x

# =======================
# Plot defaults
# =======================
plot_font = "Computer Modern"
default(
    fontfamily = plot_font,
    linewidth = 1.5,
    framestyle = :box,
    label = nothing,
    grid = false,
    dpi = 300,
    legend = :topright,
    xlim = (0, 8.0),
    minorticks = false,
    alpha = 0.8
)

# =======================
# Colors
# =======================
divide=8
colors = [ get(ColorSchemes.brg,(8)/divide),
           get(ColorSchemes.brg,(4.5)/divide),
           get(ColorSchemes.brg,(1)/divide),
         ]


# =======================
# Paths, resíduos e açúcares
# =======================
path = "./../analyses"

# Lista de resíduos desejados
residues = [5, 31, 32]

conc_1 = "0.5"
conc_2 = "1.5"
conc_3 = "2.0"

# Lista de açúcares
sugars = ["xylose"]


# =======================
# Função para plotar 12 painéis (4 resíduos x 3 açúcares)
# =======================
function plot_residues_sugars(residues, sugars)
    plt = plot(layout = (length(sugars), length(residues)))

    concentrations = [conc_1, conc_2, conc_3]
    for (col, sugar) in enumerate(sugars)
        for (row, resid) in enumerate(residues)
            sp = (row - 1) * length(sugars) + col
            for (i, conc) in enumerate(concentrations)
                label = "$(sugar) $(conc) M"
                pdb_path = "./../pdbs/solvated_xylose0.5M_cbm.pdb"
                json_path = "$path/results-$(sugar)_$(conc)M_cbm_last80ns.json"
                color = colors[i]
                try
                    atoms = readPDB(pdb_path)
                    R = load(json_path)
                    sel = PDBTools.select(atoms, "protein and sidechain and resnum $resid")
                    group = SoluteGroup(sel)
                    patches = contributions(R, group)
                    patches_cn = contributions(R, group; type=:coordination_number)
                    if isempty(sel)
                        title_resid = "res$(resid)"
                    else
                        title_resid = try
                            one = PDBTools.oneletter(sel[1].resname)
                            string(one, sel[1].resnum)
                        catch
                            string(sel[1].resname, sel[1].resnum)
                        end
                    end
                    plot!(
                        plt,
                        R.d,
                        ma(patches),
                        color = color,
                        label = label,
                        xaxis = ("r/Å",0:2:8.0),
                        ylabel = L"\mathsf{g_{ps} \ (r)}",
                        ylim = (-0.01, 0.26),
                        yticks = 0:0.05:0.25,
                        title = title_resid,
                        titlefontsize = 10,
                        subplot = sp
                    )
                    p_right = twinx(plt[sp])
                    plot!(
                        p_right,
                        R.d,
                        patches_cn,
                        color = color,
                        linestyle = :dash,
                        label = false,
                        ylabel = "CN",
                        ylim = (-0.14, 3.64),
                        yticks = 0:0.7:3.5,
                        grid = false,
                        box = :on
                    )
                catch err
                    @warn "Erro ao processar $label, resíduo $resid: $err"
                end
            end
        end
    end
    plot!(
        plt,
        size = (900, 200),
        #top_margin = -0.2Plots.Measures.cm,
        bottom_margin = 0.5Plots.Measures.cm,
        left_margin = 1Plots.Measures.cm,
        right_margin = 0.6Plots.Measures.cm
    )
    savefig("./../figures/plot_fig3.svg")
end

# =======================
# MAIN
# =======================
plot_residues_sugars(residues, sugars)
#plot!(size=(900,200))
