import Pkg
Pkg.activate(".")
Pkg.instantiate()

using MolSimToolkit.PackmolInputCreator

cosolvent = ARGS[1]
conc = parse(Float64, ARGS[2])
replica = ARGS[3]
protein = ARGS[4]
dir_path = "./../$cosolvent/$conc/$replica/$protein"
if !isdir(dir_path)
    mkpath(dir_path; mode=0o777)
end

if cosolvent == "glucose"
density_table = [ # Fraction vs. Density (g/mL) - 298.15 K, ref: https://doi.org/10.1007/s10953-006-9100-7 - Glucose
                 0.0000000  0.99705 
                 0.0000449  0.99717 
                 0.0001090  0.99744 
                 0.0002151  0.99785 
                 0.0004253  0.99864 
                 0.0006490  0.99949 
                 0.0008631  1.00030 
                 0.0010670  1.00105 
                 0.0012780  1.00184 
                 0.0015000  1.00267 
                 0.0016970  1.00340 
                 0.0019190  1.00422 
                 0.0021374  1.00502 
                 0.0032124  1.00894 
                 0.0044292  1.01333 
                 0.0053277  1.01654 
                 0.0064350  1.02043 
                 0.0085708  1.02773 
                 0.0107420  1.03503 
                 0.0127970  1.04169 
                 0.0133310  1.04339 
                 0.0150710  1.04887 
                 0.0167050  1.05393 
                 0.0189840  1.06089
                 0.0219480  1.06930 
                 0.0265730  1.08227 
                 0.0299870  1.09161 
                 0.0333010  1.09978 
                 0.0333060  1.09990 
                 0.0363150  1.10748 
                 0.0398470  1.11557 
                 0.0435090  1.12413 
                 0.0503530  1.13879 
                 0.0531760  1.14468 
                 0.0603720  1.15888
                 0.0651700  1.17600 ### ref: https://doi.org/10.1021/je960168t
                 0.0702600  1.18569
                 0.0750400  1.19459
                 0.0800100  1.20340
                 0.0850300  1.21196
                 0.0902300  1.22048
                 0.0949600  1.22792
                 1.0000000  1.54000 ### ref: https://pubchem.ncbi.nlm.nih.gov/compound/D-glucose
                ]

elseif cosolvent == "mannose"
            
density_table = [ 
                   0.0000000  0.99705 ### ref: https://doi.org/10.1007/s10953-006-9100-7
                   0.0000435  0.99721
                   0.0001070  0.99746
                   0.0002320  0.99792
                   0.0004330  0.99869
                   0.0006500  0.99950
                   0.0008470  1.00025
                   0.0010770  1.00112
                   0.0012670  1.00183
                   0.0014970  1.00269
                   0.0016960  1.00343
                   0.0019070  1.00422
                   0.0021290  1.00505
                   0.0031787  1.00891
                   0.0043089  1.01301
                   0.0052984  1.01654
                   0.0063895  1.02037
                   0.0085118  1.02771
                   0.0127810  1.04188
                   0.0137690  1.04511
                   0.0155780  1.05083
                   0.0173770  1.05635
                   0.0190840  1.06153
                   0.0197080  1.06342
                   0.0217570  1.06941
                   0.0226650  1.07207
                   0.0243460  1.07680
                   0.0248680  1.07828
                   0.0259190  1.08128
                   0.0274760  1.08553
                   0.0344250  1.10358
                   1.0000000  1.54000 ### ref: Google
                  ]

elseif cosolvent == "xylose"
            
density_table = [ 
                   0.0000000  0.99705  ### ref: https://doi.org/10.1007/s10953-005-2751-y
                   0.0035901  1.00772  # 0.2000    1.00772    # molality vs density at 298.15 K                   
                   0.0071546  1.01797  # 0.4000    1.01797                            
                   0.0106936  1.02782  # 0.6000    1.02782                            
                   0.0142075  1.03730  # 0.8000    1.03730                            
                   0.0176965  1.04640  # 1.0000    1.04640                            
                   0.0211609  1.05520  # 1.2000    1.05520                            
                   1.0000000  1.52500  # mas a 20 graus celcius, https://pubchem.ncbi.nlm.nih.gov/compound/D-Xylose#section=Density
                  ]

end

# Directory of test files
dir = "./../pdbs"
# Construction of system data structure
system = SolutionBoxUSC(
    solute_pdbfile = "$dir/$protein.pdb",
    solvent_pdbfile = "$dir/tip3p.pdb",
    cossolvent_pdbfile = "$dir/$(cosolvent).pdb",
    density_table = density_table,
    concentration_units = "x", # molar fraction
#    concentration_units = "mol/L", # molarity
    solute_molar_mass = nothing, # optional
    solvent_molar_mass = nothing, # optional
    cossolvent_molar_mass = nothing, # optional
)

convert_density_table!(system, "mol/L")

# Create the box.inp file

write_packmol_input(
    system;
    concentration = conc,
    box_sides = occursin("cd", protein) ? [120.0, 120.0, 120.0] : [60.0, 60.0, 60.0],
    input = "./../$cosolvent/$conc/$replica/$protein/box_naive.inp",
    output = "solvated.pdb"
)

### Create a new box.inp file to subtract 22 water molecules and include 22 potassium ions
function rewrite_box(conc, replica, num_to_replace)
    # Nome dos arquivos
    input_file = "./../$cosolvent/$conc/$replica/$protein/box_naive.inp"
    output_file = "./../$cosolvent/$conc/$replica/$protein/box.inp"

    # Ler o arquivo original
    lines = readlines(input_file)

    for i in eachindex(lines)
        # Procura pela linha que define a estrutura do solvente
        if occursin("/tip3p.pdb", lines[i])
            # A linha com o número de moléculas geralmente é a próxima (i+1)
            line_with_number = lines[i+1]
            
            # Extrai o número atual
            parts = split(line_with_number) # Divide a linha em palavras, ex: ["", "number", "27221"]
            original_count = parse(Int, parts[end]) # Converte o último elemento para inteiro
            
	    println(original_count)

            # Calcula o novo número
            new_count = original_count - num_to_replace
            
            # Substitui a linha antiga pela nova no array de linhas
            lines[i+1] = "  number $(new_count)"
            println("Número de 'tip3p.pdb' alterado de $original_count para $new_count.")
            break # Interrompe o loop pois já encontrou e modificou o que precisava
        end
    end

    # --- 3. Encontrar a estrutura de cosolvent.pdb para usar como modelo ---
    cosolvent_structure_block = ""
    for i in eachindex(lines)
        if occursin("/$(cosolvent).pdb", lines[i])
            cosolvent_structure_block = join(lines[i:min(i+3, length(lines))], "\n")
            break
        end
    end

    if isempty(cosolvent_structure_block)
        println("AVISO: Estrutura modelo '$(cosolvent).pdb' não foi encontrada. O novo bloco não será adicionado.")
    else
        # --- 4. Criar o novo bloco de k.pdb com o número correto ---
        # Primeiro, substitui o nome do arquivo
        k_structure_temp = replace(cosolvent_structure_block, "$(cosolvent).pdb" => "k.pdb")

        # Agora, modifica o número de moléculas nesse novo bloco
        k_lines = split(k_structure_temp, '\n') # Divide o bloco em linhas
        k_lines[2] = "  number $(num_to_replace)"  # Altera a segunda linha (a que contém "number")
        k_structure_final = join(k_lines, '\n') # Junta as linhas de volta em um bloco

        # --- 5. Inserir o novo bloco de k.pdb no final ---
        insert_index = findlast(x -> occursin("end structure", x), lines)

        if insert_index !== nothing
            insert!(lines, insert_index + 2, k_structure_final)
        else
            push!(lines, k_structure_final)
        end
    end

    # Salvar o novo arquivo
    write(output_file, join(lines, "\n"))

    println("Arquivo '$output_file' criado com sucesso.")
end

rewrite_box(ARGS[2], ARGS[3], 22::Int64)

### Creating a new topology file based in the previous model ###
function count_molecules(cosolvent, conc, replica, protein)
    # Nome do arquivo
    input_file = "./../$cosolvent/$conc/$replica/$protein/box.inp"

    # Ler o arquivo
    lines = readlines(input_file)

    # Inicializar contadores
    tip3p, sugar = 0, 0

    # Percorrer as linhas para encontrar os números de moléculas
    for i in eachindex(lines)
        if occursin("structure", lines[i])
            if occursin("tip3p.pdb", lines[i])
                match_obj = match(r"number\s+(\d+)", lines[i+1])
                tip3p = match_obj !== nothing ? parse(Int, match_obj.captures[1]) : 0
            elseif occursin("$(cosolvent).pdb", lines[i])
                match_obj = match(r"number\s+(\d+)", lines[i+1])
                sugar = match_obj !== nothing ? parse(Int, match_obj.captures[1]) : 0
            end
        end
    end

    # Retornar os valores
    return tip3p, sugar
end

nwat, ncos = count_molecules(cosolvent, conc, replica, protein)


  filem = open("./../ff/ff_protein/$protein/topol_$(cosolvent).top", "r")
  file  = open("./../$cosolvent/$conc/$replica/$protein/topol.top", "w" )

    for line in eachline(filem)
      if occursin("NCOS",line)
         println(file,replace(line,"NCOS" => "$ncos",count = 1))
      elseif occursin("NWAT",line)
        println(file,replace(line,"NWAT" => "$nwat",count = 1))
      else
        println(file,line)
      end
    end

  close(filem)
  close(file)
