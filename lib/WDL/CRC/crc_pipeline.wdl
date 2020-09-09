
# WORKFLOW DEFINITION
workflow crc_pipeline {
  String? custom_active_gene_script
  String genome
  String sample_name
  File sub_peaks_file

  #File enhancer_file
  #File chrom_path 

  call active_gene {
    input:
      custom_active_gene_script = custom_active_gene_script,
      sub_peaks_file = sub_peaks_file,
      genome = genome,
      sample_name = sample_name
  }

  # Run the validation 
  #call crc {
  #  input:
  #    enhancer_file = enhancer_file,
  #    genome = genome,
  #    chrom_path = chrom_path,
  #    sample_name = sample_name,
  #    sub_peaks_file = sub_peaks_file,
  #}

  # Outputs that will be retained when execution is complete
  output {
    File active_gene_file = active_gene.active_gene_file
    #Array[File] crc_results = crc.results
  }
}

# TASK DEFINITIONS


task active_gene {
  # Command parameters
  String? custom_active_gene_script
  File sub_peaks_file
  String genome
  String sample_name

  String? pyScript = if defined(custom_active_gene_script) then custom_active_gene_script else "/opt/ngsperl/lib/Chipseq/activeGene.py"

  Int? machine_mem_gb
  Int? disk_space_gb

  command {
    python ${pyScript} -g ${genome} -i ${sub_peaks_file} -o ${sample_name}
  }
  runtime {
    memory: select_first([machine_mem_gb, 10]) + " GB"
    cpu: "1"
    docker: "shengqh/bioinfo:novartis"
    disks: "local-disk " + select_first([disk_space_gb, 20]) + " HDD"
  }
  output {
    File active_gene_file = "${sample_name}.TSS_ACTIVE_-1000_1000.txt"
  }
}

task crc {
  # Command parameters
  File enhancer_file
  String genome
  String sample_name
  File chrom_path 
  File sub_peaks_file

  Int? machine_mem_gb
  Int? disk_space_gb

  command {
    crc -e ${enhancer_file} \\
        -g ${genome} \\
        -n ${sample_name} \\
        -c ${chrom_path} \\
        -s ${sub_peaks_file} \\
        -o .
  }
  runtime {
    memory: select_first([machine_mem_gb, 10]) + " GB"
    cpu: "1"
    docker: "shengqh/bioinfo:novartis"
    disks: "local-disk " + select_first([disk_space_gb, 20]) + " HDD"
  }
  output {
    Array[File] results = ["${sample_name}_EDGE_TABLE.txt"]
  }
}
