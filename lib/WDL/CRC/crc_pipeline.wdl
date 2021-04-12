
# WORKFLOW DEFINITION
workflow crc_pipeline {
  String? custom_active_gene_script
  String genome
  String sample_name
  
  File sub_peaks_file
  String sub_peaks_file_prefix

  String rose2_genome

  File bam_file
  File bam_file_index

  File enhancer_file
  File active_gene_file
  
  File chrom_file 
  String chrom_file_folder

  call active_gene {
    input:
      custom_active_gene_script = custom_active_gene_script,
      sub_peaks_file = sub_peaks_file,
      genome = genome,
      sample_name = sample_name
  }

  call rose2 {
    input :
      bam_file = bam_file,
      bam_file_index = bam_file_index,
      genome = rose2_genome,
      sub_peaks_file = sub_peaks_file,
      sub_peaks_file_prefix = sub_peaks_file_prefix
  }

  call crc {
    input:
      enhancer_file = rose2.enhancer_file,
      genome = genome,
      chrom_file = chrom_file,
      chrom_file_folder = chrom_file_folder,
      sample_name = sample_name,
      sub_peaks_file = sub_peaks_file,
      active_gene_file = active_gene.active_gene_file
  }

  # Outputs that will be retained when execution is complete
  output {
    #File active_gene_file = active_gene.active_gene_file
    #File enhancer_file = rose2.enhancer_file
    Array[File] crc_results = crc.results
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
    python3 ${pyScript} -g ${genome} -i ${sub_peaks_file} -o ${sample_name}
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

task rose2 {
  # Command parameters
  File bam_file
  File bam_file_index
  String genome
  File sub_peaks_file
  String sub_peaks_file_prefix

  Int? machine_mem_gb
  Int? disk_space_gb

  command {
    ROSE2 -g ${genome} -o . -i ${sub_peaks_file} -r ${bam_file}
  }
  runtime {
    memory: select_first([machine_mem_gb, 10]) + " GB"
    cpu: "1"
    docker: "shengqh/bioinfo:novartis"
    disks: "local-disk " + select_first([disk_space_gb, 20]) + " HDD"
  }
  output {
    File enhancer_file = "${sub_peaks_file_prefix}_SuperEnhancers_ENHANCER_TO_GENE.txt"
  }
}

task crc {
  # Command parameters
  File enhancer_file
  String genome
  String sample_name
  File chrom_file 
  String chrom_file_folder
  File sub_peaks_file
  File active_gene_file

  Int? machine_mem_gb
  Int? disk_space_gb

  command {
    tar -xzvf ${chrom_file}
    crc -e ${enhancer_file} -g ${genome} -n ${sample_name} -c ${chrom_file_folder} -s ${sub_peaks_file} -a ${active_gene_file} -o .
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
