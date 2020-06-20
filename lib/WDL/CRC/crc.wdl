
# WORKFLOW DEFINITION
workflow crc_pipeline {
  File enhancer_file
  String genome
  File chrom_path 
  String project_name
  File sub_peaks_file

  # Run the validation 
  call crc {
    input:
      enhancer_file = enhancer_file,
      genome = genome,
      chrom_path = chrom_path,
      project_name = project_name,
      sub_peaks_file = sub_peaks_file,
  }

  # Outputs that will be retained when execution is complete
  output {
    Array[File] crc_results = crc.results
  }
}

# TASK DEFINITIONS

task crc {
  # Command parameters
  File enhancer_file
  String genome
  String project_name
  File chrom_path 
  File sub_peaks_file

  Int? machine_mem_gb
  Int? disk_space_gb

  command {
    crc -e ${enhancer_file} \
        -g ${genome} \
        -n ${project_name} \
        -c ${chrom_path} \
        -s ${sub_peaks_file} \
        -o .
  }
  runtime {
    memory: select_first([machine_mem_gb, 10]) + " GB"
    cpu: "1"
    docker: "shengqh/bioinfo:novartis"
    disks: "local-disk " + select_first([disk_space_gb, 20]) + " HDD"
  }
  output {
    Array[File] results = ["${project_name}_EDGE_TABLE.txt"]
  }
}
