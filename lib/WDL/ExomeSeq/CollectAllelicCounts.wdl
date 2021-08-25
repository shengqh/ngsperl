version 1.0

import "/scratch/cqs/zhaos/tools/gatk/scripts/cnv_wdl/cnv_common_tasks.wdl" as CNVTasks

workflow CollectAllelicCountsWorkflow {

    input {
      ##################################
      #### required basic arguments ####
      ##################################
      File common_sites
      File? blacklist_intervals
      File tumor_bam
      File tumor_bam_idx
      File? normal_bam
      File? normal_bam_idx
      File ref_fasta_dict
      File ref_fasta_fai
      File ref_fasta
      String gatk_docker

      ##################################
      #### optional basic arguments ####
      ##################################

      File? gatk4_jar_override
      Int? preemptible_attempts
      # Use as a last resort to increase the disk given to every task in case of ill behaving data
      Int? emergency_extra_disk

      # Required if BAM/CRAM is in a requester pays bucket
      String? gcs_project_for_requester_pays

      #####################################################
      #### optional arguments for CollectAllelicCounts ####
      #####################################################
      String? minimum_base_quality
      Int? mem_gb_for_collect_allelic_counts

    }

    Int tumor_bam_size = ceil(size(tumor_bam, "GB") + size(tumor_bam_idx, "GB"))

    Int gatk4_override_size = if defined(gatk4_jar_override) then ceil(size(gatk4_jar_override, "GB")) else 0
    # This is added to every task as padding, should increase if systematically you need more disk for every call
    Int disk_pad = 20 + ceil(size(common_sites, "GB")) + gatk4_override_size + select_first([emergency_extra_disk, 0])

    Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_fasta_dict, "GB") + size(ref_fasta_fai, "GB"))

    Int collect_allelic_counts_tumor_disk = tumor_bam_size + ref_size + disk_pad
    call CNVTasks.CollectAllelicCounts as CollectAllelicCountsTumor {
        input:
            common_sites = common_sites,
            bam = tumor_bam,
            bam_idx = tumor_bam_idx,
            ref_fasta = ref_fasta,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_fai = ref_fasta_fai,
            minimum_base_quality =  minimum_base_quality,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_collect_allelic_counts,
            disk_space_gb = collect_allelic_counts_tumor_disk,
            preemptible_attempts = preemptible_attempts,
            gcs_project_for_requester_pays = gcs_project_for_requester_pays
    }

    output {
        String allelic_counts_entity_id_tumor = CollectAllelicCountsTumor.entity_id
        File allelic_counts_tumor = CollectAllelicCountsTumor.allelic_counts
    }
}
