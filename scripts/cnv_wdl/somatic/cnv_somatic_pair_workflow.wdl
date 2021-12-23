Workflow for running the GATK CNV pipeline on a matched pair. Supports both WGS and WES.
#
# Notes:
#
# - The intervals argument is required for both WGS and WES workflows and accepts formats compatible with the
#   GATK -L argument (see https://gatkforums.broadinstitute.org/gatk/discussion/11009/intervals-and-interval-lists).
#   These intervals will be padded on both sides by the amount specified by padding (default 250)
#   and split into bins of length specified by bin_length (default 1000; specify 0 to skip binning,
#   e.g., for WES).  For WGS, the intervals should simply cover the autosomal chromosomes (sex chromosomes may be
#   included, but care should be taken to 1) avoid creating panels of mixed sex, and 2) denoise case samples only
#   with panels containing only individuals of the same sex as the case samples).
#
# - Intervals can be blacklisted from coverage collection and all downstream steps by using the blacklist_intervals
#   argument, which accepts formats compatible with the GATK -XL argument
#   (see https://gatkforums.broadinstitute.org/gatk/discussion/11009/intervals-and-interval-lists).
#   This may be useful for excluding centromeric regions, etc. from analysis.  Alternatively, these regions may
#   be manually filtered from the final callset.
#
#  A reasonable blacklist for excluded intervals (-XL) can be found at:
#   hg19: gs://gatk-best-practices/somatic-b37/CNV_and_centromere_blacklist.hg19.list
#   hg38: gs://gatk-best-practices/somatic-hg38/CNV_and_centromere_blacklist.hg38liftover.list (untested)
#
# - The sites file (common_sites) should be a Picard or GATK-style interval list.  This is a list of sites
#   of known variation at which allelic counts will be collected for use in modeling minor-allele fractions.
#
# - If you opt to run FuncotateSegments (i.e. set `is_run_funcotator` to `true`), then please also ensure that you have
#       the correct value for `funcotator_ref_version`.  Treat `funcotator_ref_version` as required if
#       `is_run_funcotator` is `true`.  Valid values for `funcotator_ref_version` are `hg38` and `hg19`.
#       The latter includes GRCh37.
#
#
# - Example invocation:
#
#       java -jar cromwell.jar run cnv_somatic_pair_workflow.wdl -i my_parameters.json
#
#############

version 1.0

import "../cnv_common_tasks.wdl" as CNVTasks
import "cnv_somatic_oncotator_workflow.wdl" as CNVOncotator
import "cnv_somatic_funcotate_seg_workflow.wdl" as CNVFuncotateSegments

workflow CNVSomaticPairWorkflow {

    input {
      ##################################
      #### required basic arguments ####
      ##################################
      File common_sites
      File intervals
      File? blacklist_intervals
      File tumor_bam
      File tumor_bam_idx
      File? normal_bam
      File? normal_bam_idx
      File read_count_pon
      File ref_fasta_dict
      File ref_fasta_fai
      File ref_fasta
      String gatk_docker
      String case_name
      String control_name
      
      ##################################
      #### optional basic arguments ####
      ##################################
       # For running oncotator
      Boolean? is_run_oncotator
       # For running funcotator
      Boolean? is_run_funcotator

      File? gatk4_jar_override
      Int? preemptible_attempts
      # Use as a last resort to increase the disk given to every task in case of ill behaving data
      Int? emergency_extra_disk

      # Required if BAM/CRAM is in a requester pays bucket
      String? gcs_project_for_requester_pays

      ####################################################
      #### optional arguments for PreprocessIntervals ####
      ####################################################
      Int? padding
      Int? bin_length
      Int? mem_gb_for_preprocess_intervals

      ##############################################
      #### optional arguments for CollectCounts ####
      ##############################################
      String? collect_counts_format
      Int? mem_gb_for_collect_counts

      #####################################################
      #### optional arguments for CollectAllelicCounts ####
      #####################################################
      String? minimum_base_quality
      Int? mem_gb_for_collect_allelic_counts

      ##################################################
      #### optional arguments for DenoiseReadCounts ####
      ##################################################
      Int? number_of_eigensamples
      Int? mem_gb_for_denoise_read_counts

      ##############################################
      #### optional arguments for ModelSegments ####
      ##############################################
      Int? max_num_segments_per_chromosome
      Int? min_total_allele_count
      Int? min_total_allele_count_normal
      Float? genotyping_homozygous_log_ratio_threshold
      Float? genotyping_base_error_rate
      Float? kernel_variance_copy_ratio
      Float? kernel_variance_allele_fraction
      Float? kernel_scaling_allele_fraction
      Int? kernel_approximation_dimension
      Array[Int]+? window_sizes = [8, 16, 32, 64, 128, 256]
      Float? num_changepoints_penalty_factor
      Float? minor_allele_fraction_prior_alpha
      Int? num_samples_copy_ratio
      Int? num_burn_in_copy_ratio
      Int? num_samples_allele_fraction
      Int? num_burn_in_allele_fraction
      Float? smoothing_threshold_copy_ratio
      Float? smoothing_threshold_allele_fraction
      Int? max_num_smoothing_iterations
      Int? num_smoothing_iterations_per_fit
      Int? mem_gb_for_model_segments

      ######################################################
      #### optional arguments for CallCopyRatioSegments ####
      ######################################################
      Float? neutral_segment_copy_ratio_lower_bound
      Float? neutral_segment_copy_ratio_upper_bound
      Float? outlier_neutral_segment_copy_ratio_z_score_threshold
      Float? calling_copy_ratio_z_score_threshold
      Int? mem_gb_for_call_copy_ratio_segments

      #########################################
      #### optional arguments for plotting ####
      #########################################
      Int? minimum_contig_length
      # If maximum_copy_ratio = Infinity, the maximum copy ratio will be automatically determined
      String? maximum_copy_ratio
      Float? point_size_copy_ratio
      Float? point_size_allele_fraction
      Int? mem_gb_for_plotting

      ##########################################
      #### optional arguments for Oncotator ####
      ##########################################
      String? additional_args_for_oncotator
      String? oncotator_docker
      Int? mem_gb_for_oncotator
      Int? boot_disk_space_gb_for_oncotator

      ##################################################
      #### optional arguments for FuncotateSegments ####
      ##################################################
      String? additional_args_for_funcotator
      String? funcotator_ref_version
      Int? mem_gb_for_funcotator
      File? funcotator_transcript_selection_list
      File? funcotator_data_sources_tar_gz
      String? funcotator_transcript_selection_mode
      Array[String]? funcotator_annotation_defaults
      Array[String]? funcotator_annotation_overrides
      Array[String]? funcotator_excluded_fields
      Boolean? funcotator_is_removing_untared_datasources
      Int? funcotator_disk_space_gb
      Boolean? funcotator_use_ssd
      Int? funcotator_cpu
    }

    Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_fasta_dict, "GB") + size(ref_fasta_fai, "GB"))
    Int read_count_pon_size = ceil(size(read_count_pon, "GB"))
    Int tumor_bam_size = ceil(size(tumor_bam, "GB") + size(tumor_bam_idx, "GB"))
    Int normal_bam_size = if defined(normal_bam) then ceil(size(normal_bam, "GB") + size(normal_bam_idx, "GB")) else 0

    Int gatk4_override_size = if defined(gatk4_jar_override) then ceil(size(gatk4_jar_override, "GB")) else 0
    # This is added to every task as padding, should increase if systematically you need more disk for every call
    Int disk_pad = 20 + ceil(size(intervals, "GB")) + ceil(size(common_sites, "GB")) + gatk4_override_size + select_first([emergency_extra_disk, 0])

    File final_normal_bam = select_first([normal_bam, "null"])
    File final_normal_bam_idx = select_first([normal_bam_idx, "null"])

    Int preprocess_intervals_disk = ref_size + disk_pad
    call CNVTasks.PreprocessIntervals {
        input:
            intervals = intervals,
            blacklist_intervals = blacklist_intervals,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            padding = padding,
            bin_length = bin_length,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_preprocess_intervals,
            disk_space_gb = preprocess_intervals_disk,
            preemptible_attempts = preemptible_attempts
    }

    Int collect_counts_tumor_disk = tumor_bam_size + ceil(size(PreprocessIntervals.preprocessed_intervals, "GB")) + disk_pad
    call CNVTasks.CollectCounts as CollectCountsTumor {
        input:
            intervals = PreprocessIntervals.preprocessed_intervals,
            bam = tumor_bam,
            bam_idx = tumor_bam_idx,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            format = collect_counts_format,
            enable_indexing = false,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = "us.gcr.io/broad-gatk/gatk:4.0.0.0",
            mem_gb = mem_gb_for_collect_counts,
            disk_space_gb = collect_counts_tumor_disk,
            preemptible_attempts = preemptible_attempts,
            gcs_project_for_requester_pays = gcs_project_for_requester_pays
    }

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
            gcs_project_for_requester_pays = gcs_project_for_requester_pays,
            sample_id  = case_name
    }


    Int collect_counts_normal_disk = normal_bam_size + ceil(size(PreprocessIntervals.preprocessed_intervals, "GB")) + disk_pad
    if (defined(normal_bam)) {
        call CNVTasks.CollectCounts as CollectCountsNormal {
            input:
                intervals = PreprocessIntervals.preprocessed_intervals,
                bam = final_normal_bam,
                bam_idx = final_normal_bam_idx,
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                ref_fasta_dict = ref_fasta_dict,
                format = collect_counts_format,
                enable_indexing = false,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = "us.gcr.io/broad-gatk/gatk:4.0.0.0",
                mem_gb = mem_gb_for_collect_counts,
                disk_space_gb = collect_counts_normal_disk,
                preemptible_attempts = preemptible_attempts,
                gcs_project_for_requester_pays = gcs_project_for_requester_pays
        }

        Int collect_allelic_counts_normal_disk = normal_bam_size + ref_size + disk_pad
        call CNVTasks.CollectAllelicCounts as CollectAllelicCountsNormal {
            input:
                common_sites = common_sites,
                bam = final_normal_bam,
                bam_idx = final_normal_bam_idx,
                ref_fasta = ref_fasta,
                ref_fasta_dict = ref_fasta_dict,
                ref_fasta_fai = ref_fasta_fai,
                minimum_base_quality =  minimum_base_quality,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                mem_gb = mem_gb_for_collect_allelic_counts,
                disk_space_gb = collect_allelic_counts_normal_disk,
                preemptible_attempts = preemptible_attempts,
                gcs_project_for_requester_pays = gcs_project_for_requester_pays,
                sample_id = control_name
        }
}

    output {
        File preprocessed_intervals = PreprocessIntervals.preprocessed_intervals

        File read_counts_entity_id_tumor = CollectCountsTumor.entity_id
        File read_counts_tumor = CollectCountsTumor.counts
        File allelic_counts_entity_id_tumor = CollectAllelicCountsTumor.entity_id
        File allelic_counts_tumor = CollectAllelicCountsTumor.allelic_counts

        File? read_counts_entity_id_normal = CollectCountsNormal.entity_id
        File? read_counts_normal = CollectCountsNormal.counts
        File? allelic_counts_entity_id_normal = CollectAllelicCountsNormal.entity_id
        File? allelic_counts_normal = CollectAllelicCountsNormal.allelic_counts
    }
}
