task gatk4_m2_maflite_task_1 {
    #Inputs and constants defined here
    String pair_id
    String tumor_id
    String normal_id
    File M2_vcf_file
    String build_id
    String? detin_m2_filters
    String detin_m2_filters1=select_first([detin_m2_filters, "contamination,weak_evidence,slippage,clustered_events,multiallelic,base_qual,map_qual,haplotype,haplotype,haplotype,low_allele_frac,normal_artifact,n_ratio,orientation,position,strand_bias"])
    String output_disk_gb
    String boot_disk_gb = "10"
    String ram_gb = "8"
    String cpu_cores = "2"
    String? preemptible
    String preemptible1=select_first([preemptible, "2"])
  

    command {

    set -euo pipefail

    pwd
    which julia
    julia --version

    /bin/bash /opt/src/gatk4_m2_maflite.sh ${tumor_id} ${normal_id} ${pair_id} ${M2_vcf_file} ${build_id} ${detin_m2_filters1}

    ls -lath

    tar cvfz m2_maflite.tar.gz tmp1.tsv ${pair_id}.raw.tsv ${pair_id}.m2.all.maflite.tsv ${pair_id}.m2.pass.maflite.tsv ${pair_id}.m2.deTiN.maflite.tsv

    }

    output {
        File m2_pass_maflite="${pair_id}.m2.pass.maflite.tsv"
        File m2_all_maflite="${pair_id}.m2.all.maflite.tsv"
        File m2_maflite_tarball="m2_maflite.tar.gz"
        File m2_detin_filter_maflite="${pair_id}.m2.deTiN.maflite.tsv"
    }

    runtime {
        docker : "docker.io/chipstewart/gatk4_m2_maflite_task_1:1"
        memory: "${ram_gb}GB"
        cpu: "${cpu_cores}"
        disks: "local-disk ${output_disk_gb} HDD"
        bootDiskSizeGb: "${boot_disk_gb}"
        preemptible: "${preemptible1}"
    }


    meta {
        author : "Chip Stewart"
        email : "stewart@broadinstitute.org"
    }

}

workflow gatk4_m2_maflite {

    call gatk4_m2_maflite_task_1 {
    #    input: #**Define call inputs for gatk4_m2_maflite_task_1 here**
    }

    output {
        #**Define workflow outputs here. If defined, these will be the only
        #  outputs available in the Method Configuration**
    }
}
