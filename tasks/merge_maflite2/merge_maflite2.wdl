task merge_maflite2_task_1 {
    #Inputs and constants defined here
    String pair_id
    String tumor_id
    String normal_id
    File algorithm1_maflite_file
    File algorithm2_maflite_file
    File algorithm3_maflite_file
    File algorithm4_maflite_file
    File algorithm5_maflite_file
    File algorithm6_maflite_file
    #String algorithm6a_maflite_file =  select_first([algorithm6_maflite_file, "-"])
    
    String algorithm1
    String algorithm2
    String algorithm3
    String algorithm4
    String algorithm5
    String algorithm6
    #String algorithm6a =  select_first([algorithm6, "-"])

    String output_disk_gb = "50"
    String boot_disk_gb = "10"
    String ram_gb = "8"
    String cpu_cores = "1"
    command {

        set -euo pipefail

        julia --version

        /bin/bash /src/merge_maflite.sh ${pair_id} ${tumor_id} ${normal_id} ${algorithm1_maflite_file}  ${algorithm2_maflite_file} ${algorithm3_maflite_file} ${algorithm4_maflite_file} ${algorithm5_maflite_file} ${algorithm6_maflite_file} ${algorithm1} ${algorithm2} ${algorithm3} ${algorithm4} ${algorithm5} ${algorithm6}  

        tar cvfz ${pair_id}.merged.maflite.tar.gz *.tsv 

    }

    output {       
        File merged_maflite="${pair_id}.merged.maflite.tsv"
        File merged_maflite_tarball="${pair_id}.merged.maflite.tar.gz"
    }

    runtime {
        docker : "docker.io/chipstewart/merge_maflite2_task_1:2"
        memory: "${ram_gb}GB"
        cpu: "${cpu_cores}"
        disks: "local-disk ${output_disk_gb} HDD"
        bootDiskSizeGb: "${boot_disk_gb}"
        preemptible: 3
    }

    meta {
        author : "Chip Stewart"
        email : "stewart@broadinstitute.org"
    }

}

workflow merge_maflite2 {
    call merge_maflite2_task_1
}
