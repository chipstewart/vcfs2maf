task svaba_maf {

    #Inputs and constants defined here
    String pair_id
    String tumor_id
    String normal_id
    File Strelka_SNV_vcf_file
    File Strelka_INDEL_vcf_file
    String output_disk_gb
    String boot_disk_gb = "10"
    String ram_gb = "8"
    String cpu_cores = "2"
    command {
python_cmd="
import subprocess
def run(cmd):
    subprocess.check_call(cmd,shell=True)

run('ln -sT `pwd` /opt/execution')
run('ln -sT `pwd`/../inputs /opt/inputs')
run('/opt/src/algutil/monitor_start.py')

# start task-specific calls
##########################

run('julia --version')

run('/bin/bash /opt/src/svaba_maf.sh \"${pair_id}\" \"${tumor_id}\" \"${normal_id}\"  \"${svaba_INDEL_vcf_file}\"')

run('tar cvfz ${pair_id}.svaba_maf.tar.gz ${pair_id}.svaba.INDEL.tsv ${pair_id}.svaba_maflite.tsv')

#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd"
        python -c "$python_cmd"

    }

    output {
        File svaba_maf="${pair_id}.svaba_maflite.tsv"
        File svaba_tarball="${pair_id}.svaba_maf.tar.gz"
    }

    runtime {
        docker : "docker.io/chipstewart/svaba_maf:1"
        memory: "${ram_gb}GB"
        cpu: "${cpu_cores}"
        disks: "local-disk ${output_disk_gb} HDD"
        bootDiskSizeGb: "${boot_disk_gb}"
        preemptible: 0
    }


    meta {
        author : "Chip Stewart"
        email : "stewart@broadinstitute.org"
    }

}

workflow svaba_maf_workflow {
    call svaba_maf
}
