task svaba_maflite {

    #Inputs and constants defined here
    String pair_id
    String tumor_id
    String normal_id
    String tumor_bam
    String normal_bam
    File svaba_INDEL_vcf_file
    String MAX_NORMAL_ALT_COUNT
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

run('/bin/bash /opt/src/svaba_maflite.sh \"${pair_id}\" \"${tumor_id}\" \"${normal_id}\"  \"${tumor_bam}\" \"${normal_bam}\" \"${svaba_INDEL_vcf_file}\" ${MAX_NORMAL_ALT_COUNT} ')

run('tar cvfz ${pair_id}.SvABA_maf.tar.gz ${pair_id}.SvABA.INDEL.tsv ${pair_id}.SvABA_maflite.tsv')

#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd"
        python -c "$python_cmd"

    }

    output {
        File svaba_maflite="${pair_id}.SvABA_maflite.tsv"
        File svaba_maflite_tarball="${pair_id}.SvABA_maf.tar.gz"
    }

    runtime {
        docker : "docker.io/chipstewart/svaba_maflite:1"
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

workflow svaba_maflite_workflow {
    call svaba_maflite
}
