task gatk4_m2_maflite {

    #Inputs and constants defined here
    String pair_id
    String tumor_id
    String normal_id
    File M2_vcf_file
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

run('/bin/bash /opt/src/gatk4_m2_maflite.sh \"${tumor_id}\" \"${normal_id}\" \"${pair_id}\" \"${M2_vcf_file}\"')

run('tar cvfz m2_maflite.tar.gz tmp1.tsv ${pair_id}.raw.tsv ${pair_id}.m2.all.maflite.tsv ${pair_id}.m2.pass.maflite.tsv')

#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd"
        python -c "$python_cmd"

    }

    output {
        File m2_pass_maflite="${pair_id}.m2.pass.maflite.tsv"
        File m2_all_maflite="${pair_id}.m2.all.maflite.tsv"
        File m2_maflite_tarball="m2_maflite.tar.gz"
    }

    runtime {
        docker : "docker.io/chipstewart/gatk2_m2_maflite:1"
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

workflow gatk4_m2_maflite_workflow {
    call gatk4_m2_maflite
}
