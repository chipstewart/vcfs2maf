task m1m2SM_maf {

    #Inputs and constants defined here
    String tumor_id
    String normal_id
    File M1M2_maflite_file
    File Snowman_vcf_file
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

run('/bin/bash /opt/src/M1M2SM_maf.sh \"${tumor_id}\" \"${normal_id}\" \"${M1M2_maflite_file}\"  \"${Snowman_vcf_file}\"')

run('tar cvfz m1m2SM_maflite.tar.gz *.tsv ')

#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd"
        python -c "$python_cmd"

    }

    output {
        File m1m2SM_maf="m1m2SM_maflite.tsv"
        File greeting_tarball="m1m2SM_maflite.tar.gz"
    }

    runtime {
        docker : "docker.io/chipstewart/m1m2_maf:1"
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

workflow m1m2SM_maf_workflow {
    call m1m2SM_maf
}
