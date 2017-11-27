task m1m2_maf {

    #Inputs and constants defined here
    String tumor_id
    String normal_id
    File M1_vcf_file
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

run('/bin/bash /opt/src/M1M2_maf.sh \"${tumor_id}\" \"${normal_id}\" \"${M1_vcf_file}\"  \"${M2_vcf_file}\"')

run('tar cvfz m1m2_maf.tar.gz tmp1.tsv tmp2.tsv m1m2_maflite.tsv')

#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd"
        python -c "$python_cmd"

    }

    output {
        File m1m2_maf="m1m2_maflite.tsv"
        File greeting_tarball="m1m2_maf.tar.gz"
    }

    runtime {
        docker : "docker.io/chipstewart/m1m2_maf:1"
        memory: "${ram_gb}GB"
        cpu: "${cpu_cores}"
        disks: "local-disk ${output_disk_gb} HDD"
        bootDiskSizeGb: "${boot_disk_gb}"
        preemptible: 4
    }


    meta {
        author : "Chip Stewart"
        email : "stewart@broadinstitute.org"
    }

}

workflow m1m2_maf_workflow {
    call m1m2_maf
}
