task merge_maflite {
    #Inputs and constants defined here
    String pair_id
    File maf_file
    Int NALGORITHM
    Int LONGINDEL
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

run('/bin/bash /opt/src/merge_maflite.sh \"${pair_id}\" \"${tumor_id}\" \"${normal_id}\" \"${algorithm1_maflite_file}\"  \"${algorithm2_maflite_file}\" \"${algorithm3_maflite_file}\" \"${algorithm4_maflite_file}\"  \"${algorithm1}\" \"${algorithm2}\" \"${algorithm3}\" \"${algorithm4}\" ')

run('tar cvfz ${pair_id}.merged.maflite.tar.gz *.tsv')

#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd"
        python -c "$python_cmd"

    }

    output {
        File merged_maflite="${pair_id}.merged.maflite.tsv"
        File merged_maflite_tarball="${pair_id}.merged.maflite.tar.gz"
    }

    runtime {
        docker : "docker.io/chipstewart/merge_maflite:1"
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

workflow merge_maflite_workflow {
    call merge_maflite
}

