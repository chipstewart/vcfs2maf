task merge_maflite {

    #Inputs and constants defined here
    String pair_id
    String tumor_id
    String normal_id
    File M1_maflite_file
    File M2_maflite_file
    File Strelka_maflite_file
    File SvABA_maflite_file
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

run('/bin/bash /opt/src/merge_maflite.sh \"${tumor_id}\" \"${normal_id}\" \"${M1_maflite_file}\"  \"${M2_maflite_file}\" \"${Strelka_maflite_file}\" \"${SvABA_maflite_file}\"' \"M1"\ \"M2"\ \"Strelka"\ \"SvABA"\ ) 

run('tar cvfz ${pair_id}_merged.maflite.tar.gz ${pair_id}.merged.maflite.tsv')

#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd"
        python -c "$python_cmd"

    }

    output {
        File merged_maflite="${pair_id}_merged.maflite.tsv"
        File merged_maflite_tarball="${pair_id}_merged.maflite.tar.gz"
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
