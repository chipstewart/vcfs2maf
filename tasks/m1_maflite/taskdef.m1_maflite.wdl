task m1_maflite {

    #Inputs and constants defined here
    String pair_id
    String tumor_id
    String normal_id
    File M1_maflite_inputfile
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

run('/bin/bash /opt/src/m1_maflite.sh \"${tumor_id}\" \"${normal_id}\" \"${pair_id}\" \"${M1_maflite_inputfile}\"')

#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd"
        python -c "$python_cmd"

    }

    output {
        File m1_maflite_fix="${pair_id}.m1_maflite.tsv"
    }

    runtime {
        docker : "docker.io/chipstewart/m1_maflite:1"
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

workflow m1_maflite_workflow {
    call m1_maflite
}
