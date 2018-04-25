task copy_maf_fields {

    #Inputs and constants defined here
    String pair_id
    File InputMAF
    String Copy_Fields
    String extension
    String output_disk_gb
    String boot_disk_gb = "10"
    String ram_gb = "2"
    String cpu_cores = "1"
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

run('/bin/bash /opt/src/copy_maf_fields.sh \"${pair_id}\" \"${InputMAF}\" \"${Copy_Fields}\" \"${extension}\" ')

#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd"
        python -c "$python_cmd"

    }

    output {
        File copy_fields_maf="${pair_id}${extension}"
    }

    runtime {
        docker : "docker.io/chipstewart/copy_maf_fields:1"
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

workflow copy_maf_fields_workflow {
    call copy_maf_fields
}
