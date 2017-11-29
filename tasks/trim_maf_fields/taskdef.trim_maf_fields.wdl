task trim_maf_fields {

    #Inputs and constants defined here
    String pair_id
    File InputMAF
    File RemoveFieldsList
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

run('/bin/bash /opt/src/trim_maf_fields.sh \"${pair_id}\" \"${InputMAF}\" \"${RemoveFieldsList}\" ')

#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd"
        python -c "$python_cmd"

    }

    output {
        File trimmed_maf="${pair_id}.annotated.trim.maf"
    }

    runtime {
        docker : "docker.io/chipstewart/trim_maf_fields:1"
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

workflow trim_maf_fields_workflow {
    call trim_maf_fields
}
