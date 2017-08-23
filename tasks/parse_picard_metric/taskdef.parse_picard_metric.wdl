task parse_picard_metric {
    #Inputs and constants defined here
    File picard_metrics_file
    String metric_field 
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

run('/bin/bash /opt/src/parse_picard_metric.sh  \"${picard_metric_file}\"  \"${metric_field}\"  ')

run('/bin/bash cat   ${metric_field}.txt  ')

#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd"
        python -c "$python_cmd"

    }

    output {
        String value = read_string(stdout())
    }

    runtime {
        docker : "docker.io/chipstewart/parse_picard_metric:1"
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

workflow parse_picard_metric_workflow {
    call parse_picard_metric
}

