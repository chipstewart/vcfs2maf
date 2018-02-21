task postFilter_contEst {
    #Inputs and constants defined here
    String pair_id
    File maf_file
    Float ContEst
    Float pThreshold
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

run('/bin/bash /opt/src/postFilter_contEst.sh \"${pair_id}\"  \"${maf_file}\"  \"${ContEst}\" \"${pThreshold}\" ')

run('tar cvfz ${pair_id}.postFilter_contEst.maf.tar.gz *.maf')

#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd"
        python -c "$python_cmd"

    }

    output {
        File annotated_all_maf="${pair_id}.postFilter_contEst.all.maf"
        File annotated_pass_maf="${pair_id}.postFilter_contEst.pass.maf"
        File consensus_maf_tarball="${pair_id}.postFilter_contEst.maf.tar.gz"
    }

    runtime {
        docker : "docker.io/chipstewart/postfilter_contest:1"
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

workflow postFilter_contEst_workflow {
    call postFilter_contEst
}


