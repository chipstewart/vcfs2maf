task rebc_consensus_maf {
    #Inputs and constants defined here
    String pair_id
    File maf_file
    Int NALGORITHM
    String ALGORITHM_VOTING_SCHEME
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

run('/bin/bash /opt/src/rebc_consensus_maf.sh \"${pair_id}\"  \"${maf_file}\"  \"${NALGORITHM}\" \"${ALGORITHM_VOTING_SCHEME}\" ')

run('tar cvfz ${pair_id}.rebc_consensus.maf.tar.gz *.maf')

#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd"
        python -c "$python_cmd"

    }

    output {
        File consensus_all_maf="${pair_id}.consensus.all.maf"
        File consensus_pass_maf="${pair_id}.consensus.pass.maf"
    }

    runtime {
        docker : "docker.io/chipstewart/rebc_consensus_maf:1"
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

workflow rebc_consensus_maf_workflow {
    call rebc_consensus_maf
}

