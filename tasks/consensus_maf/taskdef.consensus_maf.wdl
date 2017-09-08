task consensus_maf {
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

run('/bin/bash /opt/src/consensus_maf.sh \"${pair_id}\"  \"${maf_file}\"  \"${NALGORITHM}\" \"${LONGINDEL}\" ')

run('tar cvfz ${pair_id}.consensus.maf.tar.gz *.maf')

#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd"
        python -c "$python_cmd"

    }

    output {
        File consensus_maf="${pair_id}.consensus.maf"
        File consensus_maf_tarball="${pair_id}.consensus.maf.tar.gz"
    }

    runtime {
        docker : "docker.io/chipstewart/consensus_maf:1"
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

workflow consensus_maf_workflow {
    call consensus_maf
}

