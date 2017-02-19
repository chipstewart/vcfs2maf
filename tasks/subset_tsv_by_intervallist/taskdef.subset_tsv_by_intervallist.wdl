task subset_tsv_by_intervallist { 

    #Inputs and constants defined here

    String ID
    File interval_list
    File input_tsv_file
    String CHROMOSOME_FIELD
    String POSITION_FIELD

    String output_disk_gb="20"
    String boot_disk_gb = "10"
    String ram_gb = "3"
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

#run('echo -s ${ID}\"  \"-i ${input_tsv_file}\"  \"-l ${interval_list}\" \"-c ${CHROMOSOME_FIELD}\"    \"-p ${POSITION_FIELD}\"   \"-f subset\"')
#run('python /opt/src/subset_tsv_by_intervallist.py  \"-s ${ID}\"  \"-i ${input_tsv_file}\" \"-l ${interval_list}\" \"-c ${CHROMOSOME_FIELD}\"    \"-p ${POSITION_FIELD}\"   \"-f subset\"')

run('echo python /opt/src/subset_tsv_by_intervallist.py  -s ${ID}  -i ${input_tsv_file} -l ${interval_list} -c ${CHROMOSOME_FIELD}   -p ${POSITION_FIELD} -f .subset')

run('python /opt/src/subset_tsv_by_intervallist.py -s ${ID} -i ${input_tsv_file} -l ${interval_list} -c ${CHROMOSOME_FIELD}   -p ${POSITION_FIELD}  -f .subset ')


run('tar cvfz ${ID}.subset.tar.gz stderr stdout ${ID}.subset.tsv')

#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd"
        python -c "$python_cmd"

    }

    output {
        File subset_tsv="${ID}.subset.tsv"
        File subset_tarball="${ID}.subset.tar.gz"
    }

    runtime {
        docker : "docker.io/chipstewart/subset_tsv_by_intervallist:1"
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

workflow subset_tsv_by_intervallist_workflow {
    call subset_tsv_by_intervallist
}
