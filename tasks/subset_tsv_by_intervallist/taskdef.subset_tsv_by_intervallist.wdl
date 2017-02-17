task subset_tsv_by_intervallist {

    #Inputs and constants defined here
    String salutation_input_string
    File name_input_file

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

run('python /opt/src/hello.py \"${salutation_input_string}\"  \"${name_input_file}\"')

run('tar cvfz greeting.tar.gz greeting.txt')

#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd"
        python -c "$python_cmd"

    }

    output {
        File greeting_txt="greeting.txt"
        File greeting_tarball="greeting.tar.gz"
    }

    runtime {
        docker : "docker.io/stewart/subset_tsv_by_intervallist:1"
        memory: "${ram_gb}GB"
        cpu: "${cpu_cores}"
        disks: "local-disk ${output_disk_gb} HDD"
        bootDiskSizeGb: "${boot_disk_gb}"
        preemptible: 0
    }


    meta {
        author : "Your Name"
        email : "Your.Email@Address.com"
    }

}

workflow subset_tsv_by_intervallist_workflow {
    call subset_tsv_by_intervallist
}
