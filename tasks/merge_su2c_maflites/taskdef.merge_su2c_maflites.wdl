task merge_su2c_maflites { 

    #Inputs and constants defined here

    String ID
    String TUMOR_BARCODE 
    String NORMAL_BARCODE 
    File PRODUCTION_MAF
    File M2_MAFLITE
    File STRELKA_MAFLITE
    File SVABA_MAFLITE

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

run('echo python /opt/src/merge_su2c_maflites.py  -i ${ID}  -t ${TUMOR_BARCODE} -n ${NORMAL_BARCODE} -P ${PRODUCTION_MAF}   -2 ${M2_MAFLITE} -s ${STRELKA_MAFLITE} -S ${SVABA_MAFLITE} ')

run('python /opt/src/merge_su2c_maflites.py  -i ${ID}  -t ${TUMOR_BARCODE} -n ${NORMAL_BARCODE} -P ${PRODUCTION_MAF}   -2 ${M2_MAFLITE} -s ${STRELKA_MAFLITE} -S ${SVABA_MAFLITE} ')


run('tar cvfz ${ID}.merge.maflite.tar.gz *.tsv')

#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd"
        python -c "$python_cmd"

    }

    output {
        File merge_tsv="${ID}.merge.maflite.tsv"
        File subset_tarball="${ID}.merge.maflite.tar.gz"
    }

    runtime {
        docker : "docker.io/chipstewart/merge_su2c_maflites:1"
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

workflow merge_su2c_maflites_workflow {
    call merge_su2c_maflites
}
