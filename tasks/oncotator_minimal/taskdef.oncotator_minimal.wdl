task subset_tsv_by_intervallist { 

    #Inputs and constants defined here

    String ID
    File interval_list
    File input_tsv_file
    String CHROMOSOME_FIELD
    String POSITION_FIELD

    String output_disk_gb
    String boot_disk_gb 
    String ram_gb = "3.5"
    String cpu_cores = "1"
    Int preemptible_tries
    
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


run('tar cvfz ${ID}.subset.tar.gz ${ID}.subset.tsv')

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
        preemptible:  "${preemptible_tries}"
    }


    meta {
        author : "Chip Stewart"
        email : "stewart@broadinstitute.org"
    }

}


task oncotator_minimal {

    #Inputs and constants defined here
    
    File oncoDBTarBall
    File IN
    File onco_config_file 
    File column_collapse_config_file 
    String OTHER_FLAGS = " --infer-onps  --no-multicore   "

    String id
    Int preemptible_tries
    String output_disk_gb 
    String boot_disk_gb 
    String ram_gb = "7"
    String cpu_cores = "1"

    command {
python_cmd="
import subprocess, os
def run(cmd):
    subprocess.check_call(cmd,shell=True)

run('ln -sT `pwd` /opt/execution')
run('ln -sT `pwd`/../inputs /opt/inputs')
#run('/opt/src/algutil/monitor_start.py')

# start task-specific calls
##########################
run('touch start.txt')
run('tar xvf \"${oncoDBTarBall}\" ')
pwd=os.getcwd()
from glob import glob
onco_dir = glob('onco*')
onco_db = onco_dir[0] + ' '
other_flags=\"${OTHER_FLAGS}\"
config_file='${onco_config_file}'
run('echo $PWD')
run('ln -s ${column_collapse_config_file} .')
in1='${IN}'
id1='${id}'
run ('ls -latrh')
cmd='/root/oncotator_venv/bin/oncotator -i MAFLITE --db-dir ' + pwd + '/' + onco_db + ' --default_config ' + config_file +  '  --tx-mode=EFFECT -o TCGAMAF    --collapse-number-annotations ' + other_flags + ' ' +in1 + ' ' +  id1 + '.annotated.maf hg19'
print(cmd)
run(cmd)
run ('ls -latrh')
run('cut -f1-43,49,66-67,104,110-133,143-400 \"${id}.annotated.maf\" > \"${id}.maf\"')
run('tar cvfz \"${id}.annotated.maf.gz\" \"${id}.annotated.maf\"')

#########################
# end task-specific calls
#run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd"
        python -c "$python_cmd"

    }

    output {
        File oncotator_full_maf_gz = "${id}.annotated.maf.gz"
        File oncotator_out = "${id}.maf"
    }

    runtime {
        docker: "broadinstitute/oncotator:1.9.0.0"
        memory: "${ram_gb}GB"
        cpu: "${cpu_cores}"
        disks: "local-disk ${output_disk_gb} HDD"
        bootDiskSizeGb: "${boot_disk_gb}"
        preemptible:  "${preemptible_tries}"
    }


    meta {
        author : "Chip Stewart"
        email : "stewart@broadinstitute.org"
    }

}

workflow maflite_oncotate_workflow {
  
    String pair_id
    String tumor_id
    String normal_id
    File maflite_file
    File interval_list
    File oncoDBTarBall

    Int preemptible_tries=4
    String output_disk_gb 
    String boot_disk_gb 


    call subset_tsv_by_intervallist {
      input:
        ID = pair_id,
        interval_list = interval_list,
        input_tsv_file = maflite_file,
        CHROMOSOME_FIELD = "chr",
        POSITION_FIELD = "start",
        preemptible_tries=preemptible_tries,
        output_disk_gb = output_disk_gb,
        boot_disk_gb = boot_disk_gb,
    }
        
    call oncotator_minimal {
      input:
        oncoDBTarBall = oncoDBTarBall,
        IN = subset_tsv_by_intervallist.subset_tsv,
        id = pair_id,
        output_disk_gb = output_disk_gb,
        boot_disk_gb = boot_disk_gb,
        preemptible_tries=preemptible_tries
    }

    output {
        oncotator_minimal.oncotator_out
        subset_tsv_by_intervallist.subset_tsv
    }
    
    
}