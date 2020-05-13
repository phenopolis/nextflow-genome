# copy from .command.run
# aws helper
nxf_s3_upload() {
    local pattern=$1
    local s3path=$2
    IFS=$'\n'
    for name in $(eval "ls -1d $pattern");do
      if [[ -d "$name" ]]; then
        /home/ec2-user/miniconda/bin/aws $aws_profile --region eu-central-1 s3 cp --only-show-errors --recursive --storage-class STANDARD "$name" "$s3path/$name"
      else
        /home/ec2-user/miniconda/bin/aws $aws_profile --region eu-central-1 s3 cp --only-show-errors --storage-class STANDARD "$name" "$s3path/$name"
      fi
    done
    unset IFS
}

nxf_s3_retry() {
    local all=("$@")
    local max_attempts=5
    local timeout=10
    local attempt=0
    local exitCode=0
    local action
    local source=$2
    local target=$3
    while (( $attempt < $max_attempts ))
    do
      tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
      err_log=${tmp_dir}/err.log
      eval "${all[@]}" 2>${err_log}
      errcode=$?
      echo "error code is ${errcode}"
      err=$(<${err_log})
      rm -rf ${tmp_dir}
      #err=$("$@" 3>&1 1>&2 2>&3 | tee >(cat - >&2))
      if [[ $errcode == 0 ]]
      then
          #echo $err
          break
      elif [[ $err = download* ]]
      then
          bucket=$(echo $source | cut -d'/' -f 3)
          key=$(echo $source | cut -d'/' -f4-)
          action=download
          # refer to https://stackoverflow.com/questions/45068864/resuming-interrupted-s3-download-with-awscli
          if [ -f $target ]
          then
              echo "retrying---"
              echo $source
              echo $bucket
              echo $key
              size=$(stat --printf="%s" $target)
              /home/ec2-user/miniconda/bin/aws s3api get-object\
                  $aws_profile \
                  --bucket $bucket \
                  --key $key \
                  --range "bytes=$size-" \
                  /dev/fd/3 3>>$target
          fi
      else
          echo "Error in download. Error code: $errcode, Error message: $err. Will retry"
      fi
      sleep $timeout
      attempt=$(( attempt + 1 ))
      timeout=$(( timeout * 2 ))
    done
}

nxf_s3_download() {
    echo "in download"
    local source=$1
    local target=$2
    local all=("$@")
    echo "${all[@]}"
    local file_name=$(basename ${source})
    #local is_dir=$(/home/ec2-user/miniconda/bin/aws $aws_profile --region eu-central-1 s3 ls $source | grep -F "PRE ${file_name}/" -c)
    cmd="/home/ec2-user/miniconda/bin/aws $aws_profile --region eu-central-1 s3 ls $source | grep -F \"PRE ${file_name}/\" -c"
    echo $cmd
    local is_dir=$(eval $cmd)
    [[ $source = */ ]] && is_dir=1
    echo $is_dir
    if [[ $is_dir == 1 ]]; then
        cmd="/home/ec2-user/miniconda/bin/aws $aws_profile --region eu-central-1 s3 sync --only-show-errors ${all[@]}"
        echo $cmd
        eval $cmd
    elif [ ! -f $target ]; then
        bucket=$(echo $source | cut -d'/' -f 3)
        key=$(echo $source | cut -d'/' -f4-)

        echo $aws_profile
        echo $bucket
        echo $key
        cmd="/home/ec2-user/miniconda/bin/aws s3api get-object $aws_profile --bucket $bucket --key $key $target"
        echo $cmd
        eval "$cmd"

       echo "==worked=="
    fi
}

nxf_parallel() {
    IFS=$'\n'
    local cmd=("$@")
    local cpus=$(nproc 2>/dev/null || < /proc/cpuinfo grep '^process' -c)
    local max=$(if (( cpus>16 )); then echo 16; else echo $cpus; fi)
    local i=0
    local pid=()
    (
    set +u
    while ((i<${#cmd[@]})); do
        echo "command: ${cmd[$i]}"
        local copy=()
        for x in "${pid[@]}"; do
          [[ -e /proc/$x ]] && copy+=($x) 
        done
        pid=("${copy[@]}")

        if ((${#pid[@]}>=$max)); then 
          sleep 1 
        else 
          eval "${cmd[$i]}" &
          pid+=($!)
          ((i+=1))
        fi 
    done
    ((${#pid[@]}>0)) && wait ${pid[@]}
    )
    unset IFS
}
