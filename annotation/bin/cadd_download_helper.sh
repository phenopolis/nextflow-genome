overlap() {
    # overlap, returns 1
    # else, returns 0
    local bed_start=$1
    local bed_end=$2
    local remote_file=$3
    regex="^([0-9]+)-([0-9]+)[^0-9]"
    if [[ ${remote_file} =~ $regex ]]; then
        remote_start=${BASH_REMATCH[1]}
        remote_end=${BASH_REMATCH[2]}

        separate_total=$(( ${bed_end}-${bed_start}+${remote_end}-${remote_start} ))
        max_end=$(( ${bed_end} > ${remote_end} ? ${bed_end} : ${remote_end} ))
        min_start=$(( ${bed_start} < ${remote_start} ? ${bed_start} : ${remote_start} ))
        merged_total=$(( ${max_end} - ${min_start} ))
        if [ "${separate_total}" -lt "${merged_total}" ]; then
            echo 0
        else
            echo 1
        fi
    else
        echo 0
    fi
}

perm ()
{
    outarray=( "$@" )

    # The algorithm used is the Fisher-Yates Shuffle
    # (https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle),
    # also known as the Knuth Shuffle.

    # Loop down through 'outarray', swapping the item at the current index
    # with a random item chosen from the array up to (and including) that
    # index
    local idx rand_idx tmp
    for ((idx=$#-1; idx>0 ; idx--)) ; do
        rand_idx=$(( RANDOM % (idx+1) ))
        # Swap if the randomly chosen item is not the current item
        if (( rand_idx != idx )) ; then
            tmp=${outarray[idx]}
            outarray[idx]=${outarray[rand_idx]}
            outarray[rand_idx]=$tmp
        fi
    done
}

