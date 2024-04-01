cd /storage/climatestor/PleioCEP/data/WRF_MedCyclones
for dir in postproc/Long*/RAINNC; do
    echo $dir
    cd $dir
    arr=(*)
    arr_len=${#arr[@]}
    echo ${arr_len}
    for i in $(seq 1 $((arr_len - 1))); do
        j=$((i - 1))
        date1=$(date -d "${arr[j]:11:10}" +%s)
        date2=$(date -d "${arr[i]:11:10}" +%s)
        difference=$((date2 - date1))
        days=$((difference / 86400))
        if [[ $days != 13 && $days != 14 && $days != 15 ]]; then
            echo ${arr[j]} ${arr[i]}
            echo $days
        fi
    done
    cd ../../..
done

