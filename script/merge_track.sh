#! /bin/bash


if [ $# -ne 5 ]
then 
    printf "\n\nUsage: %s [input path] [file format] [output path] [start index] [end index]\n\n\n" $0
    exit
fi

ipath=$1
ffmt=$2
opath=$3
start_index=$4
end_index=$5

printf "\n"
printf "input path: %s\n"  $ipath
printf "output path: %s\n" $opath
printf "file format: %s\n" $ffmt
printf "index from '%d' to '%d'\n" $start_index $end_index
printf "\n"

for (( pid=start_index; pid <= end_index; pid +=1 ))
do
    result_file=$(printf "%s/P%07d.track" $opath $pid)
    if [ -e $result_file ]
    then
        rm $result_file
    fi

    fid=1
    filename=$(printf $ipath"/"$ffmt $pid $fid)
    while [ -e $filename ]
    do
        cat $filename >> $result_file
        fid=$(($fid + 1))
        filename=$(printf $ffmt $pid $fid)
    done
done
