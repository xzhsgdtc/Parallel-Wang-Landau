#! /bin/bash

if [ $# -n 1 ]
then
    echo "Usage: ./script [proc id]"
    exit
fi


proc=P00000$(printf "%02d" $1)


count=1
flag=true

while $flag
do
    if [ -e "${proc}.histogram_data_00000$(printf "%02d" $count).dat" ]
    then
        count=$(($count+1))
    else
        count=$(($count-1))
        flag=false
    fi
done

echo ${proc}.histogram_data_00000$(printf "%02d" $count).dat

