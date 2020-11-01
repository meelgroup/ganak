#!/bin/bash
# set -x

tout_arjun=100
tout_ganak=5000
maximum_component_cache_size=2000
input_file=$1
rm -f independent*
rm -f out*
rm -f newfile*
rm -f input*
rm -rf *.out
run_ganak=true
hash_range=1
while [[ $run_ganak == "true" ]]; do
    echo "c Trying to run Ganak on $input_file  with  timeout: ${tout_ganak}"
    # echo "c Command: ../bin/doalarm $tout_ganak ../bin/ganak -m $hash_range -maxdec 5000000 500 $input_file > output"
    `../bin/doalarm $tout_ganak ../bin/ganak -m $hash_range -cs $maximum_component_cache_size -maxdec 5000000 500 $input_file > output` > /dev/null 2>&1
    hash_error=`grep "ERROR: We need to change the hash range" output`
    if [[ $hash_error == *"ERROR: We need to change the hash range"* ]]; then
        let hash_range=2*hash_range
        continue
    fi
    run_ganak=false
done
useIS=`grep "Terminating solver because the number of decisions" output`
if [[ $useIS == *"Terminating solver because the number of decisions"* ]]; then
    echo "c Ganak's initial run is not sucessful, we will have to run arjun"
    rm -rf independent_support
    touch independent_support
    # echo "c Command: ../bin/doalarm $tout_arjun ../bin/arjun $input_file > independent_support"
    `../bin/doalarm $tout_arjun ../bin/arjun $input_file > independent_support` > /dev/null 2>&1
    touch newfile
    grep -v "^c ind" $input_file > newfile
    grep "^vp" independent_support >> newfile
    sed -i 's/vp/s c ind/g' newfile
   # mv newfile $input_file
    echo "c arjun's run is successful! Run ganak again."
    run_ganak=true
    input_file=newfile
fi
hash_range=1
while [[ $run_ganak == "true" ]]; do
    echo "c Trying to run Ganak on $input_file  with  timeout: ${tout_ganak}"
    # echo "c Command: ../bin/doalarm $tout_ganak ../bin/ganak -m $hash_range $input_file > output"
    `../bin/doalarm $tout_ganak ../bin/ganak -cs $maximum_component_cache_size -m $hash_range $input_file > output` > /dev/null 2>&1
    # cat output
    hash_error=`grep "ERROR: We need to change the hash range" output`
    if [[ $hash_error == *"ERROR: We need to change the hash range"* ]]; then
        let hash_range=2*hash_range
        continue
    fi
    run_ganak=false
done
sed -i 's/s pmc/s mc/g' output
found=`grep "s mc" output`
if [[ $found == *"s mc"* ]]; then
    cat output
    rm -rf independent_support
    rm -rf newfile
    rm -rf output
    exit 0
fi
rm -rf independent_support
rm -rf newfile
rm -rf output
echo "c Ganak did NOT work"
