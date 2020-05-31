#!/bin/bash
# set -x

tout_be=100
tout_ganak=0
input_file=$1
rm -f proj*
rm -f out*
rm -f newfile*
rm -f input*

echo "c Trying to run Ganak on $input_file  with  timeout: ${tout_ganak}"
echo "../bin/ganak -t $tout_ganak $input_file > output"
`../bin/ganak $input_file > output` > /dev/null 2>&1
sed -i 's/s pmc/s mc/g' output
found=`grep "s mc" output`
cat output
if [[ $found == *"s mc"* ]]; then
    cat output
    exit 0
fi
echo "c Ganak did NOT work"

# echo "c Getting indep support, timeout: ${tout_be}"
# grep -v "^c" - | sed "s/pcnf/cnf/" | grep -v "^vp" > inputfile
# rm -f projection
# touch projection
# `./bin/doalarm ${tout_be_relaxed} ./bin/b_plus_e -B=projection -cpu-lim=${tout_be} inputfile` > /dev/null 2>&1
# sed -i "s/V/vp/" projection
# found=`grep "vp .* 0$" projection`

# echo "c found is: $found"
# if [[ $found == *"vp"* ]]; then
#     echo "c OK, B+E succeeded"
#     cp inputfile newfile2
#     cat projection >> newfile2
# else
#     echo "c WARNING B+E did NOT succeed"
#     cp inputfile newfile2
# fi

# sed 's/vp/c ind/g' newfile2 | sed 's/ pcnf/ cnf/g' > newfile3

# echo "c Trying a short run for appromc, timeout: ${tout_shortapproxmc}"

# `./bin/doalarm ${tout_shortapproxmc} ./bin/approxmc --delta 0.2 newfile3 > output` > /dev/null 2>&1
# found=`grep "s mc" output`
# if [[ $found == *"s mc"* ]]; then
#     cat output
#     exit 0
# fi

# echo "c Running ApproxMC"
# ./bin/approxmc --delta 0.2 newfile3
# exit 0