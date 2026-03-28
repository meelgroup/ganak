#!/bin/bash

rm -rf CMake* src cmake* ganak* sharp* Make*
cmake -DCMAKE_BUILD_TYPE=Release -DSTATICCOMPILE=ON ..
make -j$(nproc)
strip ganak

# Get SHA1 hashes and build versioned filename
HASHES=$(echo "" | ./ganak 2>&1)
GANAK_SHA=$(echo "$HASHES" | grep "Ganak SHA1" | grep -oP '[0-9a-f]{40}' | cut -c1-8)
ARJUN_SHA=$(echo "$HASHES" | grep "Arjun SHA1" | grep -oP '[0-9a-f]{40}' | cut -c1-8)
CMS_SHA=$(echo "$HASHES" | grep "CMS SHA1" | grep -oP '[0-9a-f]{40}' | cut -c1-8)
APPROXMC_SHA=$(echo "$HASHES" | grep "ApproxMC SHA1" | grep -oP '[0-9a-f]{40}' | cut -c1-8)

DEST="ganak_${GANAK_SHA}_${ARJUN_SHA}_${APPROXMC_SHA}_${CMS_SHA}"
cp ganak "$DEST"
echo "Copied to $DEST"
