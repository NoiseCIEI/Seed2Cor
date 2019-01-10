#!/bin/bash

function usage {
    echo "usage: $0 SEED_dir [out_dir]"
    echo "  SEED_dir: directory of SEED files"
    echo "  out_dir: absolute path to directory to output list of stations"
    exit 1
}

if [ $# -eq 0 ]; then
    usage
fi

seeddir=$1
if [ $# -eq 1 ]; then
    lstdir=.
else
    lstdir=$2
fi

if [ -e station.lst ]; then
    rm -f station.lst
fi
cd $seeddir
rm s1.lst
ls *seed >> seed_fname.lst
for sname in `cat seed_fname.lst`; do
    rdseed -S -f $sname
    cat rdseed.stations >> s1.lst
done
sort -k1,1 -u s1.lst > s2.lst
awk '{print $1,$4,$3}' s2.lst > station.lst
rm s1.lst s2.lst seed_fname.lst
cp station.lst ${lstdir}
