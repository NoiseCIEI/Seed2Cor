#!/bin/bash

function usage {
    echo "usage: $0 SEED_dir [out]"
    echo "  SEED_dir: ABSOLUTE path to directory of SEED files"
    echo "  out: ABSOLUTE path to where to output list of SEEDs"
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

if [ -e seed.lst ]; then
    rm -f seed.lst
fi
cd $seeddir
ls *seed >> seed_fname.lst
awk 'BEGIN {FS="."} ; {print $2,$3}' seed_fname.lst >monday.lst
awk 'BEGIN {FS="."} ; {print $1}' seed_fname.lst  > year.lst
awk 'BEGIN {FS="_"} ; {print $2}' year.lst  > year2.lst

cat monday.lst |  sed s/'JAN'/'1'/ | sed s/'FEB'/'2'/ | sed s/'MAR'/'3'/ | sed s/'APR'/'4'/ | sed s/'MAY'/'5'/ | sed s/'JUN'/'6'/ | sed s/'JUL'/'7'/ | sed s/'AUG'/'8'/ | sed s/'SEP'/'9'/ | sed s/'OCT'/'10'/ | sed s/'NOV'/'11'/ | sed s/'DEC'/'12'/ > monday2.lst
for fname in `cat seed_fname.lst`; do
    if [ "${fname: -1}" == "/" ]; then
        echo "${seeddir}${fname}" >> seed_fname2.lst
    else
        echo "${seeddir}/${fname}" >> seed_fname2.lst
    fi
done
paste seed_fname2.lst year2.lst monday2.lst > seed1.lst
awk '{print $1, $2, $3, $4 }' seed1.lst >seed.lst
rm year.lst year2.lst monday.lst monday2.lst seed_fname.lst seed_fname2.lst seed1.lst
mv seed.lst $lstdir
