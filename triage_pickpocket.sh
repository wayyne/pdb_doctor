#!/bin/bash

basecamp=/home/wayyne/wrk/pdb_doctor

python ${basecamp}/extract_sequence.py .

awk -v basecamp=${basecamp} ' 
NR > 1 {
  print " python "basecamp"/fold_and_transfer.py",$1,$2,substr($1,1,8)"D.pdb"
}
' incomplete.tsv | bash

mkdir input
mv *_D.pdb input/.

time python ${basecamp}/complete_pdb.py --mode partial --input input --output output

#update contact labels with sticky fingers
time python ${basecamp}/sticky_fingers.py
