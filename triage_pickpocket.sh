#!/bin/bash

basecamp=/home/wayyne/wrk/pdb_doctor

#get list of complete and incomplete pdbs
python ${basecamp}/extract_sequence.py . _pp.pdb && \

#run the incomplete pdbs through fold and transfer
awk -v basecamp=${basecamp} ' 
NR > 1 {
  print " python "basecamp"/fold_and_transfer.py",$1,$2,substr($1,1,5)"triaged.pdb"
}
' incomplete.tsv | bash && \

#copy the complete pdbs (do not need to fold or transfer)
awk -v basecamp=${basecamp} ' 
NR > 1 {
  print "cp "$1,substr($1,1,5)"triaged.pdb"
}
' complete.tsv | bash && \

#collect the 'triaged' pdbs
mkdir -p incomplete_triaged_pdbs && mv *_triaged.pdb incomplete_triaged_pdbs/. && \

#run the 'triaded' pdbs throug complete PDB to add missing atoms and optimize newly added side-chains
time python ${basecamp}/complete_pdb.py --mode partial --input incomplete_triaged_pdbs --output _triaged_pdbs && \

#confirm/refine contact labels with sticky fingers
time python ${basecamp}/sticky_fingers.py _triaged_pdbs && \

#oraganize the complete triaged (aka chimeric) PDBs
cd _triaged_pdbs
mkdir series_C && mv *_C.pdb series_C/.

#housekeeping
rm -rf *_A.pdb *_B.pdb
mv ../*.tsv .
mv ../*.log .
mv ../sticky*txt .

#return to parent directory
cd -

#finish housekeeping
rm -rf incomplete_triaged_pdbs pdb_cache fitted_folded.pdb
