#!/bin/bash

basecamp=/home/wayyne/wrk/pdb_doctor/scripts

#strip ligands from all PDBs
for pdb in *_pp.pdb; do
  awk '$1 == "ATOM" {print $0}' $pdb > ${pdb::4}_pp_stripped.pdb
done

#get list of complete and incomplete pdbs
python ${basecamp}/extract_sequence.py . _pp_stripped.pdb && \

#run the incomplete pdbs through fold and transfer
awk -v basecamp=${basecamp} ' 
NR > 1 {
  print " python "basecamp"/fold_and_transfer.py",$1,$2,substr($1,1,5)"triaged.pdb"
}
' incomplete.tsv | bash && \

#copy the complete pdbs w/o HETATMs (do not need to fold or transfer)
awk -v basecamp=${basecamp} ' 
NR > 1 {
  print "grep -v 'HETATM'",$1,">",substr($1,1,5)"triaged.pdb"
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

mkdir triage_log
mv *.txt triage_log
mv *.log triage_log
mv *.dat triage_log

#return to parent directory
cd -

#finish housekeeping
rm -rf incomplete_triaged_pdbs pdb_cache fitted_folded.pdb *_pp_stripped.pdb


