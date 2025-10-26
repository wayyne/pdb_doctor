#!/bin/bash

basecamp=/home/wayyne/lab/pdbdr/scripts

#get list of complete and incomplete pdbs
python ${basecamp}/extract_sequence.py . _pp.pdb && \

# run the incomplete pdbs through fill and transfer
awk -v basecamp="${basecamp}" '
NR > 1 {
  partial = $1
  chain   = $2
  pdbid   = substr($1, 1, 4)           # e.g., "6o83" from "6o83_pp.pdb"
  outpdb  = substr($1, 1, 5) "triaged.pdb"  # e.g., "6o83_triaged.pdb"

  # build the command; quote vars to be safe
  printf "python \"%s/fill_and_transfer.py\" \"%s\" %s \"%s\" --chain %s\n",
         basecamp, partial, pdbid, outpdb, chain
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

#run the 'triaged' pdbs through complete PDB to add missing atoms and optimize newly added side-chains
time python ${basecamp}/complete_pdb.py --mode partial --input incomplete_triaged_pdbs --output _triaged_pdbs && \

#confirm/refine contact labels with sticky fingers
time python ${basecamp}/sticky_fingers.py . _triaged_pdbs && \

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
rm -rf incomplete_triaged_pdbs pdb_cache fitted_folded.pdb
