
# conda install -c bioconda blast

psiblast -query fon-PB_O.fasta \
         -db ~/SharedFiles/Wasi/BigDatasets/Databases/nr/nr \
         -num_iterations 3 \
         -out_ascii_pssm fon-PB_O.pssm \
         -num_threads 64
# 