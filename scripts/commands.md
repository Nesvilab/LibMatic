# LibMatic pipeline
## DIA data processing
```
cd DIA/
pipeline.pyz pipeline_parameters_DIA.txt program_paths.txt --options 'memory = 100.0 GiB' 'parallelism = 18'
```
Program arguments:
1. Parameter file
2. File with program paths
3. Memory usage
4. CPU cores to be used
## DDA data processing
```
cd DDA/
pipeline.pyz pipeline_parameters_DDA.txt program_paths.txt --options 'memory = 100.0 GiB' 'parallelism = 18' --DDA
```
Program arguments:
Parameter file
1. File with program paths
2. Memory usage
3. CPU cores to be used
4. Data is DDA

## For uncombined DDA/DIA library
```
LIBGEN_dir=DDA_lib
datadir=./workdir_DDA
fasta_file=./td-up000005640.fasta
```
## For combined (DDA+DIA) library
```
LIBGEN_dir=DIA_DDA_lib
DIA_datadir=./workdir_DIA
DDA_datadir=./workdir_DDA
fasta_file=./td-up000005640.fasta
```

## Retention time alignment
```
RTalign.py ${LIBGEN_dir}/RT_align ${datadir}/iproph none /app/teog/tpp5/bin/
RTalign.py ${LIBGEN_dir}/RT_align ${DIA_datadir}/iproph ${DDA_datadir}/iproph /app/teog/tpp5/bin/
```

Program arguments:
1. Directory to store the RT aligned files
2. Directory of pep.xml files
3. Directory of pep.xml files for DDA (input “none” if not aligning with DDA)
4. Directory of TPP binaries

## Combining DIA and DDA pep.xml files
```
combine_prot_xml.py ${LIBGEN_dir}/combined_prots ${DIA_datadir}/iproph ${DDA_datadir}/iproph /app/teog/tpp5/bin/
```
## Generating spectral library
Run for combined DIA and DDA library
```
PATH=".:$PATH"
python3 bin/\
gen_con_spec_lib.py \
  ${fasta_file} \
  ${LIBGEN_dir}/RT_align/RT_aligned \
  ${LIBGEN_dir}/combined_prots/interact.prot.xml \
  ${LIBGEN_dir}/DIA_lib true &
python3 bin/\
gen_con_spec_lib.py \
  ${fasta_file} \
  ${LIBGEN_dir}/RT_align/RT_aligned_DDA \
  ${LIBGEN_dir}/combined_prots/interact.prot.xml \
  ${LIBGEN_dir}/DDA_lib true &
wait
```
For uncombined DIA/DDA library
```
PATH=".:$PATH" python3 bin/gen_con_spec_lib.py \
  ${fasta_file} \
  ${LIBGEN_dir}/RT_align/RT_aligned \
  ${datadir}/ProteinProphet/interact.prot.xml \
  ${LIBGEN_dir} true
```
Program arguments:
1. Directory of TPP binaries
2. Path of FASTA file
3. Data directory
4. RT aligned files
5. Output files directory



## Combining DIA and DDA generated library
```
# python3 -m pip -- install matplotlib_venn
(cd ${LIBGEN_dir} && 
python3 ../bin/\
combine_DIA_DDA_libs.py \
  ./DIA_lib/con_lib.tsv \
  ./DDA_lib/con_lib.tsv
)
```
Program arguments:
1. DIA library
2. DDA library
