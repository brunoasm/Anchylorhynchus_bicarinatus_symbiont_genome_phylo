
#!/bin/bash

set +e
source /home/bdemedeiros/miniconda3/etc/profile.d/conda.sh
conda activate gtdbtk-2.1.1

gtdbtk classify_wf --min_perc_aa 20 --full_tree --extension fasta --cpus 20 --pplacer_cpus 20 --genome_dir meta_assembly --out_dir gtdbtk_out
