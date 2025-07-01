source activate ArchR

Rscript "./archr_singlesample_V1.R" \
    --fragment  ./fragments_corrected_dedup_count.tsv.gz \
    --spname LD_50W_LWF_1 \
    --threads 1 \
    --prefix LD_50W_LWF_1 \
    --species chicken \
    --outdir ./  \
    --resolution 0.8 \
    --dim 8 \
    --scRNA no \
    --testMethod wilcoxon \
    --peakcells 75 \
    --startPoint 5 --archrdir ./ --fdr 0.651391556953224 
  