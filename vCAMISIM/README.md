# vCAMISIM 

Scripts for recreating viral benchmark results using genome simulations from CAMISIM. 

## Example run 
```bash
python CAMISIM_viral_benchmark_cmd.py <clusterfile> <pooled_gsa> <directoryout>
python CAMISIM_viral_benchmark_cmd.py vamb_clusters.tsv pooled_gsa.txt cami_benchmark
```
clusterfile - works with VAMB and metaBAT2 cluster files 
pooled_gsa - Basically all `anonymous_gsa.fasta.gz` from each sample in the simulation concatenated (keep only one header!)

Of course this benchmark will not properly without specifiying viral contig annotation as done in the `CAMISIM_viral_benchmark_cmd.py`. 
```bash 
virus_annotations/seeker.txt
virus_annotations/virsorter2/final-viral-score.tsv
virus_annotations/virfinder.txt
virus_annotations/viralverify/pooled_contigs.2000.tmp_result_table.csv
virus_annotations/DVF/pooled_contigs.2000.tmp.fna_gt2000bp_dvfpred.txt
contig_annotations/all.hmmVOG.tbl
contig_annotations/all.hmmMiComplete105.tbl
```
The script can easily be modified to only include predictions for some virus predictors (just remove them as input in the code).

## Output files
`benchmark` - Binning performance directory. Most insightful file is  the GenomeRecallStats.tsv file that summmarise the recovery of each binnde genome.
`viral_benchmark_contig` -  Prediction performance on contig level only for the listed tools above
`viral_benchmark_bin` - Prediction performance on genome level (including the RF model)


## Where can I find the Viral predictors used here? 

- Seeker (https://github.com/gussow/seeker)
- Virsorter2 (https://github.com/jiarong/VirSorter2)
- Virfinder (https://github.com/jessieren/VirFinder)
- DeepVirfinder (https://github.com/jessieren/DeepVirFinder)
- Viralverify (https://github.com/ablab/viralVerify)
 
