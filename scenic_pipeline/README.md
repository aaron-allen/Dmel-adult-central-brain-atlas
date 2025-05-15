# **Pipeline to run SCENIC analysis (pySCENIC)**

<br>

Authors: Lucy Garner (scenic_multi and scenic) and Devika Agarwal (scenic)
\
Modified by: Aaron M. Allen

<br>

Only pipeline_meta.py should be run directly as follows:

```
python pipeline_meta.py --cluster-queue-manager=slurm --cluster-queue=long make full
```



<br>
<br>




# Or ...

## Run in 4 steps

I have found that the way this pipeline was initially packaged is too slow as it uses a separate job (and 1 CPU core) to run and monitor each instance of the 100 scenic runs. On our cluster this adds significan delays. So I've split this up into 4 parts. The 2nd step (running all the 100 instances) is started with the wrapper pipeline running on the head node and submitting all the individual scenic runs to the cluster.


### Step 1: filtering

This is quick and can be run locally or on the head node if in a pinch. When complete I move the log files, as to not write over it.

```
# Run:
nohup nice -10 python pipeline_meta.py --cluster-queue-manager=slurm make gene_filtering dummy &> local_pipeline.log < /dev/null &
mv pipeline.log pipeline.log.step1 && mv local_pipeline.log local_pipeline.log.step1
```



### Step 2: 100 runs of scenic

```
# Run:
nohup nice -10 python pipeline_meta.py --no-cluster make pyscenic_run &> local_pipeline.log < /dev/null &
mv pipeline.log pipeline.log.step2 && mv local_pipeline.log local_pipeline.log.step2

# Check:
ll results.dir/*/*/*/adjacencies.tsv | wc -l
ll results.dir/*/*/*/adjacencies_wCor.csv | wc -l
ll results.dir/*/*/*/reg.csv | wc -l
ll results.dir/*/*/*/aucell.csv | wc -l
```



### Step 3: Generate regulons etc

```
# Run:
nohup nice -10 python pipeline_meta.py --cluster-queue-manager=slurm --cluster-queue=long make generate_regulons regulon_binarization rss_zscore aucell_heatmap &> local_pipeline.log < /dev/null &
mv pipeline.log pipeline.log.step3 && mv local_pipeline.log local_pipeline.log.step3

# Check:
ll results.dir/*/*/*/regulons.csv | wc -l
ll results.dir/*/*/*/aucell_thresholds.csv | wc -l
ll results.dir/*/*/*/aucell_zscores.csv | wc -l
ll results.dir/*/*/*/plots.dir/aucell_heatmap.png | wc -l
```


### Step 4: Aggregation

```		
# Run:
nohup nice -10 python pipeline_meta.py --cluster-queue-manager=slurm --cluster-queue=long make high_confidence_regulons aggregated_aucell aggregated_regulon_binarization aggregated_rss_zscore aggregated_aucell_heatmap &> local_pipeline.log < /dev/null &
mv pipeline.log pipeline.log.step4 && mv local_pipeline.log local_pipeline.log.step4

# Check:
ll results.dir/aggregated.dir/*/*/high_confidence_regulons.csv | wc -l
ll results.dir/aggregated.dir/*/*/aucell.csv | wc -l
ll results.dir/aggregated.dir/*/*/aucell_thresholds.csv | wc -l
ll results.dir/aggregated.dir/*/*/aucell_zscores.csv | wc -l
ll results.dir/aggregated.dir/*/*/plots.dir/aucell_heatmap.png | wc -l
```