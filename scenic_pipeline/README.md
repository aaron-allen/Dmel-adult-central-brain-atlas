**Pipeline to run SCENIC analysis (pySCENIC)**

Authors: Lucy Garner (scenic_multi and scenic) and Devika Agarwal (scenic)

Only pipeline_meta.py should be run directly as follows:

`python pipeline_meta.py --cluster-queue-manager=slurm --cluster-queue=long make full`

`nohup nice -10 python pipeline_meta.py --cluster-queue-manager=slurm --cluster-queue=long make full &> pipeline_extra.log < /dev/null &`
`nohup nice -10 python pipeline_meta.py --no-cluster make full &> local_pipeline.log < /dev/null &`










# run in 4 steps

Run:	`nohup nice -10 python pipeline_meta.py --cluster-queue-manager=slurm make gene_filtering dummy &> local_pipeline.log < /dev/null &`
Run:	`mv pipeline.log pipeline.log.step1 && mv local_pipeline.log local_pipeline.log.step1`



Run:	`nohup nice -10 python pipeline_meta.py --no-cluster make pyscenic_run &> local_pipeline.log < /dev/null &`        -  current pid 1748743
Check:	`ll results.dir/*/*/*/adjacencies.tsv | wc -l`
		`ll results.dir/*/*/*/adjacencies_wCor.csv | wc -l`
		`ll results.dir/*/*/*/reg.csv | wc -l`
		`ll results.dir/*/*/*/aucell.csv | wc -l`
Run:	`mv pipeline.log pipeline.log.step2 && mv local_pipeline.log local_pipeline.log.step2`



Run:	`nohup nice -10 python pipeline_meta.py --cluster-queue-manager=slurm --cluster-queue=long make generate_regulons regulon_binarization rss_zscore aucell_heatmap &> local_pipeline.log < /dev/null &`
Check:	`ll results.dir/*/*/*/regulons.csv | wc -l`
		`ll results.dir/*/*/*/aucell_thresholds.csv | wc -l`
		`ll results.dir/*/*/*/aucell_zscores.csv | wc -l`
		`ll results.dir/*/*/*/plots.dir/aucell_heatmap.png | wc -l`
Run:	`mv pipeline.log pipeline.log.step3 && mv local_pipeline.log local_pipeline.log.step3`
		

		
Run:	`nohup nice -10 python pipeline_meta.py --cluster-queue-manager=slurm --cluster-queue=long make high_confidence_regulons aggregated_aucell aggregated_regulon_binarization aggregated_rss_zscore aggregated_aucell_heatmap &> local_pipeline.log < /dev/null &`
Check:	`ll results.dir/aggregated.dir/*/*/high_confidence_regulons.csv | wc -l`
		`ll results.dir/aggregated.dir/*/*/aucell.csv | wc -l`
		`ll results.dir/aggregated.dir/*/*/aucell_thresholds.csv | wc -l`
		`ll results.dir/aggregated.dir/*/*/aucell_zscores.csv | wc -l`
		`ll results.dir/aggregated.dir/*/*/plots.dir/aucell_heatmap.png | wc -l`
Run:	`mv pipeline.log pipeline.log.step4 && mv local_pipeline.log local_pipeline.log.step4`

















# other re-runs

nohup nice -10 python pipeline_meta.py --cluster-queue-manager=slurm --cluster-queue=long make generate_regulons &> local_pipeline.log < /dev/null &
nohup nice -10 python pipeline_meta.py --cluster-queue-manager=slurm --cluster-queue=long make regulon_binarization rss_zscore aucell_heatmap &> local_pipeline.log < /dev/null &
nohup nice -10 python pipeline_meta.py --cluster-queue-manager=slurm --cluster-queue=long make aucell_heatmap &> local_pipeline.log < /dev/null &



# agg
nohup nice -10 python pipeline_meta.py --cluster-queue-manager=slurm --cluster-queue=long make high_confidence_regulons aggregated_aucell aggregated_regulon_binarization aggregated_rss_zscore aggregated_aucell_heatmap &> local_pipeline.log < /dev/null &












