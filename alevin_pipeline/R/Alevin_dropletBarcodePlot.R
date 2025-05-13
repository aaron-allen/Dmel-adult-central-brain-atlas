# funtions to plot standard barcode rank plot and density plot using dropletutils and royryk cutoff.
hist_bins <- 50

# Pick a cutoff on count as per https://github.com/COMBINE-lab/salmon/issues/362#issuecomment-490160480

pick_roryk_cutoff = function(bcs){
  bcs_hist = hist(log10(bcs), plot=FALSE, n=hist_bins)
  mids = bcs_hist$mids
  vals = bcs_hist$count
  wdensity = vals * (10^mids) / sum(vals * (10^mids))
  baseline <- median(wdensity)
  
  # Find highest density in upper half of barcode distribution
  
  peak <- which(wdensity == max(wdensity[((length(wdensity)+1)/2):length(wdensity)]))
  
  # Cutoff is the point before the peak at which density falls below 2X baseline
  
  10^mids[max(which(wdensity[1:peak] < (1.5*baseline)))]
}

# Plot densities 

barcode_density_plot = function(bcs, roryk_cutoff, knee, inflection,emptydrops_retain, name = 'no name') {
  bcs_hist = hist(log10(bcs), plot=FALSE, n=hist_bins)
  counts = bcs_hist$count
  mids = bcs_hist$mids
  y = counts * (10^mids) / sum(counts * (10^mids))
  qplot(y, 10^mids) + geom_point() + theme_bw() + ggtitle(name) + ylab('Count') + xlab ('Density') +
    geom_hline(aes(yintercept = roryk_cutoff, color = paste('roryk_cutoff =', length(which(bcs > roryk_cutoff)), 'cells'))) + 
    geom_hline(aes(yintercept = inflection, color = paste('dropletutils_inflection =', length(which(bcs > inflection)), 'cells'))) +
    geom_hline(aes(yintercept = knee, color = paste('dropletutils_knee =', length(which(bcs > knee)), 'cells'))) +
    geom_hline(aes(yintercept = emptydrops_retain, color = paste('emptydrops_retain =', length(which(bcs >emptydrops_retain)), 'cells'))) +
    scale_y_continuous(trans='log10') + theme(axis.title.y=element_blank(), legend.text = element_text(size=8)) + labs(color='Thresholds') + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(size=12))
}  

# Plot a more standard barcode rank plot

barcode_rank_plot <- function(br.out, roryk_total_cutoff, knee, inflection,emptydrops_retain, name='no name'){
  ggplot(data.frame(br.out), aes(x=rank, y=total)) + geom_line() +scale_x_continuous(trans='log10',breaks = c(1,100,1000,10000,50000,100000,max(br.out$rank))) +
    scale_y_continuous(trans='log10', breaks = c(1,10,100,1000,10000,100000)) + theme_bw() +  
    geom_hline(aes(yintercept = knee, color = 'dropletutils_knee')) + 
    geom_hline(aes(yintercept = inflection, color = 'dropletutils_inflection')) +
    geom_hline(aes(yintercept = roryk_total_cutoff, color = 'roryk_cutoff')) +
    geom_hline(aes(yintercept = emptydrops_retain, color = 'emptydrops_retain')) +
    ggtitle(name) + ylab('Total UMI Count') + xlab('Rank') + theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(size=12))
}





