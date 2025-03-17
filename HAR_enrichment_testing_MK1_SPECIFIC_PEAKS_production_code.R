########
##HAR enrichment testing for Senevirathne 2025
########

#find . -name "*narrowPeak" > narrows.txt

narrow_files <- readLines("narrows.txt")

library(data.table)
narrow_dat <- lapply(narrow_files, fread)

narrow_set <- list()
for (x in 1:length(narrow_dat)) {
	curr_dat <- narrow_dat[[x]]
	curr_dist <- curr_dat$V3 - curr_dat$V2
	curr_frame <- data.frame(file = narrow_files[x], curr_dist)
	narrow_set[[x]] <- curr_frame
}

##convert to bed files for future use.
for (x in 1:length(narrow_dat)) {
    fwrite(narrow_dat[[x]][,1:3], gsub(".narrowPeak", ".bed", narrow_files[x]), sep = "\t", row.names = F, col.names = F)
    print(narrow_files[x])
}

big_collapse <- do.call(rbind, narrow_set)
big_collapse$TIME <- "E53"
big_collapse$TIME[grepl("E57", big_collapse$file)] <- "E57"
big_collapse$TIME[grepl("E67", big_collapse$file)] <- "E67"
big_collapse$TIME[grepl("E72", big_collapse$file)] <- "E72"

pad_bed_slop <- function(filename, padding) {
    outname <- gsub(".bed", paste("_padded", padding, ".bed", sep=""), filename)
	##awk_command:
	#awk '{center=int(($2+$3)/2); print $1"\t"center"\t"(center+1)"\t"$4}' $1 > $2
	##get centres of peaks
    system(paste0("./awk_command.sh ", filename, " ", gsub(".bed", "_centres.bed", filename)))
	##expand peaks to padding
    system(paste0("bedtools slop -i ", gsub(".bed", "_centres.bed", filename), " -g hg38.chrom.clean.sizes -b ", padding, " > ", outname))
    return(outname)
}

for (TIME in unique(big_collapse$TIME)) {
    curr_files <- unique(big_collapse$file[big_collapse$TIME == TIME])
    print(curr_files)
    print(TIME)
}

##now, make files where I pool everything together....
for (TIME in unique(big_collapse$TIME)) {
    curr_files <- unique(big_collapse$file[big_collapse$TIME == TIME])
    curr_data <- lapply(curr_files, fread)

    ##identify cell-type specific peaks which are mutually-exclusive from all other cell-types
    try(dir.create(paste0(TIME, "_celltype_specific")))
    outfile_set <- list()
    for (x in 1:length(curr_files)) {
        curr_bool <- rep(TRUE, length(curr_files))
        curr_bool[x] <- FALSE
        all_others <- do.call('rbind', curr_data[curr_bool])
        ##write to tmpfs /thing
        try(unlink("/thing/temp.bed"))
        fwrite(all_others[, 1:3], "/thing/temp.bed", sep = "\t", row.names = F, col.names = F)
        curr_file <- fread(curr_files[x])
        try(unlink("/thing/temp_target.bed"))
        fwrite(curr_file[, 1:3], "/thing/temp_target.bed", sep = "\t", row.names = F, col.names = F)
        outfile <- paste0(TIME, "_celltype_specific/", gsub(".narrowPeak", "_specific.bed", basename(curr_files[x])))
        system(paste0("bedtools intersect -v -a /thing/temp_target.bed ", " -b /thing/temp.bed > ", outfile))
        print(outfile)
        outfile_set[[toString(x)]] <- outfile 
    }

    print("done this timepoint")
    print(TIME)
    ##for each timepoint, pool all the cell-type specific peaks together
    ##then generate randomized peaksets matching the median # of peaks and size.
    celltype_res <- do.call('rbind', lapply(outfile_set, fread))
    median_dist <- median(celltype_res$V3 - celltype_res$V2)
    ##round to the nearest 50
    median_dist <- ceiling(median_dist/50)*50

    fwrite(celltype_res[,1:3], paste0(TIME, "_celltype_specific/", TIME, "_celltype_specific_pooled_peaks.bed"), sep = "\t", row.names = F, col.names = F)
    ##then I'll need to pad, then merge...
    curr_peak_size <- median_dist
    pad_file <- pad_bed_slop(paste0(TIME, "_celltype_specific/", TIME, "_celltype_specific_pooled_peaks.bed"), curr_peak_size)
    system(paste0("bedtools sort -i ", pad_file, " | bedtools merge -i - > ", gsub(".bed", "_sort_merge.bed", pad_file)))
    print(TIME)
}


gen_feature_set_bedtools <- function(feature, chromosomefile, number, setnum=NULL, average_size =NULL) {
	
	featureset <- read.table(feature)
	
    if(is.null(average_size)) {
    average_size <- mean(as.numeric(featureset[,3]) - as.numeric(featureset[,2]))
	average_size <- ceiling(average_size)
    }
	
	library(parallel)
	cl <- makeCluster(20)
	clusterExport(cl, c("chromosomefile", "featureset", "average_size", "number"), environment())

    if(is.null(setnum)) {
setnum <- dim(featureset)[1]
		
    }
    clusterExport(cl, "setnum", environment())
		outdir <- unlist(strsplit(feature, "/"))[length(unlist(strsplit(feature, "/")))]
		outdir <- gsub(".bed", "", outdir)
		try(dir.create(outdir))
		parLapply(cl, 1:number, function(x) system(paste0("bedtools random -l ", average_size, " -n ", setnum, " -g ", chromosomefile, " > ", outdir, "/", paste0("random_set_padding_", average_size, "_size_", setnum, "_", x, ".bed"))))

	stopCluster(cl)
	gc()
}

fileset <- system("ls */*specific.bed", intern = T)
##get the median # of peaks at each timepoint.
for (age in c("E53", "E57", "E67", "E72")) {
    curr_sets <- fileset[grepl(age, fileset)]
    curr_sizes <- unlist(lapply(curr_sets, function(x) system(paste0("wc ", x, " | awk '{print $1}'"), intern = T)))
    curr_sizes <- as.numeric(curr_sizes)
    print(paste0(age, " - ", median(curr_sizes)))
}

##generate 10K randomized peaksets, matching the median # of peaks and peaksize.

gen_feature_set_bedtools("E53_celltype_specific_pooled_peaks_padded250_sort_merge.bed", "../hg38.chrom.clean.sizes", 10000, 11675, 500)
gen_feature_set_bedtools("E57_celltype_specific_pooled_peaks_padded200_sort_merge.bed", "../hg38.chrom.clean.sizes", 10000, 8288, 400)
gen_feature_set_bedtools("E67_celltype_specific_pooled_peaks_padded200_sort_merge.bed", "../hg38.chrom.clean.sizes", 10000, 5834, 400)
gen_feature_set_bedtools("E72_celltype_specific_pooled_peaks_padded200_sort_merge.bed", "../hg38.chrom.clean.sizes", 10000, 5672, 400)

calculate_target_intersections_cleaned_FEB2025 <- function(target_bed_set, region_bed, actual_bp=F) {
	
	target_lines <- readLines(target_bed_set)
	
    ##if it's a table, where I organize sample names cleanly
	if(!identical(grep("\t", target_lines), integer(0))) {
		target_set <- read.table(target_bed_set, sep="\t")
	target_files <- as.character(target_set[,2])
	target_names <- as.character(target_set[,1])
	}else{
		target_files <- target_lines
		target_names <- target_lines
	}
	
	targets_list <- list()
	targets_perbp <- list()
	size_list <- list()
	
	counter <- 1
	library(data.table)

    quick_sum <- function(target_bed) {
        #awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $1
        return(as.numeric(system(paste0("./quick_sum.sh ", target_bed), intern = T)))
    }
    quick_sum(target_files[1])

    intersect_count <- function(target_bed, region_bed) {
        return(system(paste0("bedtools intersect -u -a ", target_bed, " -b ", region_bed, " | wc | awk '{print $1}' "), intern = T))
    }

##
do_intersect <- function(target_bed, region_bed, actual_bp) {
	
	##Let's not generate any files here...
	
	bed_command <- paste0("bedtools intersect -a ", target_bed, " -b ", region_bed, " -c ")
	
	if(actual_bp) {
		bed_command <- paste0("bedtools intersect -a ", target_bed, " -b ", region_bed, " -wo ")
	}
	
	library(data.table)
    bed_output <- as.data.frame(fread(bed_command))
	
	intersect_count <- sum(as.numeric(as.character(bed_output[, dim(bed_output)[2]])))
	
	curr_size <- quick_sum(target_bed)
	
	system("echo '1' >> par_counter.txt")
	
	out_list <- list()
	out_list[[1]] <- intersect_count
	out_list[[2]] <- intersect_count/curr_size
	out_list[[3]] <- curr_size
	return(out_list)
}

library(parallel)
cl <- makeCluster(20)

clusterExport(cl, ls(), environment())
clusterExport(cl, "quick_sum", environment())

system("rm par_counter.txt")
fin_data <- parLapply(cl, target_files, function(x) do_intersect(x, region_bed, actual_bp))

size_list <- lapply(fin_data, function(x) x[[3]])
targets_list <- lapply(fin_data, function(x) x[[1]])
targets_perbp <- lapply(fin_data, function(x) x[[2]])

out_frame <- data.frame(name = target_names, target = target_files, bp_size = unlist(size_list), intersect_count = unlist(targets_list), per_bp = unlist(targets_perbp))

saveRDS(out_frame, gsub(".txt", "_intersected.rds", target_bed_set))
if(!actual_bp) {

write.csv(out_frame, gsub(".txt", paste0("_", gsub(".bed", "", region_bed), "_intersected.csv"), target_bed_set), row.names=F)
}else{
write.csv(out_frame, gsub(".txt", paste0("_", gsub(".bed", "", region_bed), "_intersected_BP_overlaps.csv"), target_bed_set), row.names=F)
}

}

##actual_bp means # of HAR bases intersected...

source("intersect.R")
region_bed <- "HARS_hg38_clean.bed"
target_bed_set <- "E53_specific_back.txt"
calculate_target_intersections_cleaned_FEB2025(target_bed_set, region_bed, actual_bp = F)
target_bed_set <- "E57_specific_back.txt"
calculate_target_intersections_cleaned_FEB2025(target_bed_set, region_bed, actual_bp = F)
target_bed_set <- "E67_specific_back.txt"
calculate_target_intersections_cleaned_FEB2025(target_bed_set, region_bed, actual_bp = F)
target_bed_set <- "E72_specific_back.txt"
calculate_target_intersections_cleaned_FEB2025(target_bed_set, region_bed, actual_bp = F)

##now, the target sets as well
region_bed <- "HARS_hg38_clean.bed"
calculate_target_intersections_cleaned_FEB2025("E53_specific_peaks.txt", region_bed, actual_bp = F)
calculate_target_intersections_cleaned_FEB2025("E57_specific_peaks.txt", region_bed, actual_bp = F)
calculate_target_intersections_cleaned_FEB2025("E67_specific_peaks.txt", region_bed, actual_bp = F)
calculate_target_intersections_cleaned_FEB2025("E72_specific_peaks.txt", region_bed, actual_bp = F)

##now, fit a beta-binomial distribution for each timepoint-specific background set, to use for significance testing.


find_HAR_distribution("E53_specific_back_HARS_hg38_clean_intersected.csv")
find_HAR_distribution("E57_specific_back_HARS_hg38_clean_intersected.csv")
find_HAR_distribution("E67_specific_back_HARS_hg38_clean_intersected.csv")
find_HAR_distribution("E72_specific_back_HARS_hg38_clean_intersected.csv")

find_HAR_distribution <- function(back_file) {
	
	curr_back <- read.csv(back_file)

	curr_rand_ints <- curr_back$per_bp
a <- curr_rand_ints[which(curr_rand_ints != 0)]
library(fitdistrplus)
fit_b <- fitdist(a, "beta")
beta_test <- bootdist(fit_b, niter=1e3)
summary(beta_test)

out_frame <- data.frame(file = back_file, shape1_median = median(beta_test$estim$shape1), shape2_median = median(beta_test$estim$shape2))

write.csv(out_frame, gsub(".csv", "_beta_params.csv", back_file), row.names=F)
}


test_sig_wrapper <- function(csv, back_file, outname, param_file, remove_zero = F) {
	
	int_data <- read.csv(csv)
	colours <- int_data$colour
	
	target_labels <- as.character(int_data$name)
	bp_values <- as.numeric(int_data$per_bp)
	
	if(!identical(grep(".csv", back_file), integer(0))) {
		print("csv")
		test_sig(bp_values, back_file, param_file, justlist=F, target_labels, outname, remove_zero = remove_zero, colours = colours)

	}else{
		test_sig(bp_values, back_file, param_file, justlist=T, target_labels, outname, remove_zero = remove_zero)
	}
	
}

test_sig <- function(targ_count, back_file, param_file, justlist=F, targetlabels, outname, plot_ecdf = T, remove_zero = F, colours) {
	
	give_standards <- function(vect, single=NULL) {
		if(is.null(single)) {
		z_score_set <- sapply(vect, function(x) ((mean(vect) - x) / sd(vect)))
		return(z_score_set)
	}else{
		z_score_single <- (mean(vect) - single) / (sd(vect))
		return(z_score_single)
	}
	
	}
	
	target_int <- targ_count
	target_int_list <- target_int ##Unnecessary when doing single value.
	
	if(!justlist) {
	
	curr_back <- read.csv(back_file)

	curr_rand_ints <- curr_back$per_bp

}else{
	curr_rand_ints <- as.numeric(as.character(readLines(back_file)))
}

	testing_list <- list()
	if(remove_zero) {
		print("For the purposes of gamma and lnnorm, remove zeros")
		curr_rand_ints <- curr_rand_ints[curr_rand_ints != 0]
	}
	
	all_standard <- give_standards(unlist(curr_rand_ints))
	ecdf_back <- ecdf(all_standard)
	ecdf_other <- ecdf(unlist(curr_rand_ints))

##read in pre-calculated beta distribution parameters	
param_data <- read.csv(param_file)
beta_1 <- param_data$shape1_median
beta_2 <- param_data$shape2_median

	if(length(target_int) > 1) {
		print("got a list then")
		
		target_int_list <- target_int
		
		pval_set <- list()
		ecdf_pval_set <- list()
		gamma_val_set <- list()
		plnorm_set <- list()
		beta_set <- list()
		for(x in 1:length(target_int_list)) {
			target_int <- target_int_list[x]
			
			targ_standard <- give_standards(unlist(curr_rand_ints), unlist(target_int))
			stand_p <- pnorm(-abs(targ_standard))
			target_label <- targetlabels[x]
			pval_set[[toString(target_label)]] <- stand_p
			ecdf_pval_set[[toString(target_label)]] <- 1 - ecdf_other(target_int)
			curr_beta <- pbeta(target_int, shape1 = beta_1, shape2 = beta_2, lower.tail=F) 
        	beta_set[[toString(target_label)]] <- curr_beta
			
		}
		
		adj_pval <- p.adjust(unlist(pval_set), method="BH")
		
		adj_ecdf_pval <- p.adjust(unlist(ecdf_pval_set), method="BH")
		
		beta_pval <- unlist(beta_set)
		beta_adj_pval <- p.adjust(beta_pval, method="BH")
		
		target_int <- target_int_list

	}else{
	
	all_standard <- give_standards(unlist(curr_rand_ints))
	targ_standard <- give_standards(unlist(curr_rand_ints), unlist(target_int))
	stand_p <- pnorm(-abs(targ_standard))
}
	
Corner_text <- function(text, location="topright"){
legend(location,legend=text, bty ="n", pch=NA) 
}

	pdf(paste0(outname, "_HAR_enrichments_10k.pdf"), width=12, height=12)
	
	xrange <- range(c(unlist(curr_rand_ints), unlist(target_int_list)))
	
	if(max(xrange) < 1e-3) {
		xrange <- c(xrange[1] - xrange[1] %% 1e-6, xrange[2] + xrange[2] %% 1e-6)
	}else{
	
	xrange <- c(xrange[1] - xrange[1] %% 1, xrange[2] + xrange[2] %% 1)
	
}
	curr_hist <- hist(sort(unlist(curr_rand_ints)), type="l", main="", xlab="HAR per bp of sequence", xlim=xrange, xaxt="n", col="grey")
	
	abline(v= unlist(target_int), col="red")
	
	if(max(xrange) < 1e-3) {
		axis(1, at=seq(xrange[1], xrange[2], by=1e-6), labels= seq(xrange[1], xrange[2], by=1e-6))#cut(seq(xrange[1], xrange[2]), 10))

	}else{
	axis(1, at=seq(xrange[1], xrange[2], by=1), labels= seq(xrange[1], xrange[2], by=1))#cut(seq(xrange[1], xrange[2]), 10))
}
	
	if(length(target_int) > 1) {
		colour_set <- rainbow(length(target_int))[sample(length(target_int))]
		colour_set <- colours
		if(!plot_ecdf) {
		super_text <- sapply(1:length(target_int_list), function(x) paste0(targetlabels[x], " : ", target_int_list[x], " : adj pval : ", adj_pval[x]))
	}else{
		super_text <- sapply(1:length(target_int_list), function(x) paste0(targetlabels[x], " : ", target_int_list[x], " : ECDF adj pval : ", round(adj_ecdf_pval[x], 3)))
}
		sapply(1:length(target_int_list), function(x) abline(v=target_int_list[x], col=colour_set[x]))
	}else{
	Corner_text(paste0("pval=", round(stand_p, 5), "\n", "target = ", target_int, "\n", "n = ", length(curr_rand_ints)))
	
}

	#text(target_int_list, rep(max(curr_hist$counts) - 10, length(targetlabels)), targetlabels, srt = 90, adj = 1)

	dev.off()
		
	p_val_set <- data.frame(targetlabels, target_vals = targ_count, pval= unlist(pval_set), ecdf_pval = unlist(ecdf_pval_set), padj=adj_pval, ecdf_padj = adj_ecdf_pval, 
	beta_pval, beta_adj_pval)#, exp_pval, exp_adj_pval)
	
    p_val_set$back <- mean(curr_rand_ints)
    p_val_set$FOLD <- log2(p_val_set$target_vals / p_val_set$back)

	write.csv(p_val_set, paste0(outname, "_HAR_ENRICH_PVALS_DISTRIBUTIONS.csv"))
	
}

test_sig_wrapper("E53_specific_peaks_HARS_hg38_clean_intersected.csv", "E53_specific_back_HARS_hg38_clean_intersected.csv", param_file = "E53_specific_back_HARS_hg38_clean_intersected_beta_params.csv", "HAR_intersect_E53_CELLTYPESPECIFIC")

test_sig_wrapper("E57_specific_peaks_HARS_hg38_clean_intersected.csv", "E57_specific_back_HARS_hg38_clean_intersected.csv", param_file = "E57_specific_back_HARS_hg38_clean_intersected_beta_params.csv", "HAR_intersect_E57_CELLTYPESPECIFIC")

test_sig_wrapper("E67_specific_peaks_HARS_hg38_clean_intersected.csv", "E67_specific_back_HARS_hg38_clean_intersected.csv", param_file = "E67_specific_back_HARS_hg38_clean_intersected_beta_params.csv", "HAR_intersect_E67_CELLTYPESPECIFIC")

test_sig_wrapper("E72_specific_peaks_HARS_hg38_clean_intersected.csv", "E72_specific_back_HARS_hg38_clean_intersected.csv", param_file = "E72_specific_back_HARS_hg38_clean_intersected_beta_params.csv", "HAR_intersect_E72_CELLTYPESPECIFIC")


##combine these all together to perform a broader p-adjustment.
small_csv <- lapply(system('ls *DISTRIBUTIONS*.csv', intern = T), read.csv)
small_collapse <- do.call(rbind, small_csv)
small_collapse$beta_adj2 <- p.adjust(small_collapse$beta_pval, method="BH")

small_collapse$TIME <- unlist(lapply(small_collapse$targetlabels, function(x) unlist(strsplit(x, "_"))[1]))
small_collapse$CELL <- basename(small_collapse$targetlabels)
small_collapse$CELL <- gsub("_peaks_specific.bed", "", small_collapse$CELL)
write.csv(small_collapse, "all_peakset_significance_results_BADJ.csv", row.names = F)

