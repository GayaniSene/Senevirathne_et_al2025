########
##HAR enrichment testing, comparing diffferent cell types at individual timepoints
##February 2025
########


dirset <- readLines("specific_dirs.txt")
startdir <- getwd()

genome <- "../hg38.chrom.clean.sizes"
target <- "../HARS_hg38_clean.bed"

library(parallel)

##read in the previously-generated cell-type-specific peak files

for (dir in dirset) {
    setwd(dir)
    all_spec_files <- system('ls *specific.bed', intern = T)
    pad_file <- system("ls *merge.bed", intern = T)
    sample_num <- 10000
    all_frames <- list()
    summary_frames <- list()
    cl <- makeCluster(22)

    already_done_flag <- FALSE
    if(file.exists("all_frames_10K.rds")) {
        all_frames <- readRDS("all_frames_10K.rds")
        already_done_flag <- TRUE
    }

    for (bed in all_spec_files) {

        target_count <- system(paste0("bedtools intersect -u -a ", target, " -b ", bed, " | wc | awk '{print $1}'"), intern = T)
        
        if(!already_done_flag) {

        clusterExport(cl, c("bed", "pad_file", "sample_num", "genome", "target"), environment())
        sample_counts <- unlist(parLapply(cl, 1:sample_num, 
            function(x) system(paste0("bedtools shuffle -incl ", pad_file, " -excl ", bed, " -g ", genome, " -i ", bed, " | bedtools intersect -u -a ", target, " -b - | wc | awk '{print $1}'"), intern = T)))
        curr_frame <- data.frame(bed = bed, as.numeric(target_count), as.numeric(sample_counts))
        all_frames[[bed]] <- curr_frame
        }else{
            sample_counts <- all_frames[[bed]][,3]
        }

        summary_frame <- data.frame(bed = bed, target = as.numeric(target_count), mean = mean(as.numeric(sample_counts)), sd = sd(as.numeric(sample_counts)))
        summary_frame$FDR <- length(which(sample_counts > as.numeric(target_count))) / sample_num
        summary_frames[[bed]] <- summary_frame
        print(bed)
    }

    if(!already_done_flag){
    saveRDS(all_frames, "all_frames_10K.rds")
    }

    summary_frame_collapse <- do.call("rbind", summary_frames)
    write.csv(summary_frame_collapse, paste0(unlist(strsplit(dir, "_"))[1], "_summarized_celltype_enrichments_10K_mk2.csv"), row.names = F)
    print(dir)
    setwd(startdir)
}

