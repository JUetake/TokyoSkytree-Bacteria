library(dada2)

#Original fastq files are available from https://figshare.com/account/home#/projects/60497
path <- "PATH/TO/FOLDER of Air_ref samples"
fns <- list.files(path)
fns

fastqs <- fns[grepl(".fastq$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
fnRs <- fastqs[grepl("_R2", fastqs)] # Just the reverse read files
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

plotQualityProfile(fnFs[[1]])
plotQualityProfile(fnFs[[2]])
plotQualityProfile(fnRs[[1]])
plotQualityProfile(fnRs[[2]])

# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Filter
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]),
                    truncLen=c(250,210), 
                    maxN=0, maxEE=c(1,3), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, verbose=TRUE)
}

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs.lrn <- dada(derepFs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errF <- dadaFs.lrn[[1]]$err_out
dadaRs.lrn <- dada(derepRs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errR <- dadaRs.lrn[[1]]$err_out
plotErrors(dadaFs.lrn[[1]], nominalQ=TRUE)

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])
seqtab <- makeSequenceTable(mergers[names(mergers)])
dim(seqtab)
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)
saveRDS(seqtab.nochim, "Skytree_air_potential_sources.rds") 

#Merge previously sequenced control sequeces with Air sample for quality control
control <- readRDS("Previous_sequences.rds")
Air_ref <- readRDS("Skytree_air_potential_sources.rds")
Skytree <- mergeSequenceTables(control, Air_ref)

#Customized silva 128 sequence file is also available from https://figshare.com/account/home#/projects/60497
taxa <- assignTaxonomy(Skytree, "silva_nr_v128_train_set_custom.fa")

Skytree <- t(Skytree)
Skytree <- cbind('#OTUID' = rownames(Skytree), Skytree)#Add '#OTUID' to the header (required for convert to biom)
write.table(Skytree, "ASV_table.txt" , sep=",")
write.table(taxa, "ASV_tax.txt" , sep=",")
