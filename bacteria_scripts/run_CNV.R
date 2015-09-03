require(cn.mops)

window_length = 1000
base_dir = "/lustre/scratch108/bacteria/jl11/pairs_analysis/haem_flu_pairs"
pairs_file = "haem_pairs_Rin.txt"

pairs <- read.table(file=paste(base_dir, pairs_file,sep="/"),header=T,sep="\t")

BAMfiles <- paste("/lustre/scratch108/bacteria/jl11/cnvs/H_influenzae", 
      pairs$CSF_lane,paste(pairs$CSF_lane,"mapping.bam",sep="."),sep="/")

BAMfiles <- c(BAMfiles, paste("/lustre/scratch108/bacteria/jl11/cnvs/H_influenzae", 
                                  pairs$blood_lane,paste(pairs$blood_lane,"mapping.bam",sep="."),sep="/"))

bamDataRanges <- getReadCountsFromBAM(BAMfiles, 
                    sampleNames=c(paste(pairs$Sample_ID,"csf",sep="_"),paste(pairs$Sample_ID,"blood",sep="_")),mode="paired",WL=window_length)
cnv_result <- haplocn.mops(bamDataRanges)
int_cnv_result <- calcIntegerCopyNumbers(cnv_result)

saveRDS(cnv_result,file="cnv_result.RData")
saveRDS(int_cnv_result,file="int_cnv_result.RData")

counts <- matrix(ncol = nrow(pairs)*2, nrow = length(cnvr(int_cnv_result)))
counts <- as.data.frame(counts)
colnames(counts) <- c(paste(pairs$Sample_ID,"csf",sep="_"),paste(pairs$Sample_ID,"blood",sep="_"))

overlaps <- matrix(ncol = nrow(pairs), nrow = length(cnvr(int_cnv_result)))
overlaps <- as.data.frame(overlaps)
colnames(overlaps) <- pairs$Sample_ID

num_mistmatched = 0
mismatches = rep(0,length(cnvr(int_cnv_result)))
for (i in seq(1,nrow(pairs))) # reverted to using a loop as data structures are complex...
{
  overlaps[,i] <- factor(mcols(cnvr(int_cnv_result))[,eval(paste("X",pairs$Sample_ID,"_csf",sep="")[i])],levels = paste("CN",seq(0,8),sep="")) == 
    factor(mcols(cnvr(int_cnv_result))[,eval(paste("X",pairs$Sample_ID,"_blood",sep="")[i])],levels = paste("CN",seq(0,8),sep=""))
  mismatched_regions <- which(!(overlaps[,i]))
  if (length(mismatched_regions) != 0)
  {
    num_mistmatched <- num_mistmatched + 1
    mismatches[mismatched_regions] <- mismatches[mismatched_regions] + 1
  }
}

print(num_mistmatched)
print(mismatches)

ranges_format <- cbind(start(ranges(cnvr(int_cnv_result))), end(ranges(cnvr(int_cnv_result))), width(ranges(cnvr(int_cnv_result))))
colnames(ranges_format) <- c("start", "end", "width")

write.csv(overlaps, "cnv_mismatches.csv")
write.csv(mcols(cnvr(int_cnv_result)),file="cnv_counts.csv")
write.csv(ranges_format,"cnv_locations.csv")

#pdf(file="cnvs.pdf") # does nothing
plot(int_cnv_result,which=which(mismatches>0)) # irritatingly this calls dev.new()
#dev.off()
