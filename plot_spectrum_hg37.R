library("GenomicRanges")
library("Rsamtools")
library("MASS")

plot_spectrum_37 = function(bed,save,add_to_title=""){
  genomeFile="/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/human/GRCh37d5/genome.fa"
  mutations=as.data.frame(bed)

  mutations$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(mutations$Chrom, IRanges(as.numeric(mutations$Pos)-1, 
                                                                                       as.numeric(mutations$Pos)+1))))
  
  # 2. Annotating the mutation from the pyrimidine base
  ntcomp = c(T="A",G="C",C="G",A="T")
  mutations$sub = paste(mutations$Ref,mutations$Alt,sep=">")
  mutations$trinuc_ref_py = mutations$trinuc_ref
  for (j in 1:nrow(mutations)) {
    if (mutations$Ref[j] %in% c("A","G")) { # Purine base
      mutations$sub[j] = paste(ntcomp[mutations$Ref[j]],ntcomp[mutations$Alt[j]],sep=">")
      mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
    }
  }

  # 3. Counting subs
  freqs = table(paste(mutations$sub,paste(substr(mutations$trinuc_ref_py,1,1),substr(mutations$trinuc_ref_py,3,3),sep="-"),sep=","))

  sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
  ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
  full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
  freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec

  xstr = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")

  if(is.null(save)){dev.new(width=12,height=4)}
  #if(is.null(save)){dev.new(width=16,height=6)} #for tiff
  colvec = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
  y = freqs_full; maxy = max(y)
  h = barplot(y, las=2, col=colvec, border=NA, ylim=c(0,maxy*1.5), space=1, cex.lab=1.5, cex.axis = 1.3, cex.names=0.9, names.arg=xstr, ylab="Number of mutations")
  #h = barplot(y, las=2, col=colvec, border=NA, ylim=c(0,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab="Number of mutations", main=paste0("Total number of mutations: ",nrow(mutations)))
  for (j in 1:length(sub_vec)) {
    xpos = h[c((j-1)*16+1,j*16)]
    rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colvec[j*16])
    text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
  } 
  if(!is.null(save)){ 
    #dev.copy(tiff,save,width=16,height=6, units="in", compression="lzw", res=300); dev.off()} #for tiff image
    dev.copy(pdf,save,width=12,height=4); dev.off()}
}
