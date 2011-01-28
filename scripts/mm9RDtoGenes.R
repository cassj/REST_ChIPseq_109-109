#!/usr/local/bin/Rscript

options(stringsAsFactors = FALSE);

library(IRanges)

qw <- function(...) {
  as.character(sys.call()[-1])
}

args <- commandArgs(trailingOnly=TRUE)
filename = args[1]

rd <- get(load(filename))

library(ChIPpeakAnno)
library(biomaRt)
ensmart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")


#load mouse transcripts
data(TSS.mouse.NCBIM37)

#The chromosomes names used in teh NCBIM37 dataset are 1..Y, MT and
#various NT_***** non-standard chrs. 

                                        #remove "chr" prefix for  ChIPpeakAnno
names(rd)<-gsub("chr","",names(rd))

#change "M" to "MT" for ChIPpeakAnno
id<-which(names(rd)=="M")
if (length(id)>0){
   names(rd)[id]<-"MT"
}

#chippeakanno discards anything in values, so hang on to them
vals <- as.data.frame(rd)
vals <- vals[,grep("values.", colnames(vals), value=T)]
colnames(vals) <- gsub('values.','',colnames(vals))


# NOTE: TSS.mouse.NCBIM37 is actually *gene* start and end positions,
# not individual transcripts.

#get the most recent annotation data from ensembl
tss <- getAnnotation(ensmart, "TSS")
save(tss, file="tss.RData")



#get the nearest start
nearest.start <- annotatePeakInBatch(rd,
                                   AnnotationData=tss,
                                   PeakLocForDistance = "middle",
                                   FeatureLocForDistance = "start",
                                   output = "nearestStart"
                                   )

#and nearest end.
nearest.end <- annotatePeakInBatch(rd,
                                   AnnotationData=tss,
                                   PeakLocForDistance = "middle",
                                   FeatureLocForDistance = "end",
                                   output = "nearestStart")


#and nearest miRNA
mirnas <- getAnnotation(ensmart, "miRNA")
save(mirnas, file="mirna.RData")
nearest.mirna <- annotatePeakInBatch(rd,
                                     AnnotationData=mirnas,
                                     PeakLocForDistance="middle",
                                     FeatureLocForDistance="middle",
                                     output ="nearestStart") 




                                        #get the nearest feature plus any features that overlap the peak that is not the nearest feature
nearest.feature <- annotatePeakInBatch(rd,
                                       AnnotationData=tss,
                                       PeakLocForDistance = "middle",
                                       FeatureLocForDistance = "TSS",
                                       output = "both",
                                       multiple = TRUE
                                       )


#save nearest.feature in case need it again, ffs its slow...
save(nearest.feature, file="nearestfeature.RData")
#nearest.feature <-get(load("nearestfeature.RData"))

#and later we'll want to know if it overlaps a coding region?
exons <- getAnnotation(ensmart, "Exon")
save(exons, file="exons.RData")

# annotatePeakInBatch is bastard slow, the IRanges findOverlaps function
# works on RangedData and is much faster
overlapping.exon <- findOverlaps(rd, exons, type="any", select="arbitrary", drop=TRUE)
inds <- is.na(overlapping.exon)
overlapping.exon[inds] <- 0
overlapping.exon[!inds] <- 1

#add exon overlap data to peak data
vals <- data.frame(vals, overlapping.exon)

#get the rd rownames
rd.df <- as.data.frame(rd)
row.nms <- as.character(rd.df$names)


#the "nearest" tables have a 1:1 relationship with the peaks:

nearest.start <- as.data.frame(nearest.start)
rownames(nearest.start) <- as.character(nearest.start$peak)
nearest.start <- nearest.start[row.nms,]
nearest.start <- cbind(nearest.start, type="START", vals)

nearest.end <- as.data.frame(nearest.end)
rownames(nearest.end) <- as.character(nearest.end$peak)
nearest.end <- nearest.end[row.nms,]
nearest.end <- cbind(nearest.end,type="END", vals)

nearest.mirna <- as.data.frame(nearest.mirna)
rownames(nearest.mirna) <- as.character(nearest.mirna$peak)
nearest.mirna <- nearest.mirna[row.nms,]
nearest.mirna <- cbind(nearest.mirna, type="miRNA", vals)


#the "nearest.feature" is nearest, plus any overlapping so it has a *:1 
#relationship with the peaks

vals <- data.frame(vals, peak.id=rd.df[,"names"])
nearest.feature.df <- as.data.frame(nearest.feature)
peak.id <- gsub(' ENSMUS.*$','',nearest.feature.df[,"names"])
nearest.feature.df <- data.frame(nearest.feature.df, peak.id)
nearest.or.overlapping<-merge(nearest.feature.df, vals, by.x="peak", by.y="peak.id")

# for each gene, keep the info only for the nearest peak
nms <- nearest.or.overlapping[,"feature"]
dis <- abs(nearest.or.overlapping[,"distancetoFeature"])
nearest.peak.to.gene<-nearest.or.overlapping[order(nms, dis, decreasing=TRUE),]
nearest.peak.to.gene<-nearest.peak.to.gene[!duplicated(nearest.peak.to.gene[,"feature"]),]



# Get extra info from the ensembl ids                
filters <- c("ensembl_gene_id")
values<-unique(c(nearest.start[,"feature"], nearest.end[,"feature"], nearest.mirna[,"feature"], nearest.or.overlapping[,"feature"] ))
attributes <- c("ensembl_gene_id","mgi_symbol", "description")

annot <- getBM(filters=filters, values=values, attributes=attributes, mart=ensmart)

#ditch any that don't have a symbol or a description
no.anno <- intersect(
                     which(annot[,"mgi_symbol"]==""),
                     which(annot[,"description"]==""))
annot <- annot[-1*no.anno,]


# a few have multiple bits of annotation
annot<-cbind(annot, alt.annot="")
dups<-annot[duplicated(annot[,"ensembl_gene_id"]), "ensembl_gene_id"]

#keep the first on and add all the others as alt.annot 
for (d in dups){
  inds <- which(annot[,"ensembl_gene_id"]==d)
  this.alt.annot <- annot[inds[-1], c("mgi_symbol", "description")]
  annot[inds[1],"alt.annot"] <- paste(paste(this.alt.annot[,1], this.alt.annot[,2]), collapse="; ")
}
annot <- annot[!duplicated( annot[,"ensembl_gene_id"] ), ]
rownames(annot) <- annot[,"ensembl_gene_id"]


#add the annotation to your nearest info

nearest.start <- cbind(nearest.start,     annot[nearest.start[,"feature"],   c("mgi_symbol", "description")])
nearest.end   <- cbind(nearest.end,     annot[nearest.end[,"feature"],   c("mgi_symbol", "description")])
nearest.mirna <- cbind(nearest.mirna,   annot[nearest.mirna[,"feature"], c("mgi_symbol", "description")])

nearest.or.overlapping <-   cbind(nearest.or.overlapping,     annot[nearest.or.overlapping[,"feature"],   c("mgi_symbol", "description")])
nearest.peak.to.gene <-   cbind(nearest.peak.to.gene,     annot[nearest.peak.to.gene[,"feature"],   c("mgi_symbol", "description")])


# make a table of nearest overall
nearest.nearest <- nearest.start
smaller <- which( abs(nearest.nearest$distancetoFeature) > abs(nearest.end$distancetoFeature) )
if (length(smaller)>0) nearest.nearest[smaller,] <- nearest.end[smaller,]
smaller <- which( abs(nearest.nearest$distancetoFeature) > abs(nearest.mirna$distancetoFeature) )
if (length(smaller)>0) nearest.nearest[smaller,] <- nearest.mirna[smaller,]

#reorder
#colnms <- qw(space, start, end, width, names, peak, strand, feature, start_position, end_position, insideFeature,distancetoFeature, #shortestDistance,fromOverlappingOrNearest,overlapping.exon, type, mgi_symbol, description , Length, Summit, nTags,neg10log10pVal, FoldEnrichment, FDR )
#nearest.tss <- nearest.tss[,colnms]
#nearest.tes <- nearest.tes[,colnms]
#nearest.mirna <- nearest.mirna[,colnms]
#nearest.nearest <- nearest.nearest[,colnms]


#save the results as csv and rd?
filepath = gsub('(.*\\/).*','\\1', filename)
if(filepath==filename) dir='.'
write.csv(nearest.start, file=paste(filepath, "nearest_start.csv",sep=""), row.names=F)
write.csv(nearest.end, file=paste(filepath,"nearest_end.csv", sep=""), row.names=F)
write.csv(nearest.mirna, file=paste(filepath, "nearest_mirna.csv",sep=""), row.names=F)
write.csv(nearest.nearest, file=paste(filepath, "nearest_nearest.csv", sep=""), row.names=F)
write.csv(nearest.or.overlapping, file=paste(filepath,"nearest_or_overlapping.csv", sep=""), row.names=F)
write.csv(nearest.peak.to.gene, file=paste(filepath, "nearest_peak_to_gene.csv",sep=""), row.names=F)








