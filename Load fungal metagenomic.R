library ("vegan")

##Importing taxa files#
pool1taxa <- read.table("pool1tax.txt", sep = "\t", header = TRUE, row.names = 1)
pool2taxa <- read.table("pool2tax.txt", sep = "\t", header = TRUE, row.names = 1)
pool3taxa <- read.table("pool3tax.txt", sep = "\t", header = TRUE, row.names = 1)
pool4taxa <- read.table("pool4tax.txt", sep = "\t", header = TRUE, row.names = 1)
itstaxa <- read.table("magic_its2_q20_tax.txt", sep = "\t", header = TRUE, row.names = 1)

#Taxa Counts#
sample_counts <- c(colSums(pool1taxa),colSums(pool2taxa),colSums(pool3taxa),colSums(pool4taxa))
its_counts <- colSums(itstaxa)
its_counts <- its_counts[-grep("BLANK", names(its_counts))]
its_counts <- its_counts[-grep("Mock", names(its_counts))]
its_low_counts <- its_counts[its_counts < 100]
sample_counts_low_its <- sample_counts[names(its_low_counts)]

summary(sample_counts_low_its)

#Default taxonomy split using lowest level of taxonomy down to species#
splittaxa.its <- strsplit(rownames(pool1taxa), ";")
for(i in 1:length(splittaxa.its)){
  if(length(splittaxa.its[[i]]) > 6){
    rownames(pool1taxa)[i] <- splittaxa.its[[i]][7] #Species-level taxa
  } else {
    rownames(pool1taxa)[i] <- paste(splittaxa.its[[i]][(length(splittaxa.its[[i]]))],collapse = "_") #Deepest taxa if not species
  }
}
pool1taxa <- aggregate(pool1taxa, by=list(rownames(pool1taxa)),sum)
rownames(pool1taxa) <- pool1taxa$Group.1
pool1taxa <- pool1taxa[,2:length(colnames(pool1taxa))]

splittaxa.its <- strsplit(rownames(pool2taxa), ";")
for(i in 1:length(splittaxa.its)){
  if(length(splittaxa.its[[i]]) > 6){
    rownames(pool2taxa)[i] <- splittaxa.its[[i]][7] #Species-level taxa
  } else {
    rownames(pool2taxa)[i] <- paste(splittaxa.its[[i]][(length(splittaxa.its[[i]]))],collapse = "_") #Deepest taxa if not species
  }
}
pool2taxa <- aggregate(pool2taxa, by=list(rownames(pool2taxa)),sum)
rownames(pool2taxa) <- pool2taxa$Group.1
pool2taxa <- pool2taxa[,2:length(colnames(pool2taxa))]

splittaxa.its <- strsplit(rownames(pool3taxa), ";")
for(i in 1:length(splittaxa.its)){
  if(length(splittaxa.its[[i]]) > 6){
    rownames(pool3taxa)[i] <- splittaxa.its[[i]][7] #Species-level taxa
  } else {
    rownames(pool3taxa)[i] <- paste(splittaxa.its[[i]][(length(splittaxa.its[[i]]))],collapse = "_") #Deepest taxa if not species
  }
}
pool3taxa <- aggregate(pool3taxa, by=list(rownames(pool3taxa)),sum)
rownames(pool3taxa) <- pool3taxa$Group.1
pool3taxa <- pool3taxa[,2:length(colnames(pool3taxa))]

splittaxa.its <- strsplit(rownames(pool4taxa), ";")
for(i in 1:length(splittaxa.its)){
  if(length(splittaxa.its[[i]]) > 6){
    rownames(pool4taxa)[i] <- splittaxa.its[[i]][7] #Species-level taxa
  } else {
    rownames(pool4taxa)[i] <- paste(splittaxa.its[[i]][(length(splittaxa.its[[i]]))],collapse = "_") #Deepest taxa if not species
  }
}
pool4taxa <- aggregate(pool4taxa, by=list(rownames(pool4taxa)),sum)
rownames(pool4taxa) <- pool4taxa$Group.1
pool4taxa <- pool4taxa[,2:length(colnames(pool4taxa))]

#Combine pools into single taxa table#
pool2v1names <- setdiff(rownames(pool2taxa),rownames(pool1taxa))
pool1taxa[pool2v1names,] <- NA
pool3v1names <- setdiff(rownames(pool3taxa),rownames(pool1taxa))
pool1taxa[pool3v1names,] <- NA
pool4v1names <- setdiff(rownames(pool4taxa),rownames(pool1taxa))
pool1taxa[pool4v1names,] <- NA
pool1taxa[is.na(pool1taxa)] <- 0
pool1taxa <- pool1taxa[order(rownames(pool1taxa)),]

pool1v2names <- setdiff(rownames(pool1taxa),rownames(pool2taxa))
pool2taxa[pool1v2names,] <- NA
pool2taxa[is.na(pool2taxa)] <- 0
pool2taxa <- pool2taxa[order(rownames(pool2taxa)),]

pool1v3names <- setdiff(rownames(pool1taxa),rownames(pool3taxa))
pool3taxa[pool1v3names,] <- NA
pool3taxa[is.na(pool3taxa)] <- 0
pool3taxa <- pool3taxa[order(rownames(pool3taxa)),]

pool1v4names <- setdiff(rownames(pool1taxa),rownames(pool4taxa))
pool4taxa[pool1v4names,] <- NA
pool4taxa[is.na(pool4taxa)] <- 0
pool4taxa <- pool4taxa[order(rownames(pool4taxa)),]


taxatable <- cbind(pool1taxa,pool2taxa)
taxatable <- cbind(taxatable,pool3taxa)
taxatable <- cbind(taxatable,pool4taxa)

remove.row <- c("s__Moesziomyces antarcticus")
taxatable <- taxatable[!rownames(taxatable) %in% remove.row,]
taxatable <- taxatable[,colSums(taxatable) > 0]

sort(rowSums(taxatable)) #Most abundant taxa

fungprev <- list(1:nrow(taxatable))
for(i in 1:nrow(taxatable)){
  fungprev[i] <- sum(taxatable[i,] > 0)
}
names(fungprev) <- paste0(rownames(taxatable)) #fungprev now a list of prevalence of each taxon
fungprev <- sort(unlist(fungprev)) #List printed out in more legible format and sorted
print (fungprev)

m.its <- read.table("metadata_full_Tim.txt", sep = "\t", header = TRUE, row.names = 1)

m.its$Patient_No <- as.character(m.its$Patient_No)
m.its <- m.its[colnames(itstaxa),]
m.its$minority <- ifelse(m.its$babyrace == c("White"),m.its$babyrace,c("Minority"))
m.its$currentfeed_bf <- ifelse(m.its$currentfeed_bf == 1, yes="Yes", ifelse(m.its$currentfeed_bf == 0, yes="No", no=NA))
m.its[is.na(m.its$currentfeed_bf),c("currentfeed_bf")] <- c("UNKNOWN")
m.its$deliverymethod <- ifelse(m.its$deliverymethod == c("C"), yes="C-Section", ifelse(m.its$deliverymethod == c("V"), yes="Vaginal", no=NA))
m.its$abxmat <- as.character(m.its$abxmat)


m.its$Timecluster <- NA
m.its <- m.its[grepl("^NA", rownames(m.its))==F,] 
for(i in 1:length(m.its$Timeline_Weeks)){
  maprow <- m.its[i,c("Timeline_Weeks")]
  if(maprow < 5){
    m.its[i,c("Timecluster")] <- 1
  }
  if(maprow > 4 & maprow < 14){
    m.its[i,c("Timecluster")] <- 2
  }
  if(maprow > 13 & maprow < 25){
    m.its[i,c("Timecluster")] <- 3
  }
  if(maprow > 24 & maprow < 38){
    m.its[i,c("Timecluster")] <- 4
  }
  if(maprow > 37 & maprow < 51){
    m.its[i,c("Timecluster")] <- 5
  }
  if(maprow > 50 & maprow < 64){
    m.its[i,c("Timecluster")] <- 6
  }
  if(maprow > 63){
    m.its[i,c("Timecluster")] <- 7
  }
}


m.its$Timecluster <- as.factor(m.its$Timecluster)
m.its$matbmicat <- NA
for(i in 1:length(m.its$maternal_bmi_prepreg)){
  maprow <- m.its[i,c("maternal_bmi_prepreg")]
  if(is.na(maprow) == TRUE){
    m.its[i,c("matbmicat")] <- NA
  }
  else{
    if(maprow < 25){
      m.its[i,c("matbmicat")] <- c("Normal")
    }
    if(maprow > 25 & maprow < 30){
      m.its[i,c("matbmicat")] <- c("Overweight")
    }
    if(maprow > 30){
      m.its[i,c("matbmicat")] <- c("Obese")
    }
  }
}
#Add antifungal metadata columns#
antifungcols <- grep("antifung",colnames(m.its))
m.its$Antifung_any <- apply(m.its[,antifungcols],
                            1,
                            function(x) any(x == "Yes"))
m.its$Antifung_any_yn <- ifelse(m.its$Antifung_any == "TRUE", "Y", "N")
m.its$Antifung_any_yn[is.na(m.its$Antifung_any_yn)] <- "N"
antifungsamples <- m.its[m.its$Antifung_any_yn == "Y",antifungcols]
temp <- grep("_____",colnames(antifungsamples))
antifungsamples <- antifungsamples[,-temp]
antifungsamples$Timeline_Weeks <- m.its[rownames(antifungsamples),c("Timeline_Weeks")]
antifungsamples$Timeline_Months <- trunc(antifungsamples$Timeline_Weeks/4)
antifungmonths <- apply(antifungsamples[,c(1,3,5,7,9,11,13,15,17)],
                        1,
                        function(x) which(x == "Yes"))
antifungmonths <- rapply(antifungmonths, function(x) ifelse(x==9,24,x), how="replace")
antifungmonths <- rapply(antifungmonths, function(x) ifelse(x==8,21,x), how="replace")
antifungmonths <- rapply(antifungmonths, function(x) ifelse(x==7,18,x), how="replace")
antifungmonths <- rapply(antifungmonths, function(x) ifelse(x==6,15,x), how="replace")
antifungmonths <- rapply(antifungmonths, function(x) ifelse(x==5,12,x), how="replace")
antifungmonths <- rapply(antifungmonths, function(x) ifelse(x==4,9,x), how="replace")
antifungmonths <- rapply(antifungmonths, function(x) ifelse(x==3,6,x), how="replace")
antifungmonths <- rapply(antifungmonths, function(x) ifelse(x==2,3,x), how="replace")

for(i in 1:length(antifungmonths)){
  temprow <- antifungmonths[i]
  tempname <- names(temprow)
  temprow2 <- if(length(unlist(temprow)) > 1){temprow[[1]][2]}
  temprow1 <- temprow[[1]][1]
  antifungsamples[tempname,c("antifungmonth")] <- temprow1
  tempmatch <- antifungsamples[tempname,c("Timeline_Months","antifungmonth")]
  tempdiff <- tempmatch[,1] - tempmatch[,2]
  temppos <- if (tempdiff < 0) {FALSE} else if (tempdiff > 3) {FALSE} else if(tempdiff <= 3){TRUE}
  temppos2 <- FALSE
  if(length(unlist(temprow)) > 1){
    antifungsamples[tempname,c("antifungmonth2")] <- temprow2
    tempmatch2 <- antifungsamples[tempname,c("Timeline_Months","antifungmonth2")]
    tempdiff2 <- tempmatch[,1] - tempmatch[,2]
    temppos2 <- if (tempdiff2 < 0) {FALSE} else if (tempdiff2 > 3) {FALSE} else if(tempdiff2 < 3){TRUE}
  }
  if(temppos == "TRUE"){
    antifungsamples[tempname,c("antifungpos")] <- temppos
  }
  else {
    antifungsamples[tempname,c("antifungpos")] <- temppos2
  }
}

m.its[rownames(antifungsamples),c("Antifung_any_yn")] <- antifungsamples$antifungpos
m.its$Antifung_any_yn <- replace(m.its$Antifung_any_yn,m.its$Antifung_any_yn=="TRUE","Y")
m.its$Antifung_any_yn <- replace(m.its$Antifung_any_yn,m.its$Antifung_any_yn=="FALSE","N")

m.its[is.na(m.its$currentfeed_bf),c("currentfeed_bf")] <- c("UNKNOWN")


#adding abx data from new file#
m.add <- read.csv("CHOP_metadata_abx_rounds.txt", sep = "\t", header = TRUE, fill = TRUE, row.names = 1)
m.add <- m.add[rownames(m.its),]
m.its$Abx_any_yn <- m.add$Abx_any_yn
m.its$abx_rounds <- as.character(m.add$abx_rounds)

#Normalization and QC steps#
taxon.sub <- sweep(taxatable,1,rowSums(taxatable),'/')*100 #Basic relative abundance normalization, all samples present
#QC
summary(colSums(taxatable))[2] #Use 1st quartile as cutoff in following line
x <- summary(colSums(taxatable))[2] #x will be the low-end cutoff for minimum number of sequences in a sample
taxon.qc <- taxatable[,colSums(taxatable) > x] #eliminates samples with fewer than x reads
elimtaxon.qc <- taxatable[,colSums(taxatable) < x]
#Eliminate OTUs present in <y samples
y <- 4 #Start with presence in 5 or more samples being the cutoff
fungprev.qc <- list(1:nrow(taxon.qc))
for(i in 1:nrow(taxon.qc)){
  fungprev.qc[i] <- sum(taxon.qc[i,] > 0)
} #Loops through taxon.qc and counts prevalence of each taxon and adds it to the list
names(fungprev.qc) <- paste0(rownames(taxon.qc)) #fungprev.qc now a list of prevalence of each taxon
taxon.qc <- taxon.qc[which(fungprev.qc > y),] #taxon.qc now has only samples with >x sequences and taxa present in >y samples
sort(unlist(fungprev.qc)) #List printed out in more legible format and sorted
taxon.qc.seqs <- taxon.qc

m.qc <- m.its[colnames(taxon.qc),]
taxa.qc <- t(taxon.qc)

taxon.qc.sub <- as.data.frame(sweep(t(taxon.qc.seqs),1,rowSums(t(taxon.qc.seqs)),'/')*100) #Basic relative abundance normalization, all samples present


#CLR transformation
taxon.qc.clr <- taxon.qc
eps <- 0.5
taxon.qc.clr <- taxon.qc.clr*(1-rowSums(taxon.qc.clr==0)*eps/rowSums(taxon.qc.clr))
taxon.qc.clr[taxon.qc.clr==0] <- eps
taxon.qc.clr <- sweep(taxon.qc.clr,1,rowSums(taxon.qc.clr),'/')
ltaxon <- log(taxon.qc.clr)
taxon.qc.clr <- t(ltaxon - rowMeans(ltaxon))
taxon.qc.clr <- taxon.qc.clr[,!is.nan(colSums(taxon.qc.clr))]
m.qc.clr <- m.qc[rownames(taxon.qc.clr),]


#Median depth normalization
taxon.med <- taxon.qc
x <- median(colSums(taxon.med))
taxon.pre <- round(sweep(taxon.med,2,colSums(taxon.med),'/')*max(colSums(taxon.med)))
taxon.med <- rrarefy(t(taxon.pre),x)
m.med <- m.qc[rownames(taxon.med),]
