library(vegan)
library(GUniFrac)
library(plyr)
library(stringi)
library(ggplot2)


#ITS2 Loading and calculations#
itstaxa <- (read.table("magic_its2_q20_tax.txt", sep = "\t", fill = TRUE, header = TRUE, row.names = 1))
m.its <- read.table("metadata_full_Tim.txt", sep = "\t", header = TRUE, row.names = 1)
otu.its.all <- t(read.table("magic_its2_q20_tax.txt", sep = "\t", fill = TRUE, header = TRUE, row.names = 1))


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
m.its$abx_rounds <- m.add$abx_rounds
m.its$abx_rounds <- as.character(m.add$abx_rounds)


##new amp qc
rownames(itstaxa) <- gsub("africana","albicans",rownames(itstaxa))
splittaxa.its <- strsplit(rownames(itstaxa), ";")
itstaxa <- t(itstaxa)
for(i in 1:length(splittaxa.its)){
  if(length(splittaxa.its[[i]]) > 6){
    colnames(itstaxa)[i] <- splittaxa.its[[i]][7] #Species-level taxa
  } else {
    colnames(itstaxa)[i] <- paste(splittaxa.its[[i]][(length(splittaxa.its[[i]]))],collapse = "_") #Deepest taxa if not species
  }
}
colnames(itstaxa) <- gsub("kudriavzevii","cerevisiae",colnames(itstaxa))
itstaxa <- aggregate(t(itstaxa), by=list(rownames(t(itstaxa))),sum)
rownames(itstaxa) <- itstaxa$Group.1
itstaxa <- itstaxa[,2:length(colnames(itstaxa))]

itsblanks <- itstaxa[,grep("BLANK", colnames(itstaxa))]
itsblanks <- itsblanks[rowSums(itsblanks) > 0,]
itssynthmocks <- itstaxa[,grep("Mock", colnames(itstaxa))]
itstaxa <- itstaxa[,-grep("BLANK", colnames(itstaxa))]
itstaxa <- itstaxa[,-grep("Mock", colnames(itstaxa))]
itstaxa <- itstaxa[-grep("Synthetic", rownames(itstaxa)),]



summary(colSums(itstaxa))[2] #Use 1st quartile as cutoff in following line
x <- summary(colSums(itstaxa))[2] #x will be the low-end cutoff for minimum number of sequences in a sample
itstaxon.qc <- itstaxa[,colSums(itstaxa) > x] #eliminates samples with fewer than x reads
elimitstaxon.qc <- itstaxa[,colSums(itstaxa) < x]
#Eliminate OTUs present in <y samples
y <- 4 #Start with presence in 5 or more samples being the cutoff
fungprev.its.qc <- list(1:nrow(itstaxon.qc))
for(i in 1:nrow(itstaxon.qc)){
  fungprev.its.qc[i] <- sum(itstaxon.qc[i,] > 0)
} #Loops through itstaxon.qc and counts prevalence of each taxon and adds it to the list
names(fungprev.its.qc) <- paste0(rownames(itstaxon.qc)) #fungprev.its.qc now a list of prevalence of each taxon

fungprev <- sort(unlist(fungprev.its.qc)) #List printed out in more legible format and sorted
print (fungprev.its.qc)


itstaxon.qc <- itstaxon.qc[which(fungprev.its.qc > y),] #itstaxon.qc now has only samples with >x sequences and taxa present in >y samples
itstaxon.qc.seqs <- itstaxon.qc
m.its.qc <- m.its[colnames(itstaxon.qc),]

sort(rowSums(itstaxon.qc)) #Most abundant taxa


m.its.filt.norm.species <- m.filt.its[(rownames(taxa.its.filt.norm.species)),]
 
t.taxa.its.filt.norm.species <- t(taxa.its.filt.norm.species)
taxa.its.filt.norm.species.sub <- as.data.frame(sweep(t(t.taxa.its.filt.norm.species),1,rowSums(t(t.taxa.its.filt.norm.species)),'/')*100) #Basic relative abundance normalization, all samples present
taxa.its.filt.norm.species.sub <- t(taxa.its.filt.norm.species.sub)

#CLR transformation#
itstaxon.clr <- itstaxon.qc
eps <- 0.5
itstaxon.clr <- itstaxon.clr*(1-rowSums(itstaxon.clr==0)*eps/rowSums(itstaxon.clr))
itstaxon.clr[itstaxon.clr==0] <- eps
itstaxon.clr <- sweep(itstaxon.clr,1,rowSums(itstaxon.clr),'/')
ltaxon <- log(itstaxon.clr)
itstaxon.clr <- t(ltaxon - rowMeans(ltaxon))
itstaxon.clr <- itstaxon.clr[,!is.nan(colSums(itstaxon.clr))]
m.its.clr <- m.its.qc[rownames(itstaxon.clr),]


###median depth normalization###
taxa.med <- (itstaxon.qc)
x <- median(colSums(taxa.med))
taxa.pre <- round(sweep(taxa.med,2,colSums(taxa.med),'/')*max(colSums(taxa.med)))
taxa.med <- rrarefy(t(taxa.pre),x)
# ##added this next line because table was transposed?##
taxa.its.med <- (taxa.med)
m.its.med <- m.its.qc[rownames(taxa.its.med),]

