rm(list=ls())
require(tcltk)
source("LOST_functions_V071113.R")


# Open Dialog Box to Select Source File
# The Path will be used as the destination directory
fileName<-tclvalue(tkgetOpenFile())

if (!nchar(fileName)){
	tkmessageBox(message="No file was selected!")
}else{
	tkmessageBox(message=paste("The file selected was",fileName))
}


# Read in the selected file to get treatment information
# Data structure requires that trt assignments are pulled from a header

Data<-read.csv(file=fileName, head=FALSE, nrows=1, stringsAsFactors=FALSE)
TRTS<-levels(as.factor(t(Data)[-1,]))
  
# Opens a list that prompts for the selection of the control group
Control<-getcontrol(TRTS)

A<-strsplit(fileName,"/")[[1]]
FileName<-A[length(A)]
B<-A[1:length(A)-1]
Directory<-paste(B,collapse="/")
  
########################################################################################################################
# Begin Analysis Package
detach("package:tcltk")

require(ggplot2)
require(doBy)
require(reshape)
require(reshape2)
require(plyr)
require (nlme)
require(gplots)
require(lsmeans)

# some functions

FileInpCSV <-FileName
RefTreat<-Control

cexout=0.75
PDF.Width=16
PDF.Height=10

StatsLogOut <- gsub(".csv", "_StatsLog.txt", FileInpCSV)
csvLongOut <- gsub(".csv", "_long.csv", FileInpCSV)
aucLongOut <- gsub(".csv", "_auc.csv", FileInpCSV)
pdfOut <- gsub(".csv","_Stats.pdf",FileInpCSV)
setwd(Directory)

t <- Short2Long(FileInpCSV)

if(!is.na(RefTreat)){
  treat_vec0 <- sort(unique(t$treatid));
  treat_vec <- c(treat_vec0[treat_vec0!=RefTreat], RefTreat);
} else {
  treat_vec <- sort(unique(t$treatid))
}

t$treatid  <- factor(t$treatid, levels= treat_vec)

# --- generate data in long format --------------------
write.csv(file=csvLongOut , t, quote=FALSE, row.names=FALSE)

# --- generate AUC data --------------------
t$test <- ifelse(t$time >0,"test","hab")
t$value0<-t$value
ta <- summaryBy( value ~ aid + treatid + test, t, FUN=c(sum),keep.names=TRUE)
write.csv(file=aucLongOut , ta, quote=FALSE, row.names=FALSE)

#-------------------------------------------
# output to PDF
#
pdf(pdfOut,w=PDF.Width,h=PDF.Height)


# --- plot time course data --------------------
plot.TC <- plot_time_dose(t, title=basename(FileInpCSV), xlab="time", ylab="total distance");
print(plot.TC);

# --- plot AUC data --------------------
ta.ss <- summaryBy( value ~ treatid + test, ta, FUN=c(mean,se,length),keep.names=TRUE)
plot.AUC1 <- ggplot(ta.ss, aes(x=treatid, y=value.mean, group=treatid, fill=treatid, width=0.7 )) + 
  geom_errorbar(aes(ymin=value.mean-value.se, ymax=value.mean+value.se), width=.3, position=position_dodge(.9) ) +
  geom_bar(stat="identity") + facet_grid(. ~ test) + theme_bw() + ylab("total distance") + xlab("") +
  theme(text = element_text(size=14), axis.text.x = element_text(angle=50, vjust=0.5))

print(plot.AUC1);

# --- plot AUC boxplots data --------------------

plot.AUC2 <- ggplot(ta, aes(x=treatid, y=value, group=treatid, fill=treatid )) + 
  geom_boxplot(aes(fill = factor(treatid))) + geom_jitter(position = position_jitter(h = 0,w=0.15),size=1.5) +
  facet_grid(. ~ test) + theme_bw() + ylab("total distance") + xlab("") +
  theme(text = element_text(size=14), axis.text.x = element_text(angle=50, vjust=0.5))

print(plot.AUC2);

# One-Way Anova
#----------------------------------
#  Shapiro-Wilk test of normality
#----------------------------------

ta.test <- subset(ta, test=="test")
SW.tab <- data.frame(treatment =levels(ta.test$treatid), "p-value"=NA )
trt.levels <- levels(ta.test$treatid);
for (it in 1:length(trt.levels)){
  (q <- shapiro.test(subset(ta.test$value, ta.test$treatid==trt.levels[it])));
  SW.tab[it,2] <- q$p.value; 
} 

#export table
textplot( capture.output(SW.tab), valign="top", cex=cexout  )
title("Shapiro-Wilk test of normality")

#Homogeneity of Variances
stdout <- bartlett.test( value ~ treatid, data=ta.test)
textplot( capture.output(stdout), valign="top", cex=cexout)
title("Bartlett test of homogeneity of variances")

#If there is a violation of homogeneity of variance, then one approach is to use the 
#Welch correction (as described in your text), which adjusts the degrees of freedom to compensate for the unequal variances. 
# finally one-way anova

stdout <- oneway.test( value ~ treatid, data=ta.test)
textplot( capture.output(stdout), valign="top",  cex=cexout)
title("One-way analysis of means (not assuming equal variances)")

ta.test$treatid  <- factor(ta.test$treatid, levels= treat_vec)
mm1 <- lm ( value ~ treatid, ta.test)
stdout <- (lsmeans2 (mm1, trt.vs.ctrlk ~ treatid, adjust="BH" ))
textplot( capture.output(stdout[[2]]), valign="top", cex=cexout)
title("least-squares means. differences from control")

# -----------------------------
# longitudinal model
# -----------------------------

tt <- subset(t, test=="test")
tt$time <- as.factor(tt$time)

nestinginfo <- groupedData( value ~ as.numeric(treatid) * as.numeric(time) | aid, data= tt)
fit.cs  <- gls(value ~ factor(treatid)*factor(time), data=nestinginfo, corr=corCompSymm(, form= ~ 1 | aid), method="REML")
fit.ar1 <- gls(value ~ factor(treatid)*factor(time), data=nestinginfo, corr=corAR1(, form= ~ 1 | aid), method="REML")
fit.ar1het <- gls(value ~ factor(treatid)*factor(time), data=nestinginfo, corr=corAR1(, form= ~ 1 | aid), weights=varIdent(form = ~ 1 | time), method="REML")

stdout1 <- anova(fit.cs, fit.ar1, fit.ar1het)
stdout2 <- anova(fit.ar1het)
stdout3 <- lsmeans2(fit.ar1het, trt.vs.ctrlk ~ treatid|time,  adjust="BH" )

textplot( capture.output(stdout2), valign="top", cex=cexout)
title("two-way anova on GLS model with correlations and unequal variances")

textplot( capture.output(stdout3[[2]]), valign="top", cex=cexout)
title("least-squares means. differences from control at fixed times")

dev.off()
rm(list=ls())

