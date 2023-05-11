#!/usr/bin/env/Rscript
##########
#Get Stats from hmp Files#
#by Kyle Parker

#just give the path to the function 'F:/path/hmp.txt' , Statname "Stats.Pop1"#

print("test1")

Hmpstats<-function(hmpfilePath,Statname){
  hmpfile<-read.delim(hmpfilePath)
  numLines<-(ncol(hmpfile)-11)
  chrs<-unique(hmpfile[,3])
  snps<-nrow(hmpfile)
  
  rowN<-c("Total.SNPs","Total.Het","Total.Indel","Total.NA","Total.Percent.Het","Total.Percent.NA")
  for(n in 1:length(chrs)){
    nn<-c(paste(chrs[n],".SNPs",sep = ""),paste(chrs[n],".Het",sep = ""),paste(chrs[n],".Indel",sep = ""),paste(chrs[n],".NA",sep = ""), paste(chrs[n],"Percent.Het",sep = ""),paste(chrs[n],"Percent.NA",sep = ""))
  rowN<-c(rowN,nn)
    }
  
  output<-matrix(data=NA,nrow=((6*length(chrs))+6),ncol = numLines,dimnames = list(rowN,colnames(hmpfile)[-c(1:11)]))
  
  
  for (i in 1:numLines){
    chz<-hmpfile[,c(3,4,(11+i))]
    df<-chz[,3]
    output[1,i]<-length(df)
    output[2,i]<-length(which(df != "A" & df!="T" & df != "C" & df != "G" & df != "N" & df != "-" & df !="+"))
    output[3,i]<-length(which(df =="-" | df=="+"))
    output[4,i]<-length(which(df =="N"))
    output[5,i]<-(output[2,i])/output[1,i]
    output[6,i]<-length(which(df =="N"))/length(df)
    write.table(output,file=paste(Statname,".Stats.txt",sep = ""))
    for (m in 1:length(chrs)){
      df2<- chz[which(chz[,1]==chrs[m]),3]
      output[1+(6*m),i]<-length(df2)
      output[2+(6*m),i]<-length(which(df2 != "A" & df2!="T" & df2 != "C" & df2 != "G" & df2 != "N" & df2 != "-" & df2 !="+"))
      output[3+(6*m),i]<-length(which(df2 =="-" | df2=="+"))
      output[4+(6*m),i]<-length(which(df2=="N"))
      output[5+(6*m),i]<-(output[2+(6*m),i])/length(df2)
      output[6+(6*m),i]<-length(which(df2=="N"))/length(df2)
      write.table(output,file=paste(Statname,".Stats.txt",sep = ""))
    }
  }
  write.table(output,file=paste(Statname,".Stats.txt",sep = ""))
  
}

print("Made it to stage 2")


Hmpstats(hmpfilePath = "/scratch/user/kxparker/NIHvcf/donevcf/NIH.DT.raw.hmp.txt", Statname = "DT.raw.hmp")
