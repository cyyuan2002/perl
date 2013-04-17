foldchange<-function(datasource,groupA,groupB,foldercut=2,groupAName="groupA",groupBName="groupB",zerocheck=0){
	rowN=nrow(datasource);
	grouplengthA=length(groupA);
	grouplengthB=length(groupB);
	dataA=data.frame();
	dataB=data.frame();
	if(length(groupA)==1){
		dataA=datasource[,groupA,drop=FALSE];
		colnames(dataA)=c("SumA");
	}
	else{
		dataA=as.data.frame(rowSums(datasource[,groupA]));
		rownames(dataA)=rownames(datasource);
		colnames(dataA)=c("SumA");
	}
	
	if(length(groupB)==1){
		dataB=datasource[,groupB,drop=FALSE];
		colnames(dataB)=c("SumB");
	}
	else{
		dataB=as.data.frame(rowSums(datasource[,groupA]));
		rownames(dataB)=rownames(datasource);
		colnames(dataB)=c("SumB");
	}
	SumMerge=merge(dataA,dataB,by=0);
	rownames(SumMerge)=SumMerge[,1];
	SumMerge=SumMerge[,2:3];
	SumMerge_Nozero=SumMerge[apply(SumMerge,1,function(x) all(x != 0)),];
	zerofolds=list();
	zerocontains=0;
	if(nrow(SumMerge)>nrow(SumMerge_Nozero)){
		SumMerge=SumMerge_Nozero;
		zerocontains=1;
		if(zerocheck==1){
			zerofolds=zerochange(datasource,groupA,groupB,1,groupAName,groupBName);
		}
	}
	foldA=as.data.frame(SumMerge[,1]/SumMerge[,2]);
	rownames(foldA)=rownames(SumMerge);
	colnames(foldA)=groupAName;
	selectA=foldA[foldA[,1]>=foldercut,,drop=FALSE];
	foldB=as.data.frame(SumMerge[,2]/SumMerge[,1]);
	rownames(foldB)=rownames(SumMerge);
	colnames(foldB)=groupBName;
	selectB=foldB[foldB[,1]>=foldercut,,drop=FALSE];
	finalA=merge(selectA,datasource[,c(groupA,groupB)],by=0);
	rownames(finalA)=finalA[,1];
	finalA=finalA[,-1];
	finalA=finalA[with(finalA,order(-finalA[,1])),,];
	finalB=merge(selectB,datasource[,c(groupA,groupB)],by=0);
	rownames(finalB)=finalB[,1];
	finalB=finalB[,-1];
	finalB=finalB[with(finalB,order(-finalB[,1])),,];
	if(zerocontains==0){
		list("groupA"=finalA,"groupB"=finalB)
	}
	else{
		if(zerocheck==1){
			list("groupA"=finalA,"groupB"=finalB,"zeroA"=zerofolds$groupA,"zeroB"=zerofolds$groupB,"median"=zerofolds$median);
		}
	}
}

zerochange<-function(datasource,groupA,groupB,foldercut=1,groupAName="groupA",groupBName="groupB"){
	rowN=nrow(datasource);
	grouplengthA=length(groupA);
	grouplengthB=length(groupB);
	dataA=as.data.frame();
	dataB=as.data.frame();
	if(length(groupA)==1){
		dataA=datasource[,groupA,drop=FALSE];
		colnames(dataA)=c("SumA");
		medianGP_A=median(datasource[,groupA]);
		
	}
	else{
		dataA=as.data.frame(rowSums(datasource[,groupA]));
		rownames(dataA)=rownames(datasource);
		colnames(dataA)=c("SumA");
		medianGP_A=median(dataA);
	}
	thresholdA=medianGP_A*foldercut;
	
	if(length(groupB)==1){
		dataB=datasource[,groupB,drop=FALSE];
		colnames(dataB)=c("SumB");
		medianGP_B=median(datasource[,groupB]);
		thresholdB=medianGP_B*foldercut;
	}
	else{
		dataB=as.data.frame(rowSums(datasource[,groupA]));
		rownames(dataB)=rownames(datasource);
		colnames(dataB)=c("SumB");
		medianGP_B=median(dataB);
	}
	thresholdB=medianGP_B*foldercut;

	SumMerge=merge(dataA,dataB,by=0);
	rownames(SumMerge)=SumMerge[,1];
	SumMerge=SumMerge[,2:3];
	ZeroA=SumMerge[SumMerge[,1]==0,];
	ZeroB=SumMerge[SumMerge[,2]==0,];
	selA=ZeroB[ZeroB[,1]>=thresholdA,];
	selB=ZeroA[ZeroA[,2]>=thresholdB,];
	selA=selA[with(selA,order(-selA[,1])),drop=FALSE];
	selB=selB[with(selB,order(-selB[,2])),drop=FALSE];
	dataA=datasource[rownames(selA),groupA];
	dataB=datasource[rownames(selB),groupB];
	median=list("groupA"=medianGP_A,"groupB"=medianGP_B);
	list("groupA"=selectA,"groupB"=selectB,"median"=median);
}

exptest<-function(datasource,groupA,groupB,adjMethod="fdr",adjvalue=0.001){
	library(edgeR);
	tagwise<-exactTest(datasource,pair=c(groupA,groupB));
	padj=as.data.frame(p.adjust(tagwise$table$PValue,method=adjMethod));
	rownames(padj)=rownames(tagwise$table);
	colnames(padj)=adjMethod;
	padj=merge(tagwise$table,padj,by=0);
	rownames(padj)=padj[,1];
	padj=padj[,2:5];
	padj.flt=padj[padj[,4]<adjvalue,];
	padj.flt=padj.flt[with(padj.flt,order(padj.flt[,3])),,];
	tagwise.up=padj.flt[padj.flt[,1]>0,];
	tagwise.down=padj.flt[padj.flt[,1]<0,];
	list("genenum"=nrow(padj.flt),"up"=tagwise.up,"down"=tagwise.down);
}

testaddcount<-function(testresult,datasource,cols){
	upgenecount=merge(datasource[,c(cols)],testresult$up,by=0);
	rownames(upgenecount)=upgenecount[,1];
	upgenecount=upgenecount[,-1];
	upgenecount=upgenecount[with(upgenecount,order(upgenecount[,"PValue"])),,];
	downgenecount=merge(datasource[,c(cols)],testresult$down,by=0);
	rownames(downgenecount)=downgenecount[,1];
	downgenecount=downgenecount[,-1];
	downgenecount=downgenecount[with(downgenecount,order(downgenecount[,"PValue"])),,];
	list("up"=upgenecount,"down"=downgenecount);
}

diffcof<-function(x,y){
	i=ncol(x);
	j=ncol(y);
	sumx=0;
	sumy=0;
	if(i!=j){
		print ("column number is not the same");
	}
	else{
		for(m in 1:i){
			sumx=sumx+x[,m]*y[,m];
			sumy=sumy+(x[,m]^2+y[,m]^2);
		}
		disscore=1-((2*sumx)/sumy)
	}
	as.vector(disscore);
}

diffcof.2<-function(x,y){
	i=ncol(x);
	j=ncol(y);
	sumx=0;
	if(i!=j){
		print ("column number is not the same");
	}
	else{
		for(m in 1:i){
			a=x[,m]*y[,m];
			b=x[,m]^2+y[,m]^2;
			c=a/b;
			sumx=sumx+c;
		}
		disscore=1-((2*sumx)/i)
	}
	as.vector(disscore);
}

gostats<-function(GOList){
	require("topGO");
	require("GO.db");
	allGene=GOList$allGene;
	totolGenesnum=length(allGene);
	instGenesnum=length(allGene[allGene==1]);
	MFStat=termStat(GOList$MF.topGOdata);
	MFTotalnum=MFStat["GO:0003674",1];
	MFSignum=MFStat["GO:0003674",2];
	BPStat=termStat(GOList$BP.topGOdata);
	BPTotalnum=BPStat["GO:0008150",1];
	BPSignum=BPStat["GO:0008150",2];
	CCStat=termStat(GOList$CC.topGOdata);
	CCTotalnum=CCStat["GO:0005575",1];
	CCSignum=CCStat["GO:0005575",2];
	MFroots=as.vector(GOMFCHILDREN$"GO:0003674");
	BProots=as.vector(GOBPCHILDREN$"GO:0008150");
	CCroots=as.vector(GOCCCHILDREN$"GO:0005575");
	MFCounts=MFStat[MFroots,];
	BPCounts=BPStat[BProots,];
	CCCounts=CCStat[CCroots,];
	MFCounts=MFCounts[!is.na(MFCounts$Annotated),];
	BPCounts=BPCounts[!is.na(BPCounts$Annotated),];
	CCCounts=CCCounts[!is.na(CCCounts$Annotated),];
	MFCounts$Annopercent=round(MFCounts$Annotated/MFTotalnum*100,2);
	MFCounts$Sigpercent=round(MFCounts$Significant/MFSignum*100,2);
	BPCounts$Annopercent=round(BPCounts$Annotated/BPTotalnum*100,2);
	BPCounts$Sigpercent=round(BPCounts$Significant/BPSignum*100,2);
	CCCounts$Annopercent=round(CCCounts$Annotated/CCTotalnum*100,2);
	CCCounts$Sigpercent=round(CCCounts$Significant/CCSignum*100,2);
	MFGONames<-unlist(lapply(rownames(MFCounts),function(x){Term(GOTERM[[x]])}))
	BPGONames<-unlist(lapply(rownames(BPCounts),function(x){Term(GOTERM[[x]])}))
	CCGONames<-unlist(lapply(rownames(CCCounts),function(x){Term(GOTERM[[x]])}))
	MF.pvalue=score(GOList$MF.fisher.test)[rownames(MFCounts)];
	BP.pvalue=score(GOList$BP.fisher.test)[rownames(BPCounts)];
	CC.pvalue=score(GOList$CC.fisher.test)[rownames(CCCounts)];
	MFCounts$Pval=as.vector(MF.pvalue);
	BPCounts$Pval=as.vector(BP.pvalue);
	CCCounts$Pval=as.vector(CC.pvalue);
	MFCounts$GOName= MFGONames;
	BPCounts$GOName= BPGONames;
	CCCounts$GOName= CCGONames;
	list("MF.Counts"=MFCounts,"BP.Counts"=BPCounts,"CC.Counts"=CCCounts);
}

gomatrix<-function(GOList,bgName=""){
	values=GOList[[1]];
	if(bgName!=""){
		MF.matrix=matrix(NA,nrow=nrow(values$MF.Counts),ncol=length(GOList)+1);
		BP.matrix=matrix(NA,nrow=nrow(values$BP.Counts),ncol=length(GOList)+1);
		CC.matrix=matrix(NA,nrow=nrow(values$CC.Counts),ncol=length(GOList)+1);
	}
	else{
		MF.matrix=matrix(NA,nrow=nrow(values$MF.Counts),ncol=length(GOList));
		BP.matrix=matrix(NA,nrow=nrow(values$BP.Counts),ncol=length(GOList));
		CC.matrix=matrix(NA,nrow=nrow(values$CC.Counts),ncol=length(GOList));	
	}

	MF.pval.matrix= MF.matrix;
	BP.pval.matrix= BP.matrix;
	CC.pval.matrix= CC.matrix;
	
	for (i in 1:length(GOList)){
		values=GOList[[i]];
		MF.matrix[,i]=values$MF.Counts$Sigpercent;
		MF.pval.matrix[,i]= values$MF.Counts$Pval;
		BP.matrix[,i]=values$BP.Counts$Sigpercent;
		BP.pval.matrix[,i]= values$BP.Counts$Pval;
		CC.matrix[,i]=values$CC.Counts$Sigpercent;
		CC.pval.matrix[,i]= values$CC.Counts$Pval;
	}
	
	
	if(bgName!=""){
		j=length(GOList);
		j=j+1;
		values=GOList[[1]];
		MF.matrix[,j]=values$MF.Counts$Annopercent;
		BP.matrix[,j]=values$BP.Counts$Annopercent;
		CC.matrix[,j]=values$CC.Counts$Annopercent;
		MF.pval.matrix[,j]=rep(1,nrow(values$MF.Counts));
		BP.pval.matrix[,j]=rep(1,nrow(values$BP.Counts));
		CC.pval.matrix[,j]=rep(1,nrow(values$CC.Counts));
	}
	
	MF=as.data.frame(MF.matrix,row.names=values$MF.Counts$GOName);
	MF.pval=as.data.frame(MF.pval.matrix,row.names=values$MF.Counts$GOName);
	BP=as.data.frame(BP.matrix,row.names=values$BP.Counts$GOName);
	BP.pval=as.data.frame(BP.pval.matrix,row.names=values$BP.Counts$GOName);
	CC=as.data.frame(CC.matrix,row.names=values$CC.Counts$GOName);
	CC.pval=as.data.frame(CC.pval.matrix,row.names=values$CC.Counts$GOName);
	
	if(bgName!=""){
		names(MF)=c(names(GOList),bgName);
		names(MF.pval)=c(names(GOList),bgName);
		names(BP)=c(names(GOList),bgName);
		names(BP.pval)=c(names(GOList),bgName);
		names(CC)=c(names(GOList),bgName);
		names(CC.pval)=c(names(GOList),bgName);
	}
	else{
		names(MF)=names(GOList);
		names(MF.pval)=names(GOList);
		names(BP)=names(GOList);
		names(BP.pval)=names(GOList);
		names(CC)=names(GOList);
		names(CC.pval)=names(GOList);
	}
	
	list("MF.val"=MF,"MF.pval"=MF.pval,"BP.val"=BP,"BP.pval"=BP.pval,"CC.val"=CC,"CC.pval"=CC.pval);
	
}


plotgostats<-function(GOCounts,GOPvalue="",main="",rmZero=FALSE,sort=FALSE,isBackground=TRUE){
	library(RColorBrewer);
	if(sort==TRUE){
		GOCounts=GOCounts[with(GOCounts,order(-rowSums(GOCounts))),,drop=FALSE];
	}
	if(ncol(GOCounts)==1){
		#colors=rainbow(nrow(GOCounts));
		if(rmZero==TRUE){
			GOCounts=GOCounts[GOCounts[1]>0,,drop=FALSE];
		}
		colors=colorRampPalette(brewer.pal(9, "RdYlBu"))(nrow(GOCounts));
		opt<-par(mar=c(6, 12, 2, 5),xpd=FALSE,cex.axis=0.8);
		bp= barplot(GOCounts[,1],horiz=TRUE,beside=T,main=main,col=colors);
		text(GOCounts[,1],round(bp,1),GOCounts[,1],pos=4,cex=0.8);
		des=shortName(rownames(GOCounts));
		axis(2,at=bp,labels= des ,cex.axis=0.9,las=2);
		par(opt);
	}
	else{
		if(rmZero==TRUE){
			if(isBackground==TRUE){
				colnum=ncol(GOCounts)-1;
				if(colnum>1){
					GOCounts=GOCounts[rowSums(GOCounts[,1:colnum])>0,];
				}
				else{
					GOCounts=GOCounts[GOCounts[1]>0,,drop=FALSE];
				}
			}
			else{
				GOCounts=GOCounts[rowSums(GOCounts)>0,];
			}
		}
		opt<-par (mar=c(6,12,2,8), xpd=TRUE,cex.axis=0.01);
		freq=t(as.matrix(GOCounts));
		allcolors=colorRampPalette(brewer.pal(9, "RdYlBu"))(100);
		partsize=round(100/(ncol(GOCounts)+1));
		colors=allcolors[c(1:ncol(GOCounts))*partsize];
		bp<-barplot(freq,horiz=TRUE, beside=T,legend.text =colnames(GOCounts),axes=1,col=colors,main=main);
		maxnumber=max(GOCounts);
		maxlength=floor(max(GOCounts)/20)*20;
		axis(1,cex.axis=0.8, labels=seq(0,maxlength,by=20),at=seq(0,maxlength,by=20))
		des=shortName(rownames(GOCounts));
		colnumber=ncol(GOCounts);
		text(freq, round(bp,1),freq,pos=4,cex=0.6)
		if(is.data.frame(GOPvalue)){
			GOPvalue=GOPvalue[rownames(GOCounts),];
			pvalt=t(as.matrix(GOPvalue));
			pSig=pvalt<0.01;
			psig=(pvalt >= 0.01) & (pvalt < 0.05);
			if(length(pSig[pSig==TRUE])>0){
				text(freq[pSig],round(bp[pSig],1),"**",offset=2,pos=4,cex=0.8);
			}
			if(length(psig[psig==TRUE])>0){
				text(freq[psig],round(bp[psig],1),"*",offset=2,pos=4,cex=0.8);
			}
		}
		axis(2,at=(bp[1,]+bp[colnumber,])/2, labels=des, cex.axis = 0.8, las=2);
		par(opt);
	}
}

shortName <- function (aString, aWidth=25){
    numCh <- nchar(aString)
    aString2 <- substr(aString, 1, aWidth) # Limit names of terms to 'aWidth' characters
    aString3 <- paste(aString2, ifelse(numCh > (aWidth-5), "...", ""), sep="")
    return(aString3)
}

