devcor<-function(data1,data2){
	corvalue=data.frame();
	for(i in 1:ncol(data1)){
	for(j in 1:ncol(data2)){
			corvalue[i,j]=cor(data1[,i],data2[,j]);
		}
	}
	rownames(corvalue)=colnames(data1);
	mycolor=rainbow(ncol(data2));
	matplot(corvalue,xlab="",col=mycolor,xaxt="n",pch=3,ylab="",ylim=c(0,0.5));
	axis(1,1:nrow(corvalue),labels=rownames(corvalue),las=2);
	for(j in 1:ncol(corvalue)){lines(corvalue[,j],col=mycolor[j]);}
	corvalue;
}

rowcor<-ƒƒfunction(data1,data2,corcutoff=0.85){
	corvalues=data.frame();
	pvalues=data.frame();
	sigvalue=data.frame();
	for(i in 1:nrow(data1)){
		for(j in 1:nrow(data2)){
			cor.value=cor.test(as.numeric(as.list(data1[i,])),as.numeric(as.list(data2[j,])),method="pearson");
			corvalues[i,j]=cor.value$estimate;
			pvalues[i,j]=cor.value$p.value;
		}
	}
	rownames(corvalues)=rownames(data1);
	colnames(corvalues)=rownames(data2);
	rownames(pvalues)=rownames(data1);
	colnames(pvalues)=rownames(data2);
	sigcount=0;
	for(i in 1:nrow(corvalues)){
		rowmax=0;
		rowmin=0;
		colmax=0;
		colmin=0;
		for(j in 1:ncol(corvalues)){
			if(corvalues[i,j]>rowmax){
				rowmax=corvalues[i,j];
				colmax=j;
			}
			if(corvalues[i,j]<rowmin){
				rowmin=corvalues[i,j];
				colmin=j;
			}
		}
		if(rowmax>=corcutoff){
			#sigvalue=rbind(sigvalue,c(paste(rownames(corvalues[i,]),colnames(corvalues[,colmax]),sep="#"),corvalues[i,colmax],pvalues[i,colmax]));
			sigvalue=rbind(sigvalue,data.frame(id=paste(rownames(corvalues)[i],colnames(corvalues)[colmax],sep="#"),estimate=corvalues[i,colmax],p.value=pvalues[i,colmax]));
		}
		else{
			if(abs(rowmin)>=corcutoff){
				sigvalue=rbind(sigvalue,data.frame(id=paste(rownames(corvalues)[i],colnames(corvalues)[colmin],sep="#"),estimate=corvalues[i,colmin],p.value=pvalues[i,colmin]));
			}
		}
	}
	sigvalue=sigvalue[rev(order(abs(sigvalue$estimate))),];
	list("cor.values"=corvalues,"p.values"=pvalues,"sig.values"=sigvalue);
}

genecor<-function(data1,data2,topnum){
	##used to compare orthologous gene expression patterns
	corvalues=data.frame();
	count=1;
	#topnum=round(nrow(data1)*toppercent/100);
	for(i in 1:nrow(data1)){
		cor.value=cor.test(as.numeric(as.list(data1[i,])),as.numeric(as.list(data2[i,])),method="pearson");
		if(is.na(cor.value)){
			next;
		}
		corvalues[count,1]=rownames(data1[i,]);
		corvalues[count,2]=cor.value$estimate;
		corvalues[count,3]=cor.value$p.value;
		count=count+1;
	}
	colnames(corvalues)=c("geneid","estimate","pvalue");
	corvalues=corvalues[rev(order(abs(corvalues$estimate))),];
	corfinal=corvalues[1:topnum,c(2,3)];
	rownames(corfinal)=corvalues$geneid[1:topnum];
	corfinal;
}



goanalysis<-function(allGenes,interestGenes,GOTable,JobName="",GraphyTopNode=5,GenTableTopNode=10){
	allGeneFac<-factor(as.integer(allGenes %in% interestGenes));
	names(allGeneFac)=allGenes;
	print ("Analysing GOTerm 'MF'...");
	MF.topGOdata<-new("topGOdata",ontology="MF",allGene=allGeneFac,annot=annFUN.gene2GO,gene2GO=GOTable);
	MF.Fisher.test<-runTest(MF.topGOdata,algorithm="classic",statistic="fisher");
	MF.result<-GenTable(MF.topGOdata,classic=MF.Fisher.test,orderBy="classic",ranksOf="classic",topNodes=GenTableTopNode);
	filename=paste(JobName,"MF","pdf",sep=".");
	pdf(filename);
	showSigOfNodes(MF.topGOdata,score(MF.Fisher.test),sigForAll=TRUE,useFullNames=TRUE,firstSigNodes=GraphyTopNode,useInfo="all");
	dev.off();
	print ("Analysing GOTerm 'BP'...");
	BP.topGOdata<-new("topGOdata",ontology="BP",allGene=allGeneFac,annot=annFUN.gene2GO,gene2GO=GOTable);
	BP.Fisher.test<-runTest(BP.topGOdata,algorithm="classic",statistic="fisher");
	BP.result<-GenTable(BP.topGOdata,classic=BP.Fisher.test,orderBy="classic",ranksOf="classic",topNodes=GenTableTopNode);
	filename=paste(JobName,"BP","pdf",sep=".");
	pdf(filename);
	showSigOfNodes(BP.topGOdata,score(BP.Fisher.test),sigForAll=TRUE,useFullNames=TRUE,firstSigNodes=GraphyTopNode,useInfo="all");
	dev.off();
	print ("Analysing GOTerm 'CC'...");
	CC.topGOdata<-new("topGOdata",ontology="CC",allGene=allGeneFac,annot=annFUN.gene2GO,gene2GO=GOTable);
	CC.Fisher.test<-runTest(CC.topGOdata,algorithm="classic",statistic="fisher");
	CC.result<-GenTable(MF.topGOdata,classic=MF.Fisher.test,orderBy="classic",ranksOf="classic",topNodes=GenTableTopNode);
	filename=paste(JobName,"CC","pdf",sep=".");
	pdf(filename);
	showSigOfNodes(CC.topGOdata,score(CC.Fisher.test),sigForAll=TRUE,useFullNames=TRUE,firstSigNodes=GraphyTopNode,useInfo="all");
	dev.off();
	
	list("allGene"=allGeneFac,"MF.topGOdata"=MF.topGOdata,"MF.result"=MF.result,"MF.fisher.test"=MF.Fisher.test,"BP.topGOdata"=BP.topGOdata,"BP.result"=BP.result,"BP.fisher.test"=BP.Fisher.test,"CC.topGOdata"=CC.topGOdata,"CC.result"=CC.result,"CC.fisher.test"=CC.Fisher.test);
}