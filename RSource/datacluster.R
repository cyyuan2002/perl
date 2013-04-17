datacluster<-function(datasource,fileN="cluster",cluster2way=FALSE){
	datacolnum=ncol(datasource);
	dataS.hr=hclust(as.dist(1-cor(t(datasource),method="spearman")),method="complete",member=NULL);
	if(cluster2way){
		dataS.hc=hclust(as.dist(1-cor(datasource,method="pearson")),method="complete",member=NULL);
	}
	library("gplots");
	filename=paste(fileN,".pdf",sep="");
	drive<-pdf(file=filename);
	if(cluster2way){
			heatmap.2(as.matrix(datasource),Rowv=as.dendrogram(dataS.hr),Colv=as.dendrogram(dataS.hc),col=greenred(75),scale="row",labRow="",trace="none");
			dev.off();
			list("rowcluster"=dataS.hr,"colcluster"=dataS.hc);
	}
	else{
			#heatmap.2(as.matrix(datasource),Rowv=as.dendrogram(dataS.hr),Colv=NA,col=greenred(75),scale="row",labRow="",trace="none",dendrogram="row");
			heatmap(as.matrix(datasource),Rowv=as.dendrogram(dataS.hr),Colv=NA,col=greenred(75),scale="row");
			dev.off();
			list("rowcluster"=dataS.hr);
	}
}

kmeancluster<-function(datasource,clusternum,fileN="kcluster",center=0,drawgraph=FALSE){
	dataTemp=t(scale(t(datasource)));
	dataTemp=dataTemp[!(apply(is.nan(dataTemp),1,any)),];
	if(!is.matrix(center)){
		km<-kmeans(dataTemp,clusternum);
	}
	else{
		km<-kmeans(dataTemp,centers=center);
		clusternum=nrow(center);
	}
	label=NULL;data=NULL;hr=NULL;
	if(drawgraph){
			filename=paste(fileN,".pdf",sep="");
			drive<-pdf(file=filename);
			colnum=ceiling(clusternum/4);
			par(mfrow=c(colnum,4));
	}
	for(i in 1:clusternum){
		datasource.label=labels(km$cluster[km$cluster==i]);
		label[[i]]<-datasource.label;
		datasource.data=datasource[datasource.label,];
		data[[i]]=datasource.data;
		#datasource.hr=hclust(as.dist(1-cor(t(datasource.data),method="spearman")),method="complete",member=NULL);
		#hr[[i]]=datasource.hr;
		datasource.mean=apply(datasource.data,1,mean);
		datasource.tmean=mean(datasource.mean);
		datasource.change=datasource.data-(datasource.mean-datasource.tmean);
		datasource.median=apply(datasource.change,2,median);
		if(drawgraph){

		#library("gplots");
		#heatmap(as.matrix(datasource.data),Rowv=as.dendrogram(datasource.hr),Colv=NA,col=greenred(65),scale="row")
		#matplot(t(datasource.change),xlab="",pch=20,col="black",ylab="",main=paste(fileN,i));
		#for(j in 1:nrow(datasource.change)){lines(t(datasource.change[j,]));}
		#boxplot(datasource.change,ylab="");
			plot(datasource.median,ylab="",main=paste("Cluster",i));
			lines(datasource.median);
		}
	}
	if(drawgraph){
		dev.off();
	}
	list("data"=data,"labels"=label,"hrs"=hr,"kmeans"=km);
}

kmeanclustersim<-function(dataSource,dataSim,corcutoff=0.85,iter=50){
	centers=length(dataSim$data);
	kmeans=NULL;
	corstat=data.frame();
	medianSim=medianOfkmeans(dataSim);
	for(i in 1:iter){
		kmeans[[i]]<-kmeancluster(dataSource,centers);
		median<-medianOfkmeans(kmeans[[i]]);
		cor=rowcor(medianSim,median,corcutoff=corcutoff);
		cornum=nrow(cor$sig.values);
		totalestimate=sum(cor$sig.values$estimate);
		corstat=rbind(corstat,data.frame(cornum=cornum,estimate=totalestimate));
	}
	#corstat;
	corstat=corstat[with(corstat,order(-cornum,-estimate)),];
	linenum=rownames(corstat[1,]);
	#list("corstat"=corstat,"linenum"=linenum,"kmeans"=kmeans);
	kmeans[[as.numeric(linenum)]];
}


medianOfkmeans<-function(datasource,rowname="row"){
	clusternum=length(datasource$data);
	medians=data.frame();
	for(i in 1:clusternum){
		data=datasource$data[[i]];
		mean=apply(data,1,mean);
		tmean=mean(mean);
		change=data-(mean-tmean);
		median=apply(change,2,median);
		medians=rbind(medians,median);
	}
	rowN=paste(rowname,(1:clusternum),sep="");
	rownames(medians)=rowN;
	colnames(medians)=colnames(datasource$data[[i]]);
	medians;
}

kmeandraw<-function(datasource,fileN="kmeancluster",drawallpoints=FALSE,drawinonefile=TRUE,colnum=4){
	datanum=length(datasource$data);
	if(drawinonefile){
		filename=paste(fileN,"pdf",sep=".");
		drive<-pdf(file=filename);
		clusternum=length(datasource$data);
		rownum=ceiling(clusternum/colnum);
		par(mfrow=c(rownum,colnum));
	}
	for(i in 1:datanum){
		datasource.data=datasource$data[[i]];
		datasource.mean=apply(datasource.data,1,mean);
		datasource.tmean=mean(datasource.mean);
		datasource.change=datasource.data-(datasource.mean-datasource.tmean);
		datasource.median=apply(datasource.change,2,median);
		if(!drawinonefile){
			filename=paste(fileN,i,"pdf",sep=".");
			drive<-pdf(file=filename);
		}
		if(drawallpoints){
			matplot(t(datasource.change),xlab="",pch=20,col="grey",axes=FALSE,ylab="",main=paste("Cluster",i));
			axis(1,1:ncol(datasource.data),labels=colnames(datasource.data),las=2);
			for(j in 1:nrow(datasource.change)){lines(t(datasource.change[j,]),col="grey");}
			lines(datasource.median,col="orange",lwd="3");
		}
		else{
			plot(datasource.median,ylab="",main=paste("Cluster",i));
			lines(datasource.median);
		}
		if(!drawinonefile){
			dev.off();
		}
	}
	if(drawinonefile){
		dev.off();
	}
}

calsigcor<-function(kmdataA,kmdataB,rownameA="",rownameB="",corcutoff=0.85,negtive=FALSE){
	medianA=medianOfkmeans(kmdataA,rowname=rownameA);
	medianB=medianOfkmeans(kmdataB,rowname=rownameB);
	cor=rowcor(medianA,medianB,corcutoff=corcutoff,negtive=negtive);
	cor;
}

##function of correlations 
devcor<-function(data1,data2,species=""){
	corvalue=data.frame();
	for(i in 1:ncol(data1)){
		for(j in 1:ncol(data2)){
			corvalue[i,j]=cor(data1[,i],data2[,j]);
		}
	}
	rownames(corvalue)=colnames(data1);
	mycolor=rainbow(ncol(data2));
	matplot(corvalue,xlab="",col=mycolor,xaxt="n",pch=3,ylab="",ylim=c(0,0.6));
	axis(1,1:nrow(corvalue),labels=rownames(corvalue),las=2);
	for(j in 1:ncol(corvalue)){lines(corvalue[,j],col=mycolor[j]);}
	legend("topleft",legend=paste(species,colnames(data2),sep=" "),pch=20,col=mycolor,cex=0.8,bty="n")
	corvalue;
}

rowcor<-function(data1,data2,corcutoff=0.85,negtive=FALSE){
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
			if(negtive){
				if(corvalues[i,j]<rowmin){
					rowmin=corvalues[i,j];
					colmin=j;
				}
			}
		}
		if(rowmax>=corcutoff){
			sigvalue=rbind(sigvalue,data.frame(id=paste(rownames(corvalues)[i],colnames(corvalues)[colmax],sep="#"),estimate=corvalues[i,colmax],p.value=pvalues[i,colmax]));
		}
		else{
			if(negtive){
				if(abs(rowmin)>=corcutoff){
					sigvalue=rbind(sigvalue,data.frame(id=paste(rownames(corvalues)[i],colnames(corvalues)[colmin],sep="#"),estimate=corvalues[i,colmin],p.value=pvalues[i,colmin]));
				}
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



goanalysis<-function(allGenes,interestGenes,GOTable,JobName="",GraphyTopNode=5,pvalue=0.01,drawgraph=FALSE,writetable=TRUE,nodeSize=0){
	library(topGO);
	allGeneFac<-factor(as.integer(allGenes %in% interestGenes));
	names(allGeneFac)=allGenes;
	print ("Analysing GOTerm 'MF'...");
	if(nodeSize>0){
		MF.topGOdata<-new("topGOdata",ontology="MF",allGene=allGeneFac,annot=annFUN.gene2GO,gene2GO=GOTable);
	}
	else{
		MF.topGOdata<-new("topGOdata",ontology="MF",allGene=allGeneFac,annot=annFUN.gene2GO,gene2GO=GOTable,nodeSize=nodeSize);
	}
	MF.Fisher.test<-runTest(MF.topGOdata,algorithm="classic",statistic="fisher");
	MF.pnum=length(MF.Fisher.test@score[MF.Fisher.test@score<pvalue]);
	MF.result<-GenTable(MF.topGOdata,classic=MF.Fisher.test,orderBy="classic",ranksOf="classic",topNodes=MF.pnum);
	MF.pvalue=MF.Fisher.test@score[MF.Fisher.test@score<pvalue];
	MF.fdr=p.adjust(MF.pvalue,method="fdr");
	MF.result=fdrmerge(MF.result,MF.fdr);
	if(drawgraph){
		if(GraphyTopNode>MF.pnum){
			MF.topnode=MF.pnum;
		}
		else{
			MF.topnode=GraphyTopNode;
		}
		filename=paste(JobName,"MF","pdf",sep=".");
		pdf(filename);
		showSigOfNodes(MF.topGOdata,score(MF.Fisher.test),sigForAll=TRUE,useFullNames=TRUE,firstSigNodes=MF.topnode,useInfo="all");
		dev.off();
	}
	print ("Analysing GOTerm 'BP'...");
	if(nodeSize>0){
		BP.topGOdata<-new("topGOdata",ontology="BP",allGene=allGeneFac,annot=annFUN.gene2GO,gene2GO=GOTable);
	}
	else{
		BP.topGOdata<-new("topGOdata",ontology="BP",allGene=allGeneFac,annot=annFUN.gene2GO,gene2GO=GOTable,nodeSize=nodeSize);
	}
	BP.Fisher.test<-runTest(BP.topGOdata,algorithm="classic",statistic="fisher");
	BP.pnum=length(BP.Fisher.test@score[BP.Fisher.test@score<pvalue]);
	BP.result<-GenTable(BP.topGOdata,classic=BP.Fisher.test,orderBy="classic",ranksOf="classic",topNodes=BP.pnum);
	BP.pvalue=BP.Fisher.test@score[BP.Fisher.test@score<pvalue];
	BP.fdr=p.adjust(BP.pvalue,method="fdr")
	BP.result=fdrmerge(BP.result,BP.fdr);
	if(drawgraph){
		if(GraphyTopNode>BP.pnum){
			BP.topnode=BP.pnum;
		}
		else{
			BP.topnode=GraphyTopNode;
		}
		filename=paste(JobName,"BP","pdf",sep=".");
		pdf(filename);
		showSigOfNodes(BP.topGOdata,score(BP.Fisher.test),sigForAll=TRUE,useFullNames=TRUE,firstSigNodes=BP.topnode,useInfo="all");
		dev.off();
	}
	print ("Analysing GOTerm 'CC'...");
	if(nodeSize>0){
		CC.topGOdata<-new("topGOdata",ontology="CC",allGene=allGeneFac,annot=annFUN.gene2GO,gene2GO=GOTable);
	}
	else{
		CC.topGOdata<-new("topGOdata",ontology="CC",allGene=allGeneFac,annot=annFUN.gene2GO,gene2GO=GOTable,nodeSize=nodeSize);
	}
	CC.Fisher.test<-runTest(CC.topGOdata,algorithm="classic",statistic="fisher");
	CC.pnum=length(CC.Fisher.test@score[CC.Fisher.test@score<pvalue]);
	CC.result<-GenTable(CC.topGOdata,classic=CC.Fisher.test,orderBy="classic",ranksOf="classic",topNodes=CC.pnum);
	CC.pvalue=CC.Fisher.test@score[CC.Fisher.test@score<pvalue];
	CC.fdr=p.adjust(CC.pvalue,method="fdr")
	CC.result=fdrmerge(CC.result,CC.fdr);
	if(drawgraph){
		if(GraphyTopNode>CC.pnum){
			CC.topnode=CC.pnum;
		}
		else{
			CC.topnode=GraphyTopNode;
		}
		filename=paste(JobName,"CC","pdf",sep=".");
		pdf(filename);
		showSigOfNodes(CC.topGOdata,score(CC.Fisher.test),sigForAll=TRUE,useFullNames=TRUE,firstSigNodes=CC.topnode,useInfo="all");
		dev.off();
	}
	if(writetable){
		filename=paste(JobName,"BP","txt",sep=".");
		write.table(BP.result,file=filename,sep="\t",row.names=FALSE);
		filename=paste(JobName,"MF","txt",sep=".");
		write.table(MF.result,file=filename,sep="\t",row.names=FALSE);
		filename=paste(JobName,"CC","txt",sep=".");
		write.table(CC.result,file=filename,sep="\t",row.names=FALSE);
	}

	list("allGene"=allGeneFac,"MF.topGOdata"=MF.topGOdata,"MF.result"=MF.result,"MF.fisher.test"=MF.Fisher.test,"BP.topGOdata"=BP.topGOdata,"BP.result"=BP.result,"BP.fisher.test"=BP.Fisher.test,"CC.topGOdata"=CC.topGOdata,"CC.result"=CC.result,"CC.fisher.test"=CC.Fisher.test);
}

gocompare<-function(GOInfoA,GOInfoB){
	BP.merge=merge(GOInfoA$BP.result, GOInfoB$BP.result,by="GO.ID");
	MF.merge=merge(GOInfoA$MF.result, GOInfoB$MF.result,by="GO.ID");
	CC.merge=merge(GOInfoA$CC.result, GOInfoB$CC.result,by="GO.ID");
	list("BP"=BP.merge,"MF"=MF.merge,"CC"=CC.merge);
}

IDOverlap<-function(Ortholog,IDsets_A,IDsets_B){
	IDframeA=as.data.frame(IDsets_A,row.names=IDsets_A);
	IDframeB=as.data.frame(IDsets_B,row.names=IDsets_B);
	Merged1=merge(Ortholog,IDframeA,by="row.names");
	rownames(Merged1)=Merged1$V1;
	#Merged1;
	Merged2=merge(Merged1,IDframeB,by="row.names")[,1:2];
	colnames(Merged2)=c("ID_a","ID_b");
	Merged2;
}

fdrmerge<-function(datasource,datainfo){
	colength=ncol(datasource);
	for(i in 1:nrow(datasource)){
		datasource[i,colength+1]=datainfo[datasource[i,1]];
	}
	conames=colnames(datasource);
	conames[colength+1]="FDR";
	colnames(datasource)=conames;
	datasource;
}

exptest<-function(datasource,groupA,groupB,allGenes,GOTable,geneAnno,p.adjMethod="fdr",adj.value=0.001,writeresult=FALSE){
	library(edgeR);
	de.tagwise<-exactTest(datasource,pair=c(groupA,groupB),common.disp=FALSE);
	#FDR.num=sum(p.adjust(de.tagwise$table$p.value,method="fdr")<FDR);
	#BH.num=sum(p.adjust(de.tagwise$table$p.value,method="BH")<BH);
	topnum=sum(p.adjust(de.tagwise$table$p.value,method=p.adjMethod)<adj.value);
	top.tgw<-topTags(de.tagwise,n=topnum);
	exp.up.info=top.tgw$table[top.tgw$table$logFC>0,];
	exp.down.info=top.tgw$table[top.tgw$table$logFC<0,];
	exp.up.gname=annogene(rownames(exp.up.info),geneAnno);
	exp.down.gname=annogene(rownames(exp.down.info),geneAnno);
	if(writeresult){
		filenameU=paste(groupA,groupB,"upgene","txt",sep=".");
		write.table=write.table(exp.up.gname,file=filenameU,sep="\t",col.names=FALSE);
		filenameD=paste(groupA,groupB,"downgene","txt",sep=".")
		write.table=write.table(exp.down.gname,file=filenameD,sep="\t",col.names=FALSE);
	}
	num.up=nrow(exp.up.info);
	num.down=nrow(exp.down.info);
	jobnameU=paste(groupA,groupB,"up",sep=".");
	jobnameD=paste(groupA,groupB,"down",sep=".");
	exp.up.GO=goanalysis(allGenes,rownames(exp.up.info),GOTable,JobName=jobnameU,GraphyTopNode=10,drawgraph=writeresult,writetable=writeresult);
	exp.down.GO=goanalysis(allGenes,rownames(exp.down.info),GOTable,JobName=jobnameD,GraphyTopNode=10,drawgraph=writeresult,writetable=writeresult);
	upinfo=list("info"=exp.up.info,"gene"=exp.up.gname);
	downinfo=list("info"=exp.down.info,"gene"=exp.down.gname);
	list("toptags"=top.tgw,"ups"=upinfo,"downs"=downinfo,"ups.GO"=exp.up.GO,"downs.GO"=exp.down.GO);
}

annogene<-function(genelist,geneAnno){
	##need to be optimized
	geneinfo=data.frame();
	nlength=length(genelist);
	for(i in 1:nlength){
		if(is.na(amphi.anno[genelist[i],1])){
			geneinfo[i,1]="NA";
			geneinfo[i,2]="NA";
		}
		else{
			geneinfo[i,1]=amphi.anno[genelist[i],1];
			geneinfo[i,2]=amphi.anno[genelist[i],2];
		}
	}
	rownames(geneinfo)=genelist;
	colnames(geneinfo)=c("gname","description");
	geneinfo;
}

orthologfind<-function(datasource1,datasource2,orthinfo){
	orthrows=orthinfo[rownames(datasource1),];
	orthrows=orthrows[apply(orthrows,1,function(x) !all(is.na(x))),];
	orth1=datasource1[rownames(orthrows),];
	orth2=datasource2[as.vector(orthrows$V3),];
	list("data1"=orth1,"data2"=orth2,"orthinfo"=orthrows);
}

drawcandigene<-function(genelist,expdata,centers){
	geneexp=expdata[rownames(genelist),];
	filenum=ceiling(nrow(geneexp)/12);
	for(filecount in 1:filenum){
		filename=paste("candigene",filecount,"pdf",sep=".");
		pdf(filename);
		par(mfrow=c(3,4));
		loopS=(filecount-1)*12+1;
		loopE=filecount*12;
		if(loopE>nrow(geneexp)){
			loopE=nrow(geneexp);
		}
		for(i in loopS:loopE){
			a=as.vector(as.matrix(geneexp[i,]));
			b=centers[genelist[i,1],];
			expcor=cor(a,b,method="pearson");
			plot(a,main=rownames(geneexp[i,]));
			lines(a);
			legtext=paste("Cor=",expcor,sep="");
			legend("topleft",legend=legtext,bty="n");
		}
		dev.off()
	}
}

getgeneexpression<-function(genelist,expdata){
	##this function is used to return the data of rownames of genelist;
	geneexp=expdata[rownames(genelist),];
	geneexp=geneexp[apply(geneexp,1,function(x)!any(is.na(x))),];
	geneexp;
}

changerowname<-function(genelist,expdata){
	##this function is used to change the rownames of expdata with the value of genelist;
	newrownames=vector();
	for(i in 1:nrow(expdata)){
		newrownames[i]=as.character(genelist[rownames(expdata[i,]),1]);
	}
	rownames(expdata)=newrownames;
	expdata;
}

mergeorthologs<-function(datasource1,datasource2,orthologinfo,colx,coly){
	##this program is used to merge two dataframe by orthologous information
	idAs=orthologinfo[,colx];
	idBs=orthologinfo[,coly];
	idAs.data=datasource1[levels(idAs),];
	idAs.data=na.omit(idAs.data);
	idBs.data=datasource2[levels(idBs),];
	idBs.data=na.omit(idBs.data);
	merge1=merge(orthologinfo[,c(colx,coly)],idAs.data,by.x=colx,by.y=0);
	merge2=merge(merge1,idBs.data,by.x=2,by.y=0)
	merge2;
}

