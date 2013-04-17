drawmatrix<-function(datasource,is.norm=FALSE,margins=c(5,5),...){
	rowN=nrow(datasource);
	colN=ncol(datasource);
	widthN=round(colN*1.1);
	par(mar=c(margins[1L],1,1,margins[2L]));
	if(is.norm){
		datasource.mean=apply(datasource,1,mean);
		datasource.sd=apply(datasource,1,sd);
		datasource=(datasource-datasource.mean)/datasource.sd;
	}
	plot(0,type="n",ylab="",xlab="",axes=FALSE,ylim=c(0,rowN),xlim=c(0,widthN));
	datamin=min(datasource);
	datamax=max(datasource);
	datawidth=datamax-datamin;
	cexRow = 0.2 + 1/log10(rowN);
	cexCol = 0.2 + 1/log10(colN);
	image(1:ncol(datasource),1:nrow(datasource),t(as.matrix(datasource)),axes=FALSE,xlab="",ylab="",...);
	axis(1, 1:colN, labels = colnames(datasource), las = 2, line = -0.5, tick = 0,cex.axis = cexCol);
	axis(4, 1:rowN, labels = rownames(datasource), las = 2, line = -0.5, tick = 0,cex.axis = cexRow)
}