DrawCoverage<-function(FileName,suffix)
{
	for(i in 1:length(suffix)){
		file=paste(FileName,suffix[i],sep=".");
		value=read.table(file,sep="\t");
		xx=value[,2];
		yy=value[,3];
		n=nrow(value);
		plot (yy ~ xx, type="n", axes=1,main=suffix[i],xlab="",ylab="");
		polygon(c(xx[1], xx, xx[n]), c(min(yy), yy, min(yy)),col="red", border=NA);
	}
}