draw_coverage<-function(Chrom_File,Exp_File,style=1){
	require(ggplot2);
	require(scales);
	require(vcd);
	chrom_length=read.table(Chrom_File,sep="\t",row.name=1);
	startpos=1;
	chrom=data.frame();
	for(i in 1:nrow(chrom_length)){
		chrom[i,1]=startpos;
		startpos=startpos+chrom_length[i,1];
		chrom[i,2]=startpos;
		startpos=startpos+1;
	}
	row.names(chrom)=rownames(chrom_length);
	colnames(chrom)=c("start","end");
	chrom$name=rownames(chrom_length);
	exp_info=read.table(Exp_File,sep="\t");
	exp_converted=data.frame();
	for(i in 1:nrow(exp_info)){
		exp_converted[i,1]=exp_info[i,1];
		exp_converted[i,2]=exp_info[i,2]+chrom[as.character(exp_info[i,1]),1];
		exp_converted[i,3]=exp_info[i,3];
	}
	exp_data=exp_converted;
	colnames(exp_data)=c("chrom","site","exp");
	yrng=range(exp_data$exp);
	themes<-theme (legend.position="non",axis.line=element_line(colour="black",size=0.2),axis.ticks.x=element_blank(),axis.text.x=element_blank(), panel.grid.minor=element_blank(),panel.background=element_blank());
	if(style==1){
		q<-qplot(site,exp,data=exp_data,geom="line",xlab="",ylab="",color=factor(chrom));
		q<-q+themes;
	}
	else{
		q<-qplot(site,exp,data=exp_data,geom="line",xlab="",ylab="");
		q<-q+ geom_rect(aes(NULL,NULL,xmin=start,xmax=end,fill=name),ymin=yrng[1],ymax=yrng[2],data=chrom)+scale_fill_manual(values=alpha(gg_color_hue(14),0.5)) + geom_text(aes(x=start,y= yrng[1],label=rownames(chrom)),data = chrom,size = 3, hjust = 0, vjust = 0)+themes
	}
	q;
}

draw_region<-function(Chrom_List,Region_File){
	require(ggplot2);
	require(scales);
	require(vcd);
	chrom_length=Chrom_List;
	startpos=1;
	Chrom=data.frame();
	for(i in 1:nrow(chrom_length)){
		Chrom[i,1]=startpos;
		startpos=startpos+chrom_length[i,1];
		Chrom[i,2]=startpos;
		startpos=startpos+1;
	}
	row.names(Chrom)=rownames(chrom_length);
	colnames(Chrom)=c("start","end");
	Chrom$name=rownames(Chrom_List)
	region_info<-try(read.table(Region_File,sep="\t"));
	q=NULL;
	if(class(region_info)=="try-error"){
		xrng=range(Chrom[,1:2]);
		df<-data.frame();
		themes<-theme (legend.position="non",axis.ticks=element_blank(), axis.line=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),text=element_blank());
		q<-ggplot(df)+geom_rect(aes(NULL,NULL,xmin=start,xmax=end,fill=name),ymin=0.9,ymax=1,data=Chrom,xlab="",ylab="")+ylim(0,1)+scale_x_continuous(limits=c(xrng[1],xrng[2]),expand=c(0,0));
		q<-q+themes;
	}
	else{
		region_converted=data.frame();
		for(i in 1:nrow(region_info)){
			region_converted[i,1]= region_info[i,1];
			region_converted[i,2]= region_info[i,2]+ Chrom[as.character(region_info[i,1]),1];
			region_converted[i,3]= region_info[i,3]+ Chrom[as.character(region_info[i,1]),1];
			region_converted[i,4]= region_info[i,4]*0.6+0.3;
		}
		#colnames(region_converted)=c("chrom","start","end");
		colnames(region_converted)=c("chrom","start","end","value");
		xrng=range(Chrom[,1:2]);
		df<-data.frame();
		themes<-theme (legend.position="non",axis.ticks=element_blank(), axis.line=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),text=element_blank());
		q<-ggplot(df)+geom_rect(aes(NULL,NULL,xmin=start,xmax=end,ymin=0.1,ymax=value),data=region_converted,fill="red",xlab="",ylab="")+ylim(0,1)+scale_x_continuous(limits=c(xrng[1],xrng[2]),expand=c(0,0));
		q<-q+geom_rect(aes(NULL,NULL,xmin=start,xmax=end,fill=name,ymin=0,ymax=0.1),data=Chrom,xlab="",ylab="");
		q<-q+themes;
	}
	q;
}

draw_SV_All<-function(SV_File_List,Strain_Names,Chrom_List){
	require(grid);
	chrom_length=Chrom_List;
	startpos=1;
	Chrom=data.frame();
	for(i in 1:nrow(chrom_length)){
		Chrom[i,1]=startpos;
		startpos=startpos+chrom_length[i,1];
		Chrom[i,2]=startpos;
		startpos=startpos+1;
	}
	row.names(Chrom)=rownames(chrom_length);
	colnames(Chrom)=c("start","end");
	Chrom$name=rownames(Chrom_List)
	vplayout <- function(x, y) {viewport(layout.pos.row = x, layout.pos.col = y)};
	grid.newpage();
	row_number=nrow(SV_File_List)+2;
	pushViewport(viewport(layout = grid.layout(row_number,12)));
	for(i in 1:nrow(SV_File_List)){
		p<-draw_region(Chrom_List,as.character(SV_File_List[i,1]))
		print(p,vp=vplayout(i+1,2:12))
		grid.text(as.character(Strain_Names[i,1]), vp = vplayout(i+1,1),gp=gpar(fontsize=10))
	}
	xrng=range(Chrom[,1:2]);
	df<-data.frame();
	themes<-theme (legend.position="non",axis.ticks=element_blank(), axis.line=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),text=element_blank());
	k<-ggplot(df)+geom_text(aes(x = start, y = 0.8, label = name), data = Chrom, size = 3, hjust = 0, vjust = 0)+ylim(0,1)+scale_x_continuous(limits=c(xrng[1],xrng[2]),expand=c(0,0));
	k<-k+themes;
	print(k,vp=vplayout(row_number,2:12));
}

gg_color_hue <- function(n) {
	hues = seq(15, 375, length=n+1)
	hcl(h=hues, l=65, c=100)[1:n]
}