numbercount<-function(datasource,number){
	colnum=ncol(datasource);
	rownum=nrow(datasource);
	rescount=vector();
	for(i in 1:colnum){
		count=0;
		for(j in 1:rownum){
			if(datasource[j,i]==number){
				count=count+1;
			}
		}
		rescount[i]=count;
	}
	rescount;
}