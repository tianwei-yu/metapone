metapone <-
function(pos=NULL, neg=NULL, pa, hmdbCompMZ, pos.adductlist = c("M+H", "M+NH4", "M+Na", "M+ACN+H","M+ACN+Na", "M+2ACN+H", "2M+H", "2M+Na", "2M+ACN+H"), neg.adductlist = c("M-H","M-2H","M-2H+Na","M-2H+K", "M-2H+NH4","M-H2O-H","M-H+Cl", "M+Cl", "M+2Cl"), use.fractional.count=TRUE, match.tol.ppm=5, p.threshold=0.05, n.permu=200, fractional.count.power=0.5, max.match.count=10)
{
	concate<-function(a)
	{
		b<-a[1]
		if(length(a)>1)
		{
			for(i in 2:length(a)) b<-paste(b, a[i], sep=",")
		}
		return(b)
	}



	#cl <- makeCluster(num.nodes)
	#registerDoParallel(cl)
 
	adductlist<-list(pos.adductlist, neg.adductlist)
	dat<-list(pos, neg)

	comp.mz<-new("list")
	for(i in seq_len(length(adductlist))) comp.mz[[i]]<-hmdbCompMZ[hmdbCompMZ[,3] %in% adductlist[[i]],]

	to.remove<-rep(FALSE, 2)
	for(i in seq_len(length(dat)))
	{
		if(is.null(dat[[i]])) to.remove[i]<-TRUE
	}
	
	if(sum(to.remove)>0)
	{
		dat<-dat[-which(to.remove)]
		comp.mz<-comp.mz[-which(to.remove)]
		adductlist<-adductlist[-which(to.remove)]
	}
	nnn<-1

	for(n in seq_len(length(dat)))
	{
		this.dat<-dat[[n]]
		this.adductlist<-adductlist[[n]]
		this.comp.mz<-comp.mz[[n]]
	
		this.comp.mz[,1]<-as.character(this.comp.mz[,1])
		this.comp.mz[,3]<-as.character(this.comp.mz[,3])
		this.comp.mz[,2]<-as.numeric(this.comp.mz[,2])

		#this.matched<-foreach(nnn=1:length(this.adductlist), .combine=rbind) %dopar%
		
		
		match.fun<-function(ion, comp.mz, this.dat, match.tol.ppm)
		{
			mass.match<-function(x, known.mz, match.tol.ppm=5)
			{
				mass.matched.pos<-rep(0, length(x))
				for(i in seq_len(length(x)))
				{
					this.d<-abs((x[i]-known.mz)/x[i])
					if(sum(!is.na(this.d))>0)
					{					
						if(min(this.d,na.rm=TRUE) < match.tol.ppm/1e6) mass.matched.pos[i]<-1
					}
				}
				return(mass.matched.pos)
			}
			this.kegg<-comp.mz[comp.mz[,3] %in% ion,]
			m<-mass.match(this.dat[,1], this.kegg[,2], match.tol.ppm=match.tol.ppm)
		
			this.matched<-NULL
			if(sum(m)>0)
			{
				this.matched<-matrix(0, ncol=2+ncol(this.dat)+ncol(comp.mz), nrow=1)
		
				for(i in which(m==1))
				{
					d<-abs(this.dat[i,1]-this.kegg[,2])/this.dat[i,1]
					sel<-which(d <= (match.tol.ppm*1e-6))
					for(j in seq_len(length(sel))) this.matched<-rbind(this.matched, c(i, sel[j],this.dat[i,],this.kegg[sel[j],]))
				}
			
				this.matched<-this.matched[-1,]
				if(is.null(nrow(this.matched))) this.matched<-matrix(this.matched, ncol=2)
				this.matched<-this.matched[,seq(-1,-2)]
			}
			this.matched
		}
		z<-bplapply(this.adductlist, match.fun, comp.mz=this.comp.mz, this.dat=this.dat, match.tol.ppm=match.tol.ppm)
		this.matched<-z[[1]]
		if(length(z)>1)
		{
			for(m in 2:length(z))
			{
				this.matched<-rbind(this.matched, z[[m]])
			}
		}
		
		if(n == 1)
		{
			matched<-this.matched
		}else{
			matched<-rbind(matched, this.matched)
		}
	}

	if(!is.null(nrow(matched)))
	{
		counts<-rep(1, nrow(matched))

		prod<-unlist(matched[,1]) * unlist(matched[,2])
		uniq.prod<-unique(prod)
		uniq.pval<-rep(0, length(uniq.prod))
		for(this.prod in uniq.prod)
		{
			counts[which(prod == this.prod)]<-1/sum(prod == this.prod)
			uniq.pval[which(uniq.prod == this.prod)]<-unlist(matched[which(prod == this.prod),3])[1]
		}
		counts[which(counts > max.match.count)]<-max.match.count
			
		counts<-counts^fractional.count.power
		matched<-cbind(matched, counts)

		all.mapped<-matched[,c(5,8)]
		all.mapped<-all.mapped[which(all.mapped[,1] %in% pa[,3]),]
		sig.mapped<-matched[matched[,3]<=p.threshold,c(5,8)]
		sig.mapped<-as.data.frame(sig.mapped)
		sig.mapped<-sig.mapped[which(sig.mapped[,1] %in% pa[,3]),]
		
		pathways<-unique(pa[,2])

		rec<-matrix(NA, nrow=length(pathways),ncol=6)
		rownames(rec)<-pathways
		colnames(rec)<-c("p_value", "n_significant metabolites", "n_mapped_metabolites", "n_metabolites", "significant metabolites", "mapped_metabolites")

		rec2<-new("list")    

		for(m in seq_len(length(pathways)))
		{
			pathway<-pathways[m]
			metabolites<-pa[which(pa[,2]==pathway),3]
			
			total.draws<-sum(unlist(sig.mapped[,2]))
			white.balls<-sum(unlist(all.mapped[which(all.mapped[,1] %in% metabolites),2]))
			if(white.balls > 0)
			{
				found.white.balls<-sum(unlist(sig.mapped[which(sig.mapped[,1] %in% metabolites),2]))
			

				mapped.metabolites<-metabolites[which(metabolites %in% all.mapped)]
				sig.metabolites<-metabolites[which(metabolites %in% sig.mapped)]
				
				rec[m,1]<-1
				rec[m,2]<-found.white.balls
				rec[m,3]<-white.balls
				rec[m,4]<-sum(pa[,2] == pathway)
				rec[m,5]<-concate(sig.metabolites)
				rec[m,6]<-concate(mapped.metabolites)
				
			}
			
			this.sel<-which(matched[,5] %in% pa[pa[,2] == pathway,3] & matched[,3]<=p.threshold)
			if(length(this.sel)>1) rec2[[m]]<-matched[this.sel,]
			if(length(this.sel)==1) rec2[[m]]<-unlist(matched[this.sel,])
			if(length(this.sel)==0) rec2[[m]]<-NA		
		}
		
		## permutation test

		#rec.permu<-foreach(nnn=1:n.permu, .combine=cbind) %dopar%
		
		perm.fun<-function(fake.counter, all.mapped, pathways, pa, that.matched, uniq.pval, uniq.prod, p.threshold, prod)
		{
			fakefake<-fake.counter
			r<-rep(NA, length(pathways))
			
			that.uniq.pval<-sample(uniq.pval, length(uniq.pval),replace=FALSE)
			
			for(that.prod in uniq.prod)
			{
				that.matched[which(prod == that.prod),3]<- that.uniq.pval[which(uniq.prod == that.prod)]
			}
			that.all.mapped<-that.matched[,c(5,8)]
			that.all.mapped<-that.all.mapped[which(that.all.mapped[,1] %in% pa[,3]),]
			that.sig.mapped<-that.matched[that.matched[,3]<=p.threshold,c(5,8)]
			that.sig.mapped<-as.data.frame(that.sig.mapped)
			that.sig.mapped<-that.sig.mapped[which(that.sig.mapped[,1] %in% pa[,3]),]
			
			for(m in seq_len(length(pathways)))
			{
				pathway<-pathways[m]
				metabolites<-pa[which(pa[,2]==pathway),3]
				white.balls<-sum(unlist(all.mapped[which(all.mapped[,1] %in% metabolites),2]))
				if(white.balls > 0)
				{
					found.white.balls<-sum(unlist(that.sig.mapped[which(that.sig.mapped[,1] %in% metabolites),2]))
				}
				r[m]<-found.white.balls
			}
			r
		}
		z<-bplapply(seq_len(n.permu), perm.fun, all.mapped=all.mapped, pathways=pathways, pa=pa, that.matched=matched, uniq.pval=uniq.pval, uniq.prod=uniq.prod, p.threshold=p.threshold, prod=prod)
		
		rec.permu<-matrix(0,ncol=n.permu, nrow=length(z[[1]]))
		for(i in seq_len(n.permu)) rec.permu[,i]<-z[[i]]
		
		for(m in seq_len(nrow(rec)))
		{
			if(!is.na(rec[m,1]))
			{
				this.r<-as.numeric(rec[m,2])
				new.p<-sum(rec.permu[m,]>=this.r)/ncol(rec.permu)
				rec[m,1]<-new.p
			}
		}

	
	#stopCluster(cl)
	
	rec<-as.data.frame(rec)
	for(i in seq(1,4)) rec[,i]<-as.numeric(as.vector(rec[,i]))
	for(i in seq(5,6)) rec[,i]<-as.vector(rec[,i])

	names(rec2)<-pathways

	to.return <- new("metaponeResult",test.result=rec, mapped.features=rec2)
	
	}else{
		to.return<-NA
		message("Too few matched metabolites. Cannot continue to calculate pathway significance.")
	}
	return(to.return)
}