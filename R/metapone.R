metapone <-
  function(pos=NULL, neg=NULL, pa, hmdbCompMZ, pos.adductlist = c("M+H", "M+NH4", "M+Na", "M+ACN+H","M+ACN+Na", "M+2ACN+H", "2M+H", "2M+Na", "2M+ACN+H"), neg.adductlist = c("M-H","M-2H","M-2H+Na","M-2H+K", "M-2H+NH4","M-H2O-H","M-H+Cl", "M+Cl", "M+2Cl"),use.fractional.count=TRUE, match.tol.ppm=5, p.threshold=0.05, n.permu=200, fractional.count.power=0.5, max.match.count=10, use.fgsea = FALSE, use.meta = FALSE)
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
      min.p <- min(unlist(matched[unlist(matched[,3]!=0),3]))
      matched[unlist(matched[,3])==0,3] <- min.p/10 
      
      prod<-unlist(matched[,1]) * unlist(matched[,2])
      uniq.prod<-unique(prod)
      uniq.pval<-rep(0, length(uniq.prod))
      for(this.prod in uniq.prod)
      {
        counts[which(prod == this.prod)]<-sum(prod == this.prod)
        uniq.pval[which(uniq.prod == this.prod)]<-unlist(matched[which(prod == this.prod),3])[1]
      }
      counts[which(counts > max.match.count)]<-max.match.count
      
      counts<-(1/counts)^fractional.count.power
      matched<-cbind(matched, counts)
      
	  #### limit contribution from single feature and single metabolite
	  
	  for(this.prod in uniq.prod)
	  {
		this.pos<-which(prod == this.prod)
		this.total.count<-sum(unlist(matched[this.pos,8]))
		if(this.total.count > 2)  matched[this.pos,8] <- unlist(matched[this.pos,8]) * 2 / this.total.count
	  }
	  
	  all.matched.id<-unlist(unique(matched[,5]))
	  for(matched.id in all.matched.id)
	  {
		this.pos<-which(matched[,5] == matched.id)
		this.total.count<-sum(unlist(matched[this.pos,8]))
		if(this.total.count>1) matched[this.pos,8] <- unlist(matched[this.pos,8]) * sqrt(this.total.count)/this.total.count
	  }
	  
	  ###########################
	  
      all.mapped<-matched[,c(5,8)]
      all.mapped<-all.mapped[which(all.mapped[,1] %in% pa[,3]),]
      sig.mapped<-matched[matched[,3]<=p.threshold,c(5,8)]
      #sig.mapped<-as.data.frame(sig.mapped)
      sig.mapped<-sig.mapped[which(sig.mapped[,1] %in% pa[,3]),]
      
      pathways<-unique(pa[,2])
      
      rec<-matrix(NA, nrow=length(pathways),ncol=8)
      rownames(rec)<-pathways
      colnames(rec)<-c("p_value", "n_significant_metabolites", "n_mapped_metabolites", "n_metabolites", "significant metabolites",
                       "mapped_metabolites", "lfdr", "adjust.p")
      
      rec2<-new("list")    
      
      for(m in seq_len(length(pathways)))
      {
        pathway<-pathways[m]
        metabolites<-pa[which(pa[,2]==pathway),3]
        
        total.draws<-sum(unlist(sig.mapped[,2]))
        white.balls<-sum(unlist(all.mapped[which(all.mapped[,1] %in% metabolites),2]))
        #if(white.balls > 0)
        #{
        found.white.balls<-sum(unlist(sig.mapped[which(sig.mapped[,1] %in% metabolites),2]))
        
        
        mapped.metabolites<-metabolites[which(metabolites %in% all.mapped)]
        sig.metabolites<-metabolites[which(metabolites %in% sig.mapped)]
        
        rec[m,1]<-1
        rec[m,2]<-found.white.balls
        rec[m,3]<-white.balls
        rec[m,4]<-sum(pa[,2] == pathway)
        rec[m,5]<-concate(sig.metabolites)
        rec[m,6]<-concate(mapped.metabolites)
        
        #}
        
        this.sel<-which(matched[,5] %in% pa[pa[,2] == pathway,3] & matched[,3]<=p.threshold)
        if(length(this.sel)>1) rec2[[m]]<-matched[this.sel,]
        if(length(this.sel)==1) rec2[[m]]<-unlist(matched[this.sel,])
        if(length(this.sel)==0) rec2[[m]]<-NA
      }
      
      if(use.fgsea==FALSE)
      {
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
      }else{
        rec_fgsea<-matrix(NA, nrow=length(pathways),ncol=3)
        rownames(rec_fgsea)<-pathways
        colnames(rec_fgsea)<-c("ES", "NES", "nMoreExtreme")
        rec <- cbind(rec, rec_fgsea)
        uniq.meta <- unique(matched[,5])
        meta.imp <- rep(0, length(uniq.meta))
        names(meta.imp) <- uniq.meta
        
        if(use.meta == TRUE){
          ######################## fgsea part ########################
          ############## 1 from the prospective of meta ##############
          for (this.meta in uniq.meta) {
            this.meta <- this.meta[[1]]
            contain.feature <- as.data.frame(matched[which(matched[,5]==this.meta),])
            meta.imp[this.meta] <- sum(-log(unlist(contain.feature[,3]))*unlist(contain.feature[,8]))
          }
          
          pathway.meta <- pa[which(pa[,3] %in% uniq.meta),c(2,3)]
          uniq.pathway <- unique(pathway.meta[,1])
          pathway_fgsea <- new("list")
          for (this.pathway in uniq.pathway) {
            new.list <- new("list")
            new.list <- pathway.meta[which(pathway.meta[,1]==this.pathway),2]
            pathway_fgsea <- c(pathway_fgsea, list(new.list))
          }
          names(pathway_fgsea) <- uniq.pathway
          
          f<-fgseaSimple(pathway_fgsea, sort(meta.imp), nperm=n.permu, maxSize = 500)  # fgsea result
          
          for(m in seq_len(length(pathways))){
            if(pathways[m] %in% f$pathway){
              idx <- which(f$pathway==pathways[m])
              new.p = f$pval[idx]
              rec[m,1]<-as.numeric(new.p)
              rec[m,9]<-as.numeric(f$ES[idx])
              rec[m,10]<-as.numeric(f$NES[idx])
              rec[m,11]<-as.numeric(f$nMoreExtreme[idx])
            }
          }
          
        }else{
          
          ############## 2 from the prospective of feature ##############
          pathway.meta <- pa[which(pa[,3] %in% uniq.meta),c(2,3)]
          uniq.pathway <- unique(pathway.meta[,1])
          pathway_fgsea <- new("list")
          for (this.pathway in uniq.pathway) {
            new.list <- new("list")
            new.list <- pathway.meta[which(pathway.meta[,1]==this.pathway),2]
            new.feature.list <- c()
            
            for (this.meta in new.list){
              if (length(which(matched[,5]==this.meta))!=0){
                mz_time <- unlist(matched[which(matched[,5]==this.meta),1]) * unlist(matched[which(matched[,5]==this.meta),2])
                idx <- which(uniq.prod %in% mz_time)
                new.feature.list <- c(new.feature.list, idx) 
              }
            }
            pathway_fgsea <- c(pathway_fgsea, list(new.feature.list))
          }
          
          names(pathway_fgsea) <- uniq.pathway
          
          uniq.weighted.pval <- rep(0, length(uniq.prod))
          for(this.prod in uniq.prod){
            uniq.weighted.pval[which(uniq.prod == this.prod)]<-(unlist(matched[which(prod == this.prod),3])[1])*(unlist(matched[which(prod == this.prod),8])[1])
          }
          feature.imp <- -log(uniq.weighted.pval)
          names(feature.imp) <- c(1:length(feature.imp))
          #plot(sort(-log(feature.imp*2)))
          
          f<-fgseaSimple(pathway_fgsea, sort(-log(feature.imp*2)), nperm=n.permu, maxSize = 500)
          
          for(m in seq_len(length(pathways))){
            if(pathways[m] %in% f$pathway){
              idx <- which(f$pathway==pathways[m])
              new.p = f$pval[idx]
              rec[m,1]<-as.numeric(new.p)
              rec[m,9]<-as.numeric(f$ES[idx])
              rec[m,10]<-as.numeric(f$NES[idx])
              rec[m,11]<-as.numeric(f$nMoreExtreme[idx])
            }
          }
        }
        
        #stopCluster(cl))$lfdr
      }
      
	    rec[is.na(rec[,1]),1]<-1
	  
      sel<-which(rec[,1]<1)
      this.lfdr<-fdrtool(as.numeric(rec[sel,1]), statistic="pvalue", plot=FALSE, verbose=FALSE)$lfdr
      rec[sel,7] <- this.lfdr
      
      this.BH<-p.adjust(as.numeric(rec[sel,1]), method = 'BH')
      rec[sel,8] <- this.BH
      
      rec<-as.data.frame(rec)
      for(i in c(seq(1,4),7,8)) {rec[,i]<-as.numeric(as.vector(rec[,i]))}
      for(i in seq(5,6)) rec[,i]<-as.vector(rec[,i])
      
      
      names(rec2)<-pathways
      
      to.return <- new("metaponeResult",test.result=rec, mapped.features=rec2)
    }else{
      to.return<-NA
      message("Too few matched metabolites. Cannot continue to calculate pathway significance.")
    }
    return(to.return)
  }





