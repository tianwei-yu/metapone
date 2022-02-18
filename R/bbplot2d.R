globalVariables(c("p_value","logp", "lfdr", "n_metabolites", "n_significant_metabolites"))
bbplot2d <-
  function(res, p_thres = 0.05, sig_metab_thres=1)
  {
    res_05 <- res[res$p_value < p_thres & res[,2]>= sig_metab_thres,]
    res_05 <- res_05[complete.cases(res_05[,c(1:4,7)]),]
    res_05$p_value[which(res_05$p_value==0)] <- min(as.numeric(res_05$p_value[res_05$p_value != 0]),na.rm=T)/10
	res_05$logp <- -1*log10(as.numeric(res_05$p_value))
    
    p = ggplot(res_05,aes_string(x = 'logp',y = 'n_significant_metabolites'))+
      geom_point(aes_string(size='n_significant_metabolites', color = 'lfdr')) +
      scale_color_gradient(low = 'green', high = 'red') +
      geom_text_repel(aes(`logp`, `n_significant_metabolites`, label = rownames(res_05)),  size=4) 
      labs(color=expression('lfdr'),size='n_significant_metabolites',x="-log10(Pvalue)",y="Number of significant metabolites",
           title="Overview of Enriched Pathway")+ theme_bw()
	p
  }
  