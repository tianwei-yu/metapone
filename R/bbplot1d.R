globalVariables(c("p_value","logp", "lfdr", "n_metabolites", "n_significant_metabolites"))
bbplot1d <- 
  function(res, p_thres = 0.05, sig_metab_thres=1){
    res_05 <- res[res$p_value < p_thres & res[,2]>= sig_metab_thres,]
    res_05 <- res_05[complete.cases(res_05[,c(1:4,7)]),]
    min_p = 0
    if (sum(res_05$p_value!=0)!=0){min_p <- min(as.numeric(res_05$p_value[res_05$p_value!= 0]),na.rm=T)/10}
    res_05$p_value[which(res_05$p_value==0)] <- max(10^(-8), min_p)
    res_05$logp <- -1*log10(as.numeric(res_05$p_value))
    #idx <- order(res_05$p_value)
    res_05 <- res_05[sort(res_05$logp, index.return = TRUE)$ix,]
    res_05$name <- rownames(res_05)
    p=ggplot(res_05,aes_string(x='name',y='logp'))+
      geom_point(aes_string(color='lfdr',size='n_significant_metabolites'),alpha=0.5)+
      coord_flip()+
      scale_color_gradient(low = "green",high = "red")+ 
      labs(color=expression(lfdr),size="N_sig_metab",x = " ", y = "-log10(Pvalue)",title="Overview of Enriched Pathway")+ 
      theme_bw()
    p
  }
  