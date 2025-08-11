nice_brms_table = function(fit, version = "short") {
    
    hyp = hypothesis(fit, 
                     paste(rownames(fixef(fit)), "= 0", sep = " ")) #[-1]
    
    evis = 1/hyp$hypothesis$Evid.Ratio # bf_10
    
    post = posterior_samples(fit, "^b")
    
    if(version == "long") {
        tab = cbind(fixef(fit),
                    apply(post, 2, function(x) quantile(x, .05)),
                    apply(post, 2, function(x) quantile(x, .95)),
                    apply(post, 2, function(x) sum(x>0)/length(x)),
                    apply(post, 2, function(x) sum(x<0)/length(x)),
                    apply(post, 2, function(x) max(mean(x > 0), mean(x < 0))),                
                    hyp$hypothesis$Evid.Ratio,
                    evis)
        colnames(tab) = c(colnames(tab)[1:4], "Q5.0", "Q95.0", "pp>0", "pp<0", "ppd", "BF_01", "BF_10")
    } else {
        tab = cbind(fixef(fit),
                    apply(post, 2, function(x) max(mean(x > 0), mean(x < 0))),                
                    evis)
        colnames(tab) = c(colnames(tab)[1:4], "ppd", "BF_10")
        tab = tab[,-2]
    }
    
    tab = round(tab, 3)
    
    
    #posterior probability of direction as in  (Makowski et al., 2019).
    return(tab)
    
}

