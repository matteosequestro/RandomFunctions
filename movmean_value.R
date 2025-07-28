






movmean_value = function(x,y, winlen = 2000, confint = .95, nBoot = 100) {
    sorted_y = sort(y)

    
    t_range = rep(999, length(sorted_y) / winlen)
    t_values = t_range
    t_lower = t_range
    t_upper = t_range
    
    
    for(kk in (winlen/2+1) : (length(sorted_y) - (winlen/2+1))) {
        
        lower = kk - winlen/2
        upper = kk + winlen/2
        t_y = sorted_y[lower : upper]
        t_range[kk] = mean(t_y)

        meanci = Hmisc::smean.cl.boot(x[y %in% t_y], conf.int = confint, B = nBoot )
        t_values[kk] = meanci[1]
        t_lower[kk] = meanci[2]
        t_upper[kk] = meanci[3]
        # meanci = mean_se(x[y%in% t_y])
        # # # meanci = mean_se(tid$AisChosen[tid$probA %in% t_uncertainty])
        # # 
        # t_values[kk] = meanci$y
        # t_lower[kk] = meanci$ymin
        # t_upper[kk] = meanci$ymax
        #     print(meanci$y)
    }
    
    export = data.frame(t_range, t_values, t_lower, t_upper, 1: length(t_range))    
    export = export[!is.na(export$t_values) & export$t_values < 900, ]
    colnames(export) = c("winmean", "mean", "lower", "upper", "sample")
    
    
    return(export)
}