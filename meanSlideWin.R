
# Adapted from Hardwick et al. 2019
meanSlideWin = function(x, y, winLen = 100, timeMax = max(x), timeMin = 1, confint = .95, nBoot = 1000, wantPlot = 0) {
    
    # x           = vector of the variable to slide along
    # y           = vector of the variable to be averaged
    # timeMax     = max value of y to be covered by the window
    # timeMin     = min value of y to be covered by the window
    # nboot       = number of iterations for the bootstrapped confidence interval
    # wantPlot    = if plot it will plot the result
    
    
    # Slide the windows
    ts = c()
    for (tt in timeMin : timeMax){
        # Set the center of the window at 1
        lower   = tt - (winLen/2)               # Lower half of the window
        upper   = tt + (winLen/2)               # Upper half
        points  = which(x <= upper & x >= lower)    # points of x contained in the window of y
        sample  = y[points]                    
        
        # if the window contains points, average them and bootrstrap CIs
        if (length(sample) > 1) {
            average = mean(sample)
            ci = Hmisc::smean.cl.boot(sample, conf.int = confint, B = nBoot )
            ts = rbind(ts, c(ci[1], ci[2], ci[3]))
        } else {
            ts = rbind(ts, c(NaN, NaN, NaN))
        }

    }
    
    ts = as.data.frame(ts)
    colnames(ts) = c("mean", "lower.ci", "upper.ci")
    ts$time = 1 : nrow(ts)
    
    
    return(ts)
    

}



