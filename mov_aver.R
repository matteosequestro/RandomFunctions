mov_aver <- function(x, n = 5){
    stepn = floor(n/2)
    set = c()
    for(ii in stepn : (length(x)-stepn)) {
        set = append(set, mean(x[(ii-stepn):(ii+stepn)]))
    }
    return(set)
}