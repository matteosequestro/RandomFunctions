boot_corr = function(x, y, method = "pearson", n_boot = 5000, seed = 123) {
    set.seed(seed)
    
    # Number of participants
    n = length(x)
    
    # Function to compute correlation on resampled indices
    boot_cor_single = function(indices) {
        cor(x[indices], y[indices], method = method)
    }
    
    # Run bootstrap
    boot_results = replicate(n_boot, boot_cor_single(sample(1:n, replace = TRUE)))
    
    # Compute bootstrap statistics
    ci = quantile(boot_results, probs = c(0.025, 0.975))
    mean_cor = mean(boot_results)
    
    # Return results as a list
    list(mean_cor = mean_cor, ci_lower = ci[1], ci_upper = ci[2], boot_distribution = boot_results)
}
