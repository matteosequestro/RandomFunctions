
function [cluster_full_exp, cluster_sum] = find_clusters(full_ind_estimates, parnames,cfg)


% Number of parameters
npars = size(full_ind_estimates, 3);

%% Find clusters
%(this could probably be done in a easier way but it's the same thing as
% what i was doing in the other function and non c'ho voglia to find
% another solution rn)
cluster_n       = 1;
for pred = 1 : npars
    clear cluster_full
    beta_series         = squeeze(full_ind_estimates(:,:, pred));                    % unpack estimaets for this predictor
    null_reference      = zeros(size(beta_series ));                              % vector of zeros as reference (you compare the parameter as the null effect, i.e., zero)
    [~, ~, ~, stats]    = ttest(beta_series, null_reference , cfg.perm_alpha);  % One sample t values against 0
    tresh               = tinv(cfg.perm_alpha/2, stats.df(1));                             % take the treshold (quantile at that alpha value for those degrees of freedom)
    tt                  = [(stats.tstat)', (1:length(stats.tstat))'];                         % reshape the ts
    clusters_serie      = tt(tt(:,1) < tresh | tt(:,1) > -tresh, :);              % find the values over the treshold

    % Extract clusters
    if ~isempty(clusters_serie)                                              %if you find the cluster...
        % cluster_struc   = {};
        this_cluster    = [];

        for sample = 1:length(clusters_serie(:,2))                          %then separate the clusters based on temporal contiguity in the following way:
            if sample == 1                                                 %if it's the first sample then put it on the list
                this_cluster = [this_cluster; clusters_serie(sample,:)];
            elseif clusters_serie(sample,2) == clusters_serie(sample-1, 2) + 1 %if the next sample is at the successive time point, then put it on the same list (i.e. the list for that cluster)
                this_cluster = [this_cluster; clusters_serie(sample,:)];
            elseif clusters_serie(sample,2) ~= clusters_serie(sample-1, 2) + 1 %if the next sample is not at the successive time point then form a new list
                cluster_full{cluster_n} = this_cluster;
                this_cluster = [[]; clusters_serie(sample,:)];
                cluster_n = cluster_n + 1;
            end
        end
        cluster_full{cluster_n} = this_cluster;
    else
        cluster_full = {};
    end

    % Delete small clusters
    cluster_full(cellfun(@(x) length(x) < cfg.minlength, cluster_full)) = [];

    % save info for each cluster
    % Save info for each cluster
    cluster_sum(pred,:).pred        = parnames{pred};
    cluster_sum(pred,:).length      = cellfun(@(x) size(x, 1), cluster_full);
    cluster_sum(pred,:).first       = cellfun(@(x) x(1, 2), cluster_full);
    cluster_sum(pred,:).last        = cellfun(@(x) x(end, 2), cluster_full);
    cluster_sum(pred,:).mass        = cellfun(@(x) sum(x(:,1)), cluster_full);
    cluster_sum(pred,:).abs_mass    = cellfun(@(x) sum(abs(x(:,1))), cluster_full);


    if isempty(cluster_full); cluster_full = {NaN}; end
    cluster_full_exp(pred, :) = {parnames{pred}, cluster_full};
end

end