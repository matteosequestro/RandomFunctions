
% this require the check_clusters_lme function
function [modelout, cluster_mass, obs_clusters_sum] = run_clust_perm_lme(data, formula, variablesToPermute, idvar, cfg)
% ----------------------------------------------------------------------------------------------------------------------------------------------
% Performs cluster based correction on parameter estimates from linear mixed
% effect model on a time-series
%
% INPUTS:
% data                  : full dataset having at least a column for ID, 
%                           one column for each predictor in your model, 
%                           and one column with the time series in a cell. 
%                           Each row is a trial.
% formula               : the model formula in Wilkinson notation (i.e., R-like
%                           notation, e.g., 'y ~ x1 * x2 + (1 | id)'. 
%                           It must be a string.
% VariablesToPermute    : the variables to permute in each iteration to
%                           compute the non parametric (monte carlo) 
%                           p-value. Ideally they would be each variable 
%                           in your model. The variable should be the 
%                           column name in the data table in a character 
%                           array (e.g., ["condition1", "condition2"].
% idvar                 : the variable containing the id of participants in
%                           the data table. It should be the name of the
%                           column in the table given as a character (e.g.,
%                           "id").
% cfg.niter             : number of permutations.
% cfg.wantplot_perm     : 1 if you want to plot time-series for each
%                           parameter with significant (blue) and
%                           non-significant (gray) clusters as shaded
%                           reptangles. 0 if you don't want such a plot.
% cfg.perm_alpha        : permutation alpha. Montecarlo p-values below such
%                           a value will be considered as significant.
% cfg.want_parallel_fit : 1 will fi the models across different time-points
%                           in parallel (much faster but computationally
%                           more expensive); 0 will do one point at the
%                           time.
%
% OUTPUTS:
% modelout              : ouput of the model on true dataset (i.e., not 
%                           permuted. You will get a "pars" field
%                           (parameter estimate, lower and upper confidence 
%                           interval, tvalues), a "criteria" field with
%                           log-likelihood of the model, as well as AIC,
%                           BIC and deviance criteria (for each sample 
%                           point; and an "rsqrd" field with ordinary and
%                           adjusted R-squared
% obs_clusters_sum      : output of the cluster based permutation. For each
%                           coefficient in the model and each cluster 
%                           detected on the true dataset (i.e., not     
%                           permutedyou), get a monte-carlo
%                           p-value, first sample in the cluster, length
%                           of the cluster and t-mass.
% ----------------------------------------------------------------------------------------------------------------------------------------------

cfg.wantplot_fit = 0;

% Findclusters in the observed data

% [modelout, cluster_beta_exp, cluster_t_exp, cluster_sum]
[modelout, ~ ,obs_clusters_t, obs_clusters_sum] = check_clusters_lme(data, formula, cfg);
npar = height(modelout.pars.estimates); % number of parameters (or coefficients)
tslen = width(modelout.pars.estimates); % length of the time-series

% ID list
ids = unique(data.(idvar));

% Extract the dependent variable from formula
y = formula(1:find(formula == '~')-1);
y = erase(y, ' ');

%% Permutations 
cluster_mass = cell(size(obs_clusters_t)); % you will have the cluster masses for each permutation here
cluster_mass(:,1) = obs_clusters_t(:,1);
for rep = 1:cfg.niter

    % Display the current iteration, one 10th at the time
    if ismember(rep, (0:cfg.niter/10:cfg.niter)); disp([num2str(rep), '/', num2str(cfg.niter)]); end 

    % Permute variables within subjects
    perm_data = data;
    for ii = 1 : length(ids)
        thisid_ind = strcmp(data.(idvar), ids(ii));
        tperm = data(thisid_ind, :);
        for vv = 1 : length(variablesToPermute)
            tperm(:,variablesToPermute(vv)) = tperm(randperm(height(tperm)),variablesToPermute(vv));
            perm_data(thisid_ind,:) = tperm;
        end
    end


    % Find clusters in the permuted data
    [~, ~, per_clusters_t] =check_clusters_lme(perm_data, formula, cfg);

    % Compute t-masses for each predictor in the "permutation model"
    for pred = 1 : height(per_clusters_t)
        t_masses_per = cellfun(@(x) sum(x(:,1)), per_clusters_t{pred,2});
        if isnan(t_masses_per); t_masses_per =0; end % if you don't find clusters add a 0
        
        try
            cluster_mass{cellfun(@(x) strcmp(x, per_clusters_t{pred,1}), cluster_mass(:,1)),2} = ...
                [cluster_mass{cellfun(@(x) strcmp(x, per_clusters_t{pred,1}), cluster_mass(:,1)),2}, ...
                t_masses_per];
        catch
            cluster_mass(cellfun(@(x) strcmp(x, per_clusters_t{pred,1}), cluster_mass(:,1)),2) = ...
                {t_masses_per};
        end
    end

end


%% Compute  monte-carlo p values
for pp = 1 : height(obs_clusters_t)
    tpred = obs_clusters_t{pp, 2};
    tperm_pred = cluster_mass{pp,2};
    obsmass = cellfun(@(x) sum(abs(x(:,1))),  tpred); % masses for each found cluster for this predictor
    obs_clusters_sum(pp).pval = arrayfun(@(x) mean(abs(tperm_pred) > x), obsmass);
    % tps = cell(1, length(tpred));
    %
    % for mm = 1:length(tpred)
    %     obsmass = cellfun(@(x) sum(abs(x(:,1))),  tpred); % masses for each found cluster for this predictor
    %     tps{mm} = arrayfun(@(x) mean(abs(tperm_pred) > x), obsmass);
    % 
    %     % obsmass = sum(abs(tpred{mm}(:,1)));
    %     % mcp = mean(abs(tperm_pred) > obsmass); %compute the p as the percentage of |permuted t masses| greater than the observed |tmass|
    %     % tps{mm} = mcp;
    % end
    % obs_clusters_sum(pp).pval = tps;
end


%% Plots
if cfg.wantplot_perm

    time = 1:tslen;
    for par = 1 : npar
        figure
        patch([time, fliplr(time)], [modelout.pars.lower_ci(par,:), fliplr(modelout.pars.upper_ci(par,:))], ...
            [0.5 0.5 0.5],'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off');
        hold on
        plot(modelout.pars.estimates(par,:), 'Color', [0 0 0], 'LineWidth', 2)

        xlim([1, max(time)]);
        ylabel('parameter estimate (a.u.)')
        xlabel('sample')
        title(modelout.pars.namesout{par})

        yline(0, '--', 'HandleVisibility', 'off')


        for this_cluster = 1 : length(obs_clusters_sum(par).first)
            hold on
            clFirst = obs_clusters_sum(par).first(this_cluster);
            clLength = obs_clusters_sum(par).length(this_cluster);
            ylims = get(gca, 'Ylim');

            if obs_clusters_sum(par).pval(this_cluster) >= cfg.perm_alpha
                rectColor = [.5 .5 .5];
            else
                rectColor =[0.3569 0.6078 0.8353];
            end
            rectangle('Position', [clFirst, ylims(1), clLength, sum(abs(ylims))], ...
                'FaceColor', [rectColor .10], 'EdgeColor',[rectColor .1], 'FaceAlpha', .2)
        end

        xlim([1,    tslen])
        ax = gca; % Get current axis handle
        ax.XTickLabel = ax.XTick / cfg.fs; % Divide by 10
        xlabel("time (s)")
        box on
    end


end


end




