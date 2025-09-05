
% this require the check_clusters_lme function
function [modelout, cluster_mass, obs_clusters_sum, clust_par_means] = run_clust_perm_lme2(data, formula, cfg)
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

% Fit the model and find clusters in the observed data
[modelout, obs_clusters_t, obs_clusters_sum] = time_series_lme(data, formula, cfg);
npar        = height(modelout.pars.estimates); % number of parameters (or coefficients)
tslen       = width(modelout.pars.estimates); % length of the time-series
parnames    = modelout.pars.namesout;



% Unpack individual parameter estimates
pars_series = modelout.full_ind_estimates;

% Extract the dependent variable from formula
y = formula(1:find(formula == '~')-1);
y = erase(y, ' ');

%% Permutations
cluster_mass = cell(size(obs_clusters_t)); % you will have the cluster masses for each permutation here
cluster_mass(:,1) = obs_clusters_t(:,1);
for rep = 1 : cfg.niter

    %%% Series that will be permuted
    perm_serie = pars_series;

    %%% Randomize for which participants to switch labels and then switch labels
    rswitches                   = logical(randi([0,1], size(perm_serie, 1), 1));
    perm_serie(rswitches,:,:)   = 0;

    % Find the clusters in the permuted set
    [clusters_per, ~] = find_clusters(perm_serie, parnames, cfg);

    % Compute t-masses for each predictor in the "permutation model"
    for pred = 1 : height(clusters_per)
        t_masses_per = cellfun(@(x) sum(x(:,1)), clusters_per{pred,2});
        if isnan(t_masses_per); t_masses_per = 0; end % if you don't find clusters add a 0

        try
            cluster_mass{cellfun(@(x) strcmp(x, clusters_per{pred,1}), cluster_mass(:,1)),2} = ...
                [cluster_mass{cellfun(@(x) strcmp(x, clusters_per{pred,1}), cluster_mass(:,1)),2}, ...
                t_masses_per];
        catch
            cluster_mass(cellfun(@(x) strcmp(x, clusters_per{pred,1}), cluster_mass(:,1)),2) = ...
                {t_masses_per};
        end
    end

end


%% Compute  monte-carlo p values
for pp = 1 : height(obs_clusters_t)
    tpred                       = obs_clusters_t{pp, 2};
    tperm_pred                  = cluster_mass{pp,2};
    obsmass                     = cellfun(@(x) sum(abs(x(:,1))),  tpred); % masses for each found cluster for this predictor
    obs_clusters_sum(pp).pval   = arrayfun(@(x) mean(abs(tperm_pred) > x), obsmass);
end


%% Plots
if cfg.wantplot_perm

    % Remove the underscores from the parameter names otherwise it's ugly
    parnames_to_plot = cellfun(@(x) replace(x, '_', ' '), parnames, 'UniformOutput', false);

    % Time vector
    time = 1 : tslen;

    % Number of columns and rows for the layout
    ncols_tiles = floor(npar / 2);
    nrows_tiles = npar - ncols_tiles;

    % Prepare the layout (keep it tight)
    tld             = tiledlayout(nrows_tiles, ncols_tiles);
    tld.TileSpacing = 'tight';
    tld.Padding     = 'tight';

    % Common labels
    ylabel(tld, 'parameter estimate (a.u.)');
    xlabel(tld, 'time (s)');

    % Plot each parameter
    for par = 1 : npar
        nexttile

        % Compute stats to plot
        this_parameter      = pars_series(:,:,par);
        grand_average       = mean(this_parameter, 1);
        standard_error      = std(this_parameter, [], 1)  ./ sqrt(height(this_parameter));  % standard error for the confidence interval
        ts                  = tinv([0.025  0.975], height(this_parameter)-1);               % t-score for the confidence interval
        lower_bound         = grand_average + ts(1) .* standard_error;                      % lower bound of confidence interval
        upper_bound         = grand_average + ts(2) .* standard_error;                      % upper bound of confidence interval

        % The shade (intervals)
        patch([time, fliplr(time)], [lower_bound, fliplr(upper_bound)], ...
            [0.5 0.5 0.5],'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off');
        hold on

        % Grand average
        plot(grand_average, 'Color', [0 0 0], 'LineWidth', 2)
        if isfield(cfg, "stim_onset_time")
            xline(cfg.stim_onset_time, '--')
        end

        if isfield(cfg, "plot_group_pars")
            if cfg.plot_group_pars
                plot(modelout.pars.estimates(par,:), 'r--', 'LineWidth', 1.5);
            end
        end

        % Put a title with the parameter name
        title(parnames_to_plot{par})

        % Draw an horizontal line on y = 0 as a reference
        yline(0, '--', 'HandleVisibility', 'off')

        % Draw a rectangle for each cluster, the color depends on the
        % significance
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


        % Cut the x axis to exclude parts of signal you don't want
        if isfield(cfg, "stim_onset_time") && isfield(cfg, "baseline_time")
            xlim([cfg.stim_onset_time - cfg.baseline_time, max(time)]);
        else
            xlim([1, max(time)]);
        end

        % Adjust the time axis
        ax = gca; % Get current axis handle
        if isfield(cfg, "stim_onset_time")
            ax.XTickLabel = (ax.XTick - cfg.stim_onset_time) / cfg.fs; % Divide by sample frequency
        else
            ax.XTickLabel = ax.XTick / cfg.fs; % Divide by sample frequency
        end

        % Draw a box cause I like it
        box on

    end % for par = 1 : par


end % cfg.wantplot_perm


%% Compute mean parameter estimates for each individual and significant cluster
if nargout > 3 % Do this only if it's required by the output to save time

    % Make parameters name suitable to be put as column names in the output
    parnames = cellfun(@(x) erase(x, '('), parnames, 'UniformOutput', false);
    parnames = cellfun(@(x) erase(x, ')'), parnames, 'UniformOutput', false);
    parnames = cellfun(@(x) replace(x, ':', 'X'), parnames, 'UniformOutput', false);

    % Loop through parameters
    out_all_clust = []; % preallocate array
    variable_names = [];
    for par = 1 : npar
        
        % Unpack current parameter
        this_parameter      = pars_series(:,:,par);
        
        % Loop through significant clusters and compute mean pars
        significant_clusters    = find(obs_clusters_sum(par).pval < cfg.perm_alpha);
        out_par_means           = NaN(height(this_parameter), length(significant_clusters));
        for tc = 1 : length(significant_clusters)
            this_cluster            = significant_clusters(tc);
            clFirst                 = obs_clusters_sum(par).first(this_cluster);
            clLength                = obs_clusters_sum(par).length(this_cluster);
            end_cluster = clFirst  + clLength;
            if end_cluster > width(this_parameter)
                end_cluster = width(this_parameter);
            end
            out_par_means(:,tc)     = mean(this_parameter(:, clFirst : end_cluster), 2);
            variable_names          = [variable_names, switchText([parnames{par}, '_cluster', num2str(significant_clusters(tc))])];
        end
        % Compute cluster means for this parameter
        out_all_clust = [out_all_clust, out_par_means];
    end 

    % Convert in nicer table for output
    clust_par_means = array2table(out_all_clust, 'VariableNames', variable_names);
    clust_par_means.id = modelout.ids;
    clust_par_means = clust_par_means(:, [end, 1:end-1]); % just because i like to have the id column as first
end


end




