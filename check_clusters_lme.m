
function [modelout, cluster_full_beta, cluster_full_t, cluster_sum, cluster_beta_exp, cluster_t_exp] = check_clusters_lme(data, formula, cfg)
% ----------------------------------------------------------------------------------------------------------------------------------------------
% Find clusters of effect from parameter estimates in a linear mixed
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
% cfg.wantplot_iter     : 1 if you want plots from the function. You will
%                           get time-series of parameter estimates with
%                           confidence intervals; time-series of
%                           information criteria (AIC and BIC); and
%                           time-series of R-squared (ordinary and
%                           adjusted)
% cfg.perm_alpha        : permutation alpha. Montecarlo p-values below such
%                           a value will be considered as significant.
% cfg.minlength         : minimum length (in samples) to consider a cluster 
%                           Detected clusters shorter than this will not be
%                           exported.
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
% cluster_sum           : one row for each parameter. You get length,
%                           position of first sample, length, position of
%                           last sample and beta-mass (sum of beta
%                           coefficients) for each cluster.
% cluster_beta_exp      : time series of beta coefficients for each cluster
% cluster_t_exp         : time series of t-values coefficients for each cluster
% ----------------------------------------------------------------------------------------------------------------------------------------------

% Extract the dependent variable from formula
y = formula(1:find(formula == '~')-1);
y = erase(y, ' ');

% Length of time sereis
tslen = length(data.(y){1});

% Preallocate arrays for information criteria and log likelihod to export
modelout.criteria.loglik = zeros(1, tslen);
modelout.criteria.AIC = modelout.criteria.loglik;
modelout.criteria.BIC = modelout.criteria.loglik;
modelout.criteria.deviance = modelout.criteria.loglik;

% Preallocate arrays for R-squared to export
modelout.rsqrd.ordinary = modelout.criteria.loglik;
modelout.rsqrd.adjusted = modelout.criteria.loglik;

% Backup of the dependent variable
data.y2 = data.(y);

% Run the models for each sample point
for tt = 1 : tslen
    data.(y) = cellfun(@(x) x(tt), data.y2);
    rm = fitlme(data, formula); % the model

    % Preallocate vectors for parameter estimates and ci (only once)
    if tt == 1
        npars = height(rm.fixedEffects);
        
        modelout.pars.estimates = zeros(npars, tslen);
        modelout.pars.lower_ci  = modelout.pars.estimates;
        modelout.pars.upper_ci  = modelout.pars.estimates;
        modelout.pars.t_stat    = modelout.pars.estimates;
    end
    
    % Export parameters
    modelout.pars.estimates(:,tt) = rm.fixedEffects;
    modelout.pars.lower_ci(:,tt) = rm.Coefficients(:, end-1);
    modelout.pars.upper_ci(:,tt) = rm.Coefficients(:, end);
    modelout.pars.t_stat(:,tt) = rm.Coefficients(:, 4);

    % Export information criteria and log likelihod
    modelout.criteria.loglik(:,tt) = rm.LogLikelihood;
    modelout.criteria.AIC(:,tt) = rm.ModelCriterion{1,1};
    modelout.criteria.BIC(:,tt) = rm.ModelCriterion{1,2};
    modelout.criteria.deviance(:,tt) = rm.ModelCriterion{1,3};

    % Export R squared
    modelout.rsqrd.ordinary(:,tt) = rm.Rsquared.Ordinary;
    modelout.rsqrd.adjusted(:,tt) = rm.Rsquared.Adjusted;

end

% Save predictors names
modelout.pars.namesout = rm.CoefficientNames;

% Plot stuff
if cfg.wantplot_fit
    % -------------------------------------------------------------------------------------
    % Plot parameter estimates
    % -------------------------------------------------------------------------------------
    % Provide default color if not provided as argument. Give an error if
    % they are not enough to cover all coefficients
    if ~isfield(cfg, 'colors')
        colors = ['#C0C0C0';'#808080'; '#000000';'#FFA500'; '#A52A2A';'#800000'];
        if height(colors) < npars
            error(['You need to provide ', num2str(npars), ' colors, but you only have ', num2str(height(colors)) ' by default. You need to provide colors manually' ]);
        end
    else
        if height(colors) < npars
            error(['You need to provide ', num2str(npars), ' colors, but you only have ', num2str(height(colors))  ]);
        end
    end

    % I couldn't make the function I use for plotting CIs to accept
    % hexadecimal colors, so this is converting them in RGB
    colors_rgb = cell(height(colors), 1);
    for col = 1 : height(colors)
        colors_rgb(col)  = {sscanf(colors(col, 2:end),'%2x%2x%2x',[1 3])/255};
    end
    
    figure
    % Plot each condition
    time = 1:width(modelout.pars.estimates);
    for rr = 1  : height(modelout.pars.estimates)
        patch([time, fliplr(time)], [modelout.pars.lower_ci(rr,:), fliplr(modelout.pars.upper_ci(rr,:))], ...
                colors_rgb{rr},'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off'); 
        hold on
        plot(modelout.pars.estimates(rr,:), 'Color', colors_rgb{rr}, 'LineWidth', 2, 'DisplayName',modelout.pars.namesout{rr})
    end
    xlim([1, max(time)]);
    legend('Location', 'north')
        ylabel('parameter estimate (a.u.)')
    xlabel('sample')
    title('parameter estimates')
    yline(0, '--', 'HandleVisibility', 'off')

    % -------------------------------------------------------------------------------------
    % Plot R squared
    % -------------------------------------------------------------------------------------
    figure
    plot(modelout.rsqrd.ordinary)
    hold on
    plot(modelout.rsqrd.adjusted)
    legend({'ordinary', 'adjusted'})
    title('R²')
    ylabel('R²')
    xlabel('sample')
    yline(0, '--', 'HandleVisibility', 'off')
    
    % -------------------------------------------------------------------------------------
    % Plot information criteria (AIC and BIC)
    % -------------------------------------------------------------------------------------
    figure
    plot(modelout.criteria.AIC)
    hold on
    plot(modelout.criteria.BIC)
    legend({'AIC', 'BIC'})
    title('Information Criteria')
    ylabel('criterion')
    xlabel('sample')
end


%% Find clusters
%(this could probably be done in a easier way but it's the same thing as
% what i was doing in the other function and non c'ho voglia to find 
% another solution rn)        
cluster_n = 1;
for pred = 1 : npars
    clear cluster_full_beta cluster_full_t
    beta_series = [modelout.pars.estimates(pred, :)', (1 : length(modelout.pars.estimates(pred, :)))'];
    t_series = [modelout.pars.t_stat(pred, :)', (1 : length(modelout.pars.t_stat(pred, :)))'];
    tupper_ci = modelout.pars.upper_ci(pred,:);
    tlower_ci = modelout.pars.lower_ci(pred,:);
    
    % Clusters here are based on the confidence interval including or not
    % zero (for the moment I use the default 95%CI, this could be changed
    % in the function to make it more flexible).
    clusters_beta_serie = beta_series(tlower_ci > 0 | tupper_ci <0,:);
    clusters_t_serie = t_series(tlower_ci > 0 | tupper_ci <0,:);

    if ~isempty(clusters_beta_serie)                                              %if you find the cluster...
        this_beta_cluster = [];
        this_t_cluster = [];


        for sample = 1:length(clusters_beta_serie(:,2))                          %then separate the clusters based on temporal contiguity in the following way:
            if sample == 1                                                 %if it's the first sample then put it on the list
                this_beta_cluster = [this_beta_cluster; clusters_beta_serie(sample,:)];
                this_t_cluster = [this_t_cluster; clusters_t_serie(sample,:)];

            elseif clusters_beta_serie(sample,2) == clusters_beta_serie(sample-1, 2) + 1 %if the next sample is at the successive time point, then put it on the same list (i.e. the list for that cluster)
                this_beta_cluster = [this_beta_cluster; clusters_beta_serie(sample,:)];
                this_t_cluster = [this_t_cluster; clusters_t_serie(sample,:)];

            elseif clusters_beta_serie(sample,2) ~= clusters_beta_serie(sample-1, 2) + 1 %if the next sample is not at the successive time point then form a new list
                cluster_full_beta{cluster_n} = this_beta_cluster;
                cluster_full_t{cluster_n} = this_t_cluster;

                this_beta_cluster = [[]; clusters_beta_serie(sample,:)];
                this_t_cluster = [[]; clusters_t_serie(sample,:)];
                cluster_n = cluster_n + 1;
            end
        end
        cluster_full_beta{cluster_n} = this_beta_cluster;
        cluster_full_t{cluster_n} = this_t_cluster;
    else
        cluster_full_beta = {};
        cluster_full_t = {};
    end

    % Delete clusters smaller than minlength.
    cluster_full_beta= cluster_full_beta(cellfun(@(x) length(x) >= cfg.minlength,  cluster_full_beta));
    cluster_full_t= cluster_full_t(cellfun(@(x) length(x) >= cfg.minlength,  cluster_full_t));

    % Save info for each cluster
    cluster_sum(pred,:).pred =  modelout.pars.namesout{pred};
    cluster_sum(pred,:).length  = cellfun(@(x) size(x, 1), cluster_full_beta);

    cluster_sum(pred,:).first   = cellfun(@(x) x(1, 2), cluster_full_beta);
    cluster_sum(pred,:).last    = cellfun(@(x) x(end, 2), cluster_full_beta);
    cluster_sum(pred,:).mass   = cellfun(@(x) sum(x(:,1)), cluster_full_beta);

    if isempty(cluster_full_beta); cluster_full_beta = {NaN}; end
    if isempty(cluster_full_t); cluster_full_t = {NaN}; end
    cluster_beta_exp(pred, :) = {modelout.pars.namesout{pred}, cluster_full_beta};
    cluster_t_exp(pred, :) = {modelout.pars.namesout{pred}, cluster_full_t};
end


end



