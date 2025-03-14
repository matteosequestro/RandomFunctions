% this require the check_clusters function
function [cluster_mass, cluster_mass_split, obs_clusters_sum] = run_clust_perm(iter, seriea, serieb, test, alpha, minlegth, want_plot)

%%% iter: number of the iterations
%%% seriea: first matrix subj * samplepoint
%%% serieb: second matrix subj * samplepoint
%%% test: test for the cluster detection, default is paired t-test ('ttest')
%%% alpha: for the moment I use the same alpha for the cluster detection and the permutation, this could be separated in the future
%%% minlength: minimum length for clusters. Default is 1
%%% want_plot: if one will produce the plot with grand averages, confidence intervals and all detected clusters between time series. Default is 0

%set default values
try test; catch test = 'ttest'; end
try alpha; catch alpha = 0.05; end
try minlength; catch minlength = 1; end
try want_plot; catch want_plot = 0; end


% Temporary errors untill i adjust the function
% if strcmp(test, 'ttest2'); error('Maybe you can permute with two samples t tests but you need to adjust the function'); end



%check clusters in the observed data
[~, clusters, obs_clusters_sum] = check_clusters(seriea, serieb, test, alpha, minlegth, 0, 0,  [], [], []);


%if there's no cluster then interrupt the function
if isempty(clusters)
    warning('no cluster detected')


    if want_plot & strcmp(test, 'ttest0')

        [~, ~, ci, ~] = ttest(seriea);
        time = 1:length(seriea);
        patch([time, fliplr(time)],[ci(1,:), fliplr(ci(2,:))],  [0.2588 0.2863 0.2863]  ,'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off')
        hold on
        plot(mean(seriea),'Color', [0.2588 0.2863 0.2863], 'LineWidth', 2)

%         if ~isnan(ylims); ylim(ylims); end
        xlim([1, width(seriea)]);

    end
    %if there are clusters then run the permutations
else
    %%%%%%%%%%%%%%%%%% Permutations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cluster_mass = [];                                                  %you will have the cluster masses for each permutation here
                           
    for rep = 1:iter
        %%% If t-test between subjects
        if strcmp(test, 'ttest2')
            
            %%% Join and permute the dataset
            serie = [seriea; serieb];
            serie = serie(randperm(height(serie)), :, :);
            
            %%% Split the permuted dataset to mantain the original proportion
            serieap = serie(1 : height(seriea), :, :);
            seriebp = serie((height(seriea)+1):end, :, :);
            
            %%% Check clusters
            [~, clusters_per] = check_clusters(serieap(:,:), seriebp(:,:), test, alpha, minlegth, 0);
        
            
        %%% If within subjects (here I include tests against 0)
        else
            %%% Build a structure subject x time point x condition (this is actually not a necessary passage)
            serie = cat(3, seriea, serieb);

            %%% Series that will be permuted
            perm_serie = serie;                                            

            %%% Randomize for which participants to switch labels and then switch labels
            randomed = randi(2, 1, size(serie,1));                          
            perm_serie(find(randomed == 2), :, [1 2]) = perm_serie(find(randomed == 2), :, [2 1]);  
            
            %%% Check clusters
            [~, clusters_per] = check_clusters(perm_serie(:,:,1), perm_serie(:,:,2), test, alpha, minlegth, 0); 
        end

        if ~isempty(clusters_per)                                       %if you find clusters ...
            t_masses_per = [];

            for this_cluster_per = 1:size(clusters_per,2)
                t_masses_per = [t_masses_per, sum(clusters_per{this_cluster_per}(:,1))];   %then compute the cluster mass for each cluster
            end

            %%% APPROACH 1: consider only the biggest permuted cluster
            %             [~, max_cluster] = max(abs(t_masses_per));                  %and take the biggest cluster (check in future if you want the absolute t value)
            %             cluster_mass = [cluster_mass, sum(clusters_per{max_cluster}(:,1))] ; %and join to the list


            %%% APPROACH 2: consider all the permuted clusters
            cluster_mass = [cluster_mass, t_masses_per];                % Consider all the clusters, not just the biggest

        else
            cluster_mass = [cluster_mass, 0];                           %if you find no cluster the cluster mass is zero
        end
        if ismember(rep, [0:iter/10:iter]); disp([num2str(rep), '/', num2str(iter)]); end %this line just tells you how many iterations have been done, one 10th at the time
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %Split the permuted t masses in a distribution greater or equal to 0 and another one lower or equal to 0 (0 is in common)
    %     cluster_mass_split = {cluster_mass(cluster_mass <= 0), cluster_mass(cluster_mass >= 0)};


    %compute the monte-carlo p values
    for i = 1:size(obs_clusters_sum, 2)                                     %for each t mass (for each cluster)

        %%% MY ORIGINAL APPROACH: compute negative observed tmasses with
        %%% negative part of the permuted distribution and positive
        %%% observed tmasses with the positive part of the distribution.
        %%% In future you may want to change the way you consider/compute
        %%% pvalues if you want to use this (i.e. maybe use an alpha of
        %%% 0.025 since it should still be one tailed).
        %             if obs_clusters_sum(i).tmass > 0                                          %if the t mass is greater than 0
        %                 p = sum(cluster_mass_split{2} >= obs_clusters_sum(i).tmass) / length(cluster_mass_split{2});  %the compute the p as the percentage of permuted t masses on the right (greater)
        %                 obs_clusters_sum(i).pval = p;
        %             elseif obs_clusters_sum(i).tmass < 0                                      %if the observed t mass is lower than 0
        %                 p = sum(cluster_mass_split{1} < obs_clusters_sum(i).tmass) / length(cluster_mass_split{1}); %compute the p as the percentage of permute t masses on the left (lower)
        %                 obs_clusters_sum(i).pval = p;
        %             end



        %%% CURRENT APPRAOCH (I actually prefer this): compare the absolute
        %%% value of the observed tmass with the absolute value of the
        %%% permuted distribution
        p = sum(abs(cluster_mass) > abs(obs_clusters_sum(i).tmass)) / length(cluster_mass); %compute the p as the percentage of |permuted t masses| greater than the observed |tmass|
        obs_clusters_sum(i).pval = p;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if want_plot & strcmp(test, 'ttest0')

        [~, ~, ci, ~] = ttest(seriea);
        time = 1:length(seriea);
        patch([time, fliplr(time)],[ci(1,:), fliplr(ci(2,:))],  [0.2588 0.2863 0.2863]  ,'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off')
        hold on
        plot(mean(seriea),'Color', [0.2588 0.2863 0.2863], 'LineWidth', 2)
        hold on



        for this_cluster = 1 : size(obs_clusters_sum, 2)
            hold on
            clFirst = obs_clusters_sum(1, this_cluster).first;
            clLength = obs_clusters_sum(1, this_cluster).length;
            ylims = get(gca, 'Ylim');

            if obs_clusters_sum(1, this_cluster).pval >= alpha
                rectColor = [.5 .5 .5];
            else

                rectColor =[0.3569 0.6078 0.8353];
            end

            rectangle('Position', [clFirst, ylims(1), clLength, sum(abs(ylims))], ...
                'FaceColor', [rectColor .10], 'EdgeColor',[rectColor .1])

        end
    end

%     if ~isnan(ylims); ylim(ylims); end
    xlim([1, width(seriea)]);


    %     if want_plot & strcmp(test, 'ttest0')
    %
    %         [~, ~, ci, ~] = ttest(seriea);
    %         time = 1:length(seriea);
    %         patch([time, fliplr(time)],[ci(1,:), fliplr(ci(2,:))],  [0.2588 0.2863 0.2863]  ,'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off')
    %         hold on
    %         plot(mean(seriea),'Color', [0.2588 0.2863 0.2863], 'LineWidth', 2)
    %     end



end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for clust = 1:size(obs_clusters_sum, 2)
    disp(['Cluster ' num2str(clust) ' pval: ' num2str(obs_clusters_sum(clust).pval)]);
end

end




