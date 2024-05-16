 %% useful for cluster script, check the clusters between the 2 time series
 function [tvals, cluster_full, cluster_sum] = check_clusters(seriea, serieb, test, alpha, minlegth, want_plot, want_reptangle, ylims, linewidth, colors)
    
    %%% seriea: first matrix subj * samplepoint
    %%% serieb: second matrix subj * samplepoint
    %%% test: test for cluster detection (ttest is default, ttest0 and ttest2 are supported but may need more work on the plotting part)
    %%% alpha: alpha for cluster detection, default is 0.05
    %%% minlength: minimum length for clusters. Default is 1
    %%% want_plot: if one will produce the plot with grand averages, confidence intervals and all detected clusters between time series. Default is 1
    %%% ylims. y axis limits for the plot
    %%% linewidth: line width for the plot. Default is 2
    %%% colors: horizontal array with two exadecimal colors for the plots (must be strings beginning with #)



    

    % Set defaults (this could be done in a better way in future)
    try test; catch test = 'ttest'; end
    try alpha; catch alpha = 0.05; end
    try want_plot; catch want_plot = 1; end
    try ylims; catch ylims = NaN; end
    try minlength; catch minlength = 1; end
    try linewidth; catch linewidth = 2; end


     % Look for colors or set default (default should be gray and orange)
    try 
        color1 = convertStringsToChars(colors(1));
        color2 = convertStringsToChars(colors(2));
        
        %%% I couldn't make the function I use for plotting CIs to accept
        %%% hexadecimal colors, so this is converting them in RGB
        color1 = sscanf(color1(2:end),'%2x%2x%2x',[1 3])/255;
        color2 = sscanf(color2(2:end),'%2x%2x%2x',[1 3])/255;
    catch
        color1 = [0.2588 0.2863 0.2863];    %#414949
        color2 = [1 164/255 79/255];        
    end
    
    if strcmp(test, 'ttest')
        [~, ~, ~, stats] = ttest(seriea, serieb, alpha);                               % Paired t t values
    
    elseif strcmp(test, 'ttest2')
        [~, ~, ~, stats] = ttest2(seriea, serieb, alpha);                              % 2 samples t values for the series
    
    elseif strcmp(test, 'ttest0')
        serieb = zeros(size(seriea));
        [~, ~, ~, stats] = ttest(seriea, serieb, alpha);                               % One sample t values against 0
    
    else error('Test non valid')
    end
    
   

    %% Find clusters (this could probably be done in a easier way)
    tresh = tinv(alpha/2, stats.df(1));                                     %take the treshold (quantile at that alpha value for those degrees of freedom)
    t = [(stats.tstat)', (1:length(stats.tstat))'];                         %reshape the ts
    
    clusters_serie = t(find(t(:,1) < tresh | t(:,1) > -tresh), :);          %find the values over the treshold
    
   if ~isempty(clusters_serie)                                              %if you find the cluster...
        cluster_struc = {};
        this_cluster = [];
        cluster_n = 1;
        
        for sample = 1:length(clusters_serie(:,2))                          %then separate the clusters based on temporal contiguity in the following way:
            if sample == 1;                                                 %if it's the first sample then put it on the list
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



%    %%%%%%%%%%%%%%%%%%%%%%%%%  Not using this for the moment. If uncommented will merge clusters closer than 'mergetreshold' from each other
%    mergethreshold = 10;
%    clusterSet = [];
%    for win = 1 : size(cluster_full, 2)
%        if win == 1
%            clusterSet = [clusterSet, cluster_full(1, win)];
%        else
%            if (cluster_full{1, win}(1,2) - clusterSet{1, end}(end,2)) <= mergethreshold
%                newWindow = {[clusterSet{1, end}; cluster_full{1, win}]};
%                clusterSet{1, end} = [];
%                clusterSet = [clusterSet, newWindow];
%            elseif (cluster_full{1, win}(1, 2) - clusterSet{1, end}(end, 2)) > mergethreshold;
%                clusterSet = [clusterSet, cluster_full(1, win)]; end
%        end
%    end
% 
%    try
%        clusterSet(cellfun(@isempty,clusterSet)) = [];
%    end
%    cluster_full = clusterSet;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


    
   % Delete clusters smaller than minlength. This could probably be done in
   % a single line, so I may update it in future
   toDelete = zeros(1, size(cluster_full, 2));
   for clust = 1:size(cluster_full, 2)
       if size(cluster_full{clust}, 1) < minlegth;
           toDelete(clust) = 1;
       end
   end
   cluster_full(:, find(toDelete == 1)) = [];
   clear clust



   % save info for each cluster
   cluster_sum = [];
   for clus = 1:size(cluster_full, 2)
       cluster_sum(clus).length     = size(cluster_full{clus}, 1);
       cluster_sum(clus).first      = cluster_full{clus}(1,2);
       cluster_sum(clus).last       = cluster_full{clus}(end,2);
       cluster_sum(clus).tmass      = sum(cluster_full{clus}(:,1));
   end


   %% Here you plot the clusters if you want
   if want_plot
%        f = figure;
       average_seriea = nanmean(seriea);


       %%% Plot if you have one sample t values against 0
       if strcmp(test, 'ttest0')
           [~, ~, ci, ~] = ttest(seriea);
           
           time = 1:length(seriea);
            
           patch([time, fliplr(time)],[ci(1,:), fliplr(ci(2,:))],  color1  ,'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off')
           hold on 
           plot(average_seriea,'Color', color1, 'LineWidth', linewidth)
           hold on

%            patch([time, fliplr(time)], [ci(1,:), fliplr(ci(2,:))], 'red','FaceAlpha',0.1, 'EdgeColor','none')
%            hold on 
%            plot(average_seriea, 'r')
%            hold on

           if ~isnan(ylims); ylim(ylims); end
           xlim([1, width(seriea)]);


       %%% Plot if you have one sample paired t values or 2 samples t values
       else
           average_serieb = nanmean(serieb); 
            
           %%% Compute standard error for paired t tests
           if strcmp(test, 'ttest')

               average_series = (seriea + serieb) ./ 2;

               grandAverage = nanmean([seriea; serieb]);
        

               adjFactor =  grandAverage - average_series;
        
               seriea_adj = seriea + adjFactor;
               serieb_adj = serieb + adjFactor;

               a = nanmean(seriea_adj);
               b = nanmean(serieb_adj);
               SEMa = nanstd(seriea_adj) / sqrt(height(seriea_adj));
               SEMb = nanstd(serieb_adj) / sqrt(height(serieb_adj));


     
                %%% This will plot the standard error
%                CIa_low = average_seriea- SEMa;
%                CIa_up = average_seriea + SEMa;
%             
%                CIb_low = average_serieb -  SEMb;
%                CIb_up = average_serieb +  SEMb;


               %%% This will plot the confidence intervals
               ts = tinv([0.025  0.975], height(seriea)-1);      % T-Score

               CIa_low = average_seriea + ts(1) * SEMa;
               CIa_up = average_seriea + ts(2) * SEMa;
            
               CIb_low = average_serieb + ts(1) * SEMb;
               CIb_up = average_serieb + ts(2) * SEMb;


           %%%  Confidence intervals for 2 samples t test
           elseif strcmp(test, 'ttest2')
               SEMa = nanstd(seriea)/sqrt(height(average_seriea));
               SEMb = nanstd(serieb)/sqrt(height(average_serieb));

%                [~, ~, cia, ~] = ttest(seriea);
%                [~, ~, cib, ~] = ttest(serieb);
%                 CIa_low = cia(1, :);
%                 CIa_up = cia(2, :);
%         
%                 CIb_low = cib(1, :);
%                 CIb_up = cia(2, :);

           ts = tinv([0.025  0.975], height(seriea)-1);      % T-Score
%     
           CIa_low = average_seriea + ts(1) * SEMa;
           CIa_up = average_seriea + ts(2) * SEMa;
        
           CIb_low = average_serieb + ts(1) * SEMb;
           CIb_up = average_serieb + ts(2) * SEMb;

%            CIa_low = average_seriea - SEMa;
%            CIa_up = average_seriea + SEMa;
% 
%            CIb_low = average_serieb - SEMb;
%            CIb_up = average_serieb + SEMb;


           %%% If you made a mistake
           else 
               error('Still Invalid Test')
           end
           
           %%% Common plotting plottinh
           time = 1:length(average_seriea);
    
           patch([time, fliplr(time)], [CIa_low, fliplr(CIa_up)], color1,'FaceAlpha',0.3, 'EdgeColor','none', 'HandleVisibility', 'off')
           hold on 
           plot(average_seriea, 'Color', color1, 'LineWidth', linewidth)
           hold on
           patch([time, fliplr(time)], [CIb_low, fliplr(CIb_up)],  color2  ,'FaceAlpha',0.3, 'EdgeColor','none', 'HandleVisibility', 'off')
           hold on 
           plot(average_serieb,'Color', color2, 'LineWidth', linewidth)

           % If you specified limits for the y axis as an argument of the function, the use them
           if ~isnan(ylims); ylim(ylims); end
           xlim([1, width(average_seriea)]);
       end


       %% This creates a shaded rectangle for each cluster
       if want_reptangle
           if height(cluster_sum) > 0
               for this_cluster = 1:size(cluster_sum, 2)
                   hold on
                   clFirst = cluster_sum(1, this_cluster).first;
                   clLength = cluster_sum(1, this_cluster).length;
                   ylims = get(gca, 'Ylim');
                   rectangle('Position', [clFirst, ylims(1), clLength, sum(abs(ylims))], ...
                       'FaceColor',[.5 .5 .5 .10], 'EdgeColor',[.5 .5 .5 .1])

               end
           end
       end


%%% This will create a shade between the two lines for each cluster (instead of the rectangle: uncommon but it may look nice)
%        for this_cluster = 1:size(cluster_full, 2)
%            if height(cluster_full{this_cluster}) > 1
%                hold on
%                x2 = cluster_full{this_cluster}(:,2);
%                curve1 = average_seriea(:,cluster_full{this_cluster}(:,2))';
%                curve2 = average_serieb(:,cluster_full{this_cluster}(:,2))';
%         
%                shade(x2, curve1, x2, curve2, 'FillType',[1 2; 2 1]);
%            end
%        end
%%%

   end


   tvals = stats.tstat;
end