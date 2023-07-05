function [cluster, dist_sum_clust] = timeClustercorrection(a,b,npermut,ntp, paired, alphaThresh, tail, pooledVar)
% performs cluster-based correction for time series.
% a = matrix with Nsubject X timepoints
% b = matrix with Nsubject X timepoints, which can be condition 2 or group 2. If you don't put it, it is assumed that you want to perform "a" against zero.
% npermut = select number of permutations. Default 500.
% ntp = how many time-points per cluster (minimum)?
% paired (1, 0) default 1;
% alphaThresh = your INITIAL threshold for your matrix
% tail = 'both' or 'greater' or 'less'
% pooledVar = 1 or 0;


% check input and assign values to missing variables
if ~exist('b', 'var')           || isempty(b)           ; b = zeros(size(a))    ; end
if ndims(a) ~= 2 || ndims(b) ~= 2; error('matrices should be 2-dimensional')    ; end
if ~exist('npermut', 'var')     || isempty(npermut)     ; npermut     =  500    ; end
if ~exist('ntp', 'var')         || isempty(ntp)         ; ntp         =  2      ; end
if ~exist('paired', 'var')      || isempty(paired)      ; paired      =  1      ; end
if ~exist('alphaThresh', 'var') || isempty(alphaThresh) ; alphaThresh =  0.05   ; end
if ~exist('tail', 'var')        || isempty(tail)        ; tail        =  'both' ; end
if ~exist('pooledVar', 'var') || isempty(pooledVar)     ; pooledVar   =  'equal'; end

fprintf('\n-----------------NB!!!-------------------\n')
fprintf('Running a cluster-based correction with\n')
fprintf('paried (1) unpaired (0) = %.02f\n', paired)
fprintf('tail = %s\n', tail)
fprintf('alpha threshold = %.02f\n', alphaThresh)
fprintf('min time points = %d\n', ntp)
pause(5)

cluster = struct;

% identify cluster in the real data
[cluster] = timeClusterIdentifier(a,b,1,ntp, paired, alphaThresh, tail, pooledVar);


if cluster(1).tp_start ~= 0 % if at least one cluster was found
    
    % initialise random number generator0
    seed = sum(100*clock);
    rand('twister',seed);
    randn('state',seed);    
    
    % create subjects' split
    if paired
    ons  = ones(1,round(size(a,1)/2));
    zers = zeros(1, size(a,1)-length(ons));
    onszers = [ons zers];
    
    else
        gr1length = size(a,1);
        gr2length = size(b,1);
    end
    
    
    if paired && length(onszers) == size(a,1) && length(onszers) == size(b,1) || ~paired && gr1length == size(a,1) && gr2length == size(b,1)
        
        dist_sum_clust = [];
        for ii = 1:npermut
            clc
            fprintf('permutation N = %d\n',ii)
            
            a_perm = a;
            b_perm = b;
            
            % permute subjects
            if paired
                mypermut = logical(onszers(randperm(length(onszers))));
                a_perm(mypermut,:) = b(mypermut,:);
                b_perm(mypermut,:) = a(mypermut,:);
            else
                merge_groups = cat(1,a,b);
                merge_groups = merge_groups(randperm(size(merge_groups,1)),:);
                
                a_perm = merge_groups(1:gr1length,:);
                b_perm = merge_groups(gr1length+1:gr1length+gr2length,:);  
            end
            [cluster_perm] = timeClusterIdentifier(a_perm,b_perm,0,ntp, paired, alphaThresh, tail, pooledVar);
            dist_sum_clust = [dist_sum_clust [cluster_perm.sum]];
            
        end
    end
    
    % check this
    for cls = 1:size(cluster,2)
        
        
        cluster(cls).pval          = sum(dist_sum_clust >= cluster(cls).sum) / length(dist_sum_clust);
        
        
    end
    
    signClusters = [cluster.pval] < 0.05;
    
    if sum(signClusters) > 0
    for cls = find(signClusters)
       fprintf('Cluster %d is significant after correction, p = %.02d\n', cls, cluster(cls).pval) 
    end
    else
       fprintf('No cluster survived correction\n') 
    end
else
    cluster(1).pval = 1;
end
end