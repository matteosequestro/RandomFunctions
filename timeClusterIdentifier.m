function [cluster] = timeClusterIdentifier(a,b,info,ntp, paired, alphaThresh, tail, pooledVar)
% a = matrix with Nsubject X timepoints
% b = matrix with Nsubject X timepoints, which can be condition 2 or group 2. If you don't put it, it is assumed that you want to perform "a" against zero.
% info = prints information while running
% ntp = how many time-points per cluster (minimum)?
% paired (1, 0) default 1;
% alphaThresh = your INITIAL threshold for your matrix
% tail = 'both' or 'greater' or 'less'
% pooledVar = 1 or 0;


% check input and assign values to missing variables
if ndims(a) ~= 2 || ndims(b) ~= 2; error('matrices should be 2-dimensional')    ; end
if ~exist('info', 'var')        || isempty(info)        ; info        =  0      ; end
if ~exist('ntp', 'var')         || isempty(ntp)         ; ntp         =  2      ; end
if ~exist('paired', 'var')      || isempty(paired)      ; paired      =  1      ; end
if ~exist('alphaThresh', 'var') || isempty(alphaThresh) ; alphaThresh =  0.05   ; end
if ~exist('tail', 'var')        || isempty(tail)        ; tail        =  'both' ; end
if ~exist('pooledVar', 'var') || isempty(pooledVar)     ; pooledVar   =  'equal'; end

cluster = struct;
    
    % paired or unpaired ttest
    if paired
        [thresh.H,thresh.P,thresh.CI,thresh.STATS] = ttest(a,b, 'alpha', alphaThresh, 'tail', tail, 'dim', 1);
    else
        [thresh.H,thresh.P,thresh.CI,thresh.STATS] = ttest2(a,b, 'alpha', alphaThresh, 'tail', tail, 'dim', 1);
    end
    
    % if all the values are the same it gives NaN. Change nans with 0
    thresh.H(isnan(thresh.H)) = 0;
    
    % extract custers
    scanvector = 1;
    position   = 1;
    cl_end     = [];
    cl_start   = [];
    temp_pos   = 0;
    
    while scanvector
        temp_pos  = temp_pos + find(diff(thresh.H(position:end)) == -1, 1,'first');
        
        if ~isempty(temp_pos)
            cl_end    = [cl_end   temp_pos];
            
            cl_start  = [cl_start find(diff(thresh.H(1:cl_end(end))) == 1, 1,'last')+1;];
            
            if ~isempty(cl_end) && isempty(cl_start) &&  thresh.H(1) == 1% this happens if the first tp is already significant
               cl_start  = 1;                
            end
            

        else
            if  (~isempty(find(diff(thresh.H(position:end)) == 1, 1,'last')) && isempty(cl_start)) || (~isempty(cl_start) && ~isempty(find(diff(thresh.H(position:end)) == 1, 1,'last')) && find(diff(thresh.H(position:end)) == 1, 1,'last') ~= cl_start(end))
                cl_end    = [cl_end   length(thresh.H)];
                cl_start  = [cl_start find(diff(thresh.H) == 1, 1,'last')+1];
                temp_pos = length(thresh.H);
            end
            if unique(thresh.H == 1) == 1
                cl_end    = [cl_end   length(thresh.H)];
                cl_start  = [cl_start 1];
                temp_pos = length(thresh.H);
            end
        end
        
        if temp_pos < length(thresh.H)
            position = temp_pos+1;
        else
            scanvector = 0;
        end
        
    end
    
    
    if ~isempty(cl_start)
        
        cl_length = cl_end-cl_start+1; % need to add one point otherwise for instance end = 9 - start = 8 = 1, while it's 2points
        
        % filter for clusters with at least the desired numbr of time points
        cl_start  = cl_start(cl_length >= ntp);
        cl_end    = cl_end(cl_length >= ntp);
        cl_length = cl_length(cl_length >= ntp);
        
        % inform on the number of identified clusters in the data
        if info
            fprintf('%d clusters identified', length(cl_length))
        end
        
        if ~isempty(cl_start) % if there are clusters surviving the minimum timepoint number
            % calculate cluster stat and save data
            for cl = 1:length(cl_length)
                cluster(cl).tp_start = cl_start(cl);
                cluster(cl).tp_end   = cl_end(cl);
                cluster(cl).sum   = sum(thresh.STATS.tstat(cl_start(cl):cl_end(cl)).^2);
            end
            
        else
            if info
                warning('no clusters identified')
            end
            cluster.tp_start = 0;
            cluster.tp_end   = 0;
            cluster.sum   = 0;
            
        end
        
        
    else
        if info
            warning('no clusters identified')
        end
        cluster.tp_start = 0;
        cluster.tp_end   = 0;
        cluster.sum   = 0;
    end
end