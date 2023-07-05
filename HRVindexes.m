
function [HRVtable] = HRVindexes(ibi, winLen, resamplingFactor, overlapProp, overSampFactor )
   

    %%% IBIs must be entered in ms

    %%% Default stuff
    try overSampFactor; catch; overSampFactor = 4;end
    try resamplingFactor; catch; resamplingFactor = 4; end

    try winLen; catch winLen = 256; end                                     % length of moving window for pwelch

    try overlapProp; catch; overlapProp = 0.5; end                          % percentage of overlappable samples in the moving window
    overlap = winLen * overlapProp;

    frequencies = [.003, .04, .15, .4];                                     % frequencies for spectral analysis

    HRVtable = array2table([]);                                             % table to export
    
    %% Time domain indexes
    time_differences = diff(ibi);                                           % You need successive time differences only for time domain indexes                                   

    hrseries = 6e4 ./ ibi;                                                  % Beat by beat Heart Rate
    HRVtable.avhr = mean(hrseries);                                         % Average Heart Rate (1/min)
    HRVtable.stdhr = std(hrseries);                                         % Standard deviation of Heart Rate
 
    HRVtable.MxDMn = max(ibi) - min(ibi);                                   % Max and Min IBI difference
    HRVtable.avnn = mean(ibi);                                              % Average RR
    HRVtable.sdnn = std(ibi);                                               % Standard deviation of RRs
    HRVtable.rmssd = sqrt(mean(time_differences .^ 2));                     % Root Mean Square of Successive Differences
    HRVtable.nn50 = sum(abs(time_differences) >= 50);                       % Number of successive time differences greater (or equal?) than 50 ms
    HRVtable.pnn50 = (HRVtable.nn50 / (length(ibi)-1)) * 100;               % Percentage of successive time differences greater (or equal?) than 50 ms
    
        


    %% Frequency domain indexes (inspired from mhrv package, Behar J. A., Rosenberg A. A. et al. (2018) ‘PhysioZoo: a novel open access platform for heart rate variability analysis of mammalian electrocardiographic data.’ Frontiers in Physiology.)
    
    %%% Reshape IBIs and define time axis
    ibi = ibi / 1000;                                                       % I need IBIs in seconds for this

    fsUni = resamplingFactor;                                               % probably I shouldn't call the variable resampling factor
    tibi = [0, cumsum(ibi(1:end-1))];                                       % cumulative sum of IBIs to define time

    ibi = detrend(ibi, 'linear');                                           % detrend IBIs 1
    ibi = ibi - mean(ibi);                                                  % detrend IBIs 2
    
    tibi_uni = tibi(1):  1/fsUni : tibi(end);                               % time axis
    
    winLenSec = winLen / fsUni;                                             % window length in seconds (for pwelch)
    samplesPerWindow = floor(winLenSec / (1/fsUni));                        % how many samples in the window
    nWindows = floor(length(tibi_uni) / samplesPerWindow);                  % ghow many windows


    %%% frequency axis
    ts = winLenSec / (samplesPerWindow-1);                                  % sampling time
    freqRes  = 1 / (samplesPerWindow * ts);                                 % frequency resolution
    freqRes  = freqRes / overSampFactor; 

    fAxis = (freqRes : freqRes : frequencies(end))';                        % frequency axis

    %%% interpolate IBIs
    interpibi = interp1(tibi, ibi, tibi_uni,  'spline');
    

    %%% define hamming window and overlappable samples
    window = hamming(samplesPerWindow);
    welch_overlap_samples = floor(samplesPerWindow * .5);


    %%% run welch periodogram
    [pxx,~] = pwelch(interpibi', window  , welch_overlap_samples, fAxis, fsUni);
    

    %%% convert into ms²
    ms2 = pxx * 1e6;
    

%     %% Do a very ugly plot
%     plot(fAxis, ms2)
%     xlabel('Frequency (Hz)')
%     ylabel('PSD (ms²/Hz)')
%     xline([0.04 0.15])
%       
%     firstHF = find(fAxis > .15, 1, 'first')
%     hold on
%     area(fAxis(firstHF:length(ms2)), ms2(firstHF:length(ms2)),  2,'EdgeColor','none','FaceColor','g')


    %%% Extract Indexes
    HRVtable.TOTpower   = trapz(ms2);

    HRVtable.VLFpower   = trapz(ms2(fAxis > frequencies(1) & fAxis < frequencies(2)));
    HRVtable.VLFmax     = fAxis(ms2 == max(ms2(fAxis > frequencies(1) & fAxis < frequencies(2))));
 
    HRVtable.LFpower    = trapz(ms2(fAxis > frequencies(2) & fAxis < frequencies(3)));
    HRVtable.LFmax      = fAxis(ms2 == max(ms2(fAxis > frequencies(2) & fAxis < frequencies(3))));
    HRVtable.LFpowNorm  = HRVtable.LFpower / HRVtable.TOTpower;

    HRVtable.HFpower    = trapz(ms2(fAxis > frequencies(3) & fAxis < frequencies(4)));
    HRVtable.HFmax      = fAxis(ms2 == max(ms2(fAxis > frequencies(3) & fAxis < frequencies(4))));
    HRVtable.HFpowNorm  = HRVtable.HFpower / HRVtable.TOTpower;
    
    HRVtable.LFoverHF = HRVtable.LFpower / HRVtable.HFpower;




    %% Geometric indexes
    % Triangular index
    binWidth = 1/128;
    h = histogram(ibi/1000, 'BinWidth', binWidth);
    HRVtable.TI = length(ibi) /  max(h.Values);                                  %Triangular Index - should be double checked
    close Figure 1
    


    % Triangular Interpolation of NN interval histogram
    %%%TINN = I'm confident that I may do this one day
    



    %% Baevsky's stress index. 
    % "In order to make SI less sensitive to slow changes in mean heart rate (which would increase the MxDMn and lower AMo), 
    % the very low frequency trend could be removed from the RR interval time series by using the smoothness priors method Tarvainen et al 2002. 
    % In addition, the square root of SI is taken to transform the tailed distribution of SI values towards normal distribution" (Kubios)
%     binWidth = 50;
%     h = histogram(ibi, 'BinWidth', binWidth); %,'Normalization','probability'
%     AMo = max(h.Values);
%     close Figure 1
%     Mo = median(ibi);
%     SI = (AMo * 100) / (2 * Mo * HRVtable.MxDMn); 
%     
%     clear binWidth h AMo Mo
    

    
     %% HR deceleration and acceleration capacity 
     %%% I'm confident that I'll be able to do this one day


end