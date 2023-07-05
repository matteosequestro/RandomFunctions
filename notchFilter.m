% Humbly adapted from https://dsp.stackexchange.com/questions/1088/filtering-50hz-using-a-notch-filter-in-matlab
function [filtered] = notchFilter(signal, notchFreq, sampFreq)
    
    fs = sampFreq;          % sampling rate
    f0 = notchFreq;         % notch frequency
    fn = fs/2;              % Nyquist frequency
    freqRatio = f0/fn;      % ratio of notch freq. to Nyquist freq.
    
    notchWidth = 0.1;       % width of the notch
    
    % Compute zeros
    notchZeros = [exp( sqrt(-1)*pi*freqRatio ), exp( -sqrt(-1)*pi*freqRatio )];
    
    % Compute poles
    notchPoles = (1-notchWidth) * notchZeros;
    
    b = poly( notchZeros ); %  Get moving average filter coefficients
    a = poly( notchPoles ); %  Get autoregressive filter coefficients
    
    
    % filter signal 
    filtered = filter(b, a, signal);

end