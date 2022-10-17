function [fh,fl] = get_REPAC_freq(dataset)

s = size(dataset);
fs = 500; %sampling frequency in Hz
dur = s(2); 
t = 0:1/fs:dur/fs-1/fs;

phaseFreqRange  = [ 4    12 ]; % Real EEG data
ampFreqRange    = [ 13    45 ]; % Real EEG data
numPhaseFreqs   = 10;

hasRateRange     = [ 0.1 10 ];
numHasRates      = 15;
plotType         = 0; % 0=no plot. 1=MVL. 2=MI.
normalizeColor   = 1;

for m = 1:s(1)
    [fL_temp,outputLFO,LFO] = get_L(dataset(m,:));
    fl(m,1:s(1)) = fL_temp*ones(1,s(1));
    for p = 1:s(1)
        fh(m,p) = get_H(fL_temp,outputLFO,LFO,dataset(p,:));
    end
end



end

function [fL_estim,outputLFO,LFO] = get_L(s1)
    s = size(dataset);
    fs = 500; %sampling frequency in Hz
    dur = s(2); 
    t = 0:1/fs:dur/fs-1/fs;

    phaseFreqRange  = [ 4    12 ]; % Real EEG data
    ampFreqRange    = [ 13    45 ]; % Real EEG data
    numPhaseFreqs   = 10;

    hasRateRange     = [ 0.1 10 ];
    numHasRates      = 15;
    plotType         = 0; % 0=no plot. 1=MVL. 2=MI.
    normalizeColor   = 1;
    
    EEG.data = s1;
    EEG.srate = fs;
    EEG.nbchan   = size(EEG.data,1);

    EEG = pac_scanLfoPhaseFreq(EEG, phaseFreqRange, numPhaseFreqs, ampFreqRange, hasRateRange, numHasRates, plotType, normalizeColor);           
    EEG.pacScan.phaseFreqList = logspace(log10(phaseFreqRange(1)), log10(phaseFreqRange(2)), numPhaseFreqs); % added by Giulia on 2020/02/14




    %% Select LFO: take the range above threshold
    % MIN = min( EEG.pacScan.lfoPhaseVect(1,:,1) );
    % MAX = EEG.pacScan.lfoPhaseVect(d1,d2,d3); %max( EEG.pacScan.lfoPhaseVect(1,:,1) );
    % COEFF = 0.1; % Synthetic data (COEFF = 0.5; % Real data)
    % THRESH = MAX - COEFF*( MAX-MIN );
    % ind = find( EEG.pacScan.lfoPhaseVect(d1,:,d3) > THRESH ); % USING MVL

    MIN = min( squeeze( mean(EEG.pacScan.lfoPhaseVect(1,:,:),3) ) );
    [MAX, d3] = max( [squeeze( mean(EEG.pacScan.lfoPhaseVect(1,:,:),3) )] );
    COEFF = 0.2;
    THRESH = MAX - COEFF*( MAX-MIN );
    ind = find( squeeze( mean(EEG.pacScan.lfoPhaseVect(1,:,:),3) ) > THRESH ); % USING MVL

    LFO = [EEG.pacScan.phaseFreqList( ind(1) ), EEG.pacScan.phaseFreqList( ind(end)+1 )]; % IPOTESI: punti connessi
    LFO(2) = min([LFO(2) phaseFreqRange(2)]); % should be useless
    fLcb = LFO(1) + ( LFO(2) - LFO(1) ) / 2;   % centre of band
    LFObw = LFO(2) - LFO(1);                        % bandwidth


    LFO



    %% Select has: take the highest value in the LFO range
    % vectLenMax = max(EEG.pacScan.lfoPhaseVect(:)); % with 1 channel
    % % miMax      = max(EEG.pacScan.lfoPhaseMi(:));
    % [d1, d2, d3]    = ind2sub(size(EEG.pacScan.lfoPhaseVect), find(EEG.pacScan.lfoPhaseVect == vectLenMax));
    hasRateList   = logspace(log10(hasRateRange(1)),   log10(hasRateRange(2)),   numHasRates);
    % has = hasRateList(d3); % 0.1 (or 0.3) should be generally good!
    has = max( hasRateList(ind) );



    
    %% Estimate low carrier frequency
    [outputLFO, ~, ~] = FFTbasedBP(EEG.data, LFO(1), LFO(2), fs, 0);
    xx1 = hilbert(outputLFO);
    yy1 = unwrap(double(angle(xx1)));
    [~, slope1, ~] = get_linear_trend(yy1, 0);
    fL_estim = (fs*slope1)/(2*pi)

end

function fH_estim = get_H(fL_estim,outputLFO,LFO,s2)

    phaseFreqRange  = [ 4    12 ]; % Real EEG data
    ampFreqRange    = [ 13    45 ]; % Real EEG data
    numPhaseFreqs   = 10;
    
    dur = length(s2); 
    
    hasRateRange     = [ 0.1 10 ];
    numHasRates      = 15;
    plotType         = 0; % 0=no plot. 1=MVL. 2=MI.
    normalizeColor   = 1;
    fs = 500; %sampling frequency in Hz
    EEG.data = s2;
    EEG.srate = fs;
    %% Select HFO: data-driven, based on comb frequencies
    fLcb = fL_estim;
    identify_HFO % Modified by Giulia on 2018/06/05

    if HFO(1)<0, HFO(1) = 0; end %check
    %fHcb = HFO(1) + ( HFO(2) - HFO(1) ) / 2; % centre of band
    HFObw = HFO(2) - HFO(1);                       % bandwidth            
    [outputHFO, ~, f] = FFTbasedBP(EEG.data, HFO(1), HFO(2), fs, 0);



    %% Estimate high carrier frequency
    xx2 = hilbert(outputHFO);
    yy2 = unwrap(double(angle(xx2)));
    [lintrend, slope2, offset] = get_linear_trend(yy2, 0);
    fH_estim = (fs*slope2)/(2*pi)
    fHcb = fH_estim;


end