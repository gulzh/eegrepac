clc
clear all
close all
%set path environment
addpath('data')
addpath('eeglab2020_0')
addpath('pac_Onslow/pac_code_best')
addpath('REPAC')
addpath('PACT-master')

%% COMPUTE THE PAC MEASURES (slow method)
addpath('data')
load REPEAT
load SWITCH
addpath('pac_Onslow/pac_code_best')

%REPEAT CONDITION
s = size(REPEAT);

MM_repeat = zeros(s(1),s(1),s(3));
fL = zeros(s(1),s(1),s(3));
fH = zeros(s(1),s(1),s(3));

for m = 1:s(1)
    for p = 1:s(1)
        for n_epoch = 1:s(3)
            signalPAC = REPEAT(p,:,n_epoch);
            signalMOD = REPEAT(m,:,n_epoch);

            fs = 500; %sampling frequency in Hz
            dur       = length(signalPAC); 
            t         = 0:1/fs:dur/fs-1/fs;
            fL_band   = [4 12]; % low frequency range
            fH_band   = [13 50]; % high frequency range
            
            [m_estim, fL_estim, fH_estim, freqvec_ph, freqvec_amp, pacmat] = compute_MI_no_plot(signalPAC, signalMOD, 500, t, fL_band, fH_band, 0);
            %no statistical test is performed
            MM_repeat(m,p,n_epoch) = m_estim;
            fL(m,p,n_epoch) = fL_estim; % keep track of the center band high frequency
            fH(m,p,n_epoch) = fH_estim; % keep track of the center band low frequency
            
        end
    end
end

%SWITCH CONDITION
s = size(SWITCH);

MM_switch = zeros(s(1),s(1),s(3));
fL = zeros(s(1),s(1),s(3));
fH = zeros(s(1),s(1),s(3));

for m = 1:s(1)
    for p = 1:s(1)
        for n_epoch = 1:s(3)
            signalPAC = SWITCH(p,:,n_epoch);
            signalMOD = SWITCH(m,:,n_epoch);

            fs = 500; %sampling frequency in Hz
            dur       = length(signalPAC); 
            t         = 0:1/fs:dur/fs-1/fs;
            fL_band   = [4 12]; % low frequency range
            fH_band   = [13 50]; % high frequency range
                        
            [m_estim, fL_estim, fH_estim, freqvec_ph, freqvec_amp, pacmat] = compute_MI_no_plot(signalPAC, signalMOD, 500, t, fL_band, fH_band, 0);
            %no statistical test is performed
            MM_switch(m,p,n_epoch) = m_estim;
            fL(m,p,n_epoch) = fL_estim; % keep track of the center band high frequency
            fH(m,p,n_epoch) = fH_estim; % keep track of the center band low frequency
            
        end
    end
end

%% faster method to compute pac measures (in the standard method way)

addpath('data')
load REPEAT
load SWITCH
addpath('pac_Onslow/pac_code_best')



[REPEAT_pac_mvl,REPEAT_pac_klmi,REPEAT_fh_mvl,REPEAT_fl_mvl,REPEAT_fh_klmi,REPEAT_fl_klmi] = MY_PAC_parallel(REPEAT);
[SWITCH_pac_mvl,SWITCH_pac_klmi,SWITCH_fh_mvl,SWITCH_fl_mvl,SWITCH_fh_klmi,SWITCH_fl_klmi] = MY_PAC_parallel(SWITCH);

%% USE PROFESSOR CISOTTO'S METHOD TO OBTAIN FREQUENCY BANDS OF INTEREST
addpath('REPAC')
parclust = parcluster;
s = size(REPEAT);
parfor (epoch = 1:s(3), parclust.NumWorkers)
    [REPEAT_fh(:,:,epoch),REPEAT_fl(:,:,epoch)] = get_REPAC_freq(REPEAT(:,:,epoch));
end

s = size(SWITCH);
parfor (epoch = 1:s(3), parclust.NumWorkers)
    [SWITCH_fh(:,:,epoch),SWITCH_fl(:,:,epoch)] = get_REPAC_freq(SWITCH(:,:,epoch));
end

%% COMPUTE THE PAC USING THE REPAC FREQUENCIES

addpath('data')
load REPEAT
load SWITCH
addpath('pac_Onslow/pac_code_best')

[REPEAT_repac_mvl,REPEAT_repac_klmi] = MY_PAC_parallel_mod(REPEAT,REPEAT_fl,REPEAT_fh);
[SWITCH_repac_mvl,SWITCH_repac_klmi] = MY_PAC_parallel_mod(SWITCH,SWITCH_fl,SWITCH_fh);

%% mean values repac

for i = 1:64
    for j = 1:64
        %REPEAT_repac_mvl_mean(i,j) = mean(REPEAT_repac_mvl(i,j,:));
        %SWITCH_repac_mvl_mean(i,j) = mean(SWITCH_repac_mvl(i,j,:));
        %REPEAT_repac_klmi_mean(i,j) = mean(REPEAT_repac_klmi(i,j,:));
        %SWITCH_repac_klmi_mean(i,j) = mean(SWITCH_repac_klmi(i,j,:));
        REPEAT_mvl(i,j) = mean(REPEAT_pac_mvl(i,j,:));
        REPEAT_klmi(i,j) = mean(REPEAT_pac_klmi(i,j,:));
        SWITCH_mvl(i,j) = mean(SWITCH_pac_mvl(i,j,:));
        SWITCH_klmi(i,j) = mean(SWITCH_pac_klmi(i,j,:));
        
        REPEAT_pac_fl_mvl(i,j) = mean(REPEAT_fl_mvl(i,j,:));
        REPEAT_pac_fh_mvl(i,j) = mean(REPEAT_fh_mvl(i,j,:));
        REPEAT_pac_fl_klmi(i,j) = mean(REPEAT_fl_klmi(i,j,:));
        REPEAT_pac_fh_klmi(i,j) = mean(REPEAT_fh_klmi(i,j,:));
        
        SWITCH_pac_fl_mvl(i,j) = mean(SWITCH_fl_mvl(i,j,:));
        SWITCH_pac_fh_mvl(i,j) = mean(SWITCH_fh_mvl(i,j,:));
        SWITCH_pac_fl_klmi(i,j) = mean(SWITCH_fl_klmi(i,j,:));
        SWITCH_pac_fh_klmi(i,j) = mean(SWITCH_fh_klmi(i,j,:));
        
        %REPEAT_repac_fl(i,j) = mean(REPEAT_fl(i,j,:));
        %REPEAT_repac_fh(i,j) = mean(REPEAT_fh(i,j,:));
        
        %SWITCH_repac_fl(i,j) = mean(SWITCH_fl(i,j,:));
        %SWITCH_repac_fh(i,j) = mean(SWITCH_fh(i,j,:));
    end
end

%% CREATE DATA STRUCTURES
%% PAC
PAC.repeat.mvl.pacval = REPEAT_pac_mvl;
PAC.repeat.mvl.fl = REPEAT_fl_mvl;
PAC.repeat.mvl.fh = REPEAT_fh_mvl;
PAC.repeat.klmi.pacval = REPEAT_pac_klmi;
PAC.repeat.klmi.fl = REPEAT_fl_klmi;
PAC.repeat.klmi.fh = REPEAT_fh_klmi;
PAC.repeat.mvl.condition = 'REPEAT'; PAC.repeat.klmi.condition = 'REPEAT';
PAC.repeat.mvl.measure = 'mean vector length MI'; PAC.repeat.klmi.measure = 'Kullback-Leibler MI';

PAC.repeat.mvl.meanpac = REPEAT_mvl;
PAC.repeat.klmi.meanpac = REPEAT_klmi;
PAC.repeat.mvl.meanfl = REPEAT_pac_fl_mvl;
PAC.repeat.mvl.meanfh = REPEAT_pac_fh_mvl;
PAC.repeat.klmi.meanfl = REPEAT_pac_fl_klmi;
PAC.repeat.klmi.meanfh = REPEAT_pac_fh_klmi;


PAC.switch.mvl.pacval = SWITCH_pac_mvl;
PAC.switch.mvl.fl = SWITCH_fl_mvl;
PAC.switch.mvl.fh = SWITCH_fh_mvl;
PAC.switch.klmi.pacval = SWITCH_pac_klmi;
PAC.switch.klmi.fl = SWITCH_fl_klmi;
PAC.switch.klmi.fh = SWITCH_fh_klmi;
PAC.switch.mvl.condition = 'SWITCH'; PAC.switch.klmi.condition = 'SWITCH';
PAC.switch.mvl.measure = 'mean vector length MI'; PAC.switch.klmi.measure = 'Kullback-Leibler MI';

PAC.switch.mvl.meanpac = SWITCH_mvl;
PAC.switch.klmi.meanpac = SWITCH_klmi;
PAC.switch.mvl.meanfl = SWITCH_pac_fl_mvl;
PAC.switch.mvl.meanfh = SWITCH_pac_fh_mvl;
PAC.switch.klmi.meanfl = SWITCH_pac_fl_klmi;
PAC.switch.klmi.meanfh = SWITCH_pac_fh_klmi;

save('PAC.mat','PAC')

%% REPAC

REPAC.repeat.mvl.pacval = REPEAT_repac_mvl;
REPAC.repeat.mvl.fl = REPEAT_fl;
REPAC.repeat.mvl.fh = REPEAT_fh;
REPAC.repeat.klmi.pacval = REPEAT_repac_klmi;
REPAC.repeat.klmi.fl = REPEAT_fl;
REPAC.repeat.klmi.fh = REPEAT_fh;

REPAC.repeat.mvl.meanpac = REPEAT_repac_mvl_mean;
REPAC.repeat.klmi.meanpac = REPEAT_repac_klmi_mean;
REPAC.repeat.mvl.meanfl = REPEAT_repac_fl;
REPAC.repeat.mvl.meanfh = REPEAT_repac_fh;
REPAC.repeat.klmi.meanfl = REPEAT_repac_fl;
REPAC.repeat.klmi.meanfh = REPEAT_repac_fh;
REPAC.repeat.mvl.condition = 'REPEAT'; REPAC.repeat.klmi.condition = 'REPEAT';
REPAC.repeat.mvl.measure = 'mean vector length MI'; REPAC.repeat.klmi.measure = 'Kullback-Leibler MI';



REPAC.switch.mvl.pacval = SWITCH_repac_mvl;
REPAC.switch.mvl.fl = SWITCH_fl;
REPAC.switch.mvl.fh = SWITCH_fh;
REPAC.switch.klmi.pacval = SWITCH_repac_klmi;
REPAC.switch.klmi.fl = SWITCH_fl;
REPAC.switch.klmi.fh = SWITCH_fh;

REPAC.switch.mvl.meanpac = SWITCH_repac_mvl_mean;
REPAC.switch.klmi.meanpac = SWITCH_repac_klmi_mean;
REPAC.switch.mvl.meanfl = SWITCH_repac_fl;
REPAC.switch.mvl.meanfh = SWITCH_repac_fh;
REPAC.switch.klmi.meanfl = SWITCH_repac_fl;
REPAC.switch.klmi.meanfh = SWITCH_repac_fh;
REPAC.switch.mvl.condition = 'SWITCH'; REPAC.switch.klmi.condition = 'SWITCH';
REPAC.switch.mvl.measure = 'mean vector length MI'; REPAC.switch.klmi.measure = 'Kullback-Leibler MI';


save('REPAC.mat','REPAC')

%% visualize results of the average PAC values for both conditions
% just an example for visualizing the results
signalMOD = {'Fpz'};
signalPAC = {'all'};
%signalPAC = {'all'};

MI_results(signalMOD, signalPAC,PAC.repeat.klmi);
MI_results(signalMOD, signalPAC,PAC.repeat.mvl);
MI_results(signalMOD, signalPAC,PAC.switch.klmi);
MI_results(signalMOD, signalPAC,PAC.switch.mvl);


MI_results(signalMOD, signalPAC,REPAC.repeat.klmi);
MI_results(signalMOD, signalPAC,REPAC.repeat.mvl);
MI_results(signalMOD, signalPAC,REPAC.switch.klmi);
MI_results(signalMOD, signalPAC,REPAC.switch.mvl);
