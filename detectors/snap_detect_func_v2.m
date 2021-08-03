function [det_sam, det_tim, det_amp, det_score, ncc, nlg,file_dbrms] = snap_detect_func_v2(y,kern_mat_file,min_det_score,fs,use_par)
%  
% [det_sam, det_tim, det_amp, det_score, ncc, nlg,file_dbrms]...
%             = snap_detect_func(y,kern_mat_file,min_det_score,fs,use_par)
% 
%   Where the INPUTS are:
%   y = pressure corrected waveform (unfiltered)
%       {Note that the sample rate is given in kernel file}  
% 
%  'kern_mat_file' = text string the contains name of the mat file 
%    Variable in the file are: 
%     kern = eveloped kernel function 
%     klen = length of kernel (# samples)
%     pk =  position of peak in kernel (# samples) 
%     kfs = sample rate of kernel and data 
%     flt_low & flt_high - bandpass used to make kernel 
%   
%   min_det_score is the minimum correlation coefficent needed to declare a
%   detection.  Set this low (e.g., 0.7) and then you can filter detections
%   later. 
% 
%   fs = sample rate of the data 
% 
%   use_par = 1 will run this function using MATLAB's parallel processing
%   pool 
% 
%   OUTPUTS are: 
%   det_sam - detection position in points 
%   det_tim - detection position in time (seconds relative to start of file)  
%   det_amp - amplidue of the detection uPa (after bandpass) 
%   det_score -  correlation coefficient 
%   ncc - time series of the normalize correlation coefficient (for
%   plotting purposes) 
%   nlg - time series of the lags (for plotting purposes) 
%   file_dbrms - the dB RMS of the waveform 
% 
% drbohnen@ncsu.edu 
% 
% Bohnenstiehl DR, Lillis A, Eggleston DB (2016) The Curious Acoustic 
% Behavior of Estuarine Snapping Shrimp: Temporal Patterns of Snapping 
% Shrimp Sound in Sub-Tidal Oyster Reef Habitat. PLoS ONE 11(1): 
% e0143691. https://doi.org/10.1371/journal.pone.0143691
% 

K = load(kern_mat_file);  % load the kern file 
% if fs == 96, nPoints = 100
if fs == 96000
    nPoints = 100;
elseif fs == 48000
    nPoints = 50;
end

kern=cat(2,zeros(1,100), K.kern);  % zero pad from the left by 100 points
%pk=pk+100; klen=klen+100;  % adjust peak position and length of kernel 
pk=K.pk+nPoints; klen=K.klen+nPoints;

% bandpass and envelope the waveforms 
pfilt=bandpassLocal(y,K.flt_low,K.flt_high,fs,5);  % filter 5th order butter 
penv=abs(hilbert(pfilt)); % envelope of time serires 

% do the cross correation and keep positive lags 
[cc,lg]=xcorr(penv,kern); 
a=find(lg >= 0 & lg < length(penv)-klen); 
cc=cc(a); lg=lg(a);  % only want positive lags 

% normalize the cc ---> ncc 
ncc=nan(1,length(cc)-klen);  % preallocate the normalized cc matric 
if use_par==1  % process in parallel 
    PP=parpool; 
    parfor i=1:length(cc)-klen 
    ncc(i)=cc(i)./(sqrt(sum(kern.^2))*sqrt(sum(penv(i:i+klen-1).^2))); 
    end
    delete(PP); 
else   % use_par = 0 
    for i=1:length(cc)-klen 
    ncc(i)=cc(i)./(sqrt(sum(kern.^2))*sqrt(sum(penv(i:i+klen-1).^2))); 
    end
end

% need to account for the fact that the kernel peak is at pk number of
% points into the time series 
nlg=lg(1:length(ncc)); 
nlg=nlg+pk; nlg=cat(2,1:1:pk-1,nlg);  % adjust the lags based on peak offset 
ncc=cat(2,nan(1,pk-1),ncc);  % make first pk number of points in ncc NaN  

% peak detection 
% CHANGE THIS!!! if 96, 48 play with minpeakdistance (35 for 48?)
[det_score,det_sam] = findpeaks(ncc,'MINPEAKHEIGHT',min_det_score,'MINPEAKDISTANCE',75); 

det_amp=nan(size(det_sam));  % find the amplitde of the snap 
for i=1:length(det_sam)
    det_amp(i)=nanmax(abs(pfilt(det_sam(i)-pk:det_sam(i)+pk)));  % get amplide of waveform 
end

det_tim=det_sam/fs;  % time in sec of the detections 
file_dbrms=20*log10(std(pfilt));  % RMS of the file 

% make local bandpass function
    function [d] = bandpassLocal(c,flp,fhi,Fs,n)
        % Function bandpass applies nth order butterworth filter 
     % [d]=bandpass(c,flp,fhi,Fs,n) 
     % 
     % INPUT 
     % c = input time series 
     % flp = lowpass corner frequency of filter 
     % fhi = hipass corner frequency 
     % Fs = sample rate 
     % n = filter order
     % 
     % OUTPUT 
     % d is the bandpassed waveform. 

     % Del Bohnenstiehl - NCSU 
     % drbohnen@ncsu.edu 
     % part of NCSU's soundscape tools package for MATLAB (bandpass_del) 

     if isempty(n) 
         n=5; 
     end

     fnq=0.5*Fs;  % Nyquist frequency 
     Wn=[flp/fnq fhi/fnq];    % butterworth bandpass non-dimensional frequency 
     [b,e]=butter(n,Wn); % construct the filter 
     d=filtfilt(b,e,c); % zero phase filter the data 
    end

end

