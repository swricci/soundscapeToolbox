function [meanHt,Ht]=sound_Ht(y,fs,segdur,segovlp)
% 
% sound_Ht calculates the temporal entropy within a recording. 
% 
%% Usage: [Ht,Htvec]=sound_Ht(y,fs,segdur,segovlp)
% For example: 
%  [Ht,Htvec]=sound_Ht(y,48000,30,0)
%  for data y sampled at 48kHz divides the data into (possibly multiple) 
%  30 sec long segments with 0 seconds of overlap between segments.  Then 
%  calculate Ht for each segment. 
%
%% INPUTS: 
% y - time series of pressure corrected amplitude.  
% fs - sample rate of the recording in Hz (def = 48000) 
% segdur - duration of data segment in seconds used to calculate Ht 
%         (default = all data)
% segovlp - overlap in seconds for data segment (def = 0); 
% 
%% OUTPUT: 
% meanHt - mean of the Ht values calcualted form each segment of y. 
% Ht - vector of Ht values in each of the segments  
% 
%% COMMENTS: 
% 1. The time series y is divided into multiple segments  
%    of length segdur and overlapping by segovlp
% 2. If the you do NOT want to analyze full bandwidth of the signal 
%    apply filter prior to input. Use dsmpl_bandpass.m (for very high 
%    sample rate data) or the bandpass.m functions.  See their help. 
% 3. If you are interest in calculating the spectral entropy Hf, use
%     sound_Hf, which takes in the unfilter data, but allows the user to 
%     restrict the analysis band with the flimits options [lo, hi].  
%     Often the entropy is given as H=Hf*Ht; in this case use the same 
%     [lo, hi] filter and flimits settings and identical segdur segovlp 
%     values when running sound_Hf and sound_Ht. 
% 4. The signal envelope is calculated using the Hilbert transform; 
%
% 
%% Reference: Sueur J, Pavoine S, Hamerlynck O, Duvail S (2008) 
% Rapid Acoustic Survey for Biodiversity Appraisal.
% PLoS ONE 3(12): e4065. doi:10.1371/ journal.pone.0004065 
% 
%% Del Bohnenstiehl - NCSU 
% Sept 2013; Modified 18 June 2016  
% drbohnen@ncsu.edu 
% part of NCSU's soundscape tools package for MATLAB 

%% checks 
if nargin==1 
    fs=48000; segdur=(length(y)-1)/fs; segovlp=0; 
end
if isempty(fs); fs=48000; end 
if isempty(segdur); segdur=(length(y)-1)/fs;  end 
if isempty(segovlp); segovlp=0; end 

% check that at least one segdur window is present. 
if length(y) < segdur*fs; 
    disp('input vector not segdur seconds long') 
    segdur=floor(length(y)/fs); segovlp=0; 
    disp(['setting segdur to: ' num2str(segdur)  ' seconds']) 
end

% other warnings 
if segdur < 10; 
disp('WARNING: segdur variable is perhaps unusually small, less 10 sec?\n'); 
end 

y=y-mean(y); % demean the data 

%% Break the time series up into as many segdur pieces as possible and store 
pts=floor(segdur*fs); ovlp=floor(segovlp*fs); 
y_mx=buffer(y,pts,ovlp,'nodelay');   
if y_mx(end)==0; y_mx(:,end)=[]; end  % if last colum is zero padded delete. 
[~,nseg]=size((y_mx)); 

%% now loop through and do calcuation on each column of data 
Ht=nan(1,nseg);  % reallocate 
for i=1:nseg  % for each segment 
% Ht calculations 
Xi=abs(hilbert(y_mx(:,i))); % Hilbert transform to get envelop
sumXi=sum(Xi); % Sum the envelope 
At=Xi/sumXi;  % Compute the probability mass function 
Ht(i)=-sum(At .* log2(At))/log2(length(At)); % Compute the temporal entropy 
end  % for each segment 

%% return the mean of Ht 
meanHt=mean(Ht); 

