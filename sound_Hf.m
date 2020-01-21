function [meanHf,Hf]=sound_Hf(y,fs,segdur,segovlp,specwl,specovlp,flimits)
%
% sound_Hf calculates both the spectral entropy within a sound
% recording. 
% 
%% Usage: [meanHf, Hf]=sound_H(y,fs,segdur,segovlp,specwl,specovlp,flimits)
% For example: 
%  [meanHf,Hf]=sound_Hf(y,48000,30,0, 1024, 0, [200,22000])
%  for data y sampled at 48kHz divides the data into (possibly multiple) 
%  30 sec long segments with 0 seconds of overlap between segments.  Then 
%  calculate Hf for each segment using a 1024 point STFT with 0 points 
%  overalp within the 200-20,000 Hz spectral band.    
%
%% INPUTS: 
% y - time series of pressure corrected amplitude  (unfiltered)  
% fs - sample rate of the recording in Hz (def = 96000) 
% segdur - length of data segment in seconds used to calculate Hf(def = all data)
% segovlp - overlap in seconds for data segment (default = 0); 
% specwl - window length in # of pts for spectral calcs (def = 512) 
%          specwl should be power-of-two; used as NFFT. 
% specovlp - # points of overlap in spectral calcs (def == 0); 
% flimits = vector of length 2 to select a frequency band (e.g.,[150,20000]).
% if flimits = [] then full bandwidth data are used. 
% 
%% OUTPUT: 
% meanHf - mean of Hf from each of the data segments  
% Hf - vector of spectral Hf for each individual data segment 
%
%% COMMENTS: 
% 1. The time series y is divided into multiple segdur second long windows 
% 2. Spectral calculations use a Hanning window. 
% 3. See Ht for temporal entropy calcultion 
% 4. Review literature in choosing specwl; this controls number of
%    frequency bins (i.e., spectral resolution ~ fs/specwl). 
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
    fs=48000; segdur=(length(y)-1)/fs; segovlp=0; specwl=512; specovlp=0; 
end
if isempty(fs); fs=48000; end 
if isempty(segdur); segdur=(length(y)-1)/fs; end 
if isempty(specwl); specwl=512; end 
if isempty(specovlp); specovlp=0; end 
if isempty(segovlp); segovlp=0; end 

% check that at least one segdur window is present. 
if length(y) < segdur*fs; 
    disp('input vector not segdur seconds long') 
    segdur=(length(y)-1)/fs; segovlp=0; 
    disp(['setting segdur to: ' num2str(segdur)  ' seconds']) 
end

% other warnings 
if (segdur*fs)/specwl < 10; 
disp('WARNING: spectral averaging is using less than 10 NFFT windows\n'); 
end 
if segdur < 5; 
disp('WARNING: segdur variable is perhaps unusually small,less 5 sec?\n'); 
end 

y=y-mean(y); % demean the data 

%% Breaktime series up into as many segdur pieces as possible 
pts=floor(segdur*fs); ovlp=floor(segovlp*fs); 
y_mx=buffer(y,pts,ovlp,'nodelay');   
if y_mx(end)==0; y_mx(:,end)=[]; end  %if last colum is zero padded delete. 
[~,nseg]=size((y_mx)); 


%% now loop through and do calcuation on each column of data 
Hf=nan(1,nseg);
for i=1:nseg  % for each segment 
% Hf calculations 
[s,frqs] = spectrogram(y_mx(:,i)-mean(y_mx(:,i)),hanning(specwl),...
    specovlp,[],fs); % calcualte the spectra 
s(1,:)=[]; frqs(1)=[];  %remove DC term from 
if isempty(flimits)==0 
    keep=frqs > flimits(1) & frqs < flimits(2); 
else 
    keep=frqs > 0 & frqs < fs/2; 
end
s=abs(s); S=mean(s(keep,:),2); % take frequencies of interest and average
                               % across the spectrogram to get mean
                               % spectrum S
sumS=sum(S); % Sum of mean spectrum 
Sf=S/sumS; % Compute the probability mass function for the mean spectrum 
Hf(i)= -sum(Sf .* log2(Sf))/log2(length(Sf)); %Compute the spectral entropy 

end  % for each segment 

%% report out the mean Hf value 
meanHf = mean(Hf); 
