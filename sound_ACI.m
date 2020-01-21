function  [meanACI, ACI, ACIfreq, F] = sound_ACI(y, fs, segdur, segovlp, specwl, flimits)
% 
% Acoustic Complexity Index - ACI 
%
% USE: 
% [meanACI, ACI, ACIfreq, F] = sound_ACI(y, fs,segdur, segovlp, specwl, flimits)
%
% Input:
% y = input time series in µPa (unfiltered) 
% fs=sampling frequency of wave (in Hz).
% segdur - length of data segment in seconds used to calculate ACI (def = all data)
% segovlp - overlap in seconds for data segment (default = 0); 
% specwl = FFT window length for the analysis (power of 2: e.g., 512, 1024, 2048) 
% flimits = vector of length 2 to select a frequency band (in Hz [0,FS/2]).
% 
% Output:
% meanACI - Acoustic Complexity Index (ACI) average over each data segment
% ACI - ACI within each data segment 
%        * commonly ACI is defined as the sum of these values when multiple 
%            data segments are used * 
% ACIfreq - ACI as a function of frequency F for each data segment  
% F - frequencies corresponding to rows of ACIfreq. the frequency resolution can 
% be determined by running: diff(F)
%
%  NOTES: 
% 1. Hanning window is used for spectral calculations 
% The ACI is calculated from the STFT spectrogram matrix, where temporal steps 
% are represented by rows and and frequency bins by the columns. 
% 2. steps in calculation of ACI: 
%   I. Calculate absolute difference between adjacent values of intensity
%       dk = | Ik-Ik+1| within a given frequency bin,  where k is a time step
%   II. Sum dk over the data segment within that frequency bin 
%   III. Normalize sum(d)/sum(Ik) within that frequency 
%   IV. Total ACI is then found by summing sum(d)/sum(Ik) in each frequency 
%       bin over the range of interest. 
%
% General refernce for ACI
% Pieretti, N., Farina, A., & Morri, D. (2011). A new methodology to infer 
% the singing activity of an avian community: The Acoustic Complexity Index
% (ACI). Ecological Indicators, 11(3), 868?873. 
% doi:10.1016/j.ecolind.2010.11.005
% 
% D. Bohnenstiehl & A. Lillis 
% NC State University  
% Feb. 2015; modified 18 June 2016 
% drbohnen@ncsu.edu 
% part of NCSU's soundscape tools package for MATLAB 


% define some defaults 
if nargin < 6
    error('Must define y, fs, segdur, segovlp, specwl, flimits')   
end
%

if length(y) <= specwl*2; 
   error('data series y must be at least twice the FFT specwl length')
end
% 
% check that at least one segdur window is present. 
if length(y) < segdur*fs; 
    disp('input vector not segdur seconds long') 
    segdur=(length(y)-1)/fs;  segovlp=0; 
    disp(['setting segdur to: ' num2str(segdur)  ' seconds']) 
end

%% Breaktime series up into as many segdur pieces as possible 
pts=floor(segdur*fs);ovlp=floor(segovlp*fs);
y_mx=buffer(y,pts,ovlp,'nodelay');   
if y_mx(end)==0; y_mx(:,end)=[]; end  %if last colum is zero padded delete. 
[~,nseg]=size((y_mx)); 
%disp(['number of data segments is: ' num2str(nseg) ]) 
ACI=[]; ACIfreq=[]; 

% calculate spectrogram 
for i=1:nseg
[I,F,~,~]=spectrogram(y_mx(:,i)-mean(y_mx(:,i)),hanning(specwl),0,specwl,fs);  
I(1,:)=[]; F(1)=[];  %remove DC term from STFT
I=abs(I); 
% difference in intensity in each freq bin from one time step to next
diffI=abs(diff(I,1,2));  
% sum of intensities in each freq bin
sumI=sum(I,2);  
% sum of differences between freq bins  
sumdiffI=sum(diffI,2);  
% ACI in each frequency bin 
ACIfreq(:,i)=sumdiffI./sumI;  
% find frequency bins over which to calculate ACI
if isempty(flimits) == 0 
Fl = F>flimits(1) & F<=flimits(2); 
else
    Fl=1:1:length(F); 
end
% sum those frequency bins to calcualte the ACI for that time window 
ACI(i)=sum(ACIfreq(Fl,i)); 

end 

meanACI=mean(ACI);



