function [mspec,f]=sound_MSPEC(y,fs,win,ovlp)
% 
% The function sound_MSPEC returns the mean spectra of the data
% in y using multiple spectral windows win points in length.  
% 
%% INPUT 
% y = pressure corrected time series uPa
% fs = sample rate 
% win = number of data points in each window 
%        this should be a power of 2 (making win = NFFT)
%        For example: If you want to average ~1/4 sec duration windows use 
%        2^netpow2(0.25*fs) 
% ovlp = # points of overlap in FFT calcs; must to < win  
% 
%% OUTPUT 
% mspec = mean power spectrum in uPa^2/Hz 
%         Use 10*log10(mspec) to convert to dB. 
% f = frequency bins corresponding to mspec 
% 
%% NOTES:
% 1. Spectrum is scalculated using Hanning window 
% 2. The spectrum is scaled such that the sum of the spectral bins is 
%   equal to the variance (power) in the timeseries data.  So one can 
%   sum the mspec over a range of frequencies to get an estimate of the 
%   sound pressure level in that band. 
%
%       For example to find SPL in the 100-20000 Hz band: 
%       a= find(f > 100 & f < 20000); dBspl= 10*log10(sum(mspec(a))); 
% 3.  Returned spectra is in uPa^2/Hz (NOT dB)
% 
%% Del Bohnenstiehl - NCSU 
% Modified 18 June 2016  
% drbohnen@ncsu.edu 
% part of NCSU's soundscape tools package for MATLAB 

%% Defaults 
if nargin < 4 
  disp('not enough arguemnts; NEED: y,fs,win,ovlp') 
   return 
end
if isempty(fs); fs=48000; disp('WARNING - Using 48000 Fs');  end 
if isempty(win); win= 2^14; end 
if isempty(ovlp); ovlp=0; end 


%% Break the time series up into as many segdur pieces as possible and store 
y=y-mean(y);  % demean the time series 
x=buffer(y,win,ovlp,'nodelay');   
if x(end)==0; x(:,end)=[]; end  % if last colum is zero padded delete. 
[r,Nwin]=size((x)); 
msub=repmat(mean(x,1),r,1); % calculate the mean of each column 
x=x-msub; % remove the mean of each column

%% now calculate the average spectra 
%wo=ones(floor(length(x)),1); % rectangular window 
wo=hanning(r); % hanning window 
zo=x.*repmat(wo,1,Nwin); % apply window 
nfft=2^nextpow2(win); % next power of 2, although it should be already  
Y=fft(zo,nfft,1); % two sided FFT opperating on each column 
po=2*abs(Y).^2/(nfft*sum(wo.^2)); % scale for PSD accounting for windowing 
po=po(1:ceil(nfft/2)+1,:); po(1)=0; % take first 1/2 & zero DC 
[prows,~] = size(po); % # rows in po. 
m=0:1:prows-1; f=m*fs/nfft; % define the frequency axis
mspec=mean(po,2); % average spectra 

