
function  [y,fstart_UTC, fs, metadata]=readST(WAVFILE,DIR,NSEC,STcalib)  
% 
%  [y,fstart_UTC,fs, metadata]=readST(WAVFILE,DIR,NSEC,STcalib) 
%  WAVFILE should be in single quotes 
% '67383320.170714182002.wav' 
%  DIR shold be in single quotes 
%  with or without the trailing / 
% '/Users/del/stuff/deployment1/' or '/Users/del/stuff/deployment1'
%
% This function reads in NSEC of data from 
% a ST201, 202, or 300 WAVFILE in directory DIR and return 
% a gain corrected waveform. 
%
% The function starts reading 6 seconds into the file 
% to avoid possible calinration tone. NSEC is the number of second 
% of data you want to read in must be less than the length of the 
% datafile - 6 seconds 
% 
% STcalib is a matlab table file that has variables 
% low_gain  high_gain  SN  on_date off_date  
% this table can be loaded using " load STcalibration.mat " 
% Be sure to use the latest STcalibration file the soundscape team drive. 
% 
% gain setting (high vs. low), sample rate and time in UTC are determined 
% automatically from the .log.xml file associated with each WAVFILE 
% 
%
% OUTPUT 
% y is the gain correct and demeaned waveform 
% fstart_UTC is the sample rate in UTC 
% fs is the sample rate 
% meta_data is a structure with sample rate, gain setting, gain and
% calibration that were used to correct y. 
% 
% original version 23 January 2019 
% updated  06 March 2019 


% if the directory doesn't have a trailing zero, then add it. 
if strcmp(DIR(end),'/')==0
    DIR=strcat(DIR,'/'); 
end

%% this block reads the associated .log.xml file to get info  
C = strsplit(WAVFILE, '.'); % split text string 
SN=C{1}; % serial number string from file name 
FT=C{2}; % time string from file name  -  time base unknown 
LF=strcat(DIR,SN,'.',FT,'.','log.xml'); % log file name 
logf=fileread(LF);  

% get times in UTC 
ETind1=strfind(logf, 'SamplingStopTimeUTC')+21;   
tmp=strfind(logf, '"/>'); 
ETind2=tmp(find(tmp > ETind1, 1,'first'))-1; 

ETstring=logf(ETind1:ETind2);  
ETstring=replace(ETstring,"T"," "); 
fend=datenum(ETstring);  

STind1=strfind(logf, 'SamplingStartTimeUTC')+21;    
STind2=tmp(find(tmp > STind1, 1,'first'))-1;  
STstring=logf(STind1:STind2);  
STstring=replace(STstring,"T"," "); 
fstart=datenum(STstring); 

% get the gain setting 
Gind1=strfind(logf, 'AUDIO G')+12; 
g=logf(Gind1);  

if isempty(g) 
    g='H'; 
     disp('warning: could read gain from logfile - defaulting to H; this is normal for coninuous or 100% duty cycle recordings')
else
g=g(1); 
end


if strcmp(g,'H')  
gainset=1; 
elseif strcmp(g,'L')
    gainset=0; 
else
    disp('could read gain from logfile, defaulting to H - cause unknown') 
    gainset=1;
end

% get the sample rate 
FSind1=strfind(logf, '<FS>')+4; tmp=strfind(logf, '</FS>')-2; 
FSind2=tmp(find(tmp > FSind1, 1,'first')); 
fs=str2double(logf(FSind1:FSind2)); 

%% now read the waveform 
info=audioinfo(strcat(DIR,WAVFILE)); 
if info.Duration <= NSEC+6 
   y=audioread(strcat(DIR,WAVFILE)); 
   y=y(fs*6:end);  
   fprintf('warning: the data file was short: Total length before trimming 6 seconds was: %1f\n', info.Duration);  
else 
   y=audioread(strcat(DIR,WAVFILE),[fs*6,(fs*(NSEC+6))-1]);
end 


%% now lookup the gains using the gainsetting 
a=find(STcalib.SN==str2double(SN) & datenum(STcalib.on_date) < fstart & datenum(STcalib.off_date) > fend); % row in table with valid calibration data 
if gainset==1;  gain=STcalib.high_gain(a); 
else 
    gain=STcalib.low_gain(a); 
end 
y=(y-mean(y))*(10^(gain/20)); 

metadata.sn = SN; 
metadata.gainsetting = g;
metadata.gain = gain;
metadata.calibration =10^(gain/20);
metadata.fs =fs;
fstart_UTC=fstart; 
