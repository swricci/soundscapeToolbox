function [filelist,fstart, fend]= mktableSTdir(DirIn) 

eval(['filelist=dir(''' DirIn '*.wav'')']) 
N=length(filelist); 
fs=nan(N,1); 
g=cell(1); 
fstart=nan(N,1);
fend=nan(N,1);

%% this block reads the associated .log.xml file to get info  
for i=1:N
WAVFILE=filelist(i).name; 
C = strsplit(WAVFILE, '.'); % split text string 
SN=C{1}; % serial number string from file name 
FT=C{2}; % time string from file name  -  time base unknown 
LF=strcat(DirIn,SN,'.',FT,'.','log.xml'); % log file name 
logf=fileread(LF);  

% get times in UTC 


% get times in UTC 
ETind1=strfind(logf, 'SamplingStopTimeUTC')+21;   
tmp=strfind(logf, '"/>'); 
ETind2=tmp(find(tmp > ETind1, 1,'first'))-1; 

ETstring=logf(ETind1:ETind2);  
ETstring=replace(ETstring,"T"," "); 
fend(i)=datenum(ETstring);  

STind1=strfind(logf, 'SamplingStartTimeUTC')+21;    
STind2=tmp(find(tmp > STind1, 1,'first'))-1;  
STstring=logf(STind1:STind2);  
STstring=replace(STstring,"T"," "); 
fstart(i)=datenum(STstring); 


% get the gain setting 
Gind1=strfind(logf, 'AUDIO G')+12; 
g{i}=logf(Gind1:Gind1+1); 


% get the sample rate 
FSind1=strfind(logf, '<FS>')+4; tmp=strfind(logf, '</FS>')-2; 
FSind2=tmp(find(tmp > FSind1, 1,'first')); 
fs(i)=str2double(logf(FSind1:FSind2)); 

end

%figure; plot(1:1:length(fend),(fend-fend(1)),'.')


