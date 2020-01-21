%% CALL SCRIPT FOR TOADFISH DETECTOR 
% version Sept. 14 2016 
% S. Ricci, D. Bohnenstiehl 
% NC State University 

clear; close all
%% Shannon, FYI THESE ARE THE SAVED VARIABLES 
% filename det_dt det_freq1 det_freq2 det_pks
% changed slightly for standarization 

%% file and directory settings 
allfiles=dir('/Volumes/P2/d_shannon_brown/VAsoundscapes_may2015/st7_eaglepoint/*.wav');
filelist={allfiles.name};
filelist=filelist(2892:2893);%seths point(405:2892 eaglepoint 405:2893 change 411:2899 lodges 210:1453 millpoint 412:2900 walnut 204:1448 rabbitislandeast 405:2893 littleneck 208:1452
filedir='/Volumes/P2/d_shannon_brown/VAsoundscapes_may2015/st7_eaglepoint/' ;%directory of wav files
DirOut='/Users/sbrown/Documents/MATLAB/toadfiles/harriscreek_toadanalysis_v2/st7_eaglepoint'; 
fileout='d_toadanalysis_st7_eaglepoint2.mat'; 

%% input run parameter settings
s= 10; %standard deviation of harmonic for detection kernel (Hz)
sweep = 3; % sweep applied to F1 harmonic (hz) 
F1range= [135 285];   % min and max range of the first harmonic 
                      % The actual bin centers are selected based on 
                      % the frequency resolution of the spectrogram 
ploton =1; %ploton=1 plots generated, ploton=0 plots not generated
thres = 4.75;  %detection threshold 
STnum=7;   % soundtrap number 
G= 'lo';   % gain setting 


%% set the correct sensitivity
if STnum==1; lowgain=185.4; highgain=173; end %serial: 67670026
if STnum==2; lowgain=187.3; highgain=175.3; end %serial: 67383320
if STnum==3; lowgain=187.9; highgain=175.9; end %serial: 67379211
if STnum==4; lowgain=186.1; highgain=173.7; end %serial: 67690520 %changed 6/10/15
if STnum==5; lowgain=184.7; highgain=173.1; end %serial: 67903499 %changed 6/10/15
if STnum==6; lowgain=185.4; highgain=173; end %serial: 67657752
if STnum==7; lowgain=184.4; highgain=172; end %serial: 67641354
if STnum==8; lowgain=184.9; highgain=172.5; end %serial: 67399703
if STnum==9; lowgain=182.2; highgain=169.5; end %serial: 805892123

if G == 'lo'; cal=10^(lowgain/20);end
if G == 'hi'; cal=10^(highgain/20);end


det_time=[];
calls_time=[];
det_dt=[];
det_freq1=[];
det_freq2=[];
filename=[];
det_pks=[];


for i=1:length(filelist)
    fprintf('Processing %s\n', char(filelist(i)));
[y,fs] = audioread(char(strcat(filedir,filelist(i)))); % read in wav file
    sec2skip = 3; 
    y=y(sec2skip*fs:fs*floor(length(y)/fs)); %floor to nearest whole second e.g., 60 , cut off first 3 seconds
    y=(y-mean(y))*cal;% demeaned & response correct to uPa 
    t=(1/fs)*(0:1:length(y)-1);  % time axis in seconds 

    %get start time of file in datenum format
    TMP=char(filelist(i)); %gets .wav file name
    tstarts(1,1:6)=[str2num((TMP(:,10:11))),str2num(TMP(:,12:13)),str2num(TMP(:,14:15)), str2num(TMP(:,16:17)), str2num(TMP(:,18:19)), str2num(TMP(:,20:21))]; 
    tstarts(1,1)=tstarts(1,1)+2000; %make year correct
    tstarts_dt=datenum(tstarts);
    
    if STnum==9; %st9 has more digits in serial number
    tstarts(1,1:6)=[str2num((TMP(:,11:12))),str2num(TMP(:,13:14)),str2num(TMP(:,15:16)),str2num(TMP(:,17:18)),str2num(TMP(:,19:20)),str2num(TMP(:,21:22))];
    tstarts(1,1)=tstarts(1,1)+2000;
    tstarts_dt=datenum(tstarts);
    end

 fprintf(['Start time of this file is :  ' datestr(tstarts_dt)  '\n']) 
    
%%run toadfish detector function, output is list of detection times
[DT,DF1,DF2,pks]=toadfish_detector96kHz_v2(y,F1range,s,sweep,thres,ploton); 

%totcalls(i)=length(DT);
if isempty(DT)==1
    DT=NaN;
    DF1=NaN;
    DF2=NaN;
    pks=NaN;
end

file=repmat(filelist(i),1,length(DT));


det_dt=cat(1,det_dt,datevec(tstarts_dt+(sec2skip+DT)/86400));   % need to add in the seconds we skipped 
det_freq1=cat(1,det_freq1,DF1');
det_freq2=cat(1,det_freq2,DF2');
det_pks=cat(1,det_pks,pks');
filename=cat(1,filename,file');

if STnum==9
saveas(gcf, [TMP(1:22), '.fig']); %change to 1:22 for ST 9 1:21 for all others
else 
saveas(gcf, [TMP(1:21), '.fig']); %change to 1:22 for ST 9 1:21 for all others
end
close all %S.Ricci added 91416, closes fig file
end


cd(DirOut)

eval(['save( ''' char(fileout) ''', ''filename'',''det_dt'',''det_freq1'',''det_freq2'',''det_pks'' );']) 