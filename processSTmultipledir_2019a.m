
% input and output locations and file names 
% process_ST_2019a(Site,Deployment,  DirIn, DirOut, FS, nsec. generate_jpg_plots, STcalib)
% process_ST_2019a('ESB',2, '/Volumes/P3/FKNMS_soundscape/deployment2/ESB_S9_805892123_dep2/', '/Where/I/Put/My/results_files/', 48000, 120, 0, STcalib)
%
% Site='ESB'; % text string comment with the location of the deployment
%
%
% Deployment= 2; %text comment with deployment information (e.g. study name, study period, deployment number)
% DirIn = '/Volumes/P3/FKNMS_soundscape/deployment7/ESB_S2_67383320_dep7/';  % text string with directory name where wav files are;  Needs trailing / and single quotes 
% DirOut= '/Volumes/P3/FKNMS_soundscape/ProcessST_dep7/ESB_S2_dep7/'; % directory where data matrix and figures are stored;   Needs trailing / and single quotes 
% FS= 48000;  % need to set the preallocation size 
% nsec = 120; % how much data do you want to process per file (this would almost always be 120;  we almost alway record 130 seconds 
% STcalib is the calibration table, which you can load with STcalibration.mat  % loading the master calibration table 
% plots=0; % make them for full nsec of data == 1; make for half the lenth of nsec = 2; otherwise no plots saved 
                      % note that plots will not be displayed to spped up
                      % the processing 

  clear
  
load STcalibration 
plots=1; 
load dir2process.mat    % read in a table with                     

for k=[3,11,13,20,26,33]; 
    Site=char(dir2process.Site(k));   
    DirIn =char(dir2process.DirIn(k));             
    DirOut=char(dir2process.DirOut(k)); 
    Deployment=dir2process.Deployment(k);
    FS=dir2process.FS(k);
    nsec=dir2process.nsec(k); 
    Sgate=datenum(dir2process.Sgate(k))-0.00416;  % 0.0041 days = 10 minutes 
    Egate=datenum(dir2process.Egate(k))+0.00416;
        
% generate file list, file names and diretories 
if exist('DirOut','dir') ~= 1; eval(['system(''mkdir '  DirOut ''')']); end %  make output directory if it does not exist 

[filelist,ftimes, fend]= mktableSTdir(DirIn);  

a=find(ftimes> Sgate & ftimes < Egate); 
filelist=filelist(a); 
size(a)


if FS==48000 
win= 2^14;   % number of points for mean spectrum averaging 
ovlp=0;  % overlap in spectral calcs 
poavg=nan(win/2+1,length(filelist));    % preallocate space  - Power in each freq bin and file 
rms=nan(1,length(filelist));            % preallocate space - RMS pressure in each file  
UTC=nan(1,length(filelist)); 
end

% 
if FS==96000 
win= 2^15;   % number of points for mean spectrum averaging 
ovlp=0;  % overlap in spectral calcs 
poavg=nan(win/2+1,length(filelist));    % preallocate space  - Power in each freq bin and file 
rms=nan(1,length(filelist));            % preallocate space - RMS pressure in each file  
utc=nan(1,length(filelist)); 
end


if FS==144000 
win= 2^16;   % number of points for mean spectrum averaging 
ovlp=0;  % overlap in spectral calcs 
poavg=nan(win/2+1,length(filelist));    % preallocate space  - Power in each freq bin and file 
rms=nan(1,length(filelist));            % preallocate space - RMS pressure in each file  
utc=nan(1,length(filelist)); 
end

%% loop through each wavfile in the DirIn
tic 
for i = 1:length(filelist)  % for each file in the list 

    fprintf('Processing %s\n', char(filelist(i).name)); %display file name of file being processed
    [y,fstart,fs,metadata]=readST(char(filelist(i).name),DirIn,nsec,STcalib);  
    utc(i)=fstart; 
    
% now do your calculations 
[poavg(:,i),f]=sound_MSPEC(y,fs,win,ovlp); 
rms(i)= std(y); % RMS in the time domain of this file 

% %% now the plotting 
 if plots > 0  
     if i <=20 || i > (length(filelist)-20)
 y=y(1:floor(length(y)/plots)); 
 [~,F,T,Pxx]=spectrogram(y,2^12,2^11,[],fs); %0.05sec windows, 0% overlap
 t=(0:1:(length(y)-1))*(1/fs); 
 h=figure('visible','off','Position',[500 900 1000 1300]); % figure('Position',[500 900 700 900]); 
 subplot(4,1,1); plot(t,y); set(gca,'FontSize',12); 
 xlabel('Time (sec)'); ylabel('amp (?Pa)'); grid on; xlim([0 120]);
 set(gca,'xtick',0:5:max(t));
 title(sprintf(filelist(i).name));
% %
 subplot(4,1,2); imagesc(T,F,10*log10(Pxx));set(gca,'FontSize',12); 
caxis([40,90]);  
axis xy; xlabel('Time (sec)'); 
 ylim([0 25000]);
 set(gca,'Ytick',0:5000:25000,'YtickLabel',[0,5,10,15,20,25],'ygrid','on');
 ylabel('Frequency (kHz)'); xlim([0,max(t)]); set(gca,'xtick',0:5:max(t));
 title('25000 Hz Spectrogram (colourbar = 45-90 dB)');
 
% %
subplot(4,1,3); imagesc(T,F,10*log10(Pxx));set(gca,'FontSize',12);   
caxis([40,90]); colormap jet;
axis xy; xlabel('Time (sec)');
ylim([0 3000]);
set(gca,'Ytick',0:500:3000,'YtickLabel',[0,500,1000,1500,2000,2500,3000],'ygrid','on');
ylabel('Frequency (Hz)'); xlim([0,max(t)]); set(gca,'xtick',0:5:max(t));
title('3000 Hz Spectrogram (colourbar = 45-90 dB)');
%
subplot(4,1,4); semilogx(f,10*log10(poavg(:,i))); xlim([50,27000]); 
set(gca,'FontSize',12); 
ymin=round2(nanmin(10*log10(poavg(:,i)))-2.5,5); 
ylim([ymin,ymin+65]); set(gca,'YTick',ymin:10:ymin+55);
set(gca,'XTick',[100,1000,10000],'XTickLabel',[100,1000,10000]);
xlabel('Frequency (Hz)'); ylabel('dB  (?Pa^2/Hz)'); grid on; 


% 

out_img_name=strcat(DirOut,Site, '_D',sprintf('%02.0f',Deployment),'_', datestr(utc(i),'yyyymmdd_HHMMSS'),'.png' ); 
fig_title=strcat(Site, ' D',sprintf('%02.0f',Deployment),' : ', datestr(utc(i),'yyyymmdd_HHMMSS'),' fn:',filelist(i).name );
title(fig_title); 
 saveas(h,out_img_name);
  close all 
     end
 end
 
end

toc 
eval([Site '_' sprintf('%02.0f',Deployment) '_poavg=poavg;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_rms=rms;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_utc=utc;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_filelist=filelist;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_metadata=metadata;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_f=f;']);

eval(['save ',DirOut, [Site '_' sprintf('%02.0f',Deployment) '_results.mat'] ' ' [Site '_' sprintf('%02.0f',Deployment) '_poavg']...
    ' ' [Site '_' sprintf('%02.0f',Deployment) '_rms'] ' ' [Site '_' sprintf('%02.0f',Deployment) '_utc']...
    ' ' [Site '_' sprintf('%02.0f',Deployment) '_filelist'] ' ' [Site '_' sprintf('%02.0f',Deployment) '_metadata']...
    ' ' [Site '_' sprintf('%02.0f',Deployment) '_f']]);    

clear poavg rms utc 
clear LKU* LKP* SDK* WSB* WDR* ESB* N1M* NFS* 
end


