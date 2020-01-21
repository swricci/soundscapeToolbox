function  [det_time,det_freq1,det_freq2,det_score]=toadfish_detector96kHz_v2(y,Frange,s,sweep,thres, ploton)
%
% toadfish detector - optimized for Harris Creek 
% version Sept. 14 2016 
% S. Ricci, D. Bohnenstiehl 
% NC State University 
% 
% INPUTS 
% y = demeaned pressure corrected time series 
% Frange = is a a x 2 vector with [minF1  maxF1] (range of fundamental
% frequency kernels to run 
% sweep is the change in frequency (Hz) over the duration of the call 
% s is the width or std of call in Hz 
% ploton = 1 will plot detections and spectrograms 
%
% HARDWIRED Settings 
% FFT uses 2^14 points zero padded to 2^15 to increase frequency resolution
% delta-F = 2.9287 sec;   delta-T = 341ms; 
% kernel length is 7 * delta-T = 
%
% OUTPUTS 
% det_time = time in seconds within the time seroes y 
% det_freq1 = frequency of first harmonic 
% det_freq2 = frequency of second harmonc = 2*det_freq1; 
% det_score = amplitude of the detection score 
% 


tic 

%% calculate PSD in small windows
W=2*8192;  fs=96000;   np2=nextpow2(W)+1;  % zero pad to next 
[~,F,T,Pxx]=spectrogram(y,W,floor(W*.8),2^np2,fs);  
a=find(F >=85 & F <= 620);  F=F(a); Pxx=Pxx(a,:);   % trim the frequency band 
[rowsF,~]=size(F); % number of frequencies in the spectrogram 

%% determine the harmonics for the kernels based on F1range 
F1a=F(F> Frange(1) & F < Frange(2));  F1b=F1a-sweep; F1=[F1a'; F1b']; 
F2=F1*2;  % define 2nd harmionic 
[~,nkerns]=size(F1);  % determine the number of kernels 


%% adjust the spectrogram by removing a smooth version of the mean spectra and normalizing 
PxxMod=10*log10(Pxx); 
PxxMod=PxxMod-(smooth(mean(PxxMod,2),34)*ones(1,length(T)));  
a= PxxMod > 68;   PxxMod(a)=68; b= PxxMod < -68;   PxxMod(b)=-68; 
PxxMod=PxxMod/68;  % rescale approximately -1 to 1 

cmatrix=nan(nkerns,length(T)+6);   % length of T + num rows (k) -1; 

for j=1:nkerns  % make each kernel 
%% now generate kernels 
    % upper one 
fo=F2(1,j); f1=F2(2,j); 
t=0:diff(T(1:2)):0.225; d=max(t);  % works out to be 7 * delta-T long 
k1=nan(length(F),length(t)); 
for i=1:length(t);  % for each time step 
x=F-(fo+(t(i)./d)*(f1-fo));
k1(:,i)=(1-(x.^2)/(s.^2)).*exp(-(x.^2)/(2*s^2)); 
end

% lower one 
fo=F1(1,j);  f1=F1(2,j); 
t=0:diff(T(1:2)):0.225;  d=max(t); 
k2=nan(length(F),length(t)); 
for i=1:length(t); 
x=F-(fo+(t(i)./d)*(f1-fo));
k2(:,i)=(1-(x.^2)/(s.^2)).*exp(-(x.^2)/(2*s^2)); 
end

k=k1+k2;  % add the result together to make the final c
[~,len_k]=size(k);  % how many time step in k 

%% now cross correlate with data spectrogram 
[c]=xcorr2(PxxMod,k);    
c=c(rowsF,:);  % pulls the center value that represents only a time shift 
cmatrix(j,:)=c;   % stores it as a row 

end   % repeat for each of the nkerns 

[cout,iout]=max(cmatrix(:,len_k:end));  % returns the max value in cout & witch row gave the max iout

[det_score,LOCS]=findpeaks(cout,'MinPeakHeight',thres,'MinPeakDistance',4,'MinPeakProminence',2); 
det_time=T(LOCS); 
det_freq1=F1(1,iout(LOCS));
det_freq2=F2(1,iout(LOCS));

%% if plotting is on 
if ploton==1; 
figure; ax(1)=subplot(2,1,1); imagesc(T,F,10*log10(Pxx));  axis xy; hold on; 
colormap('jet'); 
plot(T-diff(T(1:2))/2,200+cout*8,'k','LineWidth',2); grid on; 
plot(T(LOCS)-diff(T(1:2))/2,det_freq1,'ok','MarkerSize',4,'LineWidth',1);
plot(T(LOCS)-diff(T(1:2))/2,det_freq2,'ok','MarkerSize',4,'LineWidth',1);
ax(2)=subplot(2,1,2);  plot(T-diff(T(1:2))/2,cout); 
hold on; plot(T(LOCS)-diff(T(1:2))/2,det_score,'or'); grid on; 
for j=1:length(LOCS)
text(T(LOCS(j)),det_score(j),char(sprintf('%3.1f\n',det_score(j))),'FontSize',14); 
end

linkaxes([ax(2),ax(1)],'x'); 

end

toc


