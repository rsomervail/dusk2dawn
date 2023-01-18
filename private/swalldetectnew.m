function [swa_results]=swalldetectnew(datainput,orig_fs,thramp)
%% Slow Wave Detection
% Script modified on 2023.01.14 by Richard Somervail to remove Detection window popup for Dusk2Dawn plugin for EEGLAB
% Script modified on 2016.03.31 for Visual Deprivation Project (CIRS)
%
% Half-waves with duration 0.25-1.0s are used. A negative amplitude
% threshold may be defined by the user (however, a zero threshold is
% advised for most analyses).

disp(['**** Slow Wave Detection ****']);

 %% Baseline Correction

allavg=nanmean(datainput,2);
dataref=datainput-repmat(allavg,[1,size(datainput,2)]);
LoadedEEG.data=dataref; clear dataref;  

%% Filter Definition

fs=128; %sampling rate changes for decimated signal
Wp=[0.5 4.0]/(fs/2); % Filtering parameters
Ws=[0.1 10]/(fs/2); % Filtering parameters
Rp=3;
Rs=10;
[n, Wn]=cheb2ord(Wp,Ws,Rp,Rs);
[bbp,abp]=cheby2(n,Rs,Wn); % Loses no more than 3 dB in pass band and has at least 10 dB attenuation in stop band
clear pass* stop* Rp Rs W* n;

%% Detection Starts Here

clear datapoints swa_results channels datax signal dataff EEGder EEG difference ;
clear pos_index neg_index troughs peaks poscross negcross wndx b bx c cx nump maxb maxc lastpk;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([' Analyzing ', num2str(size(LoadedEEG.data,1)),' channels']);
% h = waitbar(0,'Detection Progress');
for i=1:size(LoadedEEG.data,1) % (i) is the number of the channel
%     waitbar(i/size(LoadedEEG.data,1),h,'Detection Progress')
    % Data Extraction, Resample and Filtering
    datax = squeeze(LoadedEEG.data(i,1:size(LoadedEEG.data,2),1));
	signal = resample(double(datax),fs,orig_fs);
    EEG=filtfilt(bbp, abp, signal);
    datapoints=length(EEG);
    channels(i).datalength=datapoints;
    
    %% Search Negcross and Poscross
    pos_index=zeros(length(EEG),1); % Create an empty list for positive peaks
    pos_index(find(EEG>0))=1; % Index of all positive points for EEG
    difference=diff(pos_index); % Locate the first positive and negative point (in time) for each series of consecutive points
    poscross=find(difference==1); % Return the position of all first positive points
    negcross=find(difference==-1); % Return the position of all first negative points
%     EEGder=meanfilt(diff(EEG),5); % Meanfilt is a function that uses a 5 sample moving window to smooth derivative
    EEGder=movmean(diff(EEG),5); % RS edit: replacing missing 'meanfilt' function with the built-in 'movmean' 
    pos_index=zeros(length(EEGder),1); % Repeat above procedure on smoothed signal (?)
    pos_index(find(EEGder>0))=1; % Index of all positive points above minimum threshold %%!!!!!!!!!!!!!!!!set to 0!!!!!!!!!!!!!
    difference=diff(pos_index); % Locate first positive and negative points
    peaks=find(difference==-1)+1; % Find pos ZX and neg ZX of the derivative (peaks)
    troughs=find(difference==1)+1; % Find pos ZX and neg ZX of the derivative (troughs)
    peaks(EEG(peaks)<0)=[]; % Rejects peaks below zero 
    troughs(EEG(troughs)>0)=[]; % Rejects troughs above zero
    
    %% Makes negcross and poscross same size to start
    if negcross(1)<poscross(1); 
            start=1;
        else
            start=2;
    end; 
    
    if start==2;
            poscross(1)=[];
    end;
        lastpk=NaN; % Way to look at Peak to Peak parameters if needed
        ch=i;
                
%% Wave parameters initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	channels(ch).negzx{1}=[];
	channels(ch).poszx{1}=[];
	channels(ch).wvend{1}=[];
	channels(ch).negpks{1}=[];
	channels(ch).maxnegpk{1}=[];
	channels(ch).negpkamp{1}=[];
	channels(ch).maxnegpkamp{1}=[];
	channels(ch).pospks{1}=[];
	channels(ch).maxpospk{1}=[];
	channels(ch).pospkamp{1}=[];
	channels(ch).maxpospkamp{1}=[];
	channels(ch).mxdnslp{1}=[];
	channels(ch).mxupslp{1}=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Locate Peaks
	wvi=1;
	for wndx=start:length(negcross)-1
% 	disp([num2str(wndx)]);
    wavest=negcross(wndx); % Only used for neg/pos peaks
    wavend=negcross(wndx+1); % Only used for neg/pos peaks
        
%     mxdn=abs(min(meanfilt(diff(EEG(wavest:poscross(wndx))),5)))*fs; % Matrix (27) determines instantaneous 1st segement slope on smoothed signal
%     mxup=max(meanfilt(diff(EEG(wavest:poscross(wndx))),5))*fs; % Matrix (28) determines for 2nd segement
    % RS edit: replacing missing 'meanfilt' function with the built-in 'movmean' 
    mxdn=abs(min(movmean(diff(EEG(wavest:poscross(wndx))),5)))*fs; % Matrix (27) determines instantaneous 1st segement slope on smoothed signal
    mxup=max(movmean(diff(EEG(wavest:poscross(wndx))),5))*fs; % Matrix (28) determines for 2nd segement
    
    negpeaks=troughs(troughs>wavest&troughs<wavend);
    wavepk=negpeaks(EEG(negpeaks)==min(EEG(negpeaks))); % Locate Wave Peak
    
    pospeaks=peaks(peaks>wavest&peaks<=wavend);

    if isempty(pospeaks)
        pospeaks=wavend; 
    end; %if pospeaks is empty set pospeak to pos ZX
                        
    %% period=wavend-wavest; %matrix(11) /fs
    poszx=poscross(wndx); %matrix(10)
    b=EEG(negpeaks);  
    bx=negpeaks;
    c=EEG(pospeaks);
    cx=pospeaks;
    nump=length(negpeaks); %matrix(24)

    %% Location of the max neg peak amp
    maxb=min(EEG(negpeaks)); % Max neg peak amp
    if maxb>0
        maxb=maxb(1);
    end;            
    maxbx=negpeaks(EEG(negpeaks)==maxb); % Location of the max neg peak amp

    %% Location of the max pos peak amp
    maxc=max(EEG(pospeaks)); % Max pos peak amp               
    if maxc>0
        maxc=maxc(1);
    end;            
    maxcx=pospeaks(EEG(pospeaks)==maxc); % Location of the max pos peak amp

    lastpk=maxcx;

	waveamp=abs(single(maxc))+abs(single(maxb)); %GBtest (peak2peak amplitude)
%     waveamp=abs(single(maxb)); %GBtest (neg peak amplitude) 
    wavelength=abs((single(wavest)-single(poszx))./fs); %GBtest (length)
	if wavelength>0.25 && wavelength<1.0 %GBtest (length)
        if waveamp>thramp %GBtest (p2p amplitude)
    
%% Wave parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    channels(ch).negzx{wvi}=round((single(wavest)./128).*orig_fs);
    channels(ch).poszx{wvi}=round((single(poszx)./128).*orig_fs);
    channels(ch).wvend{wvi}=round((single(wavend)./128).*orig_fs);
    channels(ch).negpks{wvi}={round((single(bx)./128).*orig_fs)};
    channels(ch).maxnegpk{wvi}=round((single(maxbx)./128).*orig_fs);
    channels(ch).negpkamp{wvi}={single(b)};
    channels(ch).maxnegpkamp{wvi}=single(maxb);
    channels(ch).pospks{wvi}={round((single(cx)./128).*orig_fs)};
    channels(ch).maxpospk{wvi}=round((single(maxcx)./128).*orig_fs);
    channels(ch).pospkamp{wvi}={single(c)};
    channels(ch).maxpospkamp{wvi}=single(maxc);
    channels(ch).mxdnslp{wvi}=single(mxdn);
    channels(ch).mxupslp{wvi}=single(mxup);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     wvi=wvi+1;

        end; % (amplitude)
    end; % (duration)
    
    clear wavest wavend poszx bx maxbx b maxb cx maxcx x maxc mxdn mxup nump negpeaks pospeaks;
    clear wavelength waveamp;
    
     end; %end wndx loop
     clear EEGder dataff signal;
end; %end channel loop

%% Save Output

swa_results.channels=channels;
% close(h)

disp(['**** Detection Completed ****']);