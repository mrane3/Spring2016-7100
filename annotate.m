

bpm_max = 360 ;
minTime = 60 / (bpm_max*2) ;

windowSize = 1024 ;
hopSize = 128 ;

[file, path] = uigetfile('../DATA_CUT_FINAL/*.*', 'Select Audio File') ;
[data, fs] = audioread([path, file]) ;

% [file, path] = uigetfile('../DATA_CUT_FINAL/*.*', 'Select Notation File') ;
% fp = fopen([path, file]) ;
% notation = textscan(fp, '%s') ;
% notation = notation{1} ;
of=cd('hpss');
[h,p]=s_hpss_IGprior(data,512,5);
cd(of);
onset_times = onset_detect(p, fs, windowSize, hopSize, minTime) ;
onset_times=onset_times/fs;
