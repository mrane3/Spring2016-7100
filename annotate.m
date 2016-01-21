

bpm_max = 360 ;
minTime = 60 / (bpm_max*2) ;

windowSize = 1024 ;
hopSize = 128 ;

[file, path] = uigetfile('../DATA_CUT_FINAL/*.*', 'Select Audio File') ;
[data, fs] = audioread([path, file]) ;

[file, path] = uigetfile('../DATA_CUT_FINAL/*.*', 'Select Notation File') ;
fp = fopen([path, file]) ;
notation = textscan(fp, '%s') ;
notation = notation{1} ;

onset_times = onset_detect(data, fs, windowSize, hopSize, minTime) ;

