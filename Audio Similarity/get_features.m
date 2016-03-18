function featureMatrix = get_features(featName)

[x,fs]=audioread('/Users/milaprane/Dropbox/wav/ben_20.wav');
x = mean(x,2) ;
windowSize = fs; 
hopSize = 1/4 * windowSize ;

[X,~,~] = spectrogram(x,windowSize,windowSize-hopSize,windowSize,fs) ;
 X = abs(X)*2/windowSize;

switch(featName)
    case 'pitchChroma'
       featureMatrix = getPitchChroma(X, fs) ;
    case 'mfcc'
       featureMatrix = FeatureSpectralMfccs(X, fs) ; 
    otherwise
        error('get_features:NoSuchFeature', 'Feature should be either ''pitchChroma'' or ''mfcc''') ;
end
% [mfcc] = FeatureSpectralMfccs(X, fs) ;
% 
% featureMatrix = [pitchChroma; mfcc] ;`