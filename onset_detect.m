function onset_times = onset_detect(x, fs, windowSize, hopSize, minTime)

ampThreshold = 0.0005 ;

windowStart = 1 ;
windowEnd = windowSize ;

numSamples = length(x) ;
f0 = [] ;
timeInSec = [] ;

noiseFlag = 1 ;

frame1 = x(windowStart:windowEnd) ;
prev2FFT = fft(frame1) ;
if(max(frame1) > ampThreshold)
    noiseFlag = 0 ;
    startTime = windowStart ;
end 
windowStart = windowStart + hopSize ;
windowEnd = windowEnd + hopSize ;

frame2 = x(windowStart:windowEnd) ;
prevFFT = fft(frame2) ;
windowStart = windowStart + hopSize ;
windowEnd = windowEnd + hopSize ;
if(max(frame2) > ampThreshold)
    noiseFlag = 0 ;
    startTime = windowStart ;
end

index = [];
deviation = [] ;

while(windowEnd < numSamples-windowSize)
    currFrame = x(windowStart:windowEnd) ;
    currFFT = fft(currFrame) ;
    if (noiseFlag)
        if(max(currFrame) > ampThreshold)
            noiseFlag = 0 ;
            startTime = windowStart ;
        end
    end        
    
    currMag = abs(currFFT) ;
    currPhase = angle(currFFT) ;
    prevMag = abs(prevFFT) ;
    prevPhase = angle(prevFFT) ;
    prev2Phase = angle(prev2FFT) ;
%     prevPhase = unwrap_phase(angle(prevFFT), angle(prev2FFT), fs, hopSize, windowSize)  ;
%     prev2Phase = unwrap_phase(angle(prevFFT), angle(prev2FFT), fs, windowSize) ;
    
    phaseChange = prevPhase - prev2Phase ;
    phaseChange(phaseChange < 0) = phaseChange(phaseChange < 0) + 2*pi ;
    expectedPhase = mod(prevPhase + phaseChange, 2*pi) ; 
%     expectedPhase = wrapTo2Pi(expectedPhase) - pi ;
    
%     deviation = [deviation, sum(abs(currMag - prevMag)) + sum(abs(currPhase - expectedPhase))] ;
%     deviation = [deviation, sum((currMag - prevMag).^2) + sum((currPhase - expectedPhase).^2)] ;
    
    deviation = [deviation, sum((currMag - prevMag).^2)] ;
    
    index = [index, windowStart] ;
    windowStart = windowStart + hopSize ;
    windowEnd = windowEnd + hopSize ;
    prev2FFT = prevFFT ;
    prevFFT = currFFT ;
end

deviation = deviation / max(abs(deviation)) ;
deviation = smooth(deviation, round(minTime/2*fs/hopSize)) ;
% deviation = smooth(deviation, round(minTime*fs/hopSize)) ;

onset_times = [] ;
% thresh = 0.4 ;
for i = 2 : length(deviation)
    try
%         thresh = 2*mean(deviation(i-round(minTime*fs/hopSize):i-1)) ;
        thresh = 2*median(deviation(i-round(minTime/2*fs/hopSize):i-1)) ;
    catch 
        thresh = 0.4 ;
    end
    if (deviation(i) > thresh)
        onset_times = [onset_times, index(i-1)] ;
    end
end
length(onset_times)
onset_times = index(deviation > 0.3) ;
[~,locs] = findpeaks(deviation) ;
onset_times = index(locs) ;
onset_times(onset_times < startTime) = [] ;
diffOnsets = diff(onset_times) ;
falseOnsets = find(diffOnsets < minTime*1.5*fs) ;
onset_times(falseOnsets) = [] ;
 deviation(falseOnsets) = [] ;

diffOnsets = diff(onset_times) ;
for i = 2 : length(diffOnsets) - 1
diffOnsets = diff(onset_times) ;
onset_times(diffOnsets < minTime*1.5*fs) = [] ;
diffOnsets = diff(onset_times) ;
onset_times(diffOnsets < minTime*2*fs) = [] ;
diffOnsets = diff(onset_times) ;
end

figure ;
plot(x) ;
hold on
plot( deviation, 'k') ;
plot(onset_times, zeros(1,length(onset_times)), 'ro') ;
title('Detected Onsets') ;
