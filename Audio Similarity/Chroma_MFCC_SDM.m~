featureMatrix_chroma = get_features('pitchChroma');
SDM_chroma = computeSelfDistMat(featureMatrix_chroma);
featureMatrix_mfcc = get_features('mfcc');
SDM_mfcc =computeSelfDistMat(featureMatrix_mfcc); 
figure(1);
imagesc(SDM_chroma);
title('SDM for Pitch Chroma');
figure(2);
imagesc(SDM_mfcc);
title('SDM for MFCC');
[nvt2c]=computeSdmNovelty(SDM_chroma,2);
[nvt8c]=computeSdmNovelty(SDM_chroma,8);
[nvt16c]=computeSdmNovelty(SDM_chroma,16);
[nvt2m]=computeSdmNovelty(SDM_mfcc,2);
[nvt8m]=computeSdmNovelty(SDM_mfcc,8);
[nvt16m]=computeSdmNovelty(SDM_mfcc,16);

[numerics, strings]=xlsread('03-Sargon-Waiting For Silence.xls', 'A1:A28');


figure(3);
subplot(3,1,1)
title('SDM-Pitch Chroma Novelty for L=2');
hold on;
plot(nvt2c);
stem(numerics, ones(1,length(numerics)), 'r') ;
subplot(3,1,2)
title('SDM-Pitch Chroma Novelty for L=8');
hold on;
plot(nvt8c);
stem(numerics, ones(1,length(numerics)), 'r') ;
subplot(3,1,3)
hold on;
title('SDM-Pitch Chroma Novelty for L=16');
plot(nvt16c);
stem(numerics, ones(1,length(numerics)), 'r') ;

figure(4);
subplot(3,1,1)
hold on;
title('SDM-MFCC Novelty for L=2');
stem(numerics, ones(1,length(numerics)), 'r') ;
plot(nvt2m);
subplot(3,1,2)
hold on;
title('SDM-MFCC Novelty for L=8');
stem(numerics, ones(1,length(numerics)), 'r') ;
plot(nvt8m);
subplot(3,1,3)
hold on;
title('SDM-MFCC Novelty for L=16');
plot(nvt16m);
stem(numerics, ones(1,length(numerics)), 'r') ;

Rc =computeLagDistMatrix(SDM_chroma);
Rm= computeLagDistMatrix(SDM_mfcc);
figure(5);
imagesc(Rc);
title('Time lag representation of Chroma SDM');

figure(6);
imagesc(Rm);
title('Time lag representation of MFCC SDM');

figure(7);
Lag_Bin_Chroma=computeBinSdm(Rc,0.1175);
colormap gray;
imagesc(Lag_Bin_Chroma);
title('Binary Pitch Chroma Lag ')

figure(8);
Lag_Bin_MFCC=computeBinSdm(Rm,0.4275);
colormap gray;
imagesc(Lag_Bin_MFCC);

figure(9);
SDM_ed_Chroma=erodeDilate(Lag_Bin_Chroma,1);
colormap gray;
imagesc(SDM_ed_Chroma);

figure(10);
SDM_ed_MFCC=erodeDilate(Lag_Bin_MFCC,1);
colormap gray;
imagesc(SDM_ed_MFCC);







 



