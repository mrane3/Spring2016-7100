%% compute self-distance matrix
% input:
%   featureMatrix: numFeatures by numSamples float matrix, feature matrix
% output:
%   SDM: numSamples by numSamples float matrix, self-distance matrix

function SDM = computeSelfDistMat(featureMatrix)

SDM = zeros(size(featureMatrix, 2)) ;

for i = 1 : size(featureMatrix, 2)
    for j = 1 : i % size(featureMatrix, 2)
        SDM(i,j) = norm(featureMatrix(:,i) - featureMatrix(:,j)) ;
        SDM(j,i) = SDM(i,j) ;
    end
end

