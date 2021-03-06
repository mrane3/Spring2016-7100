%% compute lag distance matrix 
% input:
%   SDM: numSamples by numSamples float matrix, self-distance matrix
% output:
%   R: numSamples by numSamples float matrix, lag-distance matrix
% Note: R should be a triangular matrix, xaxis = time, yaxis = lag
%       for more details, please refer to Figure 2 in the reference
%       "Paulus et al., Audio-based Music Structure Analysis, 2010"

function R = computeLagDistMatrix(SDM)

R = triu(SDM) ;
for i = 1 : size(SDM,1)
    ind = find(R(i,:)) ;
    if(~isempty(ind))
        R(i,:) = circshift(R(i,:), [0, -ind(1)+1]) ;
    end
end
