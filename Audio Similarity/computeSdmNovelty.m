%% compute novelty function from Self-distance matrix
% input:
%   SDM: float N by N matrix, self-distance matrix
%   L: int, size of the checkerboard kernel (L by L) preferably power of 2
% output:
%   nvt: float N by 1 vector, audio segmentation novelty function 

function [nvt] = computeSdmNovelty(SDM, L)

nvt = zeros(size(SDM,1), 1) ;

C = [ones(L/2), -ones(L/2); -ones(L/2), ones(L/2)] ;

% SDM = padarray(SDM, [L/2, L/2]) ;

for i = 1 : size(SDM,1)-L+1
    nvt(i) = sum(sum(C .* SDM(i:i+L-1, i:i+L-1))) ;
end

nvt = abs(nvt) ;        % Since we get just negative values, only the magnitude matters