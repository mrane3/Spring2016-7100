%% compute binary SDM matrix
% input:
%   SDM: numSamples by numSamples float matrix, self-distance matrix
%   threshold: float, constant value for thresholding the SDM
% output:
%   SDM_binary: numSamples by numSamples int matrix, binary SDM

function [SDM_binary] = computeBinSdm(SDM, threshold)

SDM_binary = SDM ;
SDM_binary(SDM_binary <= threshold) = 0 ;
SDM_binary(SDM_binary > threshold) = 1 ;

SDM_binary = logical(SDM_binary) ;