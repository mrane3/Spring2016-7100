function X=stft(x,L,M)

% STFT Short-Term Fourier Transform using a sine window.
%
% X=stft(x,L)
% X=stft(x,L,M)
%
% Inputs:
% x: 1 x T vector containing a single-channel signal
% L: length of the STFT window in samples (must be a multiple of 4)
% M: step between successive windows in samples (must be a multiple of 2,
%    a divider of L and smaller than L/2) (default: L/2)
%
% Output:
% X: (L/2+1) x N matrix containing the STFT coefficients for positive frequencies
%    with N=ceil(T/M)
%
% See also istft.

%%% Errors and warnings %%%
if nargin<2, error('Not enough input arguments.'); end
if nargin<3, M=L/2; end
[I,T]=size(x);
if I>1, error('The input signal must be a row vector.'); end
if L-4*floor(L/4), error('The window length must be a multiple of 4.'); end
if (L-M*floor(L/M)) | (M-2*floor(M/2)), error('The step size must be a multiple of 2 and a divider of the window length.'); end
if M>L/2, error('The step size must be smaller than half the window length.'); end

%%% Computing STFT coefficients %%%
% Defining sine window
win=sin((.5:L-.5)/L*pi).';
% Zero-padding
N=ceil(T/M);
x=[x,zeros(1,N*M-T)];
% Pre-processing for edges
x=[zeros(1,(L-M)/2),x,zeros(1,(L-M)/2)];
swin=zeros((N-1)*M+L,1);
for t=0:N-1,
    swin(t*M+1:t*M+L)=swin(t*M+1:t*M+L)+win.^2;
end
swin=sqrt(L*swin);
X=zeros(L/2+1,N);
for t=0:N-1,
    % Framing
    frame=x(t*M+1:t*M+L).'.*win./swin(t*M+1:t*M+L);
    % FFT
    fframe=fft(frame);
    X(:,t+1)=fframe(1:L/2+1);
end

return;