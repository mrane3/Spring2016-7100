function x=istft(X,M)

% ISTFT Inverse Short-Term Fourier Transform using a sine window.
%
% x=istft(X)
% x=istft(X,M)
%
% Inputs:
% X: (L/2+1) x N matrix containing a set of STFT coefficients for positive frequencies
% M: step between successive windows in samples (must be a multiple of 2,
%    a divider of L and smaller than L/2) (default: L/2)
%
% Output:
% x: 1 x (N*M) vector containing the inverse STFT signal
%
% If x is a signal of length T, X=stft(x,L) and y=istft(X), then x=y(1:T).
%
% See also stft.

%%% Errors and warnings %%%
if nargin<1, error('Not enough input arguments.'); end
[L,N]=size(X);
if ~(L-2*floor(L/2)), error('The number of rows of the STFT matrix must be odd.'); end
L=2*(L-1);
if nargin<2, M=L/2; end
if L-M*floor(L/M), error('The step size must be a divider of two times the number of rows of the STFT matrix minus two.'); end
if (M-2*floor(M/2)), error('The step size must be a multiple of 2.'); end
if M>L/2, error('The step size must be smaller than the number of rows of the STFT matrix minus one.'); end

%%% Computing inverse STFT signal %%%
% Defining sine window
win=sin((.5:L-.5)/L*pi);
% Pre-processing for edges
T=N*M;
swin=zeros(1,T+L-M);
for t=0:N-1,
    swin(t*M+1:t*M+L)=swin(t*M+1:t*M+L)+win.^2;
end
swin=sqrt(swin/L);
x=zeros(1,T+L-M);
for t=0:N-1,
    % IFFT
    fframe=[X(:,t+1);conj(X(L/2:-1:2,t+1))];
    frame=real(ifft(fframe));
    % Overlap-add
    x(t*M+1:t*M+L)=x(t*M+1:t*M+L)+frame.'.*win./swin(t*M+1:t*M+L);
end
% Truncation
x=x((L-M)/2+1:(L-M)/2+T);

return;