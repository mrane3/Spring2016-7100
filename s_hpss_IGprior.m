function [h,p]=s_hpss_IGprior(x,wlen,ITER) 

% SINGLE-CHANNEL harmonic and percussive component separation in musical mixtures

% Input:
% x: nsample x nchannel mixture signal
% wlen: STFT window length (4096 for 44.1 kHz and 1024 for 16 kHz sampling rate)
% ITER: number of EM iterations, (~5-10)

% Output:
% h: nsample x nchannel harmonic component
% p: nsample x nchannel percussive component

%**************************************************************************
% Copyright 2010. This work was sponsored by VERSAMUS project(http://versamus.inria.fr/)
% Author Ngoc Q. K. Duong
% This source code is distributed under the terms of the GNU General Public
% License version 3 (http://www.gnu.org/licenses/gpl.txt)

% If you find it useful, please cite the following reference:
% N. Q. K. Duong, H. Tachibana, E. Vincent, N. Ono, R. Gribonval, and S. Sagayama
% "Multichannel harmonic and percussive component separation by joint modeling of
% spatial and spectral continuity", submitted to ICASSP 2011
%**************************************************************************            

    % parameters of prior that can be changed 
    gamma=1;
    alpha_h=10;    % shape parameter for harmonic sources
    k_h=-1;       % scale para (-1 for "mean", 1 for "mode")
    alpha_p=10;    % shape parameter for percussive sources
    k_p=-1;       % scale para (-1 for "mean", 1 for "mode")
    
    M=wlen/2; X=stft(x',wlen,M);   
    F=size(X,1); N=size(X,2); T=length(x);
    
    % power spectrogram
    XX=abs(X).^2; 
    
    % initialization
    H=XX/2; 
    P=XX/2;   
    
    for iter=1:ITER        
        for f=1:F
            for t=1:N    
                % EM update for harmonic component
                W=H(f,t)/(H(f,t)+P(f,t));   
                H_hat=W*XX(f,t)*W'+H(f,t)-W*H(f,t);
                if (t==N)
                    H(f,t)=(H_hat+gamma*(alpha_h+k_h)*H(f,t-1))/(gamma*(alpha_h+k_h)+gamma+1);
                else
                    a=gamma*(alpha_h+k_h)/H(f,t+1);
                    b=gamma+1;
                    if (t==1), d=-H_hat; b=-gamma*alpha_h+gamma+1;
                    else       d=-H_hat-gamma*(alpha_h+k_h)*H(f,t-1); end
                    H(f,t)=(-b+sqrt(b^2-4*a*d))/(2*a);
                end
 
                % EM update for percussive component
                W=P(f,t)/(H(f,t)+P(f,t));   
                P_hat=W*XX(f,t)*W'+P(f,t)-W*P(f,t);
                if (f==F)
                    P(f,t)=(P_hat+gamma*(alpha_p+k_p)*P(f-1,t))/(gamma*(alpha_p+k_p)+gamma+1);
                else
                    a=gamma*(alpha_p+k_p)/P(f+1,t);
                    if (f==1), d=-P_hat; b=-gamma*alpha_p+gamma+1;
                    else       d=-P_hat-gamma*(alpha_p+k_p)*P(f-1,t); end
                    P(f,t)=(-b+sqrt(b^2-4*a*d))/(2*a);
                end   
            end
        end     
    end

    % Wiener filtering & reconstruct time-domain signals
    h=istft((H./(H+P)).*X,M); h=h(1:T); 
    p=istft((P./(H+P)).*X,M); p=p(1:T);
    




    
