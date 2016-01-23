function [h,p]=m_hpss_IGIWprior(x,wlen,ITER)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERSAMUS project: Harmonic and Percussive separation
% General MODEL: Rx(n,f)=sum_vj(n,f)Rj(n,f)
%
% SPECTRAL CONTINUITY: vj(n,f) follows the Inverse Gamma distribution 
% p(vh(n,f)|vh(n,f-1))=G/IG(alpha,(alpha+k)vh(n,f-1))
% p(vp(n,f)|vp(n-1,f))=G/IG(alpha,(alpha+k)vp(n-1,f))
%
% SPATIAL CONTINUITY
% p(Rj(n,f)|p(Rj(n-1,f))=InverseWishart(sigma*Rj(n-1,f),m)

% x: multichannel mixture
% ITER: number of EM iteration ~5-10
% wlen: STFT window length

% Ngoc Q. K. Duong, INRIA Rennes, Sep. 2010.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    I=size(x,2); J=2; lf=2; lt=2;
    
    %--------------------------------------------------------
    % parameters of prior that can be changed
    gamma1=0.5;  % trade-off parameter of spatial continuity prior
    m=5;   
    sigma=m-I;
    
    gamma2=1;   % trade-off parameter of spectral continuity prior
    alpha_h=10;    % shape parameter for harmonic sources
    k_h=-1;       % scale para (-1 for "mean", 1 for "mode")
    alpha_p=10;    % shape parameter for percussive sources
    k_p=-1;       % scale para (-1 for "mean", 1 for "mode")
    %--------------------------------------------------------
    
    M=wlen/2; for i=1:I, X(:,:,i)=stft(x(:,i)',wlen,M); end
    N=size(X,2); F=wlen/2+1; Cx=zeros(I,I,F,N);  T=length(x);
    winf=hanning(2*lf-1);   wint=hanning(2*lt-1).';
    
    % ----- Compute empirical covariance matrix from the mixture 
    for f=1:F
        for t=1:N
            indf=max(1,f-lf+1):min(wlen/2+1,f+lf-1);
            indt=max(1,t-lt+1):min(N,t+lt-1);
            nind=length(indf)*length(indt);
            wei=reshape(winf(indf-f+lf)*wint(indt-t+lt),1,nind); 
            X1=reshape(X(indf,indt,:),nind,I).';                
            Cx(:,:,f,t)=(X1.*(ones(I,1)*wei))*X1'/sum(wei);
        end
    end               
    
    % ------ Initialization
    v=zeros(F,N,J);
    R=zeros(I,I,F,N,J);
    for f=1:F,
        for t=1:N,
            temp=real(trace(Cx(:,:,f,t)));
            R(:,:,f,t,:)=repmat(Cx(:,:,f,t)/(temp+realmin)+1e-6*eye(I),[1 1 1 1 J]);
            v(f,t,:)=temp/J;
        end
    end
    
    % ----- MAP by EM
    for iter=1:ITER        
        for f=1:F
            for t=1:N
                if all(v(f,t,:)) && v(f,min(t+1,N),1) && v(min(f+1,F),t,2),
                    Rx=zeros(I,I);  for j=1:J, Rx=Rx+v(f,t,j)*R(:,:,f,t,j);  end               
                    for j=1:J
                        W=v(f,t,j)*R(:,:,f,t,j)/Rx;
                        R_hat=W*Cx(:,:,f,t)*W'+v(f,t,j)*R(:,:,f,t,j)-W*v(f,t,j)*R(:,:,f,t,j);

                        % update Rj M-step for all sources
                        A=R_hat/v(f,t,j)+gamma1*sigma*R(:,:,f,max(t-1,1),j);
                        c=gamma1*I+1;
                        if ((t==N)||(gamma1==0)), 
                            %R(:,:,f,t,j)=A/(gamma1*(m+I)+1);
                            R(:,:,f,t,j)=A/c;
                        else   
                            if (t==1), 
                                A=R_hat/v(f,t,j); %c=gamma1+1-gamma1*m; 
                            end
                            B=gamma1*sigma*inv(R(:,:,f,t+1,j));
                            tempB=sqrtm(B);
                            delta=(c^2)*eye(I)+4*tempB*A*tempB;
                            R(:,:,f,t,j)=0.5*tempB\(-c*eye(I)+sqrtm(delta))/tempB;
                        end

                        if (j==1)
                            % update vj in M-step for harmonic sources
                            if ((t==N)||(gamma2==0))
                                %v(f,t,j)=real(trace(R_hat*inv(R(:,:,f,t,j))))/(I+gamma2*(alpha_h+1+(alpha_h+k_h)*v(f,max(t-1,1),j)));
                                v(f,t,j)=(real(trace(R_hat/R(:,:,f,t,j)))+gamma2*(alpha_h+k_h)*v(f,max(t-1,1),j))/(gamma2*(alpha_h+k_h)+gamma2+I);
                            else
                                a=gamma2*(alpha_h+k_h)/v(f,t+1,j);
                                b=gamma2+I;
                                if (t==1), d=-real(trace(R_hat/R(:,:,f,t,j))); %b=-gamma2*alpha_h+gamma2+I;
                                else       d=-real(trace(R_hat/R(:,:,f,t,j)))-gamma2*(alpha_h+k_h)*v(f,t-1,j); end
                                v(f,t,j)=(-b+sqrt(b^2-4*a*d))/(2*a);
                            end
                        else
                            % update vj in M-step for percussive sources
                            if ((f==F)||(gamma2==0))
                                v(f,t,j)=(real(trace(R_hat/R(:,:,f,t,j)))+gamma2*(alpha_p+k_p)*v(max(f-1,1),t,j))/(gamma2*(alpha_p+k_p)+gamma2+I);
                                %v(f,t,j)=real(trace(R_hat*inv(R(:,:,f,t,j))))/(I+gamma2*(alpha_p+1+(alpha_p+k_p)*v(f,max(t-1,1),j)));
                            else
                                a=gamma2*(alpha_p+k_p)/v(f+1,t,j);
                                b=gamma2+I;
                                if (f==1), d=-real(trace(R_hat/R(:,:,f,t,j))); %b=-gamma2*alpha_p+gamma2+I;
                                else       d=-real(trace(R_hat/R(:,:,f,t,j)))-gamma2*(alpha_p+k_p)*v(f-1,t,j); end
                                v(f,t,j)=(-b+sqrt(b^2-4*a*d))/(2*a);
                            end
                        end
                    end
                end     
            end
        end
        % normalization after each iteration
        for f=1:F,
            for t=1:N,
                for j=1:J,
                    temp=real(trace(R(:,:,f,t,j)));
                    R(:,:,f,t,j)=R(:,:,f,t,j)/(temp+realmin);
                    v(f,t,j)=v(f,t,j)*temp;
                end
            end
        end
    end
    
    % ----- Apply Wiener filtering 
    Se=zeros(I,F,N,J);
    for f=1:F
        for t=1:N
            if any(v(f,t,:)),
                Rx=zeros(I,I); for j=1:J,  Rx=Rx+v(f,t,j)*R(:,:,f,t,j); end
                for j=1:J,  Se(:,f,t,j)=v(f,t,j)*R(:,:,f,t,j)/Rx*reshape(X(f,t,:),I,1);  end
            end
        end
    end
    % ----- Inverse time-frequency transformation & save files
    Se=permute(Se,[2,3,4,1]);   
    for i=1:I, for j=1:J, se_img(:,j,i)=istft(Se(:,:,j,i),M); end; end;
    se_img=permute(se_img(1:T,:,:),[1,3,2]);
    h=se_img(:,:,1);
    p=se_img(:,:,2);


    
