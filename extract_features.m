function feat_vector = extract_features(x, fs)
 
x = x(:);
 
H = hamming(length(x));
x = x.*H;

% Pre-emphasis Filter {H(z) = 1 – 0.98z-1}
x = filter([1, -0.98], 1, x);

len = length(x);
b = nextpow2(len);
N = 2 ^ b;
X = abs(fft(x,N));

X =log10( X.*X);

f = zeros(1,N/2) ;
for i=1:N/2
    f(i) = fs / N * i;
end

X = X(1:N/2);

% To Integrate by Mel Scale
Mf = 2595.*log(1+(f./700));

coeff = zeros(1,24);

for j=1:24
    for i = 1 : N/2
        if ((Mf(i)>300+(j-1)*150)&&(Mf(i)<600+(j-1)*150))
            if (Mf(i)<450+(j-1)*150)
                g(i) = Mf(i)-(300+150*(j-1));
            else
                g(i) = (600+150*(j-1))-Mf(i);
            end

            coeff(j) = coeff(j) + X(i)*g(i)/150;
        end
    end
end

mfcc = 10.*log10(abs(dct(coeff)));

% ZCR
zcr = sum(abs(diff(sign(x)))) / 2 ;
zcr = zcr / (length(x) / fs) ;

% Spectral Centroid
X = abs(fft(x, fs)) ;
X = X(1:fs/2) ;
f = (0 : fs/2-1)' ;
ind = min(find(f>=1000)) ;
f = f(1:ind) ;
X = X(1:ind) ;
spectral_centroid = sum(X.*f) / sum(X) ;

% Skewness
X = abs(fft(x, fs)) ;
X = X(1:fs/2) ;
mu_x    = mean(X);
std_x   = std(X);

X = X - mu_x ;
skewness = sum ((X.^3)./(std_x.^3*length(X)));

% Variance
variance = var(X) ;

% Kurtosis
mu_x    = mean(X);
std_x   = std(X);

X       = X - mu_x ;
kurt = sum ((X.^4)./(std_x.^4*length(X)));

% Spectral Roll-off
kappa   = 0.9;

afSum   = sum(X);
spectral_roll_off = find(cumsum(X(:)) >= kappa*afSum, 1); 

% Feature Vector
feat_vector = [mfcc, zcr, spectral_centroid, skewness, variance, kurt, spectral_roll_off]' ;