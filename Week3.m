%% NLMS for each sample
% 'etest2p.wav' represents dn
% 'etest2r.wav' represents xn
% 'etest4p.wav' represents dn
% 'etest4p.wav' represents xn

[d, fs1] = audioread('etest2p.wav');
[x, fs2] = audioread('etest2r.wav');

% how to determine the value of M?
M = 100; 
h = zeros(M,1);
mu = 1;
epsilon = 0.1;

e_temp = zeros(length(x)-M +1, 1);

for n = M:length(x)
    x_n = flipud(x((n-M+1):n,1)); % x[n], x[n-1], ... ,x[n-(M-1)]
    energy_x = sum(x_n .* x_n); 
    e_n = d(n,1) - h' * x_n;
    e_temp(n-M+1, 1) = e_n;
    
    h = h + (mu/(sqrt(energy_x)+epsilon)) * e_n * x_n;
end

subplot(3,1,1);
plot(x);
subplot(3,1,2);
plot(d);
subplot(3,1,3);
plot(e_temp);

filename = 'etest2_NLMS.wav';
audiowrite(filename, e_temp, fs1);

%% Adaptive filtering in the frequency domain

[d, fs1] = audioread('etest2p.wav');
[x, fs2] = audioread('etest2r.wav');

k = 8;
M = 2^k;
mu = 0.2;
lambda = 0.2;

N = length(x);
x_extended = [zeros(M,1); x; zeros(M-rem(N,M),1)];
d_extended = [zeros(M,1); d; zeros(M-rem(N,M),1)];
W = zeros(2*M, 1); 
y = zeros(2*M,1);
P = zeros(2*M,1);

e = zeros(length(x_extended),1);

for n = 1: M : length(x_extended)-M
    
    v = x_extended(n:(n+ 2*M-1),1);
    V = fft(v);
    P  = lambda*P + (1-lambda)*(V .* conj(V));
    
    Y = W .* V;
    y = ifft(Y);
    
    e(n:(n+M-1),1) = d_extended((n+M):(n+M+M-1),1) - y((M+1):2*M,1);
    e_temp = [zeros(M,1); e(n:(n+M-1),1)];
    E = fft(e_temp);
    
    Theta = E .* conj(V) ./ P;
    theta = ifft(Theta);
    theta = [theta(1:M,1); zeros(M,1)];
    Theta = fft(theta);
    W = W + mu*Theta;
    
end


subplot(3,1,1);
plot(d);
subplot(3,1,2);
plot(x);
subplot(3,1,3);
plot(e);

filename = 'etest2_NLMS_F.wav';
audiowrite(filename, e, fs1);

%% same filter on another sample

[d, fs1] = audioread('etest4p.wav');
[x, fs2] = audioread('etest4r.wav');

k = 8;
M = 2^k;
mu = 0.1;
lambda = 0.5;

N = length(x);
x_extended = [zeros(M,1); x; zeros(M-rem(N,M),1)];
d_extended = [zeros(M,1); d; zeros(M-rem(N,M),1)];
W = zeros(2*M, 1); 
y = zeros(2*M,1);
P = zeros(2*M,1);

e = zeros(length(x_extended),1);

for n = 1: M : length(x_extended)-M
    
    v = x_extended(n:(n+ 2*M-1),1);
    V = fft(v);
    P  = lambda*P + (1-lambda)*(V .* conj(V));
    
    Y = W .* V;
    y = ifft(Y);
    
    e(n:(n+M-1),1) = d_extended((n+M):(n+M+M-1),1) - y((M+1):2*M,1);
    e_temp = [zeros(M,1); e(n:(n+M-1),1)];
    E = fft(e_temp);
    
    Theta = E .* conj(V) ./ P;
    theta = ifft(Theta);
    theta = [theta(1:M,1); zeros(M,1)];
    Theta = fft(theta);
    W = W + mu*Theta;
    
end


subplot(3,1,1);
plot(d);
subplot(3,1,2);
plot(x);
subplot(3,1,3);
plot(e);

filename = 'etest4_NLMS_F.wav';
audiowrite(filename, e, fs1);