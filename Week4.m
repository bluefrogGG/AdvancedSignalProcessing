%% FIR adaptive line enhancement
%Using Week3.m

%Trouble:
% 1. doesn't work : enhanced_sin.wav cannot be heard
% 2. don't know how to get SNR of each signal

fs = 8000;
nbits = 16;

N = 80000;
ref = zeros(N,1); % this corresponds to x
pri = [1:N]'; % corresponds to d
noise = 0.2 * randn(N,1);
pri = 0.8 * sin(pi/4*pri)+ noise; 

ref(2:N,1) = pri(1:N-1,1);

audiowrite('noisysine.wav', pri, fs);

M = 100; 
h = zeros(M,1);
mu = 0.05;
epsilon = 0.1;

e_temp = zeros(length(ref)-M +1, 1);

for noise = M:length(ref)
    x_n = flipud(ref((noise-M+1):noise,1)); % x[n], x[n-1], ... ,x[n-(M-1)]
    energy_x = sum(x_n .* x_n); 
    e_n = pri(noise,1) - h' * x_n;
    e_temp(noise-M+1, 1) = e_n;
    
    h = h + (mu/(sqrt(energy_x)+epsilon)) * e_n * x_n;
end

audiowrite('enhanced_sine.wav', e_temp, fs);

subplot(2,2,1);
plot(pri);
[H1,w]=freqz(pri); subplot(2, 2, 2); plot(w,abs(H1));
subplot(2,2,3); 
plot(e_temp);
[H2,w] = freqz(e_temp); subplot(2,2,4); plot(w,abs(H2));
%How to get SNR..?

%% Chirp signal
%works fine
% Beware of integrating frequency into phase

[s,fs]= audioread('preamble10.wav');
N = length(s);
%s = s(N/2 : N); N = length(s); 
f1 = 2000; A=0.05;

chirp_end= 30000;
f_lin = zeros(N,1);
f_lin(1:chirp_end,1) = 2*pi*f1/fs/chirp_end * [0:chirp_end-1]';
f_lin(chirp_end+1 : N,1) = ones(N-chirp_end,1) * f_lin(chirp_end,1);
phase = f_lin;
for noise = 2:N
    phase(noise,1) = phase(noise-1,1)+ phase(noise,1);
end

chirp = A*sin(phase);

audiowrite('preamble_chirp.wav', chirp+s, fs);

% Now Filtering

s_chirp = s+chirp;
y = zeros(N,1); w_est = zeros(N,1);%frequency estimate
e = zeros(N,1);
a = 0; r = 0.9987; lambda  = 0.9; C = 0; D = 0; % probably the most optimal value : r = 0.9987, lambda = 0.9
w_win = zeros(3,1); % w_win(1) = w(n), w_win(2) = w(n-1), w_won(3) = w(n-2)

for n=1:N
    w_win(3) = w_win(2);
    w_win(2) = w_win(1);
    w_win(1) = s_chirp(n) - (1+r) * a * w_win(2) - r * w_win(3); % Why..?
    y(n) = w_win(1) + 2*a*w_win(2) + w_win(3); % y(n) : baseband
    e(n) = s_chirp(n) - y(n); % e(n): enhanced morse
    
    C = lambda * C + (1-lambda)*w_win(2)*(w_win(1)+w_win(3));
    D = lambda * D + (1-lambda)*2*(w_win(2)^2);
    temp = -C/D;
    if abs(temp) < 1 
        a = temp;
    end
    w_est(n) = acos(-a);
end

audiowrite('chirp_preamble.wav', e, fs);

subplot(2,2,1);
plot(s_chirp);
%[H1,w]=freqz(s_chirp); subplot(2, 2, 2); plot(w,abs(H1));
subplot(2,2,2); plot(f_lin);
subplot(2,2,3); 
plot(e);
%[H2,w] = freqz(e); subplot(2,2,4); plot(w,abs(H2));
subplot(2,2,4); plot(w_est);

%% Adaptive IIR filtering
% using eactly the same sine wave
% success
fs = 8000;
nbits = 16;

N = 80000;
x = [1:N]';
noise = 0.2 * randn(N,1);
x = 0.8 * sin(pi/4*x)+ noise;
y = zeros(N,1); w_est = zeros(N,1);%frequency estimate
e = zeros(N,1);
a = 0; r = 0.9; lambda  = 0.9; C = 0; D = 0;
w_win = zeros(3,1); % w_win(1) = w(n), w_win(2) = w(n-1), w_won(3) = w(n-2)

for n=1:N
    w_win(3) = w_win(2);
    w_win(2) = w_win(1);
    w_win(1) = x(n) - (1+r) * a * w_win(2) - r^2 * w_win(3);
    y(n) = w_win(1) + 2*a*w_win(2) + w_win(3); % y(n) : baseband
    e(n) = x(n) - y(n); % e(n): enhanced sine wave
    
    C = lambda * C + (1-lambda)*w_win(2)*(w_win(1)+w_win(3));
    D = lambda * D + (1-lambda)*2*(w_win(2)^2);
    temp = -C/D;
    if abs(temp) < 1 
        a = temp;
    end
    w_est(noise) = acos(-a);
end

audiowrite('noisysin_filtered.wav', e, fs);

subplot(2,2,1);
plot(x);
[H1,w]=freqz(x); subplot(2, 2, 2); plot(w,abs(H1));
subplot(2,2,3); 
plot(e);
[H2,w] = freqz(e); subplot(2,2,4); plot(w,abs(H2));

%% Morse code with noise
% Help => .... . .-.. .--.

fs = 6000;
dur = 500;
code = [1,0,1,0,1,0,0,0,1,0,0,0,1,0,1,1,1,0,1,0,1,0,0,0,1,0,1,1,1,0,1,1,1,0,1,0];
code = [code, code];
N = dur*length(code);
morse = zeros(1,N);
for i=1:length(code)
    for k = 1:dur
        %morse((i-1)*dur + k) = (-1)^k * code(i);
        morse((i-1)*dur + k) = code(i); 
    end
end
noise = 0.1*randn(N,1);

xs = [1:N]'; % N*1
xs_morse = xs .* morse'; 
%morse_corr = 0.3 * morse' + noise;
morse_corr = 0.3* sin(pi/4 * xs_morse) + noise;
audiowrite('morse_corrupted_help.wav', morse_corr, fs);

% now adaptive filtering

y = zeros(N,1); w_est = zeros(N,1);%frequency estimate
e = zeros(N,1);
a = 0; r = 0.99; lambda  = 0.99; C = 0; D = 0;
w_win = zeros(3,1); % w_win(1) = w(n), w_win(2) = w(n-1), w_won(3) = w(n-2)

for n=1:N
    w_win(3) = w_win(2);
    w_win(2) = w_win(1);
    w_win(1) = morse_corr(n) - (1+r) * a * w_win(2) - r^2 * w_win(3);
    y(n) = w_win(1) + 2*a*w_win(2) + w_win(3); % y(n) : baseband
    e(n) = morse_corr(n) - y(n); % e(n): enhanced morse
    
    C = lambda * C + (1-lambda)*w_win(2)*(w_win(1)+w_win(3));
    D = lambda * D + (1-lambda)*2*(w_win(2)^2);
    temp = -C/D;
    if abs(temp) < 1 
        a = temp;
    end
    w_est(n) = acos(-a);
end

audiowrite('morse_filtered.wav', e, fs);

subplot(2,2,1);
plot(morse_corr);
[H1,w]=freqz(morse_corr); subplot(2, 2, 2); plot(w,abs(H1));
subplot(2,2,3); 
plot(e);
[H2,w] = freqz(e); subplot(2,2,4); plot(w,abs(H2));

%% 
