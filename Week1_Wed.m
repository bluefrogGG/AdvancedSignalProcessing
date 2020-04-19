f=[0, 0.06, 0.1, 1];
m=[1 1 0 0];
L = 100;
h = firpm(L, f, m);
%order? 

[H,w]=freqz(h,1); % FFT
plot(w, abs(H));

%% IIR filter

a=0.7545;
g=1-a;
h=[1 -a];
[H,w]=freqz(g,h);
plot(w, abs(H));

%% reading audiofile

% read the audiofile and listen
[x, fs] = audioread('BabyElephantWalk60.wav');
sound(x, fs);
pause;
clear sound;

x_left = x(:,1);
N = length(x);
y = zeros(N,1);
alpha = 0.7535;

for n=2:N
    y(n) = alpha*y(n-1) + (1-alpha)*x(n);
end

audiowrite('ttt.wav', y, fs);
sound(y, fs);
pause;
clear sound;

%% Week 2 - FFT

%frequency domain filtering
[x, fs] = audioread('BabyElephantWalk60.wav');
x_left = x(:,1);

% low-pass filter
x1= 0.06;
x2 = 0.07;
L = 20;

f_freq = @() FFT(L, x1, x2, x_left);
f_time = @() conv_time(L, x1, x2, x_left);

t1 = timeit(f_freq)
t2 = timeit(f_time)

xs = [1 : 5];
y_t = zeros(1,5); y_f = zeros(1,5);
xs = 2 .^ xs;
xs_L = xs .* L;

for k = 1: 5
    f_freq = @() FFT(xs_L(1,k), x1, x2, x_left);
    f_time = @() conv_time(xs_L(1,k), x1, x2, x_left);

    y_t(1, k) = timeit(f_freq);
    y_f(1, k) = timeit(f_time);
end

plot(xs_L, y_t, xs_L, y_f);