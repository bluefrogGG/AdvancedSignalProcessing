%% Time-domain adaptive equalization

W = 3.5;
mu = 0.075;
h = zeros(4,1);
M = 100;
e_mat = zeros(2003, M);
BER = zeros(M,1);

for m = 1:M

for n=2:4
    h(n) = (1/2)*(1+cos(2*pi*(n-3)/W));
end
h(1) = 0;

x=rand(2000,1) - 0.5;
for n=1:2000
    if(x(n) > 0)
        x(n) = 1;
    else x(n) = -1;
    end
end

r = conv(x, h);
nu = sqrt(0.001)*randn(2003,1); % conv(x,h) sizeø° ∏¬√Á¡‹..
r = r + nu;

tap = 11;
d = zeros(2003,1);
for k = 8:2003
    d(k) = x(k-7);
end
r_tr = zeros(tap,1);
w = zeros(tap,1); w((tap+1)/2) = 1;
e = zeros(2003,1); q = zeros(2003,1); % quantization result
for n=1:2003
    for k=tap:-1:2
        r_tr(k) = r_tr(k-1);
    end
    r_tr(1) = r(n);
    
    y = w' * r_tr; % dot product
    e(n) = d(n) - y;
    w = w + mu*e(n)*r_tr; % update
    
    if(y > 0) 
        q(n) = 1;
    else q(n) = -1;
    end
    if(q(n) ~= d(n))
        BER(m) = BER(m) + 1;
    end
end

BER(m) = BER(m) / 2003;  
e_mat(:,m) = e .* e;

end

e_avg_1 = e_mat * ones(M,1) * (1/M);
BER = BER' * ones(M,1) * (1/M);
plot(log(e_avg_1));
hold;

%Refer to 'Week5_1.fig'
%BER = 0.0035 -> 0.0080 

%% DCT-domain aaptive equalization

M = 11; % tap
delay = 2+floor(M/2); 
dct = zeros(M,M);
for n=2:M
    for m=1:M
        dct(n,m)=cos((2*m-1)*(n-1)*pi/2/M);
    end
end
dct(1,:)= 1/sqrt(2);

W=3.5; mu=0.075; rep = 100; gamma = 0.95;
lambda = zeros(M,1);
e = zeros(2003, rep);
for n=1:rep
    c = zeros(4,1);
    for i =2:4
        c(i) = (1/2)*(1+cos(2*pi*(i-3)/W));
    end

    x=rand(2000,1) - 0.5;
    for j=1:2000
        if(x(j) > 0)
            x(j) = 1;
        else x(j) = -1;
        end
    end
    r = conv(x, c);
    nu = sqrt(0.001)*randn(2003,1); % conv(x,h) sizeø° ∏¬√Á¡‹..
    r = r + nu;
    d = zeros(2003,1);
    for k = (1+delay):2003
        d(k) = x(k-delay);
    end
    %initialize
    h = zeros(M,1); % DCT domain
    u = zeros(M,1);
    for m = 1:2003
        for j = M:-1:2
            u(j)= u(j-1);
        end
        u(1) = r(m); % update u
        v = dct * u; % u in DCT domain
        y = h' * v;
        e(m,n) = d(m) - y; %error, save
        lambda = gamma * lambda + (1-gamma)* v .* v;
        h = h + mu * e(m,n) * v ./ lambda;
    end
    
end

e_avg_2 = (e .* e) * ones(rep,1) * (1/rep);
plot(log(e_avg_2));


