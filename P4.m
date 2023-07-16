                            %% Sparse Blind Deconvolution
                            %% Part 4
file = load('hw7.mat');
x = file.x1;
K = 5;
L = 100;
T = 2000;
num_iter = 10;

X = fftshift(fft(x));
Y = zeros(length(X),1);
Y(abs(X)>5) = X(abs(X)>5);
y = transpose(ifft(ifftshift(Y)));

tau = rand_tau(K,L,T);
alpha = rand(K,1);
w = zeros(T,1);
for i=1:num_iter
        %% W is fixed. Updating S
    w(tau) = alpha;
    W = fftshift(fft(w));
    S = Y./W;
    s = ifft(ifftshift(S));
    s = s(1:L);
    
        %% S is fixed. Updating W
    R = corr_xy(s,y);
    thresh = 1/2 * max(R);
    [pks,locs] = findpeaks(R,'MinPeakDistance',L,'MinPeakHeight',thresh,'NPeaks',5);
    tau = locs;
    alpha = pks;
    for k=1:K
       if tau(k) < L/2
          tau(k) = tau(k) + L/2; 
       end
    end
end

disp("W(t)");
disp("alpha | tau");
disp(vpa([round(alpha,2) tau]));
figure(1);
plot(1:T,y);
xlabel('t');
ylabel('Amp');
title('x(t) with decreased noise');
grid on;
figure(2);
plot(1:L,s);
xlabel('t');
ylabel('Amp');
title('s(t) spike recovered');
grid on;

        %% Local Necessary Functions
function tau = rand_tau(K,L,T)
    tau = (T-3*L)*rand(K,1) + L;
    tau = sort(tau);
    for i=2:K
       tau_temp = tau(i);
       if tau_temp - tau(i-1) <= L
          tau(i) = tau(i-1) + L; 
       end
    end
    tau = ceil(tau);
end

function R = corr_xy(s,x)
    T = length(x);
    L = length(s);
    R = zeros(T-L,1);
    for i=1:length(R)
       t1 = i;
       t2 = i+L-1;
       x_window = x(1,t1:t2);
       R(i) = x_window * s;
    end
end