                            %% Sparse Blind Deconvolution
                            %% Part 1, Part 2
file = load('hw7.mat');
X = file.X;
x1 = file.x1;
x2 = file.x2;
L = 100;
K = 5;
[~,T] = size(X);

            %% Part 1
figure(1);
plot(1:T,x1);
xlabel('t');
ylabel('Amp');
title("x_1(t)");
grid on;
num_iteration_1 = 5;
[s1,alpha1,tau1] = SCSBD(x1,K,L,T,num_iteration_1);
figure(2);
plot(1:L,s1);
xlabel('t');
ylabel('Amp');
title('Recovered Spike Signal s_1(t)');
grid on;
disp("W_1(t)");
disp("alpha | tau");
disp(vpa([round(alpha1,4) tau1]));

            %% Part 2
figure(3);
plot(1:T,x2);
xlabel('t');
ylabel('Amp');
title("x_2(t)");
grid on;
num_iteration_2 = 5;
[s2,alpha2,tau2] = SCSBD(x2,K,L,T,num_iteration_2);
figure(4);
plot(1:L,s2);
xlabel('t');
ylabel('Amp');
title('Recovered Spike Signal s_2(t)');
grid on;
disp("W_2(t)");
disp("alpha | tau");
disp(vpa([round(alpha2,4) tau2]));




            %% Local Necessary Functions

function [s,alpha,tau] = SCSBD(x,K,L,T,num_iter)
    tau = rand_tau(K,L,T);
    alpha = rand(K,1);
    for i=1:num_iter
                %% W(t) is fixed. Updating S(t)
        Y = [];
        for j=1:K
           cntr = tau(j);
           window = x(1,cntr-L/2:cntr+L/2-1);
           Y = [Y;window];
        end
        S = (transpose(alpha) * Y)/(transpose(alpha)*alpha);
        s = transpose(S);

                %% S(t) is fixed. Updating W(t)
        R = corr_xy(s,x);
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
end

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