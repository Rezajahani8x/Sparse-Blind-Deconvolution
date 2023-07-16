                            %% Sparse Blind Deconvolution
                            %% Part 3
file = load('Data.mat');
X = file.X;
L = 100;
K = 5;
[~,T] = size(X);
N = 2;

            %% Part 3
x1 = X(1,:);
x2 = X(2,:);
figure(1);
subplot(2,1,1);
plot(1:T,x1);
xlabel('t');
ylabel('Amp');
title('Channel 1 of X');
grid on;
subplot(2,1,2);
plot(1:T,x2);
xlabel('t');
ylabel('Amp');
title('Channel 2 of X');
grid on;

num_iter = 2;
[S,alpha,tau] = MCSBD(X,K,L,T,N,num_iter);
figure(2);
subplot(2,1,1);
plot(1:L,S(1,:));
xlabel('t');
ylabel('Amp');
title('Spike of S_1(t)');
grid on;
subplot(2,1,2);
plot(1:L,S(2,:));
xlabel('t');
ylabel('Amp');
title('Spike of S_2(t)');
grid on;
disp("W_1(t)");
disp("alpha | tau");
disp(vpa([round(alpha(:,1),2) tau(:,1)]));
disp("W_2(t)");
disp("alpha | tau");
disp(vpa([round(alpha(:,2),2) tau(:,2)]));


            %% Local Necessary Functions
function [S,alpha_vals,tau_vals] = MCSBD(X,K,L,T,N,num_iter)
    A = [0.3 0.7;0.4 0.6];
    for i=1:num_iter
            %% A is fixed. Updating B
       A_psuedo_inv = inv(transpose(A)*A)*transpose(A);
       B = A_psuedo_inv * X;
       alpha_vals = zeros(K,N);
       tau_vals = zeros(K,N);
       S = zeros(N,L);
       for n=1:N
           x = B(n,:);
           [s,alpha,tau] = SCSBD(x,K,L,T,num_iter);
           S(n,:) = s;
           alpha_vals(:,n) = alpha;
           tau_vals(:,n) = tau;
       end
       B = zeros(N,T);
       for n=1:N
          s_n = S(n,:);
          for k=1:K
             tau_k = tau_vals(k,n);
             alpha_k = alpha_vals(k,n);
             b_k = alpha_k * s_n;
             time_slot = tau_k - L/2:tau_k + L/2 - 1;
             B(n,time_slot) = b_k;
          end
       end
            %% B is fixed. Updating A
       B_psuedo_inv = transpose(B) * inv(B*transpose(B));
       A = X * B_psuedo_inv;
    end
end

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
