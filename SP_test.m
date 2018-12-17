clear 
clc
N = 256;
m = 128;
K = 1:3:64;
iter = 500;
Phi = randn(m,N);
error_SP = zeros(1,length(K));
error_OMP = zeros(1,length(K));
for i = 1 : length(K)
    for j= 1 : iter
        x = zeros(N,1);
        temp = randperm(N);
        index = temp(1:K(i));
        x(index) = 1;
        y = Phi*x;
        x_SP = SP(K(i),Phi,y);
        x_OMP = omp(Phi,y,K(i));
        correct_SP = sum(x-x_SP<1e-8);
        correct_OMP = sum(x-x_OMP<1e-8);
        if correct_SP ~= N
            error_SP(i) = error_SP(i) + 1;
        end
        if correct_OMP ~= N
            error_OMP(i) = error_OMP(i) + 1;
        end
    end
    error_SP(i) = error_SP(i)/iter;
    error_OMP(i) = error_OMP(i)/iter;
end
figure(2)
plot(K,1-error_SP,'*-',K,1-error_OMP,'h-')
xlabel('Signal Sparsity:K')
ylabel('Frequency of Exact Reconstruction');
title('Reconstruction Rate(500 Realizations):m=128,N=256\n 2')
grid on
legend('Subspace Pursuit(SP)','Standard OMP')

clear 
clc
N = 256;
m = 128;
K = 1:3:64;
iter = 500;
Phi = randn(m,N);
error_SP = zeros(1,length(K));
error_OMP = zeros(1,length(K));
for i = 1 : length(K)
    for j= 1 : iter
        x = zeros(N,1);
        temp = randperm(N);
        index = temp(1:K(i));
        x(index) = randn(K(i),1);
        y = Phi*x;
        x_SP = SP(K(i),Phi,y);
        x_OMP = omp(Phi,y,K(i));
        correct_SP = sum(x-x_SP<1e-8);
        correct_OMP = sum(x-x_OMP<1e-8);
        if correct_SP ~= N
            error_SP(i) = error_SP(i) + 1;
        end
        if correct_OMP ~= N
            error_OMP(i) = error_OMP(i) + 1;
        end
    end
    error_SP(i) = error_SP(i)/iter;
    error_OMP(i) = error_OMP(i)/iter;
end
figure(1)
plot(K,1-error_SP,'*-',K,1-error_OMP,'h-')
xlabel('Signal Sparsity:K')
ylabel('Frequency of Exact Reconstruction');
title('Reconstruction Rate(500 Realizations):m=128,N=256\n 2')
grid on
legend('Subspace Pursuit(SP)','Standard OMP')