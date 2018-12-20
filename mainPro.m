M = 100;
k = 40;
ks = 3;
Tp = 40:10:100;
SNR = 20;
len = length(Tp);
NMSE_GA = zeros(1,len);    % Genie_aided
NBG_GA = zeros(1,len);
NMSE_SLS = zeros(1,len);   % StaticLS
NBG_SLS = zeros(1,len);
NMSE_omp = zeros(1,len);   % omp
NBG_omp = zeros(1,len);
NMSE_CSA = zeros(1,len);   % CS_aided
NBG_CSA = zeros(1,len); 
NMSE_MSP = zeros(1,len);   % M-SP
NBG_MSP = zeros(1,len);

for i = 1 : length(Tp)
    [NMSE_GA(i),NBG_GA(i)] = Genie_aided(M,k,ks,Tp(i),SNR);
    [NMSE_SLS(i),NBG_SLS(i)] = StaticLS(M,k,ks,Tp(i),SNR);
    [NMSE_omp(i),NBG_omp(i)] = omp_measure(M,k,ks,Tp(i),SNR);
    [NMSE_CSA(i),NBG_CSA(i)] = CS_aided(M,k,ks,Tp(i),SNR);
    [NMSE_MSP(i),NBG_MSP(i)] = M_SP_measure(M,k,ks,Tp(i),SNR);
end

figure(1)
plot(Tp,NMSE_GA,'*-',Tp,NMSE_SLS,'d-',...
    Tp,NMSE_omp,'>-',Tp,NMSE_CSA,'o-',...
    Tp,NMSE_MSP,'h-')
axis([40 100 0 0.5])
grid on
xlabel('Length of training sequence T_p')
ylabel('Normalized mean squared error(NMSE)')
legend('Genie-aided LS','Static LS','OMP','CS-aided','M-SP')
title({'NMSE versus the length of the training sequence with M=100,',...
    'k=40,ks=3,and SNR=20dB'});

figure(2)
plot(Tp,NBG_GA,'*-',Tp,NBG_SLS,'d-',...
    Tp,NBG_omp,'>-',Tp,NBG_CSA,'o-',...
    Tp,NBG_MSP,'h-')
axis([40 100 0 1])
grid on
xlabel('Length of training sequence T_p')
ylabel('Nromalized beamforming gain')
legend('Genie-aided LS','Static LS','OMP','CS-aided','M-SP')
title({'Normalized beamforming gain versus the length of the training,',...
    'sequence with M=100,k=40,ks=3,and SNR=20dB'});

Tp = 60;
SNR = 0:5:30;
len = length(SNR);
for i = 1:len
    [NMSE_GA(i),NBG_GA(i)] = Genie_aided(M,k,ks,Tp,SNR(i));
    [NMSE_SLS(i),NBG_SLS(i)] = StaticLS(M,k,ks,Tp,SNR(i));
    [NMSE_omp(i),NBG_omp(i)] = omp_measure(M,k,ks,Tp,SNR(i));
    [NMSE_CSA(i),NBG_CSA(i)] = CS_aided(M,k,ks,Tp,SNR(i));
    [NMSE_MSP(i),NBG_MSP(i)] = M_SP_measure(M,k,ks,Tp,SNR(i));
end
figure(3)
plot(SNR,NMSE_GA,'*-',SNR,NMSE_SLS,'d-',...
    SNR,NMSE_omp,'>-',SNR,NMSE_CSA,'o-',...
    SNR,NMSE_MSP,'h-')
axis([0 30 0 1])
grid on
xlabel('SNR(dB)')
ylabel('Normalized mean squared error(NMSE)')
legend('Genie-aided LS','Static LS','OMP','CS-aided','M-SP')
title({'NMSE versus SNR with M=100,k=40,ks=3,and Tp=60'})

Tp = 80;
SNR = 20;
kse = 1:1:10;
len = length(kse);
for i = 1:len
    [NMSE_omp(i),NBG_omp(i)] = omp_measure(M,k,ks,Tp,SNR,100,kse(i));
    [NMSE_CSA(i),NBG_CSA(i)] = CS_aided(M,k,ks,Tp,SNR,100,kse(i));
    [NMSE_MSP(i),NBG_MSP(i)] = M_SP_measure(M,k,ks,Tp,SNR,100,kse(i));
end
figure(4)
plot(kse,NMSE_omp,'>-',kse,NMSE_CSA,'o-',...
    kse,NMSE_MSP,'h-')
axis([1 10 0 0.7])
grid on
xlabel('Assumed size of common support ks')
ylabel('Normalized mean squared error(NMSE)')
legend('OMP','CS-aided','M-SP')
title({'NMSE versus the mismatched parameter k_s^e with M=100,k=40,ks=3,and SNR=20'})