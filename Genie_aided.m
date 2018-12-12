% Genie-aided LS:
% The support is assumed to be known.Using the informatino,the pilot is
% designed to train the subspace the support spans,and LS estimation is
% used.This ideal scheme provides a lower bound for MSE.

function [NMSE, NBG, Tp] = Genie_aided()

% NMSE :       Normailized mean squared error
% NBG  :       Normalized beamforming gain

M = 100;    % ��վ��������
k = 40;     % ϡ����
ks = 3;
SNR = 20;   % �����[dB]
Tp = 40:10:100;   % ��Ƶ����
len = length(Tp);
Pshi = DFTM(M);
s = zeros(M,1);
temp = randperm(M);
index = temp(1:k);
s(index) = mcrandn(k,1);
h = Pshi*s;
pul = 10^(SNR/10);
NMSE = zeros(len,1);
I = eye(k);     % k*k��С��λ����
for i = 1:len
    tp = Tp(i);
    XT_ = I(:,mod(1:tp,k)+1);
    y = sqrt(pul)*s(index)'*XT_+mcrandn(1,tp);
    y = y';
    s_ = (XT_*XT_')^(-1)*XT_*y/sqrt(pul);
    h_ = sum(Pshi(:,index)*diag(s_),2);
    NMSE(i) = (norm(h-h_)/norm(h))^2;
    NBG(i) = abs(h'*h_)/norm(h)/norm(h_);
end






