% Genie-aided LS:
% The support is assumed to be known.Using the informatino,the pilot is
% designed to train the subspace the support spans,and LS estimation is
% used.This ideal scheme provides a lower bound for MSE.

function [NMSE, NBG] = Genie_aided(M,k,ks,Tp,SNR,L)

% Input:
%      M    :  天线数量
%      k    :  稀疏程度
%      ks   :  反映信道变化快慢
%      Tp   :  导频数量
%      SNR  :  信噪比
%      L    :  块衰落信道数量
% Output:
%      NMSE :       Normailized mean squared error
%      NBG  :       Normalized beamforming gain

if nargin < 6
    L = 100;
end

Pshi = DFTM(M);
s = zeros(M,1);
temp = randperm(M);
index = temp(1:k);         % 第一个块的支撑区域
s(index) = mcrandn(k,1);
h = Pshi*s;
pul = dBChange(SNR,'dB2SNR');
I = eye(k);     % k*k大小单位矩阵
NMSE = 0;
NBG = 0;
realIndex = index;     % 实际信号的支撑区域

for i = 1 : L         % 信道变换
    XT_ = I(:,mod(1:Tp,k)+1);
    y = sqrt(pul)*s(realIndex)'*XT_+mcrandn(1,Tp);
    y = y';
    s_ = (XT_*XT_')^(-1)*XT_*y/sqrt(pul);
    h_ = sum(Pshi(:,realIndex)*diag(s_),2);
    NMSE = NMSE + (norm(h-h_)/norm(h))^2;
    NBG = NBG + abs(h'*h_)/norm(h)/norm(h_);
    [s,h,realIndex] = channelChange(s,realIndex,ks);  % 根据第一个块产生变换
end

NMSE = NMSE/L;
NBG = NBG/L;






