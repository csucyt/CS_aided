% Static LS:
% Only the support of the first block is know.The pilot is designed using
% the outdate support for all time and LS is used for estimation.

function [NMSE, NBG] = StaticLS(M,k,ks,Tp,SNR,L)

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
% realIndex = index;         % 当前块的支撑区域
s(index) = mcrandn(k,1);
h = Pshi*s;
pul = 10^(SNR/10);
NMSE = 0;
NBG = 0;
I = eye(k);     % k*k大小单位矩阵
XT_ = I(:,mod(1:Tp,k)+1);


for j = 1 : L
    y = sqrt(pul)*s(index)'*XT_+mcrandn(1,Tp);
    y = y';
    s_ = (XT_*XT_')^(-1)*XT_*y/sqrt(pul);
    h_ = sum(Pshi(:,index)*diag(s_),2);
    NMSE = NMSE + (norm(h-h_)/norm(h))^2;
    NBG = NBG + abs(h'*h_)/norm(h)/norm(h_);
    [s,h,~] = channelChange(s,index,ks);
end
NMSE = NMSE/L;
NBG = NBG/L;
