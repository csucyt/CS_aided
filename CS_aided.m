% Compressed Sensing-Aided:
% Exploits the observation that the channel statistic change slowly in
% time.
% Utilizing a conventional least squares approach and a CS technique
% simultaneously.

function [NMSE, NBG] = CS_aided(M,k,ks,Tp,SNR,L)

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
index = temp(1:k);       % 第一个块支撑区域生成
realIndex = index;       % 当前块的支撑区域
s(index) = mcrandn(k,1);
h = Pshi*s;
pul = 10^(SNR/10);
NMSE = 0;
NBG = 0;
Xd_ = eye(k);
index_d = realIndex;     
index_s = setdiff(1:M,index_d);
s_ = zeros(M,1);

if Tp == k
    [NMSE,NBG] = StaticLS(M,k,ks,Tp,SNR,L);
else
    Phi = randn(Tp-k,M-k);
    Phi = Phi*sqrt(Tp-k)/sqrt(trace(Phi'*Phi));
    for j = 1 : L
        yd = sqrt(pul)*s(index_d)'*Xd_+mcrandn(1,k);
        yd = yd.';
        ys = sqrt(pul)*s(index_s)'*Phi'+mcrandn(1,Tp-k);
        ys = ys.';
        sd_ = Xd_*conj(yd)/sqrt(pul);
        ss_ = omp(sqrt(pul)*conj(Phi),ys,ks);
        ss_ = conj(ss_);
        s_(index_d) = sd_;
        s_(index_s) = ss_;
        [s_,realIndex] = maxN(s_,k);
        h_ = Pshi*s_;
        NMSE = NMSE + (norm(h-h_)/norm(h))^2;
        NBG = NBG + norm(h'*h_)/norm(h)/norm(h_);
        index_d = realIndex;
        index_s = setdiff(1:M,index_d);
        [s,h,~] = channelChange(s,realIndex.',ks);
    end
    NMSE = NMSE/L;
    NBG = NBG/L;
end

