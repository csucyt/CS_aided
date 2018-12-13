% OMP:
% A random Gaussian matrix of size M*Tp is sequence of length Tp is used
% for training and LS is used as an estimation filter.

function [NMSE, NBG] = omp_measure(M,k,ks,Tp,SNR,L)

% Input:
%      M    :  ��������
%      k    :  ϡ��̶�
%      ks   :  ��ӳ�ŵ��仯����
%      Tp   :  ��Ƶ����
%      SNR  :  �����
%      L    :  ��˥���ŵ�����
% Output:
%      NMSE :       Normailized mean squared error
%      NBG  :       Normalized beamforming gain

if nargin < 6
    L = 100;
end

Pshi = DFTM(M);
s = zeros(M,1);
temp = randperm(M);
index = temp(1:k);      % ��һ�����֧������
realIndex = index;      % ��ǰ��ʵ�ʵ�֧������
s(index) = mcrandn(k,1);
h = Pshi*s;
pul = dBChange(SNR,'dB2SNR');
NMSE = 0;
NBG = 0;
XT = randn(M,Tp);
XT = XT*sqrt(Tp)/sqrt(trace(XT'*XT));  % let trace(XT'*XT) = Tp;

for i = 1 : L
    y = sqrt(pul)*XT.'*conj(Pshi)*s+mcrandn(Tp,1);
    s_ = omp(sqrt(pul)*XT.'*conj(Pshi),y,k);
    h_ = Pshi*s_;
    NMSE = NMSE + (norm(h-h_)/norm(h))^2;
    NBG = NBG + norm(h'*h_)/norm(h)/norm(h_);
    [s,h,realIndex] = channelChange(s,realIndex,ks);
end
NMSE = NMSE/L;
NBG = NBG/L;

