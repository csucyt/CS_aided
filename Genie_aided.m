% Genie-aided LS:
% The support is assumed to be known.Using the informatino,the pilot is
% designed to train the subspace the support spans,and LS estimation is
% used.This ideal scheme provides a lower bound for MSE.

function [NMSE, NBG] = Genie_aided(M,k,ks,Tp,SNR,L)

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
index = temp(1:k);         % ��һ�����֧������
s(index) = mcrandn(k,1);
h = Pshi*s;
pul = dBChange(SNR,'dB2SNR');
I = eye(k);     % k*k��С��λ����
NMSE = 0;
NBG = 0;
realIndex = index;     % ʵ���źŵ�֧������

for i = 1 : L         % �ŵ��任
    XT_ = I(:,mod(1:Tp,k)+1);
    y = sqrt(pul)*s(realIndex)'*XT_+mcrandn(1,Tp);
    y = y';
    s_ = (XT_*XT_')^(-1)*XT_*y/sqrt(pul);
    h_ = sum(Pshi(:,realIndex)*diag(s_),2);
    NMSE = NMSE + (norm(h-h_)/norm(h))^2;
    NBG = NBG + abs(h'*h_)/norm(h)/norm(h_);
    [s,h,realIndex] = channelChange(s,realIndex,ks);  % ���ݵ�һ��������任
end

NMSE = NMSE/L;
NBG = NBG/L;






