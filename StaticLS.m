% Static LS:
% Only the support of the first block is know.The pilot is designed using
% the outdate support for all time and LS is used for estimation.

function [NMSE, NBG] = StaticLS(M,k,ks,Tp,SNR,L)

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
% realIndex = index;         % ��ǰ���֧������
s(index) = mcrandn(k,1);
h = Pshi*s;
pul = 10^(SNR/10);
NMSE = 0;
NBG = 0;
I = eye(k);     % k*k��С��λ����
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
