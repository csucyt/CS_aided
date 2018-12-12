% In this model,an angular channel vector is assumed to have exactly k
% nonzero elements,while the remaining M-k elements are zero.The support is
% randomly selected such that k-ks indices are preserved from the previous
% support.

function [ns,nh,nindex,info] = channelChange(s,index,ks)
% Input:
%    s    :    �Ƕ���ϵ��
%    index:    ϡ��ϵ����λ��
%    ks   :    ϡ��仯��Χ
% Output:
%    ns   :    ��һ����ĽǶ���ϵ��
%    ks   :    ��һ������ŵ�
%    info :    �任��Ϣ

k = length(index);    % ϡ����
len = length(s);
nochange = k-ks;   

temp = randperm(k);   % û�б任��λ��
nochangeIndex = index(temp(1:nochange));

temp = [];
temp = setdiff(1:len,nochangeIndex);
tempi = randperm(len-k);
changeIndex = temp(tempi(1:ks));

ns = zeros(len,1);
ns(nochangeIndex) = mcrandn(nochange,1);
ns(changeIndex) = mcrandn(ks,1);
nindex = [nochangeIndex,changeIndex];
nh = DFTM(len)*ns;
info = [];