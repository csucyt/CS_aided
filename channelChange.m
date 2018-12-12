% In this model,an angular channel vector is assumed to have exactly k
% nonzero elements,while the remaining M-k elements are zero.The support is
% randomly selected such that k-ks indices are preserved from the previous
% support.

function [ns,nh,nindex,info] = channelChange(s,index,ks)
% Input:
%    s    :    角度域系数
%    index:    稀疏系数的位置
%    ks   :    稀疏变化范围
% Output:
%    ns   :    下一个块的角度域系数
%    ks   :    下一个块的信道
%    info :    变换信息

k = length(index);    % 稀疏性
len = length(s);
nochange = k-ks;   

temp = randperm(k);   % 没有变换的位置
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