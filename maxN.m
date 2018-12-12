function [ss,cindex] = maxN(s,num)

ss = zeros(size(s));
[~,index] = sort(abs(s),'descend');
cindex = index(1:num);
ss(cindex) = s(cindex);