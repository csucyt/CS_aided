% Static LS:
% Only the support of the first block is know.The pilot is designed using
% the outdate support for all time and LS is used for estimation.

L = 100;   % 块数量
M = 100;
k = 40;
ks = 3;
SNR = 20;   % 信噪比[dB]
Tp = 40:10:100;   % 导频数量
len = length(Tp);
Pshi = DFTM(M);
s = zeros(M,1);
temp = randperm(M);
index = temp(1:k);
s(index) = mcrandn(k,1);
h = Pshi*s;
pul = 10^(SNR/10);
NMSE = zeros(len,1);
NBG = zeros(len,1);
I = eye(k);     % k*k大小单位矩阵

for i = 1 : len
    cs = s;
    ch = h;
    cindex = index;
    tp = Tp(i);
    XT_ = I(:,mod(1:tp,k)+1);
    for j = 1 : L
        y = sqrt(pul)*cs(index)'*XT_+mcrandn(1,tp);
        y = y';
        s_ = (XT_*XT_')^(-1)*XT_*y/sqrt(pul);
        h_ = sum(Pshi(:,index)*diag(s_),2);
        NMSE(i) = NMSE(i)+(norm(ch-h_)/norm(ch))^2;
        NBG(i) = NBG(i)+abs(ch'*h_)/norm(ch)/norm(h_);
%         ks = randi(3);
        [ns,nh,nindex] = channelChange(cs,cindex,ks);
        cs = ns; ch = nh;  cindex = nindex;
    end
    NMSE(i) = NMSE(i)/L;
    NBG(i) = NBG(i)/L;
end
