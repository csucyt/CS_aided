% OMP:
% A random Gaussian matrix of size M*Tp is sequence of length Tp is used
% for training and LS is used as an estimation filter.

L = 1000;   % 块数量
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

for i = 1:len
    cs = s;
    ch = h;
    cindex = index;
    tp = Tp(i);
    XT = randn(M,tp);
    XT = XT*sqrt(tp)/sqrt(trace(XT'*XT));
    for j = 1 : L
        y = sqrt(pul)*XT.'*conj(Pshi)*cs+mcrandn(tp,1);
        s_ = omp(sqrt(pul)*XT.'*conj(Pshi),y,k);
        h_ = Pshi*s_;
        NMSE(i) = NMSE(i)+(norm(ch-h_)/norm(ch))^2;
        NBG(i) = NBG(i)+norm(ch'*h_)/norm(ch)/norm(h_);
        [cs,ch,cindex] = channelChange(cs,cindex,ks);
    end
    NMSE(i) = NMSE(i)/L;
    NBG(i) = NBG(i)/L;
end
