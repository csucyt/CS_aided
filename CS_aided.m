% Compressed Sensing-Aided:
% Exploits the observation that the channel statistic change slowly in
% time.
% Utilizing a conventional least squares approach and a CS technique
% simultaneously.

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
Xd_ = eye(k);
index_d = index;
index_s = setdiff(1:M,index_d);
s_ = zeros(M,1);
NMSE = zeros(len,1);
NBG = zeros(len,1);
cindex = index;

for i = 2:len
    cs = s;
    ch = h;
    cindex = index;
    index_d = index;
    index_s = setdiff(1:M,index_d);
    tp = Tp(i);
    Phi = randn(tp-k,M-k);
    Phi = Phi*sqrt(tp-k)/sqrt(trace(Phi'*Phi));
    for j = 1 : L
        yd = sqrt(pul)*cs(index_d)'*Xd_+mcrandn(1,k);
        yd = yd.';
        ys = sqrt(pul)*cs(index_s)'*Phi'+mcrandn(1,tp-k);
        ys = ys.';
        sd_ = Xd_*conj(yd)/sqrt(pul);
        ss_ = omp(sqrt(pul)*conj(Phi),ys,ks);
        ss_ = conj(ss_);
        s_(index_d) = sd_;
        s_(index_s) = ss_;
        [s_,cindex] = maxN(s_,k);
        h_ = Pshi*s_;
        NMSE(i) = NMSE(i) + (norm(ch-h_)/norm(ch))^2;
        NBG(i) = NBG(i) + norm(ch'*h_)/norm(ch)/norm(h_);
        index_d = cindex;
        index_s = setdiff(1:M,index_d);
        [cs,ch,cindex] = channelChange(cs,cindex.',ks);
    end
    NMSE(i) = NMSE(i)/L;
    NBG(i) = NBG(i)/L;
end
%
% h_ = Pshi*s_;
% norm(h-h_)/norm(h)
% norm(h'*h_)/norm(h)/norm(h_)
 plot(Tp,NMSE,'*-')
axis([50 100 0 1])
hold on
plot(Tp,NBG,'*-')
axis([50 100 0 1])



