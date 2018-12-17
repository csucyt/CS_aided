% Compressive Sensing With Prior Support Quality Information
% Compressive Sensing Model:
%       Y = Phi*X + N

function [X,T] = M_SP(Y,Phi,s,sc,T0,d,y)
% Input:
%      Y      : measurements
%      Phi    : measurement matrix
%      s      : chunk sparsity level(CSL)
%      sc     : the quality of the prior support T0
%      T0     : the prior support
%      d      : chunk size
%      y      : thershold
% Output:
%      X      : the estimated signal for X
%      T      : the estimated chunk support

[~,L] = size(Y);   %
[~,N] = size(Phi); % N :相当于天线数量
K = N/d;

% Initialization
i = 0;      % iteration index
Ti = [];    % 初始化chunk support
Ri = Y;     % Residue matrix

% Iteration
while true
    
    % (A)Support Merge
    Tb = maxN(Phi'*Ri,T0,sc,d,L);
    T2 = setdiff(1:K,Tb);
    Tc = maxN(Phi'*Ri,T2,s-sc,d,L);
    Ta = union(Ti,union(Tb,Tc));
    
    % (B)LS Estimation
    newTa = extendT(Ta,d);
    Z = zeros(N,L);
    Z(newTa,:) = pinv(Phi(:,newTa))*Y;
    
    % (C)Support Refinement
    Tb = maxN(Z,T0,sc,d,L);
    Tc = maxN(Z,setdiff(1:K,Tb),s-sc,d,L);
    T = union(Tb,Tc);
    
    % (D)Signal Estimation
    newT = extendT(T,d);
    X = zeros(N,L);
    X(newT,:) = pinv(Phi(:,newT))*Y;
    
    % (E)Residue
    R = Y - Phi(:,newT)*X(newT,:);
    
    % (F)Stopping Condition and Output
    if norm(R) < y
        break;
    end
    if norm(R) >= norm(Ri)
        T = Ti;
        if i ~= 0
            X = pX;
        end
        break;
    end
    Ri = R;
    Ti = T;
    i = i+1;
    pX = X;
end

end

function newT = extendT(T,d)
len = length(T);
newT = zeros(1,len*d);
index = 1;
for i = 1 : len
    for j= 1 : d
        newT(index) = (T(i)-1)*d+1 + j-1;
        index = index+1;
    end
end
end

function T = maxN(temp,T0,s,d,L)
% 找到前s个Frobenius norm最大数的Index

len = length(T0);    % 支撑集的大小
Corr = zeros(len,1);
for i = 1 : len
    Corr(i) = norm(temp(T0(i):T0(i)+d-1,1:L));
%     Corr(i) = norm(temp(T0(i):T0(i)+d-1,1:L),'fro');
end
[~,I] = sort(Corr,'descend');
T = T0(I(1:s));

end