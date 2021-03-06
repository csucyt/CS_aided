% Subspace Pursuit Algorithm
function x = SP(K,Phi,y)
% Formulation:
%      y = Phi*x;
% Input:
%      K    :  a signal sparsity level
%      Phi  :  a sampling matrix
%      y    :  the observation vector
% Output:
%      x    :  the estimated signal

% Initialization
[~,N] = size(Phi);
% T = maxNMag(Phi.'*y,K);
T = maxNMag(Phi'*y,K);
yr = resid(y,Phi(:,T));

% Iteration
iter = 1;
while true
    T_last = T;     % 表示上一次迭代的支撑集
    yr_last = yr;   % 表示上一次迭代的冗余
    %     T_ = [T_last.', maxNMag(Phi.'*yr_last,K).'];
    %     T_ = T_.';
    T_ = [T_last.', maxNMag(Phi'*yr_last,K).'];
    T_ = T_.';
    % T_ = union(T_last,maxNMag(Phi.'*yr_last,K));
    xp = pseudoInv(Phi(:,T_))*y;
    %     xp = proj(y,Phi(:,T_));
    temp = maxNMag(xp,K);
    T = T_(temp);
    yr = resid(y,Phi(:,T));
    iter = iter+1;
    if norm(yr) > norm(yr_last) || abs(norm(yr)-norm(yr_last))<1e-10 
        T = T_last;
        break;
    end
end

x = zeros(N,1);
x(T) = pseudoInv(Phi(:,T))*y;

end

function T = maxNMag(v,K)
% K indices corresponding to the largest magnitude entries in the vector v.
[~,I] = sort(abs(v),'descend');
T = I(1:K);
end

function yr = resid(y,Phi)
% the residue vector of the projection
yr = y - proj(y,Phi);
end


function yp = proj(y,Phi)
% the projection of y onto span(Phi)
yp = Phi*pseudoInv(Phi)*y;
end


function pinv_Phi = pseudoInv(Phi)
% the pseduo-inverse of the matrix Phi
% pinv_Phi = (Phi.'*Phi)^(-1)*Phi.';
pinv_Phi = pinv(Phi);
end

