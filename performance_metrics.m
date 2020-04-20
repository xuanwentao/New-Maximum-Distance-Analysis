function [E_aad,E_aid,E_sad,E_sid,E_rmse,Aest] = performance_metrics(A,Aest,mixed,abf,M,N,D,c)
%%
warning off;
AA = [1e-5*Aest;ones(1,length(Aest(1,:)))];
sest = zeros(length(Aest(1,:)),M*N);
for j=1:M*N
    r = [1e-5*mixed(:,j); 1];
%   s_fcls(:,j) = nnls(AA,r);
    sest(:,j) = lsqnonneg(AA,r);
end
CRD = corrcoef([A Aest]);%求两个矩阵的相关系数
DD = abs(CRD(c+1:2*c,1:c));% 求绝对值
perm_mtx = zeros(c,c);
aux=zeros(c,1);
for i=1:c
    [ld cd]=find(max(DD(:))==DD); 
    ld=ld(1);cd=cd(1); % in the case of more than one maximum
    perm_mtx(ld,cd)=1; 
    DD(:,cd)=aux; DD(ld,:)=aux';
end
Aest = Aest*perm_mtx;
sest = sest'*perm_mtx;
Sest = reshape(sest,[M,N,c]);
sest = sest';
E_rmse = sqrt(sum(sum(((abf-sest).*(abf-sest)).^2))/(M*N*c))

% the angle between abundances
nabf = diag(abf*abf'); 
nsest = diag(sest*sest');
ang_beta = 180/pi*acos( diag(abf*sest')./sqrt(nabf.*nsest));
E_aad = mean(ang_beta.^2)^.5

% cross entropy between abundance 丰度之间的交叉熵
E_entropy = sum(abf.*log((abf+1e-9)./(sest+1e-9))) + sum(sest.*log((sest+1e-9)./(abf+1e-9)));
E_aid = mean(E_entropy.^2)^.5

% the angle between material signatures
nA = diag(A'*A);
nAest = diag(Aest'*Aest);
ang_theta = 180/pi*acos( diag(A'*Aest)./sqrt(nA.*nAest) );
E_sad = mean(ang_theta.^2)^.5

% the spectral information divergence
pA = A./(repmat(sum(A),[length(A(:,1)) 1]));
qA = Aest./(repmat(sum(Aest),[length(A(:,1)) 1])); 
qA = abs(qA); 
SID = sum(pA.*log((pA+1e-9)./(qA+1e-9))) + sum(qA.*log((qA+1e-9)./(pA+1e-9)));
E_sid = mean(SID.^2)^.5