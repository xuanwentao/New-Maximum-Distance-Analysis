%%
clear all;
clc;
tic;
load('data5');
%%
[Np,Nb] = size(mixed);
R=randperm(Np);
R = R(1:3);
E = [];
%%
for i = 1:c
R = [R,E];
if length(E)<3 %||length(E)==3
R = R(end-2:end);
newM = mixed(R,:);
v_p2 = mixed-mixed(R(1),:);
else 
  newM = mixed(E,:);
  v_p2 = mixed-mixed(E(1),:);
end
r =rank(newM);
Gen_sol = null(newM,r);
n = Gen_sol(:,1);
[V_p_123,E_lab] = max(abs(v_p2*n)/norm(n));
if V_p_123>0.00001
    E(:,i) = E_lab;
else
    break
end
end
%%
mixed = mixed';
Aest = mixed(:,E);
Nx = 58;Ny = 58;Nb = 188;
[E_aad,E_aid,E_sad,E_sid,E_rmse,EM] = performance_metrics(A,Aest,mixed,abf,Nx,Ny,Nb,c);
Result = [E_sid,E_aid,E_rmse];
toc;