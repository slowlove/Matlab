function [As,Bs,Cs,Ds] = mor(num,den)
% Ham giam bac mo hinh
% Ham` tra' ve la tham so' cua' he trong khong gian can bang noi
sys_org=tf(num,den)
step(sys_org), grid on, hold on
[A,B,C,D]=tf2ss(num,den)
[~,d]=eig(A);
if diag(d)<0    % Kiem' tra tinh' dieu khien' duoc cua' he !
% giai he phuong trinh Lyapunov
Qc = B * B';
Qo = C' * C;
Wc = lyap(A,A',Qc);
Wo = lyap(A',A,Qo);
% Cheo' hoa' ma tran
Vc = orth(Wc)
Lc=(Vc' * Wc * Vc)^(1/2)
% cheo' hoa truc giao ma tran
W = (Vc * Lc)' * Wo * (Vc * Lc)
P = orth(W);
L = P' * W * P
%------------------------
S = Vc * Lc * P * L^(-1/2);
%------------------------
Wc2=S^(-1) * Wc * (inv(S))';
Wo2=(inv(S))' * Wo * S;
% He trong khong gian can bang noi
As=S ^(-1) * A * S;
Bs=S^(-1) * B;
Cs=C * S;
Ds=0;
%--------------------
r=max(size(L))-1;
while r>=1
    At=As(1:r,1:r)
    Bt=Bs(1:r,:)
    Ct=Cs(:,1:r)
    Dt=0;
    sys_mod = tf(ss(At,Bt,Ct,Dt))
    step(sys_mod), grid on, hold on
    r=r-1;
end
else
    error ('He khong dieu khien duoc')
end
end