clc
num=1.67e8;
den=[0.03 2.02 35 33.33];
sys_org=tf(num,den)
[A,B,C,D]=tf2ss(num,den)
[v,d]=eig(A)
%------------------------
Qc = B * B';
Qo = C' * C;
Wc = lyap(A,A',Qc);
Wo = lyap(A',A,Qo);
%------------------------
Vc = orth(Wc);
Lc=(Vc' * Wc * Vc)^(1/2);
%------------------------
W = (Vc * Lc)' * Wo * (Vc * Lc);
P = orth(W);
L = P' * W * P;
%------------------------
S = Vc * Lc * P * L^(-1/2);
%------------------------
Wc2=S^(-1) * Wc * (inv(S))';
Wo2=(inv(S))' * Wo * S;
%------------------------
As=S ^(-1) * A * S
Bs=S^(-1) * B
Cs=C * S
Ds=0;
%--------------
r=max(size(L))-1;      % bac cua' he can bang noi
At=As(1:r,1:r);
Bt=Bs(1:r,:);
Ct=Cs(:,1:r);
Dt=0;
%------------------------
sys_mod = tf(ss(At,Bt,Ct,Dt))
%------------------------
step(sys_org,sys_mod), grid on, hold on
%bode(sys3,sys2), grid on, hold on
