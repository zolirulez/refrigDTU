clearvars
syms KvV KvG TauV TauG R VGC VR DdDp1 DdDh1 DdDp2 DdDh2 DdDpR DdDhR
syms p1 p2 pR pEC h1 h2 hMT hL hG hR d1 d2 dR DQ CRV CRG DmV DmMT DmG DmL 
% High pressure valve
DDmV = 1/TauV*(-DmV + CRV*KvV*sqrt(d2*(p2 - pR)));
% Gas cooler
Dm21 = 1/R*sqrt(d1*(p1-p2));
Dd1 = 1/VGC*(DmMT-Dm21);
Dd2 = 1/VGC*(Dm21-DmV);
Joint1 = [-1 d1; DdDp1 DdDh1];
Joint2 = [-1 d2; DdDp2 DdDh2];
DPsi1 = 1/VGC*(DmMT*hMT-Dm21*h1 + DQ/2);
DPsi2 = 1/VGC*(Dm21*h1-DmV*h2 + DQ/2);
Dph1 = Joint1\[DPsi1; Dd1];
Dph2 = Joint2\[DPsi2; Dd2];
% Receiver
DdR = 1/VR*(DmV-DmL-DmG);
JointR = [-1 dR; DdDpR DdDhR];
DPsiR = 1/VR*(DmV*h2 - DmL*hL - DmG*hG);
DphR = JointR\[DPsiR; DdR];
% Receiver valve
DDmG = 1/TauG*(-DmG + CRG*KvG*sqrt(dR*(pR - pEC)));
% Augmentation
Dx = [Dph1; Dd1; Dph2; Dd2; DDmV; DphR; DdR; DDmG];
x = [p1; h1; d1; p2; h2; d2; DmV; pR; hR; dR; DmG];
u = [CRV; CRG];
d = [DmMT; hMT; DmL; hL; pEC; DQ];
% Assumptions, do we need them?
%assume([KvV TauV R VGC VR DdDp1 -DdDh1 DdDp2 -DdDh2 DdDpR -DdDhR...
%    p1 p2 pR h1 h2 hMT hL hG hR d1 d2 dR -DQ CRV DmV DmMT DmG DmL]>0)
% Jacobians
A = jacobian(Dx,x);
size(A)
B = jacobian(Dx,u)
G = jacobian(Dx,d)
