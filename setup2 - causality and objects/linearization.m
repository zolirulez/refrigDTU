clearvars
close all
syms KvV KvG TauV TauG R VGC VR VC VMT DdDp1 DdDh1 DdDp2 DdDh2 DdDpR DdDhR DdDpC DdDhC eS DdDpMT DdDhMT
syms p1 p2 pR pC h1 h2 hMT hL hG hR hC hCi hF d1 d2 dR dG dC dMT dA DQ1 DQ2 DQC DQF CRV CRG fMT fA DVA DmV DmMT DmG DmL DmCi DmCo DmLT DmF DVA
syms s s0 k cp TA1 TA2 T1 T2 TA0 DTDp1 DTDh1 DTDp2 DTDh2
% ----------------------- STATIC EQUATIONS --------------------------------
% TODO: grouping and ordering!
% ---------------------- DYNAMIC EQUATIONS --------------------------------
% High pressure valve
DDmV = 1/TauV*(-DmV + CRV*KvV*sqrt(d2*(p2 - pR)));
% Gas cooler
s = 1/(1/(s0 + k*DVA) + 1/(dA*DVA*cp));
w = dA*DVA*cp/s;
T1 = DTDp1*p1 + DTDh1*h1; % Note: constant is missing
T2 = DTDp2*p2 + DTDh2*h2; % Note: constant is missing
TA1 = 1/(w+1)*T1 + 1/(1/w+1)*TA0;
TA2 = 1/(w+1)*T2 + 1/(1/w+1)*TA1;
DQ1 = (TA1 - T1)*s;
DQ2 = (TA2 - T2)*s;
Dm21 = 1/R*sqrt(d1*(p1-p2));
Dd1 = 1/VGC*(DmMT-Dm21);
Dd2 = 1/VGC*(Dm21-DmV);
Joining1 = [-1 d1; DdDp1 DdDh1];
Joining2 = [-1 d2; DdDp2 DdDh2];
DPsi1 = 1/VGC*(DmMT*hMT-Dm21*h1 + DQ1);
DPsi2 = 1/VGC*(Dm21*h1-DmV*h2 + DQ2);
Dph1 = Joining1\[DPsi1; Dd1];
Dph2 = Joining2\[DPsi2; Dd2];
% Receiver Joint
DmL = DmF + DmCi;
% Receiver
DdR = 1/VR*(DmV-DmL-DmG);
JointR = [-1 dR; DdDpR DdDhR];
DPsiR = 1/VR*(DmV*h2 - DmL*hL - DmG*hG);
DphR = JointR\[DPsiR; DdR];
% Receiver valve: density is set as constant
DDmG = 1/TauG*(-DmG + CRG*KvG*sqrt(dG*(pR - pC)));
% Freezer BC 
DmF = DQF/(hF-hL);
% MT Joint 
DmCo = - DmG + DmMT - DmF;
hMTi = (DmG*hG + DmCo*hC + DmF*hF)/DmMT;
dMT = DdDpMT*pC + DdDhMT*hMTi; % Note: constant is missing
% Compressor
DmMT = fMT*dMT*VMT;
hMT = hMTi + DhDpMT*(p1 - pC)/eS;
% Cooler volume
DmCi = DQC/(hC-hL);
DdC = 1/VC*(DmCi-DmCo);
JoiningC = [-1 dC; DdDpC DdDhC];
DPsiC = 1/VGC*(DmCi*hCi-DmCo*hC);
DphC = JoiningC\[DPsiC; DdC];
% Fan (A refers to air)
DDVA = dA*DVA*fA;
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
B = jacobian(Dx,u);
G = jacobian(Dx,d);
% Substitutions check the values TODO
A = subs(A,{p1 p2 pR pEC},num2cell([86 84.8 38 30]*10^5));
A = subs(A,{h1 h2 hMT hL hG hR},num2cell([500 350 525 250 400 300]*10^3));
A = subs(A,{d1 d2 dR dG},num2cell([2.8 7 2.2 1]*10^2));
A = subs(A,{DmV DmMT DmG DmL},num2cell([0.321 0.321 0.123 0.198]*10^0));
A = subs(A,{CRV CRG},num2cell([0.25 0.25]*10^0));
A = subs(A,{DQ},num2cell([-74.3]*10^3));
A = subs(A,{KvV KvG},num2cell([0.8 2]*8.7841e-06));
A = subs(A,{TauV TauG},num2cell([0.1 0.1]*10^0));
A = subs(A,{R},num2cell([1.5]*10^5)); % TODO
A = subs(A,{VGC VR},num2cell([19.2 133]*10^-3));
A = subs(A,{DdDp1 DdDp2 DdDpR},num2cell([3 3 4]*10^-5));
A = subs(A,{DdDh1 DdDh2 DdDhR},num2cell([-2 -4 -1.5]*10^-3));
A = double(A);
% Normalizing with maximum deviations
% x = [p1; h1; d1; p2; h2; d2; DmV; pR; hR; dR; DmG];
T = diag(1./[20*10^5; 50*10^3; 0.5*10^2; 20*10^5; 50*10^3; 0.5*10^2; 0.2;...
    20*10^5; 50*10^3; 0.5*10^2; 0.2]);
A = T*A/T;
B = T*B;
G = T*G;
expA = sign(A).*exp(abs(A));
expA(A==0) = 0;
figure(1)
imagesc(expA,[-10 10])
colormap('jet')
colorbar
xlabel(char(x))
ylabel(['Diff' char(x)])
title('Pieceswise exponential of normalized system A, bounded by -10...10')
% Modal analysis
[~,eigA,W] = eig(A);
% Separating real and imaginary value?
W = W';
% W = [W(1:6,:); real(W(7,:)); imag(W(7,:)); W(9:end,:)];
% Is this right? should not be first the inverse, then the reordering?
% (This was just a copy paste from the exercise on LCD2)
%V = [V(:,1) real(V(:,2)) imag(V(:,2)) real(V(:,4)) imag(V(:,4)) V(:,6:end)];
figure(2)
subplot(131)
imagexpeigA = sign(imag(eigA)).*exp(abs(imag(eigA)));
imagexpeigA(imag(eigA)==0) = 0;
imagesc(imagexpeigA,[-10^0.5 10^0.5])
colormap('jet')
colorbar
xlabel('z')
ylabel('Dz')
title('imagexpeigA')
subplot(132)
realexpeigA = sign(real(eigA)).*exp(abs(real(eigA)));
realexpeigA(real(eigA)==0) = 0;
imagesc(realexpeigA,[-10 10])
colormap('jet')
colorbar
xlabel('z')
ylabel('Dz')
title('realexpeigA')
subplot(133)
imagesc(W,[-10^0.2 10^0.2])
% imagesc(V',[-1 1])
colormap('jet')
colorbar
xlabel(char(x))
ylabel('z')
title('W^T')