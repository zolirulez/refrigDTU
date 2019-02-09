clearvars
close all
syms KvV KvG TauV TauG TauMT TauA R VGC VR VC VMT DdDp1 DdDh1 DdDp2 DdDh2 DdDpR DdDhR DdDpC DdDhC eS DdDpMT DdDhMT DhDpMT
syms p1 p2 pR pC h1 h2 hMT hL hG hR hC hCi hF d1 d2 dR dG dC dMT dA DQC DQF CRV CRG fMT fA DVA DmV DmMT DmG DmL DmCi DmCo DmLT DmF DVA
syms s s0 k cp TA1 TA2 T1 T2 TA0 DTDp1 DTDh1 DTDp2 DTDh2 MxDVA
% ----------------------- STATIC EQUATIONS --------------------------------
% Gas cooler, heat transfer
s = 1/(1/(s0 + k*DVA) + 1/(dA*DVA*cp));
w = dA*DVA*cp/s;
TA1 = 1/(w+1)*T1 + 1/(1/w+1)*TA0;
TA2 = 1/(w+1)*T2 + 1/(1/w+1)*TA1;
DQ1 = (TA1 - T1)*s;
DQ2 = (TA2 - T2)*s;
% Receiver valve: density is set as constant
% Freezer BC 
DmF = DQF/(hF-hL);
% MT Joint 
DmCo = - DmG + DmMT - DmF;
hMTi = (DmG*hG + DmCo*hC + DmF*hF)/DmMT;
% Compressor
hMT = hMTi + DhDpMT*(p1 - pC)/eS;
% Cooler volume
DmCi = DQC/(hC-hL);
JoiningC = [-1 dC; DdDpC DdDhC];
DPsiC = 1/VGC*(DmCi*hCi-DmCo*hC);
% Receiver Joint
DmL = DmF + DmCi;
% Receiver
JointR = [-1 dR; DdDpR DdDhR];
DPsiR = 1/VR*(DmV*h2 - DmL*hL - DmG*hG);
% Gas cooler, fluid side
Dm21 = 1/R*sqrt(d1*(p1-p2));
Joining1 = [-1 d1; DdDp1 DdDh1];
Joining2 = [-1 d2; DdDp2 DdDh2];
DPsi1 = 1/VGC*(DmMT*hMT-Dm21*h1 + DQ1);
DPsi2 = 1/VGC*(Dm21*h1-DmV*h2 + DQ2);
% Fan (A refers to air)
% ---------------------- DYNAMIC EQUATIONS --------------------------------
% High pressure valve
DDmV = 1/TauV*(-DmV + CRV*KvV*sqrt(d2*(p2 - pR)));
% Gas cooler
Dd1 = 1/VGC*(DmMT-Dm21);
Dd2 = 1/VGC*(Dm21-DmV);
Dph1 = Joining1\[DPsi1; Dd1];
Dph2 = Joining2\[DPsi2; Dd2];
DT1 = [DTDp1 DTDh1]*Dph1;
DT2 = [DTDp2 DTDh2]*Dph2;
% Receiver Joint
% Receiver
DdR = 1/VR*(DmV-DmL-DmG);
DphR = JointR\[DPsiR; DdR];
% Receiver valve: density is set as constant
DDmG = 1/TauG*(-DmG + CRG*KvG*sqrt(dG*(pR - pC)));
% Freezer BC 
% MT Joint 
DdMT = [DdDpMT DdDhMT]*[pC; hMTi];
% Compressor
DDmMT = 1/TauMT*(-DmMT + dMT*VMT*fMT);
% Cooler volume
DdC = 1/VC*(DmCi-DmCo);
DphC = JoiningC\[DPsiC; DdC];
% Fan (A refers to air)
DDVA = 1/TauA*(-DVA + dA*MxDVA*fA);
% ------------------------ LINEARIZATION ----------------------------------
% Augmentation
Dx = [DDVA; Dph1; Dd1; DT1; Dph2; Dd2; DT2; DDmV; DphR; DdR; DDmG; DphC; DdC; DdMT; DDmMT];
x = [DVA; p1; h1; d1; T1; p2; h2; d2; T2; DmV; pR; hR; dR; DmG; pC; hC; dC; dMT; DmMT];
u = [fA; CRV; CRG; fMT];
d = [dA; TA0; DmMT; hL; hG; dG; hCi; hF; DQC; DQF];
% Jacobians
A = jacobian(Dx,x);
B = jacobian(Dx,u);
G = jacobian(Dx,d);
% Substitutions check the values TODO
A = subs(A,{p1 p2 pR pC},num2cell([86 84.8 38 30]*10^5));
A = subs(A,{h1 h2 hMT hL hG hR hC hCi hF},num2cell([450 350 525 212 430 300 460 460 460]*10^3));
A = subs(A,{d1 d2 dR dG dC dMT dA},num2cell([280 700 220 100 80 80 1.25]));
A = subs(A,{TA0 TA1 TA2 T1 T2},num2cell([30 35 40 90 50]+273.15));
A = subs(A,{DmV DmMT DmG DmL DmCi DmCo DmLT DmF},num2cell([0.321 0.321 0.123 0.198 0.151 0.151 0.046 0.046]));
A = subs(A,{DVA MxDVA},{3.33 6.66});
A = subs(A,{CRV CRG fMT fA},num2cell([0.25 0.25 0.25 0.6]));
A = subs(A,{DQC DQF},num2cell([36 11]*10^3));
A = subs(A,{KvV KvG},num2cell([0.8 2]*8.7841e-06));
A = subs(A,{TauV TauG TauMT TauA},num2cell([0.1 0.1 1 1]*10^0));
A = subs(A,{R},num2cell(1.5*10^5));
A = subs(A,{eS},num2cell(0.6));
A = subs(A,{s s0 k cp},num2cell([1000 200 900 1000]));
A = subs(A,{VGC VR VC VMT},num2cell([19.2 133 100 0.1]*10^-3));
A = subs(A,{DdDp1 DdDp2 DdDpR DdDpC DdDpMT},num2cell([3 3 4 2.7 2.4]*10^-5)); % TODO
A = subs(A,{DdDh1 DdDh2 DdDhR DdDhC DdDhMT},num2cell([-2 -4 -1.5 -42 -33]*10^-3)); % TODO
A = subs(A,{DTDp1 DTDp2},num2cell([7.2 5.6]*10^-6)); % TODO
A = subs(A,{DTDh1 DTDh2},num2cell([4.7 55]*10^-4)); % TODO
A = subs(A,{DhDpMT},num2cell(0.01));
A = double(A);
% ----------------------- POSTPROCESSING ----------------------------------
% Normalizing with maximum deviations
% x = [DVA; p1; h1; d1; T1; p2; h2; d2; T2; DmV; pR; hR; dR; DmG; pC; hC; dC; dMT; DmMT];
T = diag(1./[1; 20*10^5; 50*10^3; 50; 5; 20*10^5; 50*10^3; 50; 5; 0.2;...
    20*10^5; 50*10^3; 50; 0.2; 20*10^5; 50*10^3; 50; 50; 0.2]);
A = T*A/T;
B = T*B;
G = T*G;
expA = sign(A).*exp(abs(A));
expA(A==0) = 0;
% -------------------------- PLOTTING -------------------------------------
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