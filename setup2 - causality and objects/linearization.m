clearvars
close all
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
B = jacobian(Dx,u);
G = jacobian(Dx,d);
% Substitutions check the values TODO
A = subs(A,{p1 p2 pR pEC},num2cell([86 84.8 38 30]*10^5));
A = subs(A,{h1 h2 hMT hL hG hR},num2cell([500 350 525 250 400 300]*10^3));
A = subs(A,{d1 d2 dR},num2cell([2.8 7 4]*10^2));
A = subs(A,{DmV DmMT DmG DmL},num2cell([0.321 0.321 0.123 0.198]*10^0));
A = subs(A,{CRV CRG},num2cell([0.25 0.25]*10^0));
A = subs(A,{DQ},num2cell([-74.3]*10^3));
A = subs(A,{KvV KvG},num2cell([0.8 2]*10^0));
A = subs(A,{TauV TauG},num2cell([0.1 0.1]*10^0));
A = subs(A,{R},num2cell([1.5]*10^5)); % TODO
A = subs(A,{VGC VR},num2cell([19.2 133]*10^-3));
A = subs(A,{DdDp1 DdDp2 DdDpR},num2cell([3 3 3]*10^-5)); % TODO
A = subs(A,{DdDh1 DdDh2 DdDhR},num2cell([-2 -4 -5]*10^-3)); % TODO
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
[V,eigA] = eig(A);
% Separating real and imaginary values
% Is this right? should not be first the inverse, then the reordering?
% (This was just a copy paste from the exercise on LCD2)
V = [V(:,1) real(V(:,2)) imag(V(:,2)) real(V(:,4)) imag(V(:,4)) V(:,6:end)];
figure(2)
subplot(131)
imagexpeigA = sign(imag(eigA)).*exp(abs(imag(eigA)));
imagexpeigA(imag(eigA)==0) = 0;
imagesc(imagexpeigA,[-1*10^108.5 1*10^108.5])
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
imagesc(inv(V),[-10 10])
% imagesc(V',[-1 1])
colormap('jet')
colorbar
xlabel(char(x))
ylabel('z')
title('V^-^1 reordered')