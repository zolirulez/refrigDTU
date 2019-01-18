clearvars
joint = Joint;
DmR = 0.123;
DmLT = 0.046;
hR = 433e3;
hLT = 500e3;
Dm = [DmR; DmLT];
h = [hR; hLT];
DmMT = -0.321;
hE = 462e3;
[Dm2, h1] = joint.noAccummulation(Dm,h,DmMT,hE);
DmE = Dm2
hMTs = h1
% Second function
Dm = [Dm; DmE];
h = [h; hE];
[Dm2, h1] = joint.noAccummulation(Dm,h);
DmMT = -Dm2
hMTs = h1