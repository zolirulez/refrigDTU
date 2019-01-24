Qfreezer = 11000;
Qcooler = 36000;
% Hot conditions
Rwall = (30-25)/(Qfreezer+Qcooler);
Rcooler = (25-5)/Qcooler;
Rfreezer = (25+30)/Qfreezer;
% Cold conditions
Qin = (15+30)/Rfreezer + (15-5)/Rcooler;
HRR = 1.1;
QHR = HRR*Qin;
RHR = (45-15)/QHR;