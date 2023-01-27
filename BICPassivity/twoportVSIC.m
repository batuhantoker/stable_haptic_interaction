clc
clear
syms xd Fe Fint aFint Cf x y E Fh P Fm Cv s k aFint2 Vh Ve Rd
assume(xd, 'real')
assume(Rd, 'real')
assume(Fe, 'real')
assume(Fint, 'real')
assume(aFint, 'real')
assume(aFint2, 'real')
assume(Cf, 'real')
assume(x, 'real')
assume(y, 'real')
assume(E, 'real')
assume(Fh, 'real')
assume(P, 'real')
assume(Fm, 'real')
assume(Cv, 'real')
assume(Vh, 'real')
assume(s, 'real')
assume(k, 'real')
assume(Ve, 'real')
% Fm=(xd*Cv+Fint*P*Cv)/(1+P*Cv);
% xd=(Fe-Fint)*Cf;
% %y=E*(Fint+Fh);
% y=-Ve/s;
% x=(Fm-Fint)*P;
% aFint=(x-y)*k;
% aF=coeffs(eval(aFint),Fint);
% aFint2=aF(1)/(-aF(2)+1)
% aF2=coeffs(eval(aFint2),Fe);
% ave2=coeffs(eval(aF2(1)),Ve);
% H21=1/-aF2(2)
% H22=(-1/E+ave2)/-aF2(2)
% alpha=P*k*((Cf*Cv - Cv*P)/(Cv*P + 1) + 1);
% pretty(H21)
% pretty(H22)
% G11=-(1+H21/H22*Rd)
% G12=1/(Rd*H22);
% G21=-H21/H22;
% G22=1/H22;
% G=[G11 G12;
%     G21 G22]

eqn1 = 0==-Fm+(xd-x)*s*Cv;
eqn2 = 0==-xd+(Fe-Fint)*Cf;
eqn5 = 0==-Ve/s+E*(Fint+Fh);
eqn3 = 0==-Fint+(x+Ve/s)*k;
eqn6= 0==-Vh+Ve;
eqn4 = 0==-x+(Fm-Fint)*P;
vars = [Fe,Ve,Vh,Fh,Fint,Fm]
eqns = [eqn1,eqn2,eqn3,eqn4,eqn5,eqn6]
sol = solve(eqns,vars)


