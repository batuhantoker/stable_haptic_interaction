clc
clear
syms kc bc ic Jm bm nd dd1 dd2 dd3 t x nm
assume(nm, 'real')
assume(Jm, 'real')
assume(bm, 'real')
assume(kc, 'real')
assume(bc, 'real')
assume(ic, 'real')
assume(nd, 'real')
assume(dd1, 'real')
assume(dd2, 'real')
assume(dd3, 'real')
assume(t, 'real')
assume(x, 'real')

P=poly2sym(str2sym('[nm]'))/(poly2sym(str2sym('[Jm bm 0]')));
Rd=poly2sym(str2sym('[nd]'))/(poly2sym(str2sym('[dd1 dd2 dd3]')))
F(1)=(poly2sym(str2sym('[bc kc]')))/(poly2sym(str2sym('[t 1]')))
F(2)=(poly2sym(str2sym('[ic bc kc]')))/(poly2sym(str2sym('[t 1]')))
F(3)=(poly2sym(str2sym('[ic ic/bc 1/ic]')))/(poly2sym(str2sym('[t 1]')))
for i=1:length(F)
H(1,1,i)=P/(1+P*F(i));
H(1,2,i)=P*F(i)/(1+P*F(i));
H(2,1,i)=P*F(i)/(1+P*F(i));
H(2,2,i)=-F(i)/(1+P*F(i));
end
for  i=1:length(F)
    ThetaC(i)=(P-Rd-Rd^2*F(i))/(Rd*(1+F(i)*(P+Rd)));
    R(i)=H(1,1,i)+(H(1,2,i)*H(2,1,i)*Rd)/(1-H(2,2,i)*Rd);
    ThetaO(i)=P/R(i)-1;
end
Jm=1;
bm=0.775;
nm=1550;
nd=474;
dd1=1;
dd2=5.03;
dd3=158;
bc=0.203;
kc=17.1;
ic=0.3;
t=100;
for i=1:length(F)
    OL(i)=sym2tf(eval(ThetaO(i)))
    CL(i)=sym2tf(eval(ThetaC(i)))
end
figure
title('Open-loop')
hold
for i=1:length(F)
    bode(OL(i))
end
figure
title('Closed-loop')
hold
for i=1:length(F)
    bode(CL(i))
end   
figure
title('Closed-loop')
hold
for i=1:length(F)
    step(CL(i))
end  
figure
title('Open-loop')
hold
for i=1:length(F)
    step(OL(i))
end