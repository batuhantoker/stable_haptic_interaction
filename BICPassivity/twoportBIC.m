clc
clear
syms plant F_e Fint aFint C_f x y E Fh P m Cv s k aFint2 V_h V_e Rd B_11 B_22 B_12 B_21 yd yc H T_s err Pd A w J_m b_m D k_d d_d

assume(Rd, 'real')
assume(F_e, 'real')
assume(Fint, 'real')
assume(aFint, 'real')
assume(D, 'real')
assume(C_f, 'real')
assume(w, 'real')
assume(y, 'real')
assume(Pd, 'real')
assume(A, 'real')
assume(E, 'real')
assume(m, 'real')
assume(P, 'real')
assume(err, 'real')
assume(T_s, 'real')
assume(V_h, 'real')
assume(s, 'real')
assume(k, 'real')
assume(V_e, 'real')
assume(J_m, 'real')
assume(b_m, 'real')
assume(B_11, 'real')
assume(B_12, 'real')
assume(B_21, 'real')
assume(B_22, 'real')
assume(k_d, 'real')
assume(d_d, 'real')
s=j*w;
E=1/(m*s^2);
C_f=P+D*s;
plant=J_m*s^2+b_m*s

Rd=k_d+d_d*s;
%T_s=k;
B_11=-(T_s+err);
Pd=plant*C_f/(1+T_s*plant*(C_f+1));
B_22=0;
B_12=Rd*Pd;
B_21=-A;

H_1=[E/(1-B_11*E) B_12*E/(1-B_11*E);
    B_21*E/(1-B_11*E) B_22+B_12*B_21*E/(1-B_11*E)];


A=1;
m=1;
err=0;
pretty(real(eval(expand(H_1(1,1)))))

b=real(eval(expand(H_1(2,2))))
a=sym2poly(b)

%V_C = [-C_f C_f ; -1 1]

% X = linsolve(inv(V_C),[F_e;V_e]);
% 
% Y = linsolve(inv(H_1),[Fh;xd]);

%Fint=-F_e
% F_e=V_C(2,2)/(1-V_C(2,1));
% yc=V_C(1,1)*-F_e+V_C(1,2);
% yd=P*yc
% H=[H_1(1,1) H_1(1,2)*yd;
%     H_1(2,1) H_1(2,2)*yd]
% eqn1 = y==H_1(1,1)*Fh+H_1(1,2)*xd;
% eqn2 = F_e==H_1(1,1)*Fh+H_1(1,2)*xd;
% eqn3 = yd==V_C(1,1)*F_e+V_C(1,2)*V_e;
% eqn4 = F_e==Fint;
% eqn1 = 0 == Z(2)-X(2);
% eqn2 = F==X(1);
% eqn3 = V_e+Fh==X(1)+Z(1);
% eqn4 = 0 == P*W(1)-Z(2);
% vars = [F_e,V_e,y, Fh]
% eqns = [eqn1,eqn2,eqn3,eqn4]
% sol = solve(eqns,vars)

