clc
clear
syms  F_int s F_e b_c m_c K k_f k_c m_e b_e k_e D_t P_f I_f P_v I_v k w v_e C_m B J P_m I_m tau v_h P_t I_t Z_d F_f F_d v_f v_d C_f P F_c F_v v_m v_v C_v Z_c C F_m C_m1 C_f1 P1 F_m F_int F_e v_e v_h Z_d F_f v_f C_f P F_c F_v v_m v_v C_v Z_c


%BASIC IMPEDANCE CONTROLLER
eqn1 = v_e == -v_h;
eqn2 = F_e == Z_d*v_e;
eqn3 = F_f == F_e-F_int;
eqn4 = v_f == F_f*C_f;
eqn7 = v_m == F_m*P;
eqn8 = F_m == v_f-F_int;
eqn9 = F_int == (v_m-v_h)*Z_c;
vars = [v_e,F_e,F_f,v_f,v_m,F_int,F_m];
eqns = [eqn1,eqn2,eqn3,eqn4,eqn7,eqn8,eqn9];
sol = solve(eqns,vars)
F_int=sol.F_int
%v_h=sol.v_h
x=s;
ss=symvar(F_int);
syms(ss(:),'real')
Z_b=-coeffs(F_int,v_h)

C_f=poly2sym(str2sym('[P_t I_t]'))/(poly2sym(str2sym('[1 0]')))%+poly2sym(str2sym('[D_t 0]'))/poly2sym(str2sym('[1/tau 1]'));
C_m=poly2sym(str2sym('[P_m I_m]'))/(poly2sym(str2sym('[1 0]')));
P=1/(poly2sym(str2sym('[J B]')))
P=eval(P)
C_f=eval(C_f)
C_m=eval(C_m)

Z_c=poly2sym(str2sym('[K]'))/(poly2sym(str2sym('[1 0]')));
Z_c=eval(Z_c) % transfer function in s domain

%Inertance model rendering

Z_d=poly2sym(str2sym('[0 0]'))/(poly2sym(str2sym('[1 0]')));%poly2sym(str2sym('[0 b_e*k_e]'))/(poly2sym(str2sym('[b_e k_e]')))
%poly2sym(str2sym('[0 0]'))/(poly2sym(str2sym('[1 0]')));%
Z_d=eval(Z_d)
Z_c=eval(Z_c)
Z2=eval(Z_b)
ss=symvar(eval(Z2));
syms(ss(:),'real')
x=s;
res=Z2*eval(poly2sym(str2sym('[1 0]')))
resid=limit(res,s,0)
[num1 den1]=numden(expand(Z2))
syms s;
num=fliplr(coeffs(num1,s))
den=fliplr(coeffs(den1,s))
routh=routh_hurwitz(den);
s=j*w; % to frequency domain
syms(w,'real');
[n_11,d_11] = numden(eval(Z2)) % extracting num and den
n_11=expand(n_11*conj(eval(d_11)))
d_11=expand(d_11*conj(eval(d_11)))
bx=real(n_11)
aa1=coeffs(real(n_11),w)'
