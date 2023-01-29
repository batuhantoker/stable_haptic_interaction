clc
clear
syms  F_int t_d s tau_d x F_e Y_e omega_0 b_f a_1 b_f m_f d_1 K D_m tau_m b_c m_e b_e k_e m_c k_f k_c tau_f tau_e tau D_t alpha_1 alpha_2a alpha_8 d_0 d_1 B P_m a_0 J P_t  alpha_3 I_m I_t beta_4 alpha_4 alpha_5 alpha_6 alpha_7 sigma_1 sigma_2 sigma_3 k_c beta_1 beta_2 beta_3 P_f I_f P_v I_v k w v_e C_m B J P_m I_m tau v_h P_t I_t Z_d F_f F_d v_f v_d C_f P F_c F_v v_m v_v C_v Z_c C F_m C_m1 C_f1 alpha_1 alpha_2a alpha_8 d_0 d_1 B P_m a_0 J P_t  alpha_3 I_m I_t beta_4 alpha_4 alpha_5 alpha_6 alpha_7 sigma_1 sigma_2 sigma_3 k_c beta_1 beta_2 beta_3 P1 F_m F_int F_e v_e v_h Z_d F_f v_f C_f P F_c F_v v_m v_v C_v Z_c


%VELOCITY SOURCED IMPEDANCE CONTROLLER
eqn1 = v_e == -F_int*Y_e;
eqn2 = v_m == v_e*P
eqn3 = F_int == (v_m-v_h)*Z_c;
vars = [v_e,v_m,F_int];
eqns = [eqn1,eqn2,eqn3];
sol = solve(eqns,vars)
F_int=sol.F_int
%v_h=sol.v_h
ss=symvar(F_int);
syms(ss(:),'real')
Z_v=-coeffs(F_int,v_h)
x=s;
t_d=1%exp(-x*tau_d)
C_f=poly2sym(str2sym('[P_t I_t]'))/(poly2sym(str2sym('[1 0]')));
C_m=poly2sym(str2sym('[P_m I_m]'))/(poly2sym(str2sym('[1 0]')));
P=poly2sym(str2sym('[omega_0]'))/(poly2sym(str2sym('[1 omega_0]')))
P=eval(P)
C_f=eval(C_f)
Y_e=poly2sym(str2sym('[1 0]'))/poly2sym(str2sym('[m_e b_e 0]')) % *poly2sym(str2sym('[1-k_e/K]'))
Y_e=eval(Y_e)
Z_c=poly2sym(str2sym('[0 b_f K]'))/(poly2sym(str2sym('[1 0]')));%poly2sym(str2sym('[b_f K]'))/(poly2sym(str2sym('[1 0]')))+poly2sym(str2sym('[m_f 0]'))/poly2sym(str2sym('[1/tau_f 1]'));; %m_f b_f
Z_c=eval(Z_c)
Z_v=eval(Z_v) % transfer function in s domain
Z2=eval(Z_v)
ss=symvar(eval(Z2));
syms(ss(:),'real')
x=s;
x=s;
res=Z2*eval(poly2sym(str2sym('[1 0]')))
resid=limit(res,s,0)
s=j*w; % to frequency domain
syms(w,'real');
[n_11,d_11] = numden(eval(Z2)) % extracting num and den
n_11=expand(n_11*conj(eval(d_11)))
d_11=expand(d_11*conj(eval(d_11)))
bx=real(n_11)
aa1=coeffs(bx,w)'
[num1 den1]=numden(expand(Z2))
syms s;
den=fliplr(coeffs(den1+num1,s))
%routh=routh_hurwitz(den);
[numa dena]=numden(eval(Z2));
numa=eval(fliplr(coeffs(numa,s)));
dena=eval(fliplr(coeffs(dena,s)));
% sys=tf(numa,[dena 0]);