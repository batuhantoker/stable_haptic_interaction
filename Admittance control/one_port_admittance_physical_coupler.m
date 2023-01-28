clc
clear
syms  F_int s F_e v_c F_h m_ee E Y_e m b a_1 b_f m_f d_1 K D_m tau_m b_c m_e b_e k_e m_c k_f k_c tau_f tau_e tau D_t alpha_1 alpha_2a alpha_8 d_0 d_1 B P_m a_0 J P_t  alpha_3 I_m I_t beta_4 alpha_4 alpha_5 alpha_6 alpha_7 sigma_1 sigma_2 sigma_3 k_c beta_1 beta_2 beta_3 P_f I_f P_v I_v k w v_e C_m B J P_m I_m tau v_h P_t I_t Z_d F_f F_d v_f v_d C_f P F_c F_v v_m v_v C_v Z_c C F_m C_m1 C_f1 alpha_1 alpha_2a alpha_8 d_0 d_1 B P_m a_0 J P_t  alpha_3 I_m I_t beta_4 alpha_4 alpha_5 alpha_6 alpha_7 sigma_1 sigma_2 sigma_3 k_c beta_1 beta_2 beta_3 P1 F_m F_int F_e v_e v_h Z_d F_f v_f C_f P F_c F_v v_m v_v C_v Z_c


%CASCADED ADMITTANCE CONTROLLER
% % % % Equations
eqn3 = F_d == -F_int;
eqn4 = v_d == F_d*Y_e;
eqn5 = v_c == v_d-v_m;
eqn6 = F_c == C_m*v_c;
eqn7 = v_m == F_m*P;
eqn8 = F_m == F_c-F_int;
eqn9 = F_int == (v_m-v_h)*Z_c;
%eqn10 = F_h ==F_int-v_h/E
vars = [F_d,v_d,v_c,F_c,v_m,F_m,F_int];
eqns = [eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9];
sol = solve(eqns,vars)
F_int=sol.F_int
%F_h=sol.F_h
%v_h=sol.v_h %For two port analysis
ss=symvar(F_int); %Takes symbolic variables of solution
syms(ss(:),'real') %Set variable types as real
%Z_v=-coeffs(F_h,v_h)
Z_v=-coeffs(F_int,v_h) %Output impedance
x=s;
Y_e=(poly2sym(str2sym('[1 0]')))/poly2sym(str2sym('[m_e b_e k_e]')); %Admittance of virtual environment
C_m=poly2sym(str2sym('[P_m I_m]'))/(poly2sym(str2sym('[1 0]'))); %Motion controller
P=1/(poly2sym(str2sym('[J B]'))) %Haptic device
P=eval(P) 
Y_e=eval(Y_e)
C_m=eval(C_m)
%E=(poly2sym(str2sym('[1 0]')))/poly2sym(str2sym('[m_ee 0 0 ]'));
%E=eval(E)

Z_c=poly2sym(str2sym('[K]'))/(poly2sym(str2sym('[1 0]'))); %Physical coupler/Force sensing unit
Z_c=eval(Z_c)
Z2=eval(Z_v)
Z3=(m_e*s^3+(b_e+P_m)*s^2+(I_m+k_e)*s)/(J*m_e*s^4+s^3*(B*m_e+b_e*J+P_m*m_e)+(B*b_e+b_e*P_m+I_m*m_e+J*k_e)*s^2+(I_m*b_e+k_e*(B+P_m))*s+I_m*k_e)
ss=symvar(eval(Z2));
syms(ss(:),'real')
x=s;
x=s;
res=Z2*eval(poly2sym(str2sym('[1 0]')))
res2=1/Z3*eval(poly2sym(str2sym('[1 0]')))%check residue 
resid=limit(res,s,0)
resid2=limit(res2,s,0)
s=j*w; % to frequency domain
syms(w,'real');
[n_11,d_11] = numden(eval(Z2)) % extracting num and den
n_11=expand(n_11*conj(eval(d_11)))
d_11=expand(d_11*conj(eval(d_11)))
bx=real(n_11)
aa1=coeffs(bx,w)'