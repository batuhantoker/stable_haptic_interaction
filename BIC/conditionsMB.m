a=den
xi=a(2)*a(3)*a(4)-a(1)*a(4)^2-a(5)*a(2)^2

delta_k=K-k_e
g_2=I_t^2*b_e*k_e^2*K^2
g_4=K^2*k_e*(((B+P_t*b_e)*(P_t+1)-I_t*J)*k_e-I_t*b_e^2)-K*(B*I_t*b_e*k_e^2)
g_6=K*b_e*(B*K*b_e-J*P_t*k_e^2+b_e*(B*P_t-I_t*J)*delta_k)
%g_6=K^2*b_e^2*(B*(P_t+1)-I_t*J)+K*b_e*k_e*(b_e*(J*I_t-B*P_t)-J*P_t*k_e)

beta=g_4^2-4*g_2*g_6
beta_cond=collect(beta,J)

%% Damper
J_max=(B*(P_t+1))/I_t
syms eps
beta_new=subs(beta,J,J_max-eps)
solve(beta_new,eps)
beta_new_simp=-4*subs(g_2,k_e,0)*K*(I_m*I_t*b_e - I_t*P_m*(K ))*eps+subs(g_4,k_e,0)^2%-4*g_2*(K - I_t*P_m*b_e)*(K - 1)*(B + P_m)
I_tmax=(K*P_m-B*I_m)*(B+P_m)/(J*K*P_m^2)
J_limit=[J_min+subs(g_4,k_e,0)^2/(4*subs(g_2,k_e,0)*K*((I_m*I_t*b_e - I_t*P_m*(K - k_e)))) J_min]