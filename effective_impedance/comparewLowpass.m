J=0.0033
B=3
K=250
tau=0.001*50 % in second
k_e=300
b_e=0.1
m_e=1.5%% .1
P_t=0
I_t=1
I_m=100
P_m=0.1
omega_0=30
%eval(eqnset)
sys1=tf([ K*m_e, K*m_e*omega_0, K*b_e*omega_0, K*k_e*omega_0],[ m_e, m_e*omega_0,0, K*omega_0 0]) % low pass
sys2=tf([ J*K*m_e, B*K*m_e + K*P_m*m_e, K*P_m*b_e + I_m*K*m_e, I_m*K*b_e + K*P_m*k_e, I_m*K*k_e], [ J*m_e, B*m_e + P_m*m_e, I_m*m_e + K*m_e, K*P_m, I_m*K 0])
env=tf([b_e k_e],[1 0]);
bode(sys1,sys2,env,{0.1,100})