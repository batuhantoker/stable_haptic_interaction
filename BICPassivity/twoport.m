clc; clear;
syms B_11 B_12 B_21 B_22 C_11 C_12 C_21 C_22 th u F_int F_e V_e y E J_m b_m P I b_f k_f w m 

eqn1 = u== B_11*y+B_12*th;
eqn2 =F_int==B_21*y+B_22*th;;
eqn3 =th==C_11*F_int+C_12*F_e
eqn4 =V_e==C_21*F_int+C_22*F_e
eqn5 = solve(eqn1,th)
eqn6 = subs(eqn2,th,eqn5)
eqn7 = solve(eqn6,F_int)
eqn8 = subs(eqn3,th,eqn5)
eqn9 = subs(eqn8,F_int,eqn7)
eqn10 = subs(eqn4,F_int,eqn7)
eqn11= solve(eqn9,u)
eqn12= subs(eqn10,u,eqn11)
eqnA= solve(eqn9,u)
eqnB= solve(eqn12,V_e)
A = equationsToMatrix([eqnA,eqnB],[y,F_e])
H=[E/(1-A(1,1)*E) A(1,2)*E/(1-A(1,1)*E);
    A(2,1)*E/(1-A(1,1)*E) A(2,2)+A(1,2)*A(2,1)*E/(1-A(1,1)*E)];

%%
%         assume(J_m, 'real')
%         assume(b_m, 'real')
%         assume(b_f, 'real')
%         assume(k_f, 'real')
%         assume(P, 'real')
%         assume(I, 'real')
%         assume(w, 'real')
%         assume(w,'positive')
%         assume(E, 'real')
%         plant=1/(poly2sym(str2sym('[J_m b_m 0]')));
%         comp=poly2sym(str2sym('[b_f k_f]'))/poly2sym([1 0]);
%         cf=poly2sym(str2sym('[P I]'))/poly2sym([1 0]);
%         B_11=comp;
%         B_22=comp;
%         B_12=-comp;
%         B_21=-comp;
%         B=[B_11 B_12; B_21 B_22];
%         C_11=-plant*cf/(1+plant*cf);
%         C_12=plant*cf/(1+plant*cf);
%         C_21=-cf;
%         C_22=cf;
%         C=[C_11 C_12; C_21 C_22];
%         x=i*w
%         real(eval(B(1,1)))*real(eval(B(2,2)))-(abs(eval((B(2,1))+eval(B(1,2)))/2))^2
%         real(eval(C(1,1)))*real(eval(C(2,2)))-(abs(eval((C(2,1))+eval(C(1,2)))/2))^2
%%

% sol = solve([eqn3,eqn4,eqn5,eqn6],[F_e,V_e,u,y])
% sol = solve([eqn1,eqn2,eqn3,eqn4,eqn5,eqn6],[F_e,V_e,u,y,F_int,th])