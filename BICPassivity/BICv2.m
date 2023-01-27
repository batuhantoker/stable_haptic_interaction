clc
clear
syms P D kd dd Jm bm kf w T bf t x s
assume(Jm, 'real')
assume(bm, 'real')
assume(kd, 'real')
assume(dd, 'real')
assume(P, 'real')
assume(D, 'real')
assume(bf, 'real')
assume(kf, 'real')
assume(w, 'real')
assume(t, 'real')
assume(s, 'real')
Jm=str2sym('Jm');
imp=poly2sym(str2sym('[dd kd]'));
cf=poly2sym(str2sym('[D P]'));
Tact=poly2sym(str2sym('[bf kf]'));
robot=(poly2sym(str2sym('[Jm bm 0]')));

Tsens=poly2sym(str2sym('[kf]'));

x=j*w
envnum=Tact*(imp*cf+robot)
envden=(robot+Tact*(cf+1))*poly2sym(str2sym('[1 0]'))
envnum1=expand(eval(envnum));
envden1=expand(eval(envden));
envnum2=expand(envnum1*conj(envden1));
envden2=expand(conj(envden1)*envden1);

pretty(real(envnum2))
a=coeffs(real(eval(envnum2)),w);
d=coeffs(real(eval(envden2)),w); 
for i=1:length(a)
    b(i,1:length(coeffs(a(i),t)))=coeffs(a(i),t);
    b1(i,1:length(coeffs(d(i),t)))=coeffs(d(i),t);
end
