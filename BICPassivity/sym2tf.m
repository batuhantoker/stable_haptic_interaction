function G = sym2tf(g)
[n,m]=size(g);
for i=1:n
    for j=1:m
        [num,den]=numden(g(i,j));
        num_n=sym2poly(num);
        den_n=sym2poly(den);
        G(i,j)=tf(num_n,den_n);
    end
end