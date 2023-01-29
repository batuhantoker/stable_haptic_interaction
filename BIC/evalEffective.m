clear s; 
syms s;
x=s
[numa dena]=numden(eval(Z2));
numa=eval(fliplr(coeffs(numa,s)));
dena=eval(fliplr(coeffs(dena,s)));
%sys=tf(numa,[dena]); % maxwell
sys=tf([numa 0],[dena]); % null
%Z=tf([1 20 100],[100 0])
w = logspace(-2,2,400)
[re,im,freq] = nyquist(sys,w);
syms w
syms(w,'real');
s=j*w
%Z2=(s^2+20*s+100)/(100*s)
%Z2=eval(Z2)
im_pos_index=find(( im>0 ))
im_neg_index=find(( im<0 ))
re_pos_index=find(( re>0 ))
ES=[]
EM=[]
ED=[]
a=1
for i=1:length(freq)
    if any(im_neg_index==i)
    ES(i)=-freq(i)*im(:,:,i);
    else
    ES(i)=0;
    end
end
for i=1:length(freq)
    if any(im_pos_index==i)
    EM(i)=freq(i)^-1*im(:,:,i);
    else
    EM(i)=0;
    end
end
for i=1:length(freq)
    if any(re_pos_index==i)
    ED(i)=re(:,:,i);
    else
    ED(i)=0;
    end
end