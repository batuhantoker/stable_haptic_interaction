passivity=0;
if sign(eval(d_6s))==1
    if sign(eval(d_2s))==1
        if sign(eval(d_4s))==1
            if sign(eval(xi))==1
                passivity=1;
            else
                passivity=0;
            end
        else
            if sign(beta)==-1
                passivity=1;
            end
        end
    end
end
passivity;