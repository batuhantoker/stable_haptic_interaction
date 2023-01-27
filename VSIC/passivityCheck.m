passivity=0;
if sign(eval(a_6))==1
    if sign(eval(a_2))==1
        if sign(eval(a_4))==1
            if sign(eval(routh(:,1)))==[1 1 1 1]
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