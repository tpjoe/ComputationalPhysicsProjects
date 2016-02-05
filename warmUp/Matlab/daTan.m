function [ d ] = daTan( x, h, method, precision )
%DATAN calculate 1st derivative of atan(x) under steplength h using two method 
%   Detailed explanation goes here
if(strcmp(precision,'single'))
    x = single(x);
    h = single(h);
end


switch(method)
    case 1
        d = (atan(x+h)-atan(x))./h;
        if(strcmp(precision,'single'))
            d = single(d);
        end
    case 2
        d = (atan(x+h)-atan(x-h))./2./h;
        if(strcmp(precision,'single'))
            d = single(d);
        end
    otherwise
        disp('Wrong methods');
end



end

