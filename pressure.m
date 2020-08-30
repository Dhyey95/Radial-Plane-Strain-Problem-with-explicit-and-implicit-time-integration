function P= pressure(tn,t)
if (tn-(t/2))<= (10^(-10))
    P = sin(2*pi*tn/t);
else
    P = 0;
end
end
