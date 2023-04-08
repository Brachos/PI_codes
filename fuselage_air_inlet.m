syms a
syms b
for i=1:10
    b=sqrt(0.04192/i);
    a=b*i;
    if 2*a>=0.57912 && 2*b>=0.57912
        fprintf('a=%.2f and b=%.2f',sa,sb);
    end
end 