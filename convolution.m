close, clear, clc;

x = [1 2 3 4 5 6 7 8 9 10 ];
y = [11 12 13 14];
z2 = conv(x, y);

xlen = length(x);
ylen = length(y);

zlen = xlen + ylen -1;
z = zeros(1, zlen);

for b=1:zlen
    for a =1: b-1
        
        z(b) = z(b)+  x(a)*y(b-a);
        
      
    end  
   disp(z(b));
end
