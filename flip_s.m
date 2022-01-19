function Y = FlipS(num,i,j)
 
 Y = bitxor(num, 2^i+2^j);
 
end
