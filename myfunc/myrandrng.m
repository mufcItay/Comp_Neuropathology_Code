function [x] = myrandrng(N,a,b)
x = (b-a).*rand(N,1) + a;
end

