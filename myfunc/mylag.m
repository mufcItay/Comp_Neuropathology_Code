function [xx] = mylag(x,lag)
for i=1:length(x)
    if i<=lag
    xx(i)=NaN;
    else
    xx(i)=x(i-lag);
    end
end
end


