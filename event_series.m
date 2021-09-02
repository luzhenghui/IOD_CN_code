function [X,n,ts,te] = event_series(x,th)

a = x(1);
j = 1;
for i = 1:(length(x)-1);
    if (x(i+1)-x(i))==1
        a = [a,x(i+1)];
    else
        b{j} = a;
        j = j+1;
        clear a;
        a = x(i+1);
    end
end
b{j} = a;
        
n  = zeros(length(b),1);
ts = zeros(length(b),1);
te = zeros(length(b),1);

for i = 1:length(b)
    n(i) = length(b{i});
    ts(i) = b{i}(1);
    te(i) = b{i}(end);
end

ind1 = find(n>=th);
ind2 = find(n<th);

n(ind2)  = [];
ts(ind2) = [];
te(ind2) = [];

X = [];
for i = 1:length(ind1);
    X = [X,b{ind1(i)}];    
end




end