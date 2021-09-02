function [wa,ra,ta,w,r,t,w0,r0,s] = link_cal2(k,x,L,Ltau,w1,w2,b,v)

%%% Input
%%% k: 
%%% x: grid-time matrix data
%%%

m = size(x,1);
ntim = floor(w1/w2);

wa = zeros(m,m,ntim);
ra = zeros(m,m,ntim);
ta = zeros(m,m,ntim);
s  = zeros(m,ntim);

if strcmp(b,'1') == 1
    w  = zeros(m,m,ntim);
    r  = zeros(m,m,ntim);
    t = zeros(m,m,ntim);
elseif strcmp(b,'2') == 1
    w  = zeros(m,m,ntim);
    r  = zeros(m,m,ntim);
    t = zeros(m,m,ntim);
    
    r0 = zeros(m,m,ntim);
    w0 = zeros(m,m,ntim);
end

c = zeros(m,m,2*Ltau+1);

for j = 1:ntim
    ts = w1*k+w2*j-L;
    te = w1*k+w2*j;
%%%  w1 = 365 w2 = 30;
%%%  w1 = 12  w2 = 1;  
    y = x;
    
    tt  = ts:te;

    c(:,:,Ltau+1) = corr(y(:,tt)',y(:,tt)','type',v);
    for tau = 1:Ltau
        c(:,:,Ltau+1-tau) = corr(y(:,tt-tau)',y(:,tt)','type',v);
        c(:,:,tau+Ltau+1) = corr(y(:,tt)',y(:,tt-tau)','type',v);
    end
    s(:,j) = nanstd(y(:,tt),1,2);
    
    if strcmp(b,'1') == 1
        
        [tmp1,ind1] = max(c,[],3);
        [tmp2,ind2] = min(c,[],3);
        tmp3 = tmp1+tmp2;
        tmp1(tmp3<=0) = 0;
        ind1(tmp3<=0) = 0;
        tmp2(tmp3> 0) = 0;
        ind2(tmp3> 0) = 0;
        
        tmp4 = tmp1+tmp2;
        t(:,:,j) = ind1+ind2-Ltau-1;
        r(:,:,j) = tmp4;
        w(:,:,j) = (tmp4-mean(c,3))./std(c,1,3);
    elseif strcmp(b,'2') == 1        
        [tmp1,ind1] = max(c,[],3);
        [tmp2,ind2] = min(c,[],3);
        tmp3 = tmp1+tmp2;
        tmp1(tmp3<=0) = 0;
        ind1(tmp3<=0) = 0;
        tmp2(tmp3> 0) = 0;
        ind2(tmp3> 0) = 0;
        
        tmp4 = tmp1+tmp2;
        t(:,:,j) = ind1+ind2-Ltau-1;
        r(:,:,j) = tmp4;
        w(:,:,j) = (tmp4-mean(c,3))./std(c,1,3);
        
        r0(:,:,j) = c(:,:,Ltau+1);
        w0(:,:,j) = (r0(:,:,j)-mean(c,3))./std(c,1,3);
    end
    
    c = abs(c);
    [ra(:,:,j),I] = max(c,[],3);
    ta(:,:,j) = I-Ltau-1;
    wa(:,:,j) = (max(c,[],3)-mean(c,3))./std(c,1,3);
end

end