function [f,s,s0k_r,s0k_w]=psd_ynm(DATA,x)
% Power Spectral Density(PSD)
% provided by Yuan Naiming 
% INPUT PARAMETERS-----------------------------------------------------
% DATA:       input DATA
% ts:            sampling time interval 
% OUTPUT VARIABLES-----------------------------------------------------
% f:              frequency
% s:              power spectrum density
%-----------------------------------------------------
N=length(DATA);
m=ceil(N/3);
DATA_mean=mean(DATA);
DATA_std=std(DATA,1);
for j=0:m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for t=1:N-j
%         a_1(t)=(DATA(t)-DATA_mean)/DATA_std;
%         a_2(t)=(DATA(t+j)-DATA_mean)/DATA_std;
%     end
    t = 1:N-j;
    a_1 = (DATA(t)-DATA_mean)/DATA_std;
    a_2 = (DATA(t+j)-DATA_mean)/DATA_std;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a=a_1.*a_2;
    r(j+1)=sum(a)/(N-j);
    clear a t;
end  
s_1(1)=(r(1)+r(m+1))/(2*m)+(sum(r)-r(1)-r(m+1))/m;
% s_1(1)=(r(1)+r(m+1)+2*(sum(r)-r(1)-r(m+1)))/m;
for k=1:m-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for j=1:m-1
%         s_sa(j)=r(j+1)*cos((k*pi*j)/m);
%     end
    j = 1:m-1;
    s_sa = r(j+1).*cos((k*pi*j)/m);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s_s(k)=sum(s_sa);
    clear s_sa  j;
    s_1(k+1)=(r(1)+r(m+1)*cos(k*pi)+2*s_s(k))/m;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:m-1
    s_sb(j)=((-1)^j)*r(j+1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_sbb=sum(s_sb)/m;
s_1(m+1)=(r(1)+(-1)^m*r(m+1))/(2*m)-s_sbb;
% s_sbb=sum(s_sb);
% s_1(m+1)=(r(1)+(-1)^m*r(m+1)+2*s_sbb)/m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s(1)=0.5*s_1(1)+0.5*s_1(2);
for k=1:m-1
    s(k+1)=0.25*s_1(k)+0.5*s_1(k+1)+0.25*s_1(k+2);
end
s(m+1)=0.5*s_1(m)+0.5*s_1(m+1);

for k=1:m
    T(k)=2*m/k;
    f(k) = 1/T(k);
end

s_mean = (s(1)+s(m+1))/(2*m)+(sum(s)-s(1)-s(m+1))/m;

k = 0:m;
s0k = s_mean*bsxfun(@rdivide, 1-r(2)^2, 1+r(2)^2+2*r(2)*cos(pi*k/m) );
s0k_r = s0k*bsxfun(@rdivide,x,(2*N-m/2)/m );

s0k_w = s_mean*bsxfun(@rdivide,x,(2*N-m/2)/m );

end 
%plot(T,s(2:m+1));
