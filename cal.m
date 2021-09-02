clc;
clear;

cal01_noaa = 0;
cal02_K_v2 = 0;
cal03_shuffle = 0;
cal04 = 0;
cal05 = 0;
cal06 = 0;
cal07 = 0;
cal08 = 0;

if cal01_noaa == 1
    path = '../../myDATA/sst.day.anom.v2.nc/';
    
    lon = double(ncread([path,'sst.day.mean.1983.nc'],'lon',1,inf,8));
    lat = flip(double(ncread([path,'sst.day.mean.1983.nc'],'lat',1,inf,8)));

    Nlon = length(lon);
    Nlat = length(lat);
    
    sst  = zeros(Nlon,Nlat,365,39);
    sst(:,:,:,:) = nan;
    
    for i = 1:39
        yr = i+1981;
        
        if yr <=2015
            tmp = double(ncread(['../../myDATA/sst.day.anom.v2.nc/sst.day.mean.',num2str(yr),'.nc'],'sst',  [1,1,1],[inf,inf,inf],[8,8,1]));
        else
            tmp = double(ncread(['../../myDATA/sst.day.anom.v2.1.nc/sst.day.mean.',num2str(yr),'.nc'],'sst',[1,1,1],[inf,inf,inf],[8,8,1]));
        end
            
        tmp(abs(tmp)>100) = nan;
        
        if mod(yr,4) == 0
            sst(:,:,1:size(tmp,3)-1,i) = tmp(:,:,[1:59,61:end]);
        else
            sst(:,:,1:size(tmp,3),i) = tmp(:,:,1:end);
        end
        disp(i);
    end

    sst = flip(sst,2);
    ssta = bsxfun(@minus,sst,nanmean(sst(:,:,:,1:30),4));
    
    x = reshape(ssta((lon(:)>=40 & lon(:)<=110),abs(lat(:))<=11,:,:),[],365*39);
    x = bsxfun(@minus,x,nanmean(x,1));
    
    for k = 2:38
        disp(k);
        [wa,~,~,w,r,t] = link_cal2(k,x,365,200,365,30,'1','Pearson');
        
        save(['cal01_IOD_NOAA_',num2str(k+1982),'_v2.mat'],'w','r','t','-v7.3');      
    end

end

if cal02_K_v2 == 1
    m = 385;
    
    path = '../../myDATA/sst.day.anom.v2.nc/';
    
    lon = double(ncread([path,'sst.day.mean.1983.nc'],'lon',1,inf,8));
    lat = flip(double(ncread([path,'sst.day.mean.1983.nc'],'lat',1,inf,8)));

    
    [lon1,lat1] = meshgrid(lon(lon>=40 & lon<=110),lat(abs(lat)<=11));
    lon1 = lon1';
    lat1 = lat1';
    
    ind1 = find(lon1(:)>=90 & lon1(:)<=110 & lat1(:)<=0);
    ind2 = find(lon1(:)>=50 & lon1(:)<=70);
        
    X1 = zeros(m,m,12,37);
    T1 = zeros(m,m,12,37);
    
    X2 = zeros(m,m,12,37);
    T2 = zeros(m,m,12,37);    
    for i = 1:37
        load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2.mat'],'w','t','r');
        X1(:,:,:,i) = w;
        T1(:,:,:,i) = t;
        
        load(['cal04_shuffle_',num2str(i+1983),'_b.mat'],'w','t');
        
        X2(:,:,:,i) = w;
        T2(:,:,:,i) = t;
    end
    
    XX = X1;
    XX(abs(T1)>150) = nan;
    K = reshape(nansum(nansum(XX(ind1,ind2,:,:),1),2),[],1);
    
%     X2(isnan(XX)==1) = nan;
%     K2 = reshape(nansum(nansum(X2(ind1,ind2,:,:),1),2),[],1);

    save('cal02_K_v2_T150.mat','K');
end

if cal03_shuffle == 1
    path = '../../myDATA/sst.day.anom.v2.nc/';
    
    lon = double(ncread([path,'sst.day.mean.1983.nc'],'lon',1,inf,8));
    lat = flip(double(ncread([path,'sst.day.mean.1983.nc'],'lat',1,inf,8)));

    Nlon = length(lon);
    Nlat = length(lat);
    
    sst  = zeros(Nlon,Nlat,365,39);
    sst(:,:,:,:) = nan;
    
    for i = 1:39
        yr = i+1981;
        
        if yr <=2015
            tmp = double(ncread(['../../myDATA/sst.day.anom.v2.nc/sst.day.mean.',num2str(yr),'.nc'],'sst',  [1,1,1],[inf,inf,inf],[8,8,1]));
        else
            tmp = double(ncread(['../../myDATA/sst.day.anom.v2.1.nc/sst.day.mean.',num2str(yr),'.nc'],'sst',[1,1,1],[inf,inf,inf],[8,8,1]));
        end
            
        tmp(abs(tmp)>100) = nan;
        
        if mod(yr,4) == 0
            sst(:,:,:,i) = tmp(:,:,[1:59,61:end]);
        else
            sst(:,:,:,i) = tmp(:,:,1:end);
        end
        disp(i);
    end
    
    sst = flip(sst,2);
    ssta = bsxfun(@minus,sst,nanmean(sst(:,:,:,1:30),4));
%%%%% a
    ssta((lon(:)>=40 & lon(:)<=110),abs(lat(:))<=11,:,:) = ...
    bsxfun(@minus,ssta((lon(:)>=40 & lon(:)<=110),abs(lat(:))<=11,:,:),...
    nanmean(nanmean(ssta((lon(:)>=40 & lon(:)<=110),abs(lat(:))<=11,:,:),1),2));
%%%%    
    [lon1,lat1] = meshgrid(lon,lat);
    lon1 = lon1';
    lat1 = lat1';
    
    ind1 = find(lon1(:)>=90 & lon1(:)<=110 & lat1(:)<=0 & abs(lat1(:))<=11);
    
    N1 = length(ind1);
    
    x = reshape(ssta,[],365*39);
    x1 = x(ind1,:);
    
    clear ssta sst tmp;
    
%%
    path = '../../myDATA/sst.day.anom.v2.nc/';
    
    lon = double(ncread([path,'sst.day.mean.1983.nc'],'lon',1,inf,8));
    lat = double(ncread([path,'sst.day.mean.1983.nc'],'lat',1,inf,8));
    
    lat_ind = find(abs(lat)<=30);
    lat = lat(lat_ind);
    
    Nlon = length(lon);
    Nlat = length(lat_ind);
    
    sst  = zeros(Nlon,Nlat,365,39);
    sst(:,:,:,:) = nan;
    
    for i = 1:39
        yr = i+1981;
        
        if yr <=2015
            tmp = double(ncread(['../../myDATA/sst.day.anom.v2.nc/sst.day.mean.',num2str(yr),'.nc'],...
                'sst',[1,1,1],[inf,inf,inf],[8,8,1]));
        else
            tmp = double(ncread(['../../myDATA/sst.day.anom.v2.1.nc/sst.day.mean.',num2str(yr),'.nc'],...
                'sst',[1,1,1],[inf,inf,inf],[8,8,1]));
        end
            
        tmp(abs(tmp)>100) = nan;
        
        if mod(yr,4) == 0
            sst(:,:,:,i) = tmp(:,lat_ind,[1:59,61:end]);
        else
            sst(:,:,:,i) = tmp(:,lat_ind,1:end);
        end
        disp(i);
    end
    
    ssta = bsxfun(@minus,sst,nanmean(sst(:,:,:,1:30),4));
    ssta = bsxfun(@minus,ssta,nanmean(nanmean(ssta((lon(:)>=40 & lon(:)<=110),abs(lat(:))<=11,:,:),1),2));
    ssta(lon(:)>=90 & lon(:)<=110, lat(:)<=0 & abs(lat(:))<=11,:,:) = nan;
    
    x = reshape(ssta,[],365*39);
    x(isnan(x(:,1))==1,:) = [];
    
    clear sst ssta;
        
    for k = 2:38
        [w,r,t,~,~] = link_cal3(k,x1,x,365,200,365,30,'2','Pearson');
        save(['cal03_shuffle_',num2str(k+1982),'_2R.mat'],'w','t','-v7.3');      
    end

end

if cal04 == 1
    m = 385;
    
    path = '../../myDATA/sst.day.anom.v2.nc/';
    
    lon = double(ncread([path,'sst.day.mean.1983.nc'],'lon',1,inf,8));
    lat = flip(double(ncread([path,'sst.day.mean.1983.nc'],'lat',1,inf,8)));

    [lon1,lat1] = meshgrid(lon(lon>=40 & lon<=110),lat(abs(lat)<=11));
    lon1 = lon1';
    lat1 = lat1';
    
    ind1 = find(lon1(:)>=90 & lon1(:)<=110 & lat1(:)<=0);
    ind2 = find(lon1(:)>=50 & lon1(:)<=70);
    
    for i = 1:37
        load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2.mat'],'w','t');
        
        X(:,:,:,i) = w;
        T(:,:,:,i) = t;
    end
    
    XX = X;
    XX(abs(T)>150) = nan;
    XX = reshape(XX,m,m,[]);
    XX = XX(ind1,ind2,:);
    
    nd = zeros(size(XX,1),size(XX,3));
    
    for i = 1:length(ind1)
        for j = 1:size(XX,3)
            tmp = reshape(XX(i,:,j),[],1);
            nd(i,j) = length(find(isnan(tmp)==0));
        end
    end
    nd = reshape(nd,[],12,37);
    
    d = zeros(1000,length(ind1),12,37);
    for i = 1:37
       load(['cal03_shuffle_',num2str(i+1983),'_2R.mat'],'w');
       for j = 1:12
           for jj = 1:length(ind1)
               a = reshape(w(jj,:,j),[],1);
               a(isnan(a)==1) = [];
               for p = 1:1000
                   b = randperm(length(a),squeeze(nd(jj,j,i)));
                   d(p,jj,j,i) = nansum(a(b));
               end
           end
       end
       disp(i);
    end
    
    save('com13_v2_T150_2R.mat','d','nd');
    
end

if cal05 == 1
   load('cal02_K_v2_T150.mat');
   K1 = smooth(K,3);
   K1 = reshape(K1,12,[]);
   
   load('com13_v2_T150_2R.mat','d');
   dmin = reshape(nansum(squeeze(min(d,[],1)),1),[],1);
   th_K = nanmean(reshape(dmin,12,[]),1);

   [min_K,min_loc] = min(K1,[],1);
   min_K(16) = K1(10,16);
   min_K(24) = K1(12,24);
   min_K(25) = nan;
   min_K(26) = K1(12,26);
   min_K(28) = K1(11,28);
   min_K(29) = K1(5,29);
   min_K(31) = K1(5,31);
   min_K(37) = K1(9,37);
   
   min_loc(16) = 10;
   min_loc(24) = 12;
   min_loc(25) = nan;
   min_loc(26) = 12;
   min_loc(28) = 11;
   min_loc(29) = 5;
   min_loc(31) = 5;
   min_loc(37) = 9;

   save('min_K_v2_T150.mat','min_K','min_loc','th_K');

   
end

if cal06 == 1
    lon = double(ncread('../../myDATA/Heat_content/Monthly/CZ16_1_2000m_Temp_year_1980_month_01.nc','lon'));
    lat = double(ncread('../../myDATA/Heat_content/Monthly/CZ16_1_2000m_Temp_year_1980_month_01.nc','lat'));
    dep = double(ncread('../../myDATA/Heat_content/Monthly/CZ16_1_2000m_Temp_year_1980_month_01.nc','depth_std'));
    
    Nlon = length(lon);
    Nlat = length(lat);
    
    X = zeros(Nlon,Nlat,12,37);
    X(:,:,:,:) = nan;
    
    for i = 1:37
        for j = 1:12
            tmp = ncread(['../../myDATA/Heat_content/Monthly/CZ16_1_2000m_Temp_year_',...
                num2str(i+1983),'_month_',num2str(j,'%02d'),'.nc'],'temp');
            tmp(tmp>=900) = nan; 
%% H300 
%             for k = 1:40
%                 tmp2(k,:,:) = ((tmp(k+1,:,:)+tmp(k,:,:))/2)*(dep(k+1)-dep(k));
%             end
%             X(:,:,j,i) = squeeze(nansum(tmp2(1:18,:,:),1))/300;
%% D20
            tmp1 = interpn(dep,lon,lat,tmp,1:500,lon,lat);
            tmp1(tmp1<20) = nan;
            [~,X(:,:,j,i)] = min(tmp1,[],1);
        end
        disp(i);
    end
    
    X(X==1) = nan;
    Z = X;
    
    X = bsxfun(@minus,X,nanmean(X(:,:,:,1:30),4));
    X = reshape(X,Nlon,Nlat,12,37);
%% H300
%     Y = X;
%     Y(:,:,:,:) = nan;
%     for i = 1:Nlat
%         Y(:,i,:,:) = X(:,i,:,:)*cosd(lat(i));
%     end
%     GT = squeeze(nanmean(nanmean(Y,1),2));
%% D20
    GT = squeeze(nanmean(nanmean(X,1),2));


%     save('cal06_ohc_CZ_H300.mat','X','GT','lon','lat','-v7.3');
    save('cal06_ohc_CZ_D20.mat','X','GT','lon','lat','-v7.3');
end

if cal07 == 1
    setup_nctoolbox;

    a = ncgeodataset('../../myDATA/ERA5/surface/u10m.ERA5.198401.grib');
    lon = a.data('lon');
    lat = a.data('lat');
    
    lon = lon(1:4:end);
    lat = lat(1:4:end);
    
    Nlon = length(lon);
    Nlat = length(lat);
    
    U = zeros(Nlon,Nlat,12,37);
    U(:,:,:,:) = nan;
    
    V = zeros(Nlon,Nlat,12,37);
    V(:,:,:,:) = nan;
    
    for i = 1:37
        for j = 1:12
            b = ncgeodataset(['../../myDATA/ERA5/surface/u10m.ERA5.',num2str(i+1983),num2str(j,'%02d'),'.grib']);
            c = b.variables;
            tmp = squeeze(double(b.data(c{2,1})))';
            U(:,:,j,i) = tmp(1:4:end,1:4:end);
            
            b = ncgeodataset(['../../myDATA/ERA5/surface/v10m.ERA5.',num2str(i+1983),num2str(j,'%02d'),'.grib']); 
            c = b.variables;
            tmp = squeeze(double(b.data(c{2,1})))';            
            V(:,:,j,i) = tmp(1:4:end,1:4:end);

        end
        disp(i);
    end
    
    U = bsxfun(@minus,U,nanmean(U(:,:,:,1:30),4));
    U = reshape(U,Nlon,Nlat,12,37);

    Y = U;
    Y(:,:,:,:) = nan;
    for i = 1:length(Nlat)
        Y(:,i,:,:) = U(:,i,:,:)*cosd(lat(i));
    end
    UGT = squeeze(nanmean(nanmean(Y,1),2));
    
    V = bsxfun(@minus,V,nanmean(V(:,:,:,1:30),4));
    V = reshape(V,Nlon,Nlat,12,37);

    Y = V;
    Y(:,:,:,:) = nan;
    for i = 1:Nlat
        Y(:,i,:,:) = V(:,i,:,:)*cosd(lat(i));
    end
    VGT = squeeze(nanmean(nanmean(Y,1),2));
    
    save('cal07_10Wind.mat','U','V','UGT','VGT','lon','lat','-v7.3');

end

if cal08 == 1

    path = '../../myDATA/sst.day.anom.v2.nc/';
    
    lon = double(ncread([path,'sst.day.mean.1983.nc'],'lon',1,inf,8));
    lat = flip(double(ncread([path,'sst.day.mean.1983.nc'],'lat',1,inf,8)));

    Nlon = length(lon);
    Nlat = length(lat);
    
    sst  = zeros(Nlon,Nlat,365,38);
    sst(:,:,:,:) = nan;
    
    for i = 1:38
        yr = i+1981;
        
        if yr <=2015
            tmp = double(ncread(['../../myDATA/sst.day.anom.v2.nc/sst.day.mean.',num2str(yr),'.nc'],'sst',  [1,1,1],[inf,inf,inf],[8,8,1]));
        else
            tmp = double(ncread(['../../myDATA/sst.day.anom.v2.1.nc/sst.day.mean.',num2str(yr),'.nc'],'sst',[1,1,1],[inf,inf,inf],[8,8,1]));
        end
            
        tmp(abs(tmp)>100) = nan;
        
        if mod(yr,4) == 0
            sst(:,:,:,i) = tmp(:,:,[1:59,61:end]);
        else
            sst(:,:,:,i) = tmp(:,:,1:end);
        end
        disp(i);
    end
    
    sst = flip(sst,2);
    sst = reshape(sst,Nlon,Nlat,[]);
    y = zeros(Nlon,Nlat,12*38);
    y(:,:,:) = nan;
    for i = 1:Nlon
        for j = 1:Nlat
            y(i,j,:) = day2mon(sst(i,j,:));
        end
    end    
    y = reshape(y,Nlon,Nlat,12,[]);
    ssta = bsxfun(@minus,y,nanmean(y(:,:,:,1:30),4));    
    x = reshape(ssta((lon(:)>=40 & lon(:)<=110),abs(lat(:))<=11,:,:),[],12*38);
    
    [a,b,c,d] = eof_dan(x',10);
    
    b = reshape(b,[],10);
    ind1 = find(b(:,2)>=0);
    ind2 = find(b(:,2) <0);
    
    save('cal08.mat','a','b','c','d','ind1','ind2','-v7.3');
end
