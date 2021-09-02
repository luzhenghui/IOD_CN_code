clc;
clear;

cal01_revision_data = 0;
cal01_revision_eof = 0;
cal01_revision_filter = 0;

if cal01_revision_data == 1
    path = '../../myDATA/sst.day.anom.v2.nc/';
    
    lon = double(ncread([path,'sst.day.mean.1983.nc'],'lon',1,inf,8));
    lat = flip(double(ncread([path,'sst.day.mean.1983.nc'],'lat',1,inf,8)));

    Nlon = length(lon);
    Nlat = length(lat);
    
    sst  = zeros(Nlon,Nlat,365,40);
    sst(:,:,:,:) = nan;
    
    for i = 1:40
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
    
    x = reshape(ssta((lon(:)>=40 & lon(:)<=110),abs(lat(:))<=11,:,:),[],365*40);
    x = bsxfun(@minus,x,nanmean(x,1));
    x = bsxfun(@minus,x,nanmean(x,2));
    
    save('cal01_revision_data.mat','x','-v7.3');

end

if cal01_revision_eof == 1
    load('cal01_revision_data.mat');
    m = size(x,1);
    
    ind1 = find(isnan(x(:,1))==1);
    ind2 = find(isnan(x(:,2))==0);
    
    x(ind1,:) = [];
    
    [lam, lds, pcs, per] = eof_dan(x',length(ind2));
    
    pcs_b = zeros(size(x,2),length(ind2));    
    pcs_r = zeros(size(x,2),length(ind2));    
    
    for i = 1:length(ind2)
        if i == 1
            pcs_b(:,i) = lanczosfilter(pcs(:,i),1,1/1100,1100,'low');
            pcs_r(:,i) = pcs(:,i) - pcs_b(:,i);
        else
            pcs_b(:,i) = pcs(:,i);
            pcs_r(:,i) = pcs(:,i);
        end
    end
    
    y = (pcs_b*lds^(-1))';
    
    XB = zeros(m,14235);
    XB(:,:) = nan;
    XB(ind2,:) = y;

    y = (pcs_r*lds^(-1))';

    XR = zeros(m,14235);
    XR(:,:) = nan;
    XR(ind2,:) = y;    
    
    for k = 2:38
        
        disp(k);
        [~,~,~,w,r,t] = link_cal2(k,XB,365,200,365,30,'1','Pearson');
        save(['cal01_IOD_NOAA_',num2str(k+1982),'_v2_revision_remove_eof1_L1100d.mat'],'w','r','t','-v7.3');
        clear w r t;
        
        [~,~,~,w,r,t] = link_cal2(k,XR,365,200,365,30,'1','Pearson');
        save(['cal01_IOD_NOAA_',num2str(k+1982),'_v2_revision_remove_eof1_H1100d.mat'],'w','r','t','-v7.3');
        clear w r t;        
        
    end
end

if cal01_revision_filter == 1    
    load('cal01_revision_data.mat');
    m = size(x,1);
    
    ind1 = find(isnan(x(:,1))==1);
    ind2 = find(isnan(x(:,2))==0);
    
    x(ind1,:) = [];
    
    [lam, lds, pcs, per] = eof_dan(x',length(ind2));
    
    pcs_b = zeros(size(x,2),length(ind2));    
    pcs_r = zeros(size(x,2),length(ind2));    
    
    for i = 1:length(ind2)
        if i == 1
            pcs_b(:,i) = lanczosfilter(pcs(:,i),1,1/400,400,'low')-...
                lanczosfilter(pcs(:,i),1,1/1000,1000,'low');
            
            pcs_r(:,i) = pcs(:,i) - pcs_b(:,i);
        else
            pcs_b(:,i) = pcs(:,i);
            pcs_r(:,i) = pcs(:,i);
        end
    end
    
    y = (pcs_b*lds^(-1))';
    
    XB = zeros(m,14235);
    XB(:,:) = nan;
    XB(ind2,:) = y;

    y = (pcs_r*lds^(-1))';

    XR = zeros(m,14235);
    XR(:,:) = nan;
    XR(ind2,:) = y;    
    
    for k = 2:38
        
        disp(k);
        [~,~,~,w,r,t] = link_cal2(k,XB,365,200,365,30,'1','Pearson');
        save(['cal01_IOD_NOAA_',num2str(k+1982),'_v2_revision_filter_band_400d_1000d.mat'],'w','r','t','-v7.3');
        clear w r t;
        
        [~,~,~,w,r,t] = link_cal2(k,XR,365,200,365,30,'1','Pearson');
        save(['cal01_IOD_NOAA_',num2str(k+1982),'_v2_revision_filter_band_400d_1000d_residual.mat'],'w','r','t','-v7.3');
        clear w r t;        
        
    end
end