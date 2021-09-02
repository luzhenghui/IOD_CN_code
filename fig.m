clear;

fig1 = 0;
fig2 = 0;
fig3 = 0;
fig4 = 0;
figS1 = 0;
figS2 = 0;
figS3 = 0;
figS4 = 0;
figS5 = 0;
figS6 = 0;
figS7 = 0;
figS8 = 0;
figS9 = 0;
figS10 = 0;
figS11 = 0;
figS12 = 0;

if fig1 == 1
    path = '../../myDATA/sst.day.anom.v2.nc/';
    
    lon = double(ncread([path,'sst.day.mean.1983.nc'],'lon',1,inf,8));
    lat = flip(double(ncread([path,'sst.day.mean.1983.nc'],'lat',1,inf,8)));    
    
    load('cal08.mat');
    tmp = ind1;
    ind1 = ind2;
    ind2 = tmp;
    
    [Lon1,Lat1] = meshgrid(lon(lon(:)>=40 & lon(:)<=110),lat(abs(lat(:))<=11));
    Lon1 = Lon1';
    Lat1 = Lat1';    
    
    Lat2 = Lat1(:);
    Lon2 = Lon1(:);
    
    m = 385;

    X = zeros(m,m,12,36);
    T = zeros(m,m,12,36);
    
    for i = 1:36
        load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2.mat'],'w','t');
        
        X(:,:,:,i) = w;
        T(:,:,:,i) = t;
    end
    
    XX = X;
    XX(abs(T)>150) = nan;
    K = nanmean(nanmean(XX,4),3);
    
    K1 = nansum(K(ind1,:),1);
    K1(ind1) = nan;
    K1 = reshape(K1,35,11);
    K1(K1==0) = nan;
    
    K2 = nansum(K(:,ind2),2);
    K2(ind2) = nan;
    K2 = reshape(K2,35,11);
    K2(K2==0) = nan;
    
    figure(1)
    subplot(5,4,1:8)
    m_proj('miller','lon',[40 110],'lat',[-15 15]);
    m_grid('xtick',40:10:110,'ytick',20:-5:-20,'fontsize',10,'tickdir','out');
    hold on;
    m_coast('patch',[0.9,0.9,0.9],'edgecolor','k','linewidth',1);
    hold on;
    m_plot(Lon2(ind1),Lat2(ind1),'ro','markersize',3,'markerfacecolor','r');
    hold on;
    m_plot(Lon2(ind2),Lat2(ind2),'bo','markersize',3,'markerfacecolor','b');
    hold on;
    m_text(60,0,'A','fontsize',25,'fontweight','bold','horizontalalignment','center');
    hold on;
    m_text(100,-5,'B','fontsize',25,'fontweight','bold','horizontalalignment','center');
    hold on;
    m_line([50,50,70,70,50],[-10,10,10,-10,-10],'linewidth',2,'color','k');
    hold on;
    m_line([90,90,110,110,90],[-10,0,0,-10,-10],'linewidth',2,'color','k');
    m_text(42,8,'(a)','fontsize',10);

    [Lon1,Lat1] = meshgrid(lon(lon(:)>=40 & lon(:)<=112),lat(abs(lat(:))<=11));
    Lon1 = Lon1';
    Lat1 = Lat1';
    
    K1 = cat(1,K1,K1(end,:));
     
    subplot(5,4,[9,10,13,14]);
    m_proj('miller','lon',[40 110],'lat',[-15 15]);
    m_grid('xtick',40:10:110,'ytick',20:-5:-20,'fontsize',10,'tickdir','out');
    hold on;
    m_pcolor(Lon1,Lat1,K1);
    shading flat;
    hold on;
    m_coast('patch',[0.9,0.9,0.9],'edgecolor','k','linewidth',1);
    hold on;
    m_line([50,50,70,70,50],[-10,10,10,-10,-10],'linewidth',2,'color','k');
    caxis([-100,100]);
    m_text(42,8,'(b)','fontsize',10);

    K2 = cat(1,K2,K2(end,:));    
    
    subplot(5,4,[11,12,15,16]);
    m_proj('miller','lon',[40 110],'lat',[-15 15]);
    m_grid('xtick',40:10:110,'ytick',20:-5:-20,'fontsize',10,'tickdir','out','yaxislocation','right');
    hold on;
    m_pcolor(Lon1,Lat1,K2);
    shading flat;
    hold on;
    m_coast('patch',[0.9,0.9,0.9],'edgecolor','k','linewidth',1);    
    hold on;
    m_line([90,90,110,110,90],[-10,0,0,-10,-10],'linewidth',2,'color','k');
    caxis([-100,100]);
    m_text(42,8,'(c)','fontsize',10);
    
    h = colorbar('location','southoutside');
    set(h,'position',[0.2,0.22,0.6,0.03]);
    
    print(gcf,'-dpng', 'Figure01.png');
    print(gcf,'-depsc','Figure01.eps');
    
end

if fig2 == 1
   load('cal02_K_v2_T150.mat','K');
   load('min_K_v2_T150.mat');
   load('../../myDATA/IOD/SINTEX_DMI.csv');
   
   dmi = SINTEX_DMI(:,2);
   t = 1982:1/12:2020.999;
   dmi(t<1984)= [];
   t(t<1984)  = [];
   dmi = reshape(dmi,12,[]);
   dmi = bsxfun(@minus,dmi,nanmean(dmi,2));  
   y1 = smooth(dmi,3);
   y2 = reshape(y1,12,[]);
   th1 = 0.5;
   ind_p = find(y1>=+th1);
   ind_n = find(y1<=-th1);
   [ind_p2,n_p,ts_p,te_p] = event_series(ind_p,3);
   [ind_n2,n_n,ts_n,te_n] = event_series(ind_n,3);   
   tt = reshape(t,12,[]);
   T = 1984:1:2020;
   
   load('cal07_10Wind.mat','U','lon','lat');
   Hc = reshape(nanmean(nanmean(U(lon>=70 & lon<=100, lat>=-5 & lat<=5,:,:),1),2),[],1);    
   Hc = reshape(Hc,12,[]);
   
   scrsz = get(0,'screensize');
   figure('position',scrsz);

   subplot(3,8,1:5);    
   h = bar(1984:1:2020,min_K,'histc');
   set(h,'facecolor','k','edgecolor','w');
   alpha(0.2);
   hold on;
   plot(1984:1/12:2020.999,K(:),'k.-','markersize',10);
   hold on;
   plot(1984.5:1:2020.5,th_K,'ro-','linewidth',1,'markersize',4);
   ylabel('Total Degree (x10^3)','fontsize',10);
   xlim([1984,2022]);
   ylim([-8000,4000]);   
   hold on;
   j = 1;
   for i = 1:37
       if min_K(i)<th_K(i)
           plot(tt(12,i)-0.5,-7000,'markeredgecolor','none','marker','^','markersize',8,'markerfacecolor',[50,205,50]/255);
           hold on;
%            if j == 3 || j == 7 || j == 9 || j == 10 || j == 12 || j == 17
%                plot([tt(12,i)-0.5,tt(12,i)+0.5],[-7000,-8000],'color',[50,205,50]/255,'linestyle',':');
%            elseif j == 18
%                plot([tt(12,i)-0.5,tt(12,i)+0.5],[-7000,-8000],'color','k','linestyle','-');
%            else               
               plot([tt(12,i)-0.5,tt(12,i)+0.5],[-7000,-8000],'color',[50,205,50]/255,'linestyle','-','linewidth',1);
%            end
           j = j+1;
       end
       hold on;
   end
   set(gca,'position',[0.06,0.65,0.53,0.2],'xtick',1984:1:2021,'xticklabel',{'';'';'';'';'';'';...
       '';'';'';'';'';'';'';'';'';'';'';'';'';'';'';'';'';'';'';'';'';...
       '';'';'';'';'';'';'';'';'';''},'ytick',-8000:2000:4000,'yticklabel',{'-8';'-6';'-4';'-2';'0';'2';'4'},...
       'fontsize',10,'tickdir','in');   
   text(1985,3000,'(a)','fontsize',10);
   
   subplot(3,8,9:13);
   min_K(min_K>th_K) = nan;
   plot(t,y1,'color','k');
   xlim([1984,2022])
   ylim([-2,2]);
   ylabel('DMI (^oC)','fontsize',10);
%    xlabel('Year','fontsize',10);
   hold on;
   h = bar(T+1,min_K*0-2,'histc');
   set(h,'facecolor',[50,205,50]/255,'edgecolor','w');
   alpha(0.2);
   hold on;
   h = bar(T+1,min_K*0+2,'histc');
   set(h,'facecolor',[50,205,50]/255,'edgecolor','w');
   alpha(0.2);
   hold on;   
   tmp = zeros(37,1);
   tmp(:) = nan;
   tmp([11,15,23,27,30,36]) = min_K([11,15,23,27,30,36]);
   h = bar(T+1,tmp*0+2,'histc');   
   set(h,'facecolor','w','linestyle','none');
   hatch(h,30,'k',':',3);   
   set(h,'edgecolor','w');
   hold on;
   h = bar(T+1,tmp*0-2,'histc');   
   set(h,'facecolor','w','linestyle','none'); 
   hatch(h,30,'k',':',3);   
   hold on;
   tmp = zeros(37,1);
   tmp(:) = nan;
   tmp(37) = min_K(37);
   h = bar(T+1,tmp*0+2,'histc');   
   set(h,'facecolor','w','linestyle','none');
   hatch(h,30,'k','-',5,0.5);
   hold on;
   h = bar(T+1,tmp*0-2,'histc');   
   set(h,'facecolor','w','linestyle','none');      
   hatch(h,30,'k','-',5,0.5);
   hold on;
   plot(1984:2022,zeros(39,1)-th1,'k:',1984:2022,zeros(39,1)+th1,'k:');
   hold on;
   y = y1;
   y(setdiff(1:length(y),ind_p2)) = +th1;
   area(t,y,+th1,'facecolor','r','linestyle','none');
   hold on;
   y = y1;
   y(setdiff(1:length(y),ind_n2)) = -th1;
   area(t,y,-th1,'facecolor','b','linestyle','none');
   hold on;
   plot(t,y1,'color','k','linewidth',1);
   axis tight;
   
   set(gca,'position',[0.06,0.5,0.53,0.15],'xtick',1984:1:2021,'xticklabel',{'';'';'';'';'';'';...
       '1990';'';'';'';'';'1995';'';'';'';'';'2000';'';'';'';'';'2005';'';'';'';'';'2010';...
       '';'';'';'';'2015';'';'';'';'';'2020';''},'ytick',-1:1:1,'fontsize',10,'tickdir','in');
   text(1985,1.5,'(b)','fontsize',10);
    
   subplot(3,8,17:21);
   load('calS15_neg.mat');
   [ax,h1,h2] = plotyy(t,smooth(K,3),t,smooth(kn,3));
   set(ax(1),'ycolor','k');
   set(ax(2),'ycolor','b');
   set(h1,'linewidth',1.5,'color','k');
   set(h2,'linewidth',1.5,'color','b');
   xlim(ax(1),[1984,2022]);
   xlim(ax(2),[1984,2022]);
   ylim(ax(1),[-8000,4000]);
   ylim(ax(2),[400,4000]);
   ylabel(ax(1),'Total Degree (x10^3)','fontsize',10);   
%    ylabel(ax(2),'Number of Links','fontsize',10);
%    set(ax(2),'yaxis','left');
   h = legend(h2,'Number of Links (x10^3)','location','southwest');
   set(h,'box','on','fontsize',10);
   set(ax(1),'position',[0.06,0.25,0.53,0.25],'xtick',1984:1:2021,'xticklabel',{'';'1985';'';'';'';'';...
       '1990';'';'';'';'';'1995';'';'';'';'';'2000';'';'';'';'';'2005';'';'';'';'';'2010';...
       '';'';'';'';'2015';'';'';'';'';'2020';''},'ytick',-8000:2000:2000,'yticklabel',{'-8';'-6';'-4';'-2';'0';'2'},...
       'fontsize',10,'tickdir','in');   
   set(ax(2),'position',[0.06,0.3,0.53,0.2],'xtick',1984:1:2021,'xticklabel',{'';'1985';'';'';'';'';...
       '1990';'';'';'';'';'1995';'';'';'';'';'2000';'';'';'';'';'2005';'';'';'';'';'2010';...
       '';'';'';'';'2015';'';'';'';'';'2020';''},'ytick',400:600:3400,'yticklabel',...
       {'0.4';'1.0';'1.6';'2.2';'2.8';'3.4'},'fontsize',10,'tickdir','in');
   text(1985,3000,'(c)','fontsize',10);
   text(1995,3000,['Corr. Coef. : ',num2str(corr(smooth(kn,3),smooth(K,3)),'%.2f')],'fontsize',10,'color','k');
   
%%   
   subplot(3,8,[6:8,14:16,22:24]);    
   for i = 1:36
       if min_K(i)<th_K(i)

           tmp = abs(y2(:,i+1));
           [~,ic] = max(tmp);
           if nanmean(y2(ic,i+1),1)>0.5
               plot(Hc(12,i),nanmean(y2(ic,i+1),1),'ro','markersize',20,'markerfacecolor','r');
           elseif nanmean(y2(ic,i+1),1)<0.5 && nanmean(y2(ic,i+1),1)>-0.5
               plot(Hc(12,i),nanmean(y2(ic,i+1),1),'ko','markersize',20,'markerfacecolor','k','linewidth',2);
           else
               plot(Hc(12,i),nanmean(y2(ic,i+1),1),'bo','markersize',20,'markerfacecolor','b');
           end              
           text(Hc(12,i),nanmean(y2(ic,i+1),1),num2str(i+1984),'color','w','horizontalalignment',...
               'center','fontsize',8,'fontweight','bold')
           hold on;
       end
   end
   plot(-4:0.1:4,zeros(81,1),'k-',zeros(41,1),-2:0.1:2,'k-',-4:0.1:4,ones(81,1)*0.5,'k--',-4:0.1:4,-ones(81,1)*0.5,'k--');
   xlim([-4,4]);
   ylim([-2,2]);
   set(gca,'yaxis','right');
   axis tight;
   
   xlabel('10m U in Last Dec. Prior to Predicted Year (m/s)','fontsize',10);
   ylabel('DMI in Predicted Year (^oC)','fontsize',10);
   set(gca,'position',[0.63,0.3,0.30,0.55],'fontsize',10,'xtick',-4:1:4,'ytick',-2:0.5:2,...
       'yticklabel',{'-2.0';'-1.5';'-1.0';'-0.5';'0.0';'0.5';'1.0';'1.5';'2.0'});
   text(-3.8,1.8,'(d)','fontsize',10);
   
   print(gcf,'-dpng', 'Figure02.png');
   print(gcf,'-depsc','Figure02.eps');
end

if fig3 == 1
    lon = double(ncread('../../myDATA/Heat_content/Monthly/CZ16_1_2000m_Temp_year_1980_month_01.nc','lon'));
    lat = double(ncread('../../myDATA/Heat_content/Monthly/CZ16_1_2000m_Temp_year_1980_month_01.nc','lat'));
    
    Nlon = length(lon);
    Nlat = length(lat);
    
    [Lon,Lat] = meshgrid(lon,lat);
    Lon = Lon';
    Lat = Lat';
    
    load('cal06_ohc_CZ_H300.mat');
    %%
    load('../../Research/color/BlueWhiteOrangeRed.rgb');
    BlRe = BlueWhiteOrangeRed/255;
    
    load('min_K_v2_T150.mat');
    T = 1984:1:2019;
    t = 1984:1/12:2019.999;
    
    neg_y = [8,12,14,21,32];
    pos_y = [10,13,28,33,34,35];
    
    neg_m = min_loc([8,12,14,21,32]);
    pos_m = min_loc([10,13,28,33,34,35]);

    %%
    pos_dy = zeros(Nlon,Nlat,12);
    pos_dy(:,:,1) = nanmean(X(:,:,12,pos_y),4);
    for i = 1:11
        pos_dy(:,:,i+1) = nanmean(X(:,:,i,pos_y+1),4);
    end
    
    neg_dy = zeros(Nlon,Nlat,12);
    neg_dy(:,:,1) = nanmean(X(:,:,12,neg_y),4);
    for i = 1:11
        neg_dy(:,:,i+1) = nanmean(X(:,:,i,neg_y+1),4);
    end
    
    %%
    load('cal06_ohc_CZ_D20.mat','X');
    pos_d20 = zeros(Nlon,Nlat,12);
    pos_d20(:,:,1) = nanmean(X(:,:,12,pos_y),4);
    for i = 1:11
        pos_d20(:,:,i+1) = nanmean(X(:,:,i,pos_y+1),4);
    end
    
    neg_d20 = zeros(Nlon,Nlat,12);
    neg_d20(:,:,1) = nanmean(X(:,:,12,neg_y),4);
    for i = 1:11
        neg_d20(:,:,i+1) = nanmean(X(:,:,i,neg_y+1),4);
    end
    
    %%
    load('cal07_10Wind.mat','U','V','lat','lon');
    pos_u = zeros(length(lon),length(lat),12);
    pos_v = zeros(length(lon),length(lat),12);
    pos_u(:,:,1) = nanmean(U(:,:,12,pos_y),4);
    pos_v(:,:,1) = nanmean(V(:,:,12,pos_y),4);
    for i = 1:11
        pos_u(:,:,i+1) = nanmean(U(:,:,i,pos_y+1),4);
        pos_v(:,:,i+1) = nanmean(V(:,:,i,pos_y+1),4);
    end
    
    neg_u = zeros(length(lon),length(lat),12);
    neg_v = zeros(length(lon),length(lat),12);
    neg_u(:,:,1) = nanmean(U(:,:,12,neg_y),4);
    neg_v(:,:,1) = nanmean(V(:,:,12,neg_y),4);
    for i = 1:11
        neg_u(:,:,i+1) = nanmean(U(:,:,i,neg_y+1),4);
        neg_v(:,:,i+1) = nanmean(V(:,:,i,neg_y+1),4);        
    end
    
    [Lon1,Lat1] = meshgrid(lon,lat);
    Lon1 = Lon1';
    Lat1 = Lat1';
    %%
    
    cmin = -10;
    cmax =  10;
        
    figure(1)

    tit = {'(a) nIOD Dec(-1)';'(b) nIOD Feb(0)';'(c) nIOD Apr(0)';'(d) nIOD Jun(0)';'(e) nIOD Aug(0)';'(f) nIOD Oct(0)';...
           '(g) pIOD Dec(-1)';'(h) pIOD Feb(0)';'(i) pIOD Apr(0)';'(j) pIOD Jun(0)';'(k) pIOD Aug(0)';'(l) pIOD Oct(0)'};
    
    for j = 1:2
        if j == 1
            tmp1 = neg_d20;
            tmp3 = neg_u;
            tmp4 = neg_v;
        else
            tmp1 = pos_d20;
            tmp3 = pos_u;
            tmp4 = pos_v;
        end
        
        for i = 1:6
            subplot(4,3,i+(j-1)*6);
            h = pcolor(Lon,Lat,tmp1(:,:,(i-1)*2+1)); 
            set(h,'linestyle','none');
            caxis([cmin,cmax]);
            colormap(BlRe);
            title(tit{i+(j-1)*6},'fontsize',10);
            axis([30,120,-20,20]);
            hold on;
            quiver(Lon1(1:4:end,1:4:end),Lat1(1:4:end,1:4:end),tmp3(1:4:end,1:4:end,(i-1)*2+1),tmp4(1:4:end,1:4:end,(i-1)*2+1),1,'color','k');
            geoshow('landareas.shp','FaceColor', [0.8 0.8 0.8],'linestyle','none');
            set(gca,'ytick',-20:10:20,'yticklabel',{'20°„S';'10°„S';'0';'10°„N';'20°„N'},...
                'xtick',30:30:120,'xticklabel',{'30°„E';'60°„E';'90°„E';'120°„E'},'fontsize',8);
            
         
            
            
        end
    end
    
    h = colorbar;
    set(h,'position',[0.94,0.2,0.02,0.6]);
    
%     h = axes;
%     quiver(65,0,20,0,1,'color','k','linewidth',2);
%     set(h,'position',[0.41,0.05,0.2,0.02]);
%     axis([30,120,-20,20]);
%     axis off;
%     text(65,-35,'20m/s','fontsize',10);

    h = axes;
    quiver(65,0,20,0,1,'color','k','linewidth',2);
    set(h,'position',[0.85,0.88,0.2,0.02]);
    axis([30,120,-20,20]);
    axis off;
    text(65,-35,'20m/s','fontsize',10);
    text(72,-140,'m','fontsize',10);
    
    print(gcf,'-dpng', 'Figure03_D20.png');
    print(gcf,'-depsc','Figure03_D20.eps');    
    
end

if fig4 == 1
   load('min_K_v2_T150.mat');
   min_K(min_K>th_K) = nan;
   T = 1984:1:2020;
   min_loc(isnan(min_K)==1) = nan;
   min_loc([11,15,23,27,30]) = nan;
   
   load('../../myDATA/IOD/SINTEX_DMI.csv');
   dmi = SINTEX_DMI(:,2);
   t = 1982:1/12:2020.999;
   dmi(t<1984)= [];
   t(t<1984)  = [];
   dmi = reshape(dmi,12,[]);
   dmi = bsxfun(@minus,dmi,nanmean(dmi,2));
   th1 = 0.5;
   y1 = smooth(dmi,3);
   ind_p = find(y1>=+th1);
   ind_n = find(y1<=-th1);
   [ind_p2,n_p,ts_p,te_p] = event_series(ind_p,3);
   [ind_n2,n_n,ts_n,te_n] = event_series(ind_n,3);
   
   figure(1)
   subplot(5,1,[1,2]);
   min_K(min_K>th_K) = nan;
   plot(t,y1,'color','k');
   xlim([1984,2022])
   ylim([-2,2]);
   ylabel('DMI (^oC)','fontsize',10);
%    xlabel('Year','fontsize',10);   
   hold on;
   plot(1984:2022,zeros(39,1),'k',1984:2022,zeros(39,1)-th1,'k:',1984:2022,zeros(39,1)+th1,'k:');
   
   hold on;
   tmp = min_K;
   tmp(:) = nan;
   tmp([8,11,12,14,21,23,32,36]) = min_K([8,11,12,14,21,23,32,36]);
   h = bar(T+1,tmp*2,'histc');
   set(h,'facecolor','b','edgecolor','w');
   alpha(0.2);
   h = bar(T+1,tmp*-2,'histc');
   set(h,'facecolor','b','edgecolor','w');
   alpha(0.2);   
   hold on;
   
   tmp = min_K;
   tmp(:) = nan;
   tmp([10,13,15,27,28,30,33,34,35,37]) = min_K([10,13,15,27,28,30,33,34,35,37]);
%    tmp([10,13,15,24,25,28,30,33,34,35]) = min_K([10,13,15,24,25,28,30,33,34,35]);
   h = bar(T+1,tmp*2,'histc');
   set(h,'facecolor','r','edgecolor','w');
   alpha(0.2);
   h = bar(T+1,tmp*-2,'histc');
   set(h,'facecolor','r','edgecolor','w');
   alpha(0.2);
   
   hold on;
   tmp = min_K;
   tmp([1:10,12:14,16:22,24:26,28:29,32:35,37]) = nan;
   h = bar(T+1,tmp*2,'histc');
   set(h,'facecolor','w','linestyle','none');
   hatch(h,30,'k',':',3);
   alpha(0.2);
   h = bar(T+1,tmp*-2,'histc');
   set(h,'facecolor','w','linestyle','none');
   hatch(h,30,'k',':',3);   
   alpha(0.2);
   
   hold on;
   tmp = min_K;
   tmp(1:36) = nan;
   h = bar(T+1,tmp*2,'histc');
   set(h,'facecolor','w','linestyle','none');
   hatch(h,30,'k','-',5,0.5);   
   alpha(0.2);   
   h = bar(T+1,tmp*-2,'histc');
   set(h,'facecolor','w','linestyle','none');
   hatch(h,30,'k','-',5,0.5);   
   alpha(0.2);
   
   hold on;
   y = y1;
   y(setdiff(1:length(y),ind_p2)) = +th1;
   area(t,y,+th1,'facecolor','r','linestyle','none');
   hold on;
   y = y1;
   y(setdiff(1:length(y),ind_n2)) = -th1;
   area(t,y,-th1,'facecolor','b','linestyle','none');  
   hold on;
   plot(t,y1,'color','k','linewidth',1);
   set(gca,'xtick',1984:1:2020,'xticklabel',{'';'1985';'';'';'';'';...
       '1990';'';'';'';'';'1995';'';'';'';'';'2000';'';'';'';'';'2005';'';'';'';'';'2010';...
       '';'';'';'';'2015';'';'';'';'';'2020'},'ytick',-2:1:2,'fontsize',10);
   text(1985,1.5,'(a)','fontsize',12);
   
   subplot(5,1,[3,4,5]);
   hit = [11/15,5/6,6/9];
   fal = [6/17,3/8,3/9];

%    hit = [12/15,5/6,7/9];
%    fal = [5/17,2/7,3/10];

   c = {'k';'b';'r'};
    for i = 1:3
        plot(fal(i),hit(i),c{i},'marker','o','markersize',12,'markeredgecolor',c{i},'linestyle','none','linewidth',2);
        hold on;
    end

    hit = [1/14,0/6,1/8];
    fal = [1/2,1/1,0/1];
    for i = 1:3
        plot(fal(i),hit(i),c{i},'marker','s','markersize',12,'markeredgecolor',c{i},'linestyle','none','linewidth',2);
        hold on;
    end
    
    hit = [6/11,2/5,4/6];
    fal = [4/10,2/4,2/6];
    for i = 1:3
        plot(fal(i),hit(i),c{i},'marker','*','markersize',12,'markeredgecolor',c{i},'linestyle','none','linewidth',2);
        hold on;
    end
    
    hit = [4/16,3/6,1/10];
    fal = [5/9,1/4,4/5];
    for i = 1:3
        plot(fal(i),hit(i),c{i},'marker','^','markersize',12,'markeredgecolor',c{i},'linestyle','none','linewidth',2);
        hold on;
    end
    ylabel('Hit Rate','fontsize',12);
    xlabel('False Alarm Rate','fontsize',12);
    axis([0,1,0,1]);
    
    hold on;
    plot(0:0.1:1,zeros(11,1)+0.6,'k','linestyle','--');
    hold on;
    plot(zeros(11,1)+0.4,0:0.1:1,'k','linestyle','--');

%     hold on;
%     plot(0:0.1:1,zeros(11,1)+0.7,'k','linestyle','--');
%     hold on;
%     plot(zeros(11,1)+0.3,0:0.1:1,'k','linestyle','--');

    set(gca,'fontsize',10,'ytick',0:0.1:1,'xtick',0:0.1:1,'yticklabel',...
        {'0.0';'0.1';'0.2';'0.3';'0.4';'0.5';'0.6';'0.7';'0.8';'0.9';'1.0'},...
         'position',[0.13,0.1,0.88,0.46]);
    text(0.05,0.9,'(b)','fontsize',12);
      
    h = legend('TD_{Total} @Dec(-1)','TD_{nIOD} @Dec(-1)','TD_{pIOD} @Dec(-1)',...
           'BCC\_CSM_{Total} @Dec(-1)','BCC\_CSM_{nIOD} @Dec(-1)','BCC\_CSM_{pIOD} @Dec(-1)',...
           'ECs5_{Total} @May(0)','ECs5_{nIOD} @May(0)','ECs5_{pIOD} @May(0)',...
           'CFSv2_{Total} @Feb(0)','CFSv2_{nIOD} @Feb(0)','CFSv2_{pIOD} @Feb(0)');
    set(h,'location','eastoutside','fontsize',9);

    print(gcf,'-dpng', 'Figure04.png');
    print(gcf,'-depsc','Figure04.eps');
    
    
end

if figS1 == 1
    path = '../../myDATA/sst.day.anom.v2.nc/';
    sst_nan = ncread([path,'sst.ltm.1971-2000.nc'],'sst');
    sst_nan(abs(sst_nan)>100) = nan;
    sst_nan(isnan(sst_nan)~=1) = 0;
    
    lon = double(ncread([path,'sst.mnmean.nc'],'lon'));
    lat = double(ncread([path,'sst.mnmean.nc'],'lat'));
    
    Nlon = length(lon);
    Nlat = length(lat);
    
    lon1 = 40:1:110;
    lat1 = 10:-1:-10;
    
    [Lon1,Lat1] = meshgrid(lon1,lat1);
    Lon1 = Lon1';
    Lat1 = Lat1';    
    
    tmp = ncread([path,'sst.mnmean.nc'],'sst');
    tmp(:,:,1) = [];
    tmp(abs(tmp)>100) = nan;
    tmp = bsxfun(@plus,tmp,sst_nan(:,:,1));
    sst = reshape(tmp,Nlon,Nlat,12,[]);
    ssta = bsxfun(@minus,sst,mean(sst(:,:,:,1:29),4));
    ssta = reshape(ssta,Nlon,Nlat,[]);
    ssta = interpn(lon,lat,1:size(ssta,3),ssta,lon1',lat1',1:size(ssta,3));
   
    x = reshape(ssta,[],12*38);
    [a,b,c,d] = eof_dan(x',10);
    b = reshape(b,71,21,10,[]);
    
    xmean = nanmean(x,1);

    dmi = load('../../myDATA/IOD/SINTEX_DMI.csv');
    dmi = dmi(:,2);
    t1 = 1982:1/12:2019.999;
    
    figure(1)
    for i = 1:2
        subplot(4,6,[1,2,3,7,8,9]+(i-1)*12);
        m_proj('miller','lon',[40 110],'lat',[-15 15]);
        m_grid('xtick',40:10:110,'ytick',20:-5:-20,'fontsize',9,'tickdir','out');
        hold on;
        h = m_pcolor(Lon1,Lat1,-b(:,:,i));
        set(h,'linestyle','none');
        shading interp;
        hold on;
        m_contour(Lon1,Lat1,-b(:,:,i),'color','k','showtext','on','levellist',-0.1:0.01:0.1);
        hold on;
        m_contour(Lon1,Lat1,-b(:,:,i),'color','k','showtext','on','levellist',0.00,'linewidth',2);
        caxis([-0.04,0.04]);
        h = colorbar('location','southout');
        hold on;
        m_coast('patch',[0.9,0.9,0.9],'edgecolor','k','linewidth',1);
        hold on;
        text(-0.1,0.3,['EOF',num2str(i),' ',num2str(d(i),'%.01f'),'%'],'fontsize',10);
        if i == 1
            set(gca,'position',[0.06,0.6,0.38,0.25]);
            m_text(42,8,'(a)','fontsize',10);
        else
            set(gca,'position',[0.06,0.2,0.38,0.25]);
            m_text(42,8,'(c)','fontsize',10);            
        end
        
    end
    
    subplot(4,6,[4:6,10:12]);
    plot(t1, -c(:,1)/40,'k',t1,xmean,'r','linewidth',2);
    xlim([1980,2021]);
    ylim([-1,1]);
    set(gca,'fontsize',10,'ytick',-1:0.5:1,'xtick',1980:5:2020,'yticklabel',...
        {'-1.0';'-0.5';' 0.0';' 0.5';' 1.0'});
    text(1982,0.8,['(b) Corr. = ',num2str(corr(-c(:,1),xmean'),'%.2f')],'fontsize',10);
    set(gca,'position',[0.5,0.55,0.49,0.29]);
    legend('PC1/40','Area Mean (^oC)','location','southeast');
%     ylabel(ax(1),'PC1','fontsize',10);
%     ylabel(ax(2),'Regional Mean (^oC)','fontsize',10);    
    
    subplot(4,6,[16:18,22:24]);
    plot(t1,c(:,2)/10,'k',t1,dmi,'r','linewidth',2);
    xlim([1980,2021])
    ylim([-3,3]);
    set(gca,'fontsize',10,'ytick',-3:1:3,'xtick',1980:5:2020,...
        'yticklabel',{'-3.0';'-2.0';'-1.0';'0.0';' 1.0';' 2.0';' 3.0'});   
    text(1982,2.5,['(d) Corr. = ',num2str(corr(c(:,2),dmi),'%.2f')],'fontsize',10);
    set(gca,'position',[0.5,0.15,0.49,0.29]);
    legend('PC2/10','DMI (^oC)','location','southeast');
%     ylabel(ax(1),'PC2','fontsize',10);
%     ylabel(ax(2),'DMI (^oC)','fontsize',10);    

    print(gcf,'-dpng', 'figS1.png');
    print(gcf,'-depsc','figS1.eps');
    
    
end

if figS2 == 1
   load('cal02_K_v2_T150.mat','K');
   K = reshape(K,12,[]);
   KC = mean(K,2);
   [f,s,s0k_r,s0k_w] = psd_ynm(K(:),13);
 
   figure(1)
   subplot('position',[0.1,0.4,0.38,0.32]);
   semilogx(1./f,s(2:end),'k','linewidth',2);
   hold on;
   semilogx(1./f,flip(s0k_r(2:end)),'r--','linewidth',2);    
   xlim([1,1000]);
   ylim([0,0.15]);
   set(gca,'ytick',0:0.03:0.15,'yticklabel',{'0.00';'0.03';'0.06';'0.09';'0.12';'0.15'});
   ylabel('Spectral Density Analysis','fontsize',14);
   xlabel('Period (mo)','fontsize',14);
   text(1.3,0.137,'(a)','fontsize',14);
   
   subplot('position',[0.6,0.4,0.38,0.32]);
   plot(1:12,KC,'k.-','linewidth',2,'markersize',20);
   xlim([1,12]);
   ylim([-2000,-800]);
   set(gca,'ytick',-2000:200:-800,'xtick',1:12,'xticklabel',{'Jan';'Feb';'Mar';'Apr';'May';...
       'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'});
   ylabel('Total Degree','fontsize',14);
   text(1.5,-900,'(b)','fontsize',14);
   
%    xlabel('Month','fontsize',12);   
   h = gca;
   th = rotateticklabel(h,60);
    
    print(gcf,'-dpng', 'FigureS2.png');
    print(gcf,'-depsc','FigureS2.eps'); 
   
end

if figS3 == 1
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
    
    for i = 1:37
        load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2.mat'],'w','t','r');
        X1(:,:,:,i) = w;
        T1(:,:,:,i) = t;
    end
    
    XX = X1;
    XX(abs(T1)>150) = nan;
    
    PP = XX;
%     PP(PP>0) = nan;
    
    k1d = squeeze(nansum(PP(ind1,:,:,:),1));
    k1d(setdiff(1:m,ind2),:,:) = nan;
    k1d = reshape(k1d,35,11,[]);

    k2d = squeeze(nansum(PP(:,ind2,:,:),2));
    k2d(setdiff(1:m,ind1),:,:) = nan;
    k2d = reshape(k2d,35,11,[]);    

    PP = XX;
    PP(PP>0) = nan;    
    
    k1s = squeeze(nanmean(PP(ind1,:,:,:),1));
    k1s(setdiff(1:m,ind2),:,:) = nan;
    k1s = reshape(k1s,35,11,[]);
    
    k2s = squeeze(nanmean(PP(:,ind2,:,:),2));
    k2s(setdiff(1:m,ind1),:,:) = nan;
    k2s = reshape(k2s,35,11,[]);
    
    PP(isnan(PP)==0) = 1;
    
    k1n = squeeze(nansum(PP(ind1,:,:,:),1));
    k1n(setdiff(1:m,ind2),:,:) = nan;
    k1n = reshape(k1n,35,11,[]);
    
    k2n = squeeze(nansum(PP(:,ind2,:,:),2));
    k2n(setdiff(1:m,ind1),:,:) = nan;
    k2n = reshape(k2n,35,11,[]);
       
    
    load('cal02_K_v2_T150.mat','K');    
    load('min_K_v2_T150.mat');
    
    th = zeros(12,37);
    for i = 1:37
        th(:,i) = th_K(i);
    end
    th = reshape(th,[],1);
    
    ind_anom = [];
    ind_norm = [];
    for i = 1:length(K)-1;
        if K(i)<=th(i)
            ind_anom = [ind_anom,i];
        elseif K(i)>th(i)
            ind_norm = [ind_norm,i];
        end
    end    
    
    load('../../research/color/blRe.rgb');
    blRe = blRe/255;
    scrsz = get(0,'screensize');
    figure('position',scrsz);
    x1 = -160;
    x2 =  160;
    y1 =  -60;
    y2 =  60;
    z1 = -3;
    z2 = -2 ;
    for i = 1:9
        subplot(3,3,i);
        m_proj('miller','lon',[40 110],'lat',[-15 15]);
        m_grid('xtick',40:10:110,'ytick',20:-5:-20,'fontsize',10,'tickdir','out');
        hold on;
        m_coast('patch',[0.9,0.9,0.9],'edgecolor',[0.9,0.9,0.9],'linewidth',0.5);
        hold on;
        colormap(jet);
        if i == 1
            h = m_contourf(lon1,lat1,mean(k1d(:,:,ind_anom),3),'color','w','linewidth',2,'levellist',-160:20:0);
            hold on;
            m_contourf(lon1,lat1,mean(k2d(:,:,ind_anom),3),'color','w','linewidth',2,'levellist',-160:20:0);            
%             set(gca,'position',[0.1,0.75,0.28,0.15]);
            text(0,1.7,'Outgoing Node Degree','fontsize',15,'HorizontalAlignment','center');
            title('(a) Anomaly','fontsize',12);
            caxis([x1,x2]);
            h2 = colorbar('location','southoutside');
            set(h2,'fontsize',10);
        elseif i == 2
            m_contourf(lon1,lat1,nanmean(k1n(:,:,ind_anom),3),'color','w','linewidth',2,'levellist',0:5:70);
            hold on;
            m_contourf(lon1,lat1,nanmean(k2n(:,:,ind_anom),3),'color','w','linewidth',2,'levellist',0:5:70);         
%             set(gca,'position',[0.4,0.75,0.28,0.15]);
            text(0,1.7,'Number of Outgoing-links','fontsize',15,'HorizontalAlignment','center');
            title('(d) Anomaly','fontsize',12);
            caxis([y1,y2]);
            h2 = colorbar('location','southoutside');
            set(h2,'fontsize',10);            
        elseif i == 3
            m_contourf(lon1,lat1,nanmean(k1s(:,:,ind_anom),3),'color','w','linewidth',2,'levellist',-3:0.1:3);
            hold on;
            m_contourf(lon1,lat1,nanmean(k2s(:,:,ind_anom),3),'color','w','linewidth',2,'levellist',-3:0.1:3);            
%             set(gca,'position',[0.7,0.75,0.28,0.15]);
            text(0,1.7,'Average Outgoing-link Strength','fontsize',15,'HorizontalAlignment','center');
            title('(g) Anomaly','fontsize',12);
            caxis([z1,z2]);
            h2 = colorbar('location','southoutside');
            set(h2,'fontsize',10);            
        elseif i == 4
            m_contourf(lon1,lat1,nanmean(k1d(:,:,ind_norm),3),'color','w','linewidth',2,'levellist',-160:20:0);
            hold on;
            m_contourf(lon1,lat1,nanmean(k2d(:,:,ind_norm),3),'color','w','linewidth',2,'levellist',-160:20:0);            
%             shading interp;
            %             set(gca,'position',[0.2,0.47,0.3,0.2]);
            title('(b) Normal','fontsize',12);
            caxis([x1,x2]);
            h2 = colorbar('location','southoutside');
            set(h2,'fontsize',10);
        elseif i == 5
            m_contourf(lon1,lat1,nanmean(k1n(:,:,ind_norm),3),'color','w','linewidth',2,'levellist',0:5:70);
            hold on;
            m_contourf(lon1,lat1,nanmean(k2n(:,:,ind_norm),3),'color','w','linewidth',2,'levellist',0:5:70);           
%             shading interp;
            %             set(gca,'position',[0.2,0.47,0.3,0.2]);
            title('(e) Normal','fontsize',12);
            caxis([y1,y2]);
            h2 = colorbar('location','southoutside');
            set(h2,'fontsize',10);            
        elseif i == 6
            m_contourf(lon1,lat1,nanmean(k1s(:,:,ind_norm),3),'color','w','linewidth',2,'levellist',-3:0.1:3);
            hold on;
            m_contourf(lon1,lat1,nanmean(k2s(:,:,ind_norm),3),'color','w','linewidth',2,'levellist',-3:0.1:3);           
%             shading interp;
            %             set(gca,'position',[0.2,0.47,0.3,0.2]);
            title('(h) Normal','fontsize',12);

            caxis([z1,z2]);
            h2 = colorbar('location','southoutside');
            set(h2,'fontsize',10);            
        elseif i == 7
            m_contourf(lon1,lat1,nanmean(k1d(:,:,ind_anom),3)-nanmean(k1d(:,:,ind_norm),3),'color','w','linewidth',2,'levellist',-160:20:0);
            hold on;
            m_contourf(lon1,lat1,nanmean(k2d(:,:,ind_anom),3)-nanmean(k2d(:,:,ind_norm),3),'color','w','linewidth',2,'levellist',-160:20:0);
            
            %             set(gca,'position',[0.2,0.24,0.3,0.2]);
            title('(c) Difference','fontsize',12);
            colorbar;

            caxis([x1,x2]);
            h2 = colorbar('location','southoutside');
            set(h2,'fontsize',10);
%             
            x1 = nanmean(k1d(:,:,ind_anom),3);
            x2 = nanmean(k1d(:,:,ind_norm),3);
            n1 = length(ind_anom);
            n2 = length(ind_norm);
            s2 = (sum(bsxfun(@minus,k1d(:,:,ind_anom),x1).^2,3)+sum(bsxfun(@minus,k1d(:,:,ind_norm),x2).^2,3))/(n1+n2-2);
            test = (x1-x2)./sqrt(s2)/sqrt(1/n1+1/n2);
            ind11 = find(abs(test(:))>2.576);
            hh = m_plot(lon1(ind11),lat1(ind11),'k+','markersize',2);

            hold on;
            x1 = nanmean(k2d(:,:,ind_anom),3);
            x2 = nanmean(k2d(:,:,ind_norm),3);
            n1 = length(ind_anom);
            n2 = length(ind_norm);
            s2 = (sum(bsxfun(@minus,k2d(:,:,ind_anom),x1).^2,3)+sum(bsxfun(@minus,k2d(:,:,ind_norm),x2).^2,3))/(n1+n2-2);
            test = (x1-x2)./sqrt(s2)/sqrt(1/n1+1/n2);
            ind11 = find(abs(test(:))>2.576);
            hh = m_plot(lon1(ind11),lat1(ind11),'k+','markersize',2);            
   
        elseif i == 8
            h = m_contourf(lon1,lat1,nanmean(k1n(:,:,ind_anom),3)-nanmean(k1n(:,:,ind_norm),3),'color','w','linewidth',2,'levellist',0:5:70);
            hold on;
            h = m_contourf(lon1,lat1,nanmean(k2n(:,:,ind_anom),3)-nanmean(k2n(:,:,ind_norm),3),'color','w','linewidth',2,'levellist',0:5:70);
            

            %             set(gca,'position',[0.2,0.24,0.3,0.2]);
            title('(f) Difference','fontsize',12);
            caxis([y1,y2]);
            h2 = colorbar('location','southoutside');
            set(h2,'fontsize',10);
            
            x1 = nanmean(k1n(:,:,ind_anom),3);
            x2 = nanmean(k1n(:,:,ind_norm),3);
            n1 = length(ind_anom);
            n2 = length(ind_norm);
            s2 = (nansum(bsxfun(@minus,k1n(:,:,ind_anom),x1).^2,3)+nansum(bsxfun(@minus,k1n(:,:,ind_norm),x2).^2,3))/(n1+n2-2);
            test = (x1-x2)./sqrt(s2)/sqrt(1/n1+1/n2);
            ind11 = find(abs(test(:))>2.576);
            hh = m_plot(lon1(ind11),lat1(ind11),'k+','markersize',2);            
            hold on;
            x1 = nanmean(k2n(:,:,ind_anom),3);
            x2 = nanmean(k2n(:,:,ind_norm),3);
            n1 = length(ind_anom);
            n2 = length(ind_norm);
            s2 = (nansum(bsxfun(@minus,k2n(:,:,ind_anom),x1).^2,3)+nansum(bsxfun(@minus,k2n(:,:,ind_norm),x2).^2,3))/(n1+n2-2);
            test = (x1-x2)./sqrt(s2)/sqrt(1/n1+1/n2);
            ind11 = find(abs(test(:))>2.576);
            hh = m_plot(lon1(ind11),lat1(ind11),'k+','markersize',2);            
            
            
            
        else
            h = m_contourf(lon1,lat1,nanmean(k1s(:,:,ind_anom),3)-nanmean(k1s(:,:,ind_norm),3),'color','w','linewidth',2,'levellist',-3:0.1:3);
            hold on;
            h = m_contourf(lon1,lat1,nanmean(k2s(:,:,ind_anom),3)-nanmean(k2s(:,:,ind_norm),3),'color','w','linewidth',2,'levellist',-3:0.1:3);
            
%             set(gca,'position',[0.2,0.24,0.3,0.2]);
            title('(i) Difference','fontsize',12);
            caxis([0,0.5]);
            h2 = colorbar('location','southoutside');
            set(h2,'fontsize',10);
%             
            x1 = nanmean(k1s(:,:,ind_anom),3);
            x2 = nanmean(k1s(:,:,ind_norm),3);
            n1 = length(ind_anom);
            n2 = length(ind_norm);
            s2 = (sum(bsxfun(@minus,k1s(:,:,ind_anom),x1).^2,3)+sum(bsxfun(@minus,k1s(:,:,ind_norm),x2).^2,3))/(n1+n2-2);
            test = (x1-x2)./sqrt(s2)/sqrt(1/n1+1/n2);
            ind11 = find(abs(test(:))>2.576);
            hh = m_plot(lon1(ind11),lat1(ind11),'k+','markersize',2);            
            hold on;
            x1 = nanmean(k2s(:,:,ind_anom),3);
            x2 = nanmean(k2s(:,:,ind_norm),3);
            n1 = length(ind_anom);
            n2 = length(ind_norm);
            s2 = (sum(bsxfun(@minus,k2s(:,:,ind_anom),x1).^2,3)+sum(bsxfun(@minus,k2s(:,:,ind_norm),x2).^2,3))/(n1+n2-2);
            test = (x1-x2)./sqrt(s2)/sqrt(1/n1+1/n2);
            ind11 = find(abs(test(:))>2.576);
            hh = m_plot(lon1(ind11),lat1(ind11),'k+','markersize',2);            
        end
        hold on;
        m_coast('patch',[0.9,0.9,0.9],'edgecolor',[0.9,0.9,0.9],'linewidth',0.5);
%         m_grid('xtick',40:10:110,'ytick',20:-5:-20,'fontsize',7,'tickdir','out');
    end
    
    
   print(gcf,'-dpng', 'FigureS3.png');
   print(gcf,'-depsc','FigureS3.eps');    

end

if figS4 == 1
    path = '../../myDATA/sst.day.anom.v2.nc/';
    
    lon = double(ncread([path,'sst.day.mean.1983.nc'],'lon',1,inf,8));
    lat = flip(double(ncread([path,'sst.day.mean.1983.nc'],'lat',1,inf,8)));
    
    [lon1,lat1] = meshgrid(lon(lon>=40 & lon<=110),lat(abs(lat)<=11));
    lon1 = lon1';
    lat1 = lat1';
    
    ind1 = find(lon1(:)>=90 & lon1(:)<=110 & lat1(:)<=0);
    ind2 = find(lon1(:)>=50 & lon1(:)<=70);

    m = 385;

    X = zeros(m,m,12,37);
    T = zeros(m,m,12,37);
    
    TD = zeros(4,444);
    
    for j = 1:4
        for i = 1:37
            if j == 1
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1.mat'],'w','t');                
            elseif j == 2
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof2.mat'],'w','t');
            elseif j == 3
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof3.mat'],'w','t');
            else               
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2.mat'],'w','t');
            end
                
            X(:,:,:,i) = w;
            T(:,:,:,i) = t;
        end
        XX = X;
        XX(abs(T)>150) = nan;
        TD(j,:) = smooth(reshape(nansum(nansum(XX(ind1,ind2,:,:),1),2),[],1),3);          
    end 
        
    t = 1984:1/12:2020.999;
 
    scrsz = get(0,'screensize');
    figure('position',scrsz);    

    subplot(2,1,1,'position',[0.1,0.5,0.8,0.35]);
    plot(t,smooth(TD(4,:),3),'k',t,smooth(TD(1,:),3),'r',t,smooth(TD(2,:),3),'b',t,smooth(TD(3,:),3),'g',...
        'linewidth',2);
    xlim([1984,2021]);
    ylim([-8000,4000]);
    ylabel('Total Degree','fontsize',12);
    legend('Orignal','PC1 Removed','PC2 Removed','PC3 Removed','location','south','orientation','horizontal');
    set(gca,'xtick',1984:1:2020,'xticklabel',{'';'1985';'';'';'';'';'1990';'';'';'';'';'1995';...
        '';'';'';'';'2000';'';'';'';'';'2005';'';'';'';'';'2010';'';'';'';'';'2015';'';'';'';'';'2020';''},'fontsize',10);    
    text(1985,3000,'(a)','fontsize',12);
    
    path = '../../myDATA/sst.day.anom.v2.nc/';
    
    lon = double(ncread([path,'sst.day.mean.1983.nc'],'lon',1,inf,8));
    lat = flip(double(ncread([path,'sst.day.mean.1983.nc'],'lat',1,inf,8)));
    
    [lon1,lat1] = meshgrid(lon(lon>=40 & lon<=110),lat(abs(lat)<=11));
    lon1 = lon1';
    lat1 = lat1';
    
    ind1 = find(lon1(:)>=90 & lon1(:)<=110 & lat1(:)<=0);
    ind2 = find(lon1(:)>=50 & lon1(:)<=70);

    m = 385;

    X = zeros(m,m,12,37);
    T = zeros(m,m,12,37);
    
    TD_L = zeros(13,444);
    TD_H = zeros(13,444);
    
    for j = 1:13
        for i = 1:37
            if j == 1
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1_L100d.mat'],'w','t');                
            elseif j == 2
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1_L200d.mat'],'w','t');
            elseif j == 3
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1_L300d.mat'],'w','t');
            elseif j == 4
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1_L400d.mat'],'w','t');
            elseif j == 5
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1_L500d.mat'],'w','t'); 
            elseif j == 6
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1_L600d.mat'],'w','t');                 
            elseif j == 7
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1_L700d.mat'],'w','t'); 
            elseif j == 8
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1_L800d.mat'],'w','t');             
            elseif j == 9
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1_L900d.mat'],'w','t');             
            elseif j == 10
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1_L1000d.mat'],'w','t');
            elseif j == 11
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1_L1100d.mat'],'w','t');
            elseif j == 12
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1_L1200d.mat'],'w','t');
            else                
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2.mat'],'w','t');
            end
                
            X(:,:,:,i) = w;
            T(:,:,:,i) = t;
        end
        XX = X;
        XX(abs(T)>150) = nan;
        TD_L(j,:) = reshape(nansum(nansum(XX(ind1,ind2,:,:),1),2),[],1);
        
        for i = 1:37
            if j == 1
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1_H100d.mat'],'w','t');                
            elseif j == 2
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1_H200d.mat'],'w','t');
            elseif j == 3
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1_H300d.mat'],'w','t');
            elseif j == 4
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1_H400d.mat'],'w','t');
            elseif j == 5
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1_H500d.mat'],'w','t'); 
            elseif j == 6
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1_H600d.mat'],'w','t');                 
            elseif j == 7
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1_H700d.mat'],'w','t'); 
            elseif j == 8
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1_H800d.mat'],'w','t');            
            elseif j == 9
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1_H900d.mat'],'w','t');
            elseif j == 10
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1_H1000d.mat'],'w','t');
            elseif j == 11
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1_H1100d.mat'],'w','t');
            elseif j == 12
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_remove_eof1_H1200d.mat'],'w','t');            
            else 
                load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2.mat'],'w','t');
            end
                
                
            X(:,:,:,i) = w;
            T(:,:,:,i) = t;
        end
        XX = X;
        XX(abs(T)>150) = nan;
        TD_H(j,:) = reshape(nansum(nansum(XX(ind1,ind2,:,:),1),2),[],1);            
    end
    
    a = corr(TD_L');
    b = corr(TD_H');

    subplot(2,1,2,'position',[0.1,0.1,0.8,0.35]);
    plot(100:100:1300,a(13,:),'ro-','linewidth',2);
    hold on;
    plot(100:100:1300,b(13,:),'bo-','linewidth',2);
    legend('Low-pass','High-pass','location','southwest');    
    hold on;
    plot(100:100:1300,ones(13,1)*0.9,'k--','linewidth',2);    
    xlim([100,1200]);
    ylim([0,1]);
    set(gca,'ytick',0:0.1:1,'yticklabel',{'0.0';'0.1';'0.2';'0.3';'0.4';'0.5';...
        '0.6';'0.7';'0.8';'0.9';'1.0'},'xtick',100:100:1200,'fontsize',10);
    ylabel('Correlation','fontsize',12);
    xlabel('Filter Period (day)','fontsize',12);
    text(120,0.8,'(b)','fontsize',12);
    
    path = '../../myDATA/sst.day.anom.v2.nc/';
    
    lon = double(ncread([path,'sst.day.mean.1983.nc'],'lon',1,inf,8));
    lat = flip(double(ncread([path,'sst.day.mean.1983.nc'],'lat',1,inf,8)));
    
    [lon1,lat1] = meshgrid(lon(lon>=40 & lon<=110),lat(abs(lat)<=11));
    lon1 = lon1';
    lat1 = lat1';
    
    ind1 = find(lon1(:)>=90 & lon1(:)<=110 & lat1(:)<=0);
    ind2 = find(lon1(:)>=50 & lon1(:)<=70);

    m = 385;

    X = zeros(m,m,12,37);
    T = zeros(m,m,12,37);
    
    for i = 1:37
        load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_filter_band_400d_1000d.mat'],'w','t');        
        X(:,:,:,i) = w;
        T(:,:,:,i) = t;
    end
    XX = X;
    XX(abs(T)>150) = nan;
    TD_L = reshape(nansum(nansum(XX(ind1,ind2,:,:),1),2),[],1);
 
    X = zeros(m,m,12,37);
    T = zeros(m,m,12,37);
    
    for i = 1:37
        load(['cal01_IOD_NOAA_',num2str(i+1983),'_v2_revision_filter_band_400d_1000d_residual.mat'],'w','t');        
        X(:,:,:,i) = w;
        T(:,:,:,i) = t;
    end
    XX = X;
    XX(abs(T)>150) = nan;
    TD_H = reshape(nansum(nansum(XX(ind1,ind2,:,:),1),2),[],1);
       
    load('cal02_K_v2_T150.mat','K');    
    t = 1984:1/12:2020.999;
    
    scrsz = get(0,'screensize');
    figure('position',scrsz);    

    subplot(2,1,2,'position',[0.1,0.1,0.8,0.35]);
    plot(t,smooth(K,3),'k',t,smooth(TD_L,3),'r',t,smooth(TD_H,3),'c','linewidth',2);
    xlim([1984,2021]);
    ylabel('Total Degree','fontsize',12);
    legend('Orignal','Residual of Band-pass of PC1 Removed','Band-pass of PC1 Removed','location','south','orientation','horizontal');
    set(gca,'xtick',1984:1:2020,'xticklabel',{'';'1985';'';'';'';'';'1990';'';'';'';'';'1995';...
        '';'';'';'';'2000';'';'';'';'';'2005';'';'';'';'';'2010';'';'';'';'';'2015';'';'';'';'';'2020';''},'fontsize',10);    
    text(1985,3000,'(d)','fontsize',12);
    
    load('../../myDATA/IOD/SINTEX_DMI.csv');
    dmi = SINTEX_DMI(:,2);
    t = 1982:1/12:2020.999;
    dmi(t<1984)= [];
    t(t<1984)  = [];
    dmi = reshape(dmi,12,[]);
    dmi = reshape(bsxfun(@minus,dmi,nanmean(dmi,2)),[],1);
    th1 = 0.5;
    y1 = smooth(dmi,3);
    ind_p = find(y1>=+th1);
    ind_n = find(y1<=-th1);
    [ind_p2,n_p,ts_p,te_p] = event_series(ind_p,3);
    [ind_n2,n_n,ts_n,te_n] = event_series(ind_n,3);

    load('cal01_revision_data.mat');
    m = size(x,1);
    
    ind1 = find(isnan(x(:,1))==1);
    ind2 = find(isnan(x(:,2))==0);
    
    x(ind1,:) = [];
    
    [lam, lds, pcs, per] = eof_dan(x',length(ind2));
    
    pcs_b1 = lanczosfilter(pcs(:,1),1,1/400,400,'low')-lanczosfilter(pcs(:,1),1,1/1000,1000,'low');
   
    load('cal02_K_v2_T150.mat','K');
    t1 = 1982:1/365:2020.999;    
    t2 = 1984:1/12:2020.999;
    
    s1 = unique(floor(t(ts_p)));
    s2 = unique(floor(t(ts_n)));
    
    subplot(2,1,1,'position',[0.1,0.5,0.8,0.35]);
    for i = 1:length(s1)
    h = bar([s1(i)-1/12,s1(i)+0.5-1/12],15*ones(2,1),'histc');
    set(h,'facecolor','r','edgecolor','none');
    alpha(0.2);
    hold on;
    h = bar([s1(i)-1/12,s1(i)+0.5-1/12],-15*ones(2,1),'histc');
    set(h,'facecolor','r','edgecolor','none');
    alpha(0.2);
    hold on;
    plot(s1(i)*ones(100,1)+1-1/12,linspace(-15,15,100),'w');
    hold on;    
    end
    hold on;
    for i = 1:length(s2)
    h = bar([s2(i)-1/12,s2(i)+0.5-1/12],15*ones(2,1),'histc');
    set(h,'facecolor','b','edgecolor','none');
    alpha(0.2);
    hold on;
    h = bar([s2(i)-1/12,s2(i)+0.5-1/12],-15*ones(2,1),'histc');
    set(h,'facecolor','b','edgecolor','none');
    alpha(0.2);
    hold on;
    plot(s2(i)*ones(100,1)+1-1/12,linspace(-15,15,100),'w');
    hold on;
    end
    h1 = plot(t1,pcs(:,1),'k',t1,pcs_b1,'r','linewidth',1);
    legend(h1,'PC1','Band-pass of PC1','location','south','orientation','horizontal');
    plot(t1,pcs_b1,'r','linewidth',2);
    hold on;
    plot(t1,zeros(length(t1),1),'k--');
    xlim([1984,2021]);
    ylim([-12,12]);
    set(gca,'xtick',1984:1:2020,'xticklabel',{'';'1985';'';'';'';'';'1990';'';'';'';'';'1995';...
        '';'';'';'';'2000';'';'';'';'';'2005';'';'';'';'';'2010';'';'';'';'';'2015';'';'';'';'';'2020';''},'ytick',-15:3:15,'fontsize',10);    
    text(1985,11,'(c)','fontsize',12);

    print(gcf,'-dpng', 'FigureS4.png');
    print(gcf,'-depsc', 'FigureS4.eps');

end

if figS5 == 1
    load('../../myDATA/IOD/SINTEX_DMI.csv');
    dmi = SINTEX_DMI(:,2);
    t = 1982:1/12:2020.999;
    dmi(t<1984)= [];
    t(t<1984)  = [];
    dmi = reshape(dmi,12,[]);
    dmi = reshape(bsxfun(@minus,dmi,nanmean(dmi,2)),[],1);
    th1 = 0.5;
    y1 = smooth(dmi,3);
    ind_p = find(y1>=+th1);
    ind_n = find(y1<=-th1);
    [ind_p2,n_p,ts_p,te_p] = event_series(ind_p,3);
    [ind_n2,n_n,ts_n,te_n] = event_series(ind_n,3);

    load('cal01_revision_data.mat');
    m = size(x,1);
    
    ind1 = find(isnan(x(:,1))==1);
    ind2 = find(isnan(x(:,2))==0);
    
    x(ind1,:) = [];
    
    [lam, lds, pcs, per] = eof_dan(x',length(ind2));
    
    pcs_b1 = lanczosfilter(pcs(:,1),1,1/400,400,'low')-lanczosfilter(pcs(:,1),1,1/1000,1000,'low');
%     pcs_b2 = smooth(pcs(:,1),400,'moving')-smooth(pcs(:,1),1000,'moving');

    [s1,f1] = psd_ynm(pcs(:,1),1);
    [s2,f2] = psd_ynm(pcs_b1,1);
    [s3,f3] = psd_ynm(dmi(:),1);
    
    scrsz = get(0,'screensize');
    figure('position',scrsz);    

    subplot(2,2,1:2,'position',[0.1,0.6,0.83,0.3]);    
    semilogx(1./s2/365,f2(2:end),'r',1./s3/12,f3(2:end),'g',1./s1/365,f1(2:end),'k','linewidth',2);
    ylim([0,0.14]);
    xlim([0.01,100]);
    h = legend('Band-pass of PC1 (Daily)','DMI (Monthly)','PC1 (Daily)','location','northeast');
    ylabel('Power Spectral Density','fontsize',12);
    xlabel('Period (Year)','fontsize',12);
    set(h,'Interpreter','none');
    set(gca,'fontsize',10,'ytick',0:0.02:0.14,'yticklabel',{'0.00';'0.02';'0.04';'0.06';'0.08';...
        '0.10';'0.12';'0.14'});
    text(0.02,0.13,'(a)','fontsize',12);
     
    load('../../myDATA/IOD/SINTEX_DMI.csv');
    dmi = SINTEX_DMI(:,2);
    t = 1982:1/12:2020.999;
    dmi(t<1984)= [];
    t(t<1984)  = [];
    dmi = reshape(dmi,12,[]);
    dmi = reshape(bsxfun(@minus,dmi,nanmean(dmi,2)),[],1);
    y1 = smooth(dmi,3);

    load('cal01_revision_data.mat');
    m = size(x,1);
    
    ind1 = find(isnan(x(:,1))==1);
    ind2 = find(isnan(x(:,2))==0);
    
    x(ind1,:) = [];    
    [lam, lds, pcs, per] = eof_dan(x',length(ind2));
    
    pcs_b1 = lanczosfilter(pcs(:,1),1,1/400,400,'low')-lanczosfilter(pcs(:,1),1,1/1000,1000,'low');
    pcs_b1 = reshape(pcs_b1,365,[]);
    pcs_b1(:,1:2) = [];
    pcs1 = reshape(pcs(:,1),365,[]);
    pcs1(:,1:2) = [];
    
    dmi = reshape(dmi,12,[]);
    
    neg_all = [1992,1996,1998,2005,2010,2016];
    pos_all = [1994,1997,2006,2008,2012,2015,2017,2018,2019];
    
    neg = [1992,1996,1998,2005,2016];
    pos = [1994,1997,2012,2017,2018,2019];
     
    load('cal02_K_v2_T150.mat','K');
    t1 = 1982:1/365:2020.999;    
    t2 = 1984:1/12:2020.999;
    
    K = reshape(K,12,[]);   
    T = 1984:1:2020;
    
    neg_ind = neg-1983;
    pos_ind = pos-1983;

    p1 = [0.1,0.15,0.33,0.35];
    subplot(2,2,3,'position',p1);
    h = bar([0,1],3*ones(1,2),'histc');
    set(h,'facecolor','r','edgecolor','none');
    alpha(0.1);    
    hold on;
    h = bar([0,1],-3*ones(1,2),'histc');
    set(h,'facecolor','r','edgecolor','none');
    alpha(0.1);
    hold on;
    plot((-364:365)/365,-ones(730,1)*std(pcs_b1(:)),'r--',(-364:365)/365,ones(730,1)*std(pcs_b1(:)),'r--',...
         (-364:365)/365,ones(730,1)*std(dmi(:))*2,'g--',(-364:365)/365,-ones(730,1)*std(dmi(:))*2,'g--',...
         (-364:365)/365,-ones(730,1)*std(pcs1(:))/3,'k--',(-364:365)/365,ones(730,1)*std(pcs1(:))/3,'k--');
    hold on;
    xlim([-1,1]);
    ylim([-3,3]);
    text(-0.9,2.6,'(b)','fontsize',12);    
    [ax,hlines] = plotyyy((-364:365)/365,[mean(pcs_b1(:,pos_ind-1),2);mean(pcs_b1(:,pos_ind),2)],...
                            (-11:1:12)/12,[mean(dmi(:,pos_ind-1),2);mean(dmi(:,pos_ind),2)],...
                            (-364:365)/365,[mean(pcs1(:,pos_ind-1),2);mean(pcs1(:,pos_ind),2)],...
                            p1,{'Band-pass of PC1','DMI (^oC)','PC1'});
    xlim([-1,1]);
    ylim(ax(1),[-3,3]);
    ylim(ax(2),[-1.5,1.5]);
    ylim(ax(3),[-9,9]);
    set(ax(1),'ycolor','r','ytick',-3:3,'fontsize',10);
    set(ax(2),'ycolor','g','ytick',-1.5:0.5:1.5,'fontsize',10,'yticklabel',{'-1.5';'-1.0';'-0.5';' 0.0';' 0.5';' 1.0';' 1.5'});
    set(ax(3),'ycolor','k','ytick',-9:3:9,'fontsize',10,'yticklabel',{'-9';'-6';'-3';' 0';' 3';' 6';' 9'});
    set(hlines(1),'color','r','linewidth',2);
    set(hlines(2),'color','g','linewidth',2);
    set(hlines(3),'color','k','linewidth',2);
    ylabel(ax(1),'Band-pass of PC1','fontsize',10);
    ylabel(ax(2),'DMI (^oC)','fontsize',10);
    ylabel(ax(3),'PC1','fontsize',10);   
    xlabel(ax(1),'Year','fontsize',10);
    
    p2 = [0.57,0.15,0.33,0.35];
    subplot(2,2,4,'position',p2);    
    h = bar([0,1],3*ones(1,2),'histc');
    set(h,'facecolor','b','edgecolor','none');
    alpha(0.1);    
    hold on;
    h = bar([0,1],-3*ones(1,2),'histc');
    set(h,'facecolor','b','edgecolor','none');
    alpha(0.1);
    hold on;
    plot((-364:365)/365,-ones(730,1)*std(pcs_b1(:)),'r--',(-364:365)/365,ones(730,1)*std(pcs_b1(:)),'r--',...
         (-364:365)/365,ones(730,1)*std(dmi(:))*2,'g--',(-364:365)/365,-ones(730,1)*std(dmi(:))*2,'g--',...
         (-364:365)/365,-ones(730,1)*std(pcs1(:))/3,'k--',(-364:365)/365,ones(730,1)*std(pcs1(:))/3,'k--');
    hold on;
    text(-0.9,2.6,'(c)','fontsize',12);    
    xlim([-1,1]);
    ylim([-3,3]);    
    [ax,hlines] = plotyyy((-364:365)/365,[mean(pcs_b1(:,neg_ind-1),2);mean(pcs_b1(:,neg_ind),2)],...
                            (-11:1:12)/12,[mean(dmi(:,neg_ind-1),2);mean(dmi(:,neg_ind),2)],...
                            (-364:365)/365,[mean(pcs1(:,neg_ind-1),2);mean(pcs1(:,neg_ind),2)],...
                            p2,{'Band-pass of PC1','DMI (^oC)','PC1'});
    xlim([-1,1]);
    ylim(ax(1),[-3,3]);
    ylim(ax(2),[-1.5,1.5]);
    ylim(ax(3),[-9,9]);
    set(ax(1),'ycolor','r','ytick',-3:3,'fontsize',10);
    set(ax(2),'ycolor','g','ytick',-1.5:0.5:1.5,'fontsize',10,'yticklabel',{'-1.5';'-1.0';'-0.5';' 0.0';' 0.5';' 1.0';' 1.5'});
    set(ax(3),'ycolor','k','ytick',-9:3:9,'fontsize',10,'yticklabel',{'-9';'-6';'-3';' 0';' 3';' 6';' 9'});
    set(hlines(1),'color','r','linewidth',2);
    set(hlines(2),'color','g','linewidth',2);
    set(hlines(3),'color','k','linewidth',2);
    ylabel(ax(1),'Band-pass of PC1','fontsize',10);
    ylabel(ax(2),'DMI (^oC)','fontsize',10);
    ylabel(ax(3),'PC1','fontsize',10);   
    xlabel(ax(1),'Year','fontsize',10);
   
    print(gcf,'-dpng', 'FigureS5.png');
    print(gcf,'-depsc', 'FigureS5.eps');    
    
    
end

if figS6 == 1
    load('cal01_revision_data.mat');
    m = size(x,1);
    
    ind1 = find(isnan(x(:,1))==1);
    ind2 = find(isnan(x(:,2))==0);
    
    x(ind1,:) = [];    
    [~, ~, pcs, ~] = eof_dan(x',length(ind2));
    clear x;
    
    pcs_b1 = lanczosfilter(pcs(:,1),1,1/400,400,'low')-lanczosfilter(pcs(:,1),1,1/1000,1000,'low');
    pcs_b1 = reshape(pcs_b1,365,[]);
    pcs_b1(:,1:2) = [];    
    pcs_b1 = reshape(nanmean(reshape(pcs_b1(1:360,:),30,12,[]),1),[],1);

     
%    scrsz = get(0,'screensize');
%    figure('position',scrsz);
   for i = 1:6
        if i == 1 || i == 2
            load('cal07_10Wind.mat','U','lon','lat');
            X = U;
            a = (i-1)*3;
        elseif i == 3 || i == 4
            load('cal07_10Wind.mat','V','lon','lat');
            X = V;
            a = (i-3)*3;
        else
            load('cal06_ohc_CZ_D20.mat','X','lon','lat');
            a = (i-5)*3;
        end

        
        [LON,LAT] = meshgrid(lon,lat);
        X = reshape(X,length(lon)*length(lat),[]);        
           
        subplot(3,2,i);
        m_proj('robinson','lon',[0 360],'lat',[-90 90]);
        m_grid('xtick',0:90:360,'ytick',-90:30:90,'fontsize',8,'tickdir','out');
        hold on;
        m_coast('patch',[0.9,0.9,0.9],'edgecolor','k','linewidth',1);
        hold on;
        
        c = corr(pcs_b1((1+a):444,1),X(:,1:(444-a))');
        c = reshape(c,length(lon),length(lat)); 
        c(abs(c)<0.2)= nan;
                
        h = m_pcolor(LON,LAT,c');
        set(h,'linestyle','none');
        caxis([-1,1]);
        if i == 1
            title('U(0)','fontsize',12);
        elseif i == 2
            title('U(-3)','fontsize',12);
        elseif i == 3
            title('V(0)','fontsize',12);
        elseif i == 4
            title('V(-3)','fontsize',12);
        elseif i == 5
            title('D20(0)','fontsize',12);
        elseif i == 6
            title('D20(-3)','fontsize',12);
        end    
    end
    axes;
    axis off;
    
    h = colorbar;
    set(h,'position',[0.93,0.2,0.02,0.6],'ytick',-1:0.2:1,...
        'yticklabel',{'-1.0';'-0.8';'-0.6';'-0.4';'-0.2';' 0.0';' 0.2';' 0.4';' 0.6';' 0.8';' 1.0'});
    caxis([-1,1]);
    
    print(gcf,'-dpng', 'FigureS6.png');
    print(gcf,'-depsc','FigureS6.eps');    
    
end

if figS7 == 1

    load('../../Research/color/BlueWhiteOrangeRed.rgb');
    BlRe = BlueWhiteOrangeRed/255;

    load('min_K_v2_T150.mat');   
    min_K(min_K>th_K) = nan;
    T = 1984:1:2019;
    t = 1984:1/12:2019.999;
    
    neg_y = [8,12,14,21,32];
    pos_y = [10,13,28,33,34,35];
    
    neg_m = min_loc([8,12,14,21,32]);
    pos_m = min_loc([10,13,28,33,34,35]);
        
    load('cal07_10Wind.mat','U','lon','lat');
    
    [Lon1,Lat1] = meshgrid(lon,lat);
    Lon1 = Lon1';
    Lat1 = Lat1';
    
    U(U<0) = -1;
    U(U>0) =  1;
    
    Nlon = length(lon);
    Nlat = length(lat);
    
    tmp1 = zeros(Nlon,Nlat,length(neg_y));
    tmp2 = tmp1;
    for i = 1:length(neg_y)
        tmp1(:,:,i) = squeeze(U(:,:,neg_m(i),neg_y(i)));
        tmp2(:,:,i) = squeeze(U(:,:,      12,neg_y(i)));
    end
 
    tmp3 = zeros(Nlon,Nlat,length(pos_y));
    tmp4 = tmp3;
    for i = 1:length(pos_y)
        tmp3(:,:,i) = squeeze(U(:,:,pos_m(i),pos_y(i)));
        tmp4(:,:,i) = squeeze(U(:,:,      12,pos_y(i)));
    end
    
    Z(:,:,1) = nanmean(tmp1,3);
    Z(:,:,2) = nanmean(tmp2,3);
    Z(:,:,3) = nanmean(tmp3,3);
    Z(:,:,4) = nanmean(tmp4,3);
    Z(:,:,5) = Z(:,:,1).*Z(:,:,3);
    Z(:,:,6) = Z(:,:,2).*Z(:,:,4);
    
    cmin = -1;
    cmax =  1;
    
    figure(1)
    for i = 1:3
        subplot(3,5,[2:4]+(i-1)*5);
        h = pcolor(Lon1,Lat1,Z(:,:,i*2));
        set(h,'linestyle','none');
        caxis([cmin,cmax]);
        hold on;
        colormap(BlRe);
        axis([30,120,-20,20]);
        hold on;
        geoshow('landareas.shp','FaceColor', 'none');
        set(gca,'ytick',-20:10:20,'yticklabel',{'20°„S';'10°„S';'0';'10°„N';'20°„N'},...
            'xtick',30:10:120,'xticklabel',{'30°„E';'40°„E';'50°„E';'60°„E';'70°„E';'80°„E';'90°„E';'100°„E';'110°„E';'120°„E'},'fontsize',8);
        if i == 1
            title('(a) nIOD U@Dec','fontsize',12);
        elseif i == 2
            title('(b) pIOD U@Dec','fontsize',12);
        else
            title('(c) nIOD U@Dec X pIOD U@Dec','fontsize',12);
            hold on;
            line([70,70,100,100,70],[-5,5,5,-5,-5],'linewidth',2);
        end
    end
    h = colorbar;
    set(h,'position',[0.8,0.2,0.03,0.6],'ytick',-1:0.2:1,'yticklabel',{'-1.0';'-0.8';'-0.6';'-0.4';...
        '-0.2';' 0.0';' 0.2';' 0.4';' 0.6';' 0.8';' 1.0'});
    
    print(gcf,'-dpng', 'FigureS7.png');
    print(gcf,'-depsc','FigureS7.eps');    
 
end

if figS8 == 1
    lon = double(ncread('../../myDATA/Heat_content/Monthly/CZ16_1_2000m_Temp_year_1980_month_01.nc','lon'));
    lat = double(ncread('../../myDATA/Heat_content/Monthly/CZ16_1_2000m_Temp_year_1980_month_01.nc','lat'));
    
    Nlon = length(lon);
    Nlat = length(lat);
    
    [Lon,Lat] = meshgrid(lon,lat);
    Lon = Lon';
    Lat = Lat';
    
    load('cal06_ohc_CZ_H300.mat');

    %%
    load('../../Research/color/BlueWhiteOrangeRed.rgb');
    BlRe = BlueWhiteOrangeRed/255;
    
    load('min_K_v2_T150.mat');
    min_K(min_K>th_K) = nan;
    T = 1984:1:2019;
    t = 1984:1/12:2019.999;
   %%    
%     neg_y = [8,12,14,21,32];
%     pos_y = [10,13,28,33,34,35];
    
%     neg_m = min_loc([8,12,14,21,32]);
%     pos_m = min_loc([10,13,28,33,34,35]);
    %%
    neg_y = 14;
    pos_y = 10;

    %%
    pos_dy = zeros(Nlon,Nlat,12);
    pos_dy(:,:,1) = nanmean(X(:,:,12,pos_y),4);
    for i = 1:11
        pos_dy(:,:,i+1) = nanmean(X(:,:,i,pos_y+1),4);
    end
    
    neg_dy = zeros(Nlon,Nlat,12);
    neg_dy(:,:,1) = nanmean(X(:,:,12,neg_y),4);
    for i = 1:11
        neg_dy(:,:,i+1) = nanmean(X(:,:,i,neg_y+1),4);
    end
    
    %%
    load('cal06_ohc_CZ_D20.mat','X');
    pos_d20 = zeros(Nlon,Nlat,12);
    pos_d20(:,:,1) = nanmean(X(:,:,12,pos_y),4);
    for i = 1:11
        pos_d20(:,:,i+1) = nanmean(X(:,:,i,pos_y+1),4);
    end
    
    neg_d20 = zeros(Nlon,Nlat,12);
    neg_d20(:,:,1) = nanmean(X(:,:,12,neg_y),4);
    for i = 1:11
        neg_d20(:,:,i+1) = nanmean(X(:,:,i,neg_y+1),4);
    end
    
    %%
    load('cal07_10Wind.mat','U','V','lat','lon');
    pos_u = zeros(length(lon),length(lat),12);
    pos_v = zeros(length(lon),length(lat),12);
    pos_u(:,:,1) = nanmean(U(:,:,12,pos_y),4);
    pos_v(:,:,1) = nanmean(V(:,:,12,pos_y),4);
    for i = 1:11
        pos_u(:,:,i+1) = nanmean(U(:,:,i,pos_y+1),4);
        pos_v(:,:,i+1) = nanmean(V(:,:,i,pos_y+1),4);
    end
    
    neg_u = zeros(length(lon),length(lat),12);
    neg_v = zeros(length(lon),length(lat),12);
    neg_u(:,:,1) = nanmean(U(:,:,12,neg_y),4);
    neg_v(:,:,1) = nanmean(V(:,:,12,neg_y),4);
    for i = 1:11
        neg_u(:,:,i+1) = nanmean(U(:,:,i,neg_y+1),4);
        neg_v(:,:,i+1) = nanmean(V(:,:,i,neg_y+1),4);        
    end
    
    [Lon1,Lat1] = meshgrid(lon,lat);
    Lon1 = Lon1';
    Lat1 = Lat1';
    %%
    
    cmin = -20.0;
    cmax =  20.0;
    
    str = num2str(neg_y+1983);
    
    tmp1 = pos_dy;
    tmp2 = pos_d20;
    tmp3 = pos_u;
    tmp4 = pos_v;
    
    figure(1)

    
    tit = {['(a) ',num2str(pos_y+1983),'Dec(-1)' ];['(b) ',num2str(pos_y+1984),'Jan(0)'];...
           ['(c) ',num2str(pos_y+1984),'Feb(0)'];['(d) ',num2str(pos_y+1984),'Mar(0)'];...
           ['(e) ',num2str(pos_y+1984),'Apr(0)'];['(f) ',num2str(pos_y+1984),'May(0)'];...
           ['(g) ',num2str(pos_y+1984),'Jun(0)'];['(h) ',num2str(pos_y+1984),'Jul(0)'];...
           ['(i) ',num2str(pos_y+1984),'Aug(0)'];['(j) ',num2str(pos_y+1984),'Sep(0)'];...
           ['(k) ',num2str(pos_y+1984),'Oct(0)'];['(l) ',num2str(pos_y+1984),'Nov(0)']};
    for i = 1:12
        subplot(4,3,i);
        h = pcolor(Lon,Lat,tmp2(:,:,i));
        set(h,'linestyle','none');
        caxis([cmin,cmax]);
        hold on;
        colormap(BlRe);
        title(tit{i},'fontsize',10);
        axis([30,120,-20,20]);
        hold on;
        quiver(Lon1(1:4:end,1:4:end),Lat1(1:4:end,1:4:end),tmp3(1:4:end,1:4:end,i),tmp4(1:4:end,1:4:end,i),1,'color','k');
        geoshow('landareas.shp','FaceColor', [0.8 0.8 0.8],'linestyle','none');
        set(gca,'ytick',-20:10:20,'yticklabel',{'20°„S';'10°„S';'0';'10°„N';'20°„N'},...
            'xtick',30:30:120,'xticklabel',{'30°„E';'60°„E';'90°„E';'120°„E'},'fontsize',8);
    end    
    
    h = colorbar;
    set(h,'position',[0.94,0.2,0.02,0.6],'ytick',-20:5:20,'yticklabel',...
        {'-20';'-15';'-10';'-5';' 0';' 5';' 10';' 15';' 20'});
    text(-62,225,'pIOD','fontsize',15);
    
%     h = axes;
%     quiver(65,0,20,0,1,'color','k','linewidth',2);
%     set(h,'position',[0.41,0.05,0.2,0.02]);
%     axis([30,120,-20,20]);
%     axis off;
%     text(65,-35,'20m/s','fontsize',10);

    h = axes;
    quiver(65,0,20,0,1,'color','k','linewidth',2);
    set(h,'position',[0.85,0.88,0.2,0.02]);
    axis([30,120,-20,20]);
    axis off;
    text(65,-35,'20m/s','fontsize',10);
    text(72,-140,'m','fontsize',10);
    
    print(1,'-dpng', 'FigureS8_D20.png');
    print(1,'-depsc','FigureS8_D20.eps');    
    
    
end

if figS9 == 1
    lon = double(ncread('../../myDATA/Heat_content/Monthly/CZ16_1_2000m_Temp_year_1980_month_01.nc','lon'));
    lat = double(ncread('../../myDATA/Heat_content/Monthly/CZ16_1_2000m_Temp_year_1980_month_01.nc','lat'));
    
    Nlon = length(lon);
    Nlat = length(lat);
    
    [Lon,Lat] = meshgrid(lon,lat);
    Lon = Lon';
    Lat = Lat';
    
    load('cal06_ohc_CZ_H300.mat');

    %%
    load('../../Research/color/BlueWhiteOrangeRed.rgb');
    BlRe = BlueWhiteOrangeRed/255;
    
    load('min_K_v2_T150.mat');
    min_K(min_K>th_K) = nan;
    T = 1984:1:2019;
    t = 1984:1/12:2019.999;
   %%    
%     neg_y = [8,12,14,21,32];
%     pos_y = [10,13,28,33,34,35];
    
%     neg_m = min_loc([8,12,14,21,32]);
%     pos_m = min_loc([10,13,28,33,34,35]);
    %%
    neg_y = 14;
    pos_y = 35;

    %%
    pos_dy = zeros(Nlon,Nlat,12);
    pos_dy(:,:,1) = nanmean(X(:,:,12,pos_y),4);
    for i = 1:11
        pos_dy(:,:,i+1) = nanmean(X(:,:,i,pos_y+1),4);
    end
    
    neg_dy = zeros(Nlon,Nlat,12);
    neg_dy(:,:,1) = nanmean(X(:,:,12,neg_y),4);
    for i = 1:11
        neg_dy(:,:,i+1) = nanmean(X(:,:,i,neg_y+1),4);
    end
    
    %%
    load('cal06_ohc_CZ_D20.mat','X');
    pos_d20 = zeros(Nlon,Nlat,12);
    pos_d20(:,:,1) = nanmean(X(:,:,12,pos_y),4);
    for i = 1:11
        pos_d20(:,:,i+1) = nanmean(X(:,:,i,pos_y+1),4);
    end
    
    neg_d20 = zeros(Nlon,Nlat,12);
    neg_d20(:,:,1) = nanmean(X(:,:,12,neg_y),4);
    for i = 1:11
        neg_d20(:,:,i+1) = nanmean(X(:,:,i,neg_y+1),4);
    end
    
    %%
    load('cal07_10Wind.mat','U','V','lat','lon');
    pos_u = zeros(length(lon),length(lat),12);
    pos_v = zeros(length(lon),length(lat),12);
    pos_u(:,:,1) = nanmean(U(:,:,12,pos_y),4);
    pos_v(:,:,1) = nanmean(V(:,:,12,pos_y),4);
    for i = 1:11
        pos_u(:,:,i+1) = nanmean(U(:,:,i,pos_y+1),4);
        pos_v(:,:,i+1) = nanmean(V(:,:,i,pos_y+1),4);
    end
    
    neg_u = zeros(length(lon),length(lat),12);
    neg_v = zeros(length(lon),length(lat),12);
    neg_u(:,:,1) = nanmean(U(:,:,12,neg_y),4);
    neg_v(:,:,1) = nanmean(V(:,:,12,neg_y),4);
    for i = 1:11
        neg_u(:,:,i+1) = nanmean(U(:,:,i,neg_y+1),4);
        neg_v(:,:,i+1) = nanmean(V(:,:,i,neg_y+1),4);        
    end
    
    [Lon1,Lat1] = meshgrid(lon,lat);
    Lon1 = Lon1';
    Lat1 = Lat1';
    %%
    
    cmin = -20;
    cmax =  20;
    
    str = num2str(neg_y+1983);
    
    tmp1 = neg_dy;
    tmp2 = neg_d20;
    tmp3 = neg_u;
    tmp4 = neg_v;
    
    figure(1)

    
    tit = {['(a) ',num2str(neg_y+1983),'Dec(-1)' ];['(b) ',num2str(neg_y+1984),'Jan(0)'];...
           ['(c) ',num2str(neg_y+1984),'Feb(0)'];['(d) ',num2str(neg_y+1984),'Mar(0)'];...
           ['(e) ',num2str(neg_y+1984),'Apr(0)'];['(f) ',num2str(neg_y+1984),'May(0)'];...
           ['(g) ',num2str(neg_y+1984),'Jun(0)'];['(h) ',num2str(neg_y+1984),'Jul(0)'];...
           ['(i) ',num2str(neg_y+1984),'Aug(0)'];['(j) ',num2str(neg_y+1984),'Sep(0)'];...
           ['(k) ',num2str(neg_y+1984),'Oct(0)'];['(l) ',num2str(neg_y+1984),'Nov(0)']};
    for i = 1:12
        subplot(4,3,i);
        h = pcolor(Lon,Lat,tmp2(:,:,i));
        set(h,'linestyle','none');
        caxis([cmin,cmax]);
        hold on;
        colormap(BlRe);
        title(tit{i},'fontsize',10);
        axis([30,120,-20,20]);
        hold on;
        quiver(Lon1(1:4:end,1:4:end),Lat1(1:4:end,1:4:end),tmp3(1:4:end,1:4:end,i),tmp4(1:4:end,1:4:end,i),1,'color','k');
        geoshow('landareas.shp','FaceColor', [0.8 0.8 0.8],'linestyle','none');
        set(gca,'ytick',-20:10:20,'yticklabel',{'20°„S';'10°„S';'0';'10°„N';'20°„N'},...
            'xtick',30:30:120,'xticklabel',{'30°„E';'60°„E';'90°„E';'120°„E'},'fontsize',8);
    end    
    
    h = colorbar;
    set(h,'position',[0.94,0.2,0.02,0.6],'ytick',-20:5:20,'yticklabel',...
        {'-20';'-15';'-10';'-5';' 0';' 5';' 10';' 15';' 20'});    text(-62,225,'nIOD','fontsize',15);
    
%     h = axes;
%     quiver(65,0,20,0,1,'color','k','linewidth',2);
%     set(h,'position',[0.41,0.05,0.2,0.02]);
%     axis([30,120,-20,20]);
%     axis off;
%     text(65,-35,'20m/s','fontsize',10);

    h = axes;
    quiver(65,0,20,0,1,'color','k','linewidth',2);
    set(h,'position',[0.85,0.88,0.2,0.02]);
    axis([30,120,-20,20]);
    axis off;
    text(65,-35,'20m/s','fontsize',10);
    text(72,-140,'m','fontsize',10);
    
    print(1,'-dpng', 'FigureS9_D20.png');
    print(1,'-depsc','FigureS9_D20.eps');    
    
    
end

if figS10 == 1
    load('../../myDATA/IOD/SINTEX_DMI_20200210.csv');
    dmi = SINTEX_DMI_20200210(:,2);
    t = 1982:1/12:2019.999;
                
    dmi(t<1984)= [];
    t(t<1984)  = [];
    dmi = reshape(dmi,12,[]);
    tmp = bsxfun(@minus,dmi,nanmean(dmi,2));    
    dmi = tmp(:);    
                       
    th1 = 0.5;
    y1 = smooth(dmi,3);
    ind_p = find(y1>=+th1);
    ind_n = find(y1<=-th1);
    [ind_p2,n_p,ts_p,te_p] = event_series(ind_p,3);
    [ind_n2,n_n,ts_n,te_n] = event_series(ind_n,3);
       
    y = y1;
    plot(t,y,'k-','linewidth',0.5);
    hold on;
    load('BCC_DMI_1991_2018_iniDEC.mat');
    BCC = bsxfun(@minus,BCC_DMI,nanmean(BCC_DMI,2));    
    y2 = smooth(BCC(:),3);
    plot(1991:1/12:2018.999,y2,'m','linewidth',2);
    hold on;
    
    load('ECs5_DMI_1993-2016_iniMAY.mat');
    EC5 = bsxfun(@minus,EC5_DMI,nanmean(EC5_DMI,2));
    EC = zeros(12,24);
    EC(:,:) = nan;
    EC(6:11,:) = EC5;
    plot(1993:1/12:2016.999,EC(:),'g','linewidth',2);
    hold on;

    load('IOD_CFS2_anom_feb.mat');
    tmp = bsxfun(@minus,IOD_CFS2_anom,nanmean(IOD_CFS2_anom,2));
    CFS = zeros(12,39);
    CFS(:,:) = nan;
    CFS(2:11,:) = tmp;
    plot(1982:1/12:2020.999,CFS(:),'color',[1,0.8,0] ,'linewidth',2);
   
    legend('Observation','BCC\_CSM @Dec(-1)','ECs5 @May(0)','CFSv2 @Feb(0)',...
        'fontsize',10,'location','southwest','Orientation','horizontal');
    
    hold on;
    plot(1979:2020,ones(42,1)*+th1,'k:',1979:2020,-ones(42,1)*th1,'k:','linewidth',2);
    hold on;
    y = y1;
    y(setdiff(1:length(y),ind_p2)) = +th1;
    area(t,y,+th1,'facecolor','r','linestyle','none');
    hold on;
    y = y1;
    y(setdiff(1:length(y),ind_n2)) = -th1;
    area(t,y,-th1,'facecolor','b','linestyle','none');
    hold on;
    
    xlim([1984,2020]);
    ylim([-2.0,2.0]);
    ylabel('DMI (^oC)','fontsize',12);
    set(gca,'fontsize',12,'ytick',-2:0.5:2,'yticklabel',{'-2.0';'-1.5';'-1.0';'-0.5';'0.0';'0.5';'1.0';'1.5';'2.0'},'xtick',1980:1:2020,'xticklabel',...
        {'1980';'';'';'';'';'1985';'';'';'';'';'1990';'';'';'';'';'1995';'';'';'';'';'2000';'';'';'';'';...
         '2005';'';'';'';'';'2010';'';'';'';'';'2015';'';'';'';'';'2020'},'ycolor','k');    
    
    set(gca,'position',[0.08,0.2,0.89,0.5]);
    
    hold on;
    y = y1;
    y(setdiff(1:length(y),ind_p2)) = nan;
    h = bar(t,y*0-2,'histc');
    set(h,'facecolor','r','linestyle','none');
    alpha(h,0.1);
    h = bar(t,y*0+2,'histc');
    set(h,'facecolor','r','linestyle','none');
    alpha(h,0.1);
    hold on;
    y = y1;
    y(setdiff(1:length(y),ind_n2)) = nan;
    h = bar(t,y*0-2,'histc');
    set(h,'facecolor','b','linestyle','none');
    alpha(h,0.1);
    h = bar(t,y*0+2,'histc');
    set(h,'facecolor','b','linestyle','none');
    alpha(h,0.1);
      
    print(1,'-dpng','FigureS10.png');
    print(1,'-depsc','FigureS10.eps');

    
end

if figS11 == 1
    lon = double(ncread('../../myDATA/Heat_content/Monthly/CZ16_1_2000m_Temp_year_1980_month_01.nc','lon'));
    lat = double(ncread('../../myDATA/Heat_content/Monthly/CZ16_1_2000m_Temp_year_1980_month_01.nc','lat'));
    
    Nlon = length(lon);
    Nlat = length(lat);
    
    [Lon,Lat] = meshgrid(lon,lat);
    Lon = Lon';
    Lat = Lat';
    
    %%
    load('../../Research/color/BlueWhiteOrangeRed.rgb');
    BlRe = BlueWhiteOrangeRed/255;
    
    load('min_K_v2_T150.mat');
    min_K(min_K>th_K) = nan;
    T = 1984:1:2019;
    t = 1984:1/12:2019.999;
    
    ind = find(isnan(min_K)==0);
    
    neg_y = [11,23,36];
    pos_y = [15,27,30];   
%     norm_y = find(isnan(min_K)==1);
    norm_y = 1:36;
    
    %%
    load('cal06_ohc_CZ_D20.mat','X');
    pos_d20 = zeros(Nlon,Nlat,12,length(pos_y));
    neg_d20 = zeros(Nlon,Nlat,12,length(neg_y));
    norm_d20 = zeros(Nlon,Nlat,12,length(norm_y));    
    
    pos_d20(:,:,1,:) = X(:,:,12,pos_y);
    neg_d20(:,:,1,:) = X(:,:,12,neg_y);
    norm_d20(:,:,1,:) = X(:,:,12,norm_y);
    
    for i = 1:11
        pos_d20(:,:,i+1,:) = X(:,:,i,pos_y+1);
        neg_d20(:,:,i+1,:) = X(:,:,i,neg_y+1);
        norm_d20(:,:,i+1,:) = X(:,:,i,norm_y+1);        
    end    
    
    %%
    load('cal07_10Wind.mat','U','V','lat','lon');
    pos_u = zeros(length(lon),length(lat),12,length(pos_y));
    pos_v = zeros(length(lon),length(lat),12,length(pos_y));
    neg_u = zeros(length(lon),length(lat),12,length(neg_y));
    neg_v = zeros(length(lon),length(lat),12,length(neg_y)); 
    norm_u = zeros(length(lon),length(lat),12,length(norm_y));
    norm_v = zeros(length(lon),length(lat),12,length(norm_y));
    
    pos_u(:,:,1,:) = U(:,:,12,pos_y);
    pos_v(:,:,1,:) = V(:,:,12,pos_y);
    neg_u(:,:,1,:) = U(:,:,12,neg_y);
    neg_v(:,:,1,:) = V(:,:,12,neg_y);       
    norm_u(:,:,1,:) = U(:,:,12,norm_y);  
    norm_v(:,:,1,:) = V(:,:,12,norm_y);         
    
    for i = 1:11
        pos_u(:,:,i+1,:) = U(:,:,i,pos_y+1);
        pos_v(:,:,i+1,:) = V(:,:,i,pos_y+1);
        neg_u(:,:,i+1,:) = U(:,:,i,neg_y+1);
        neg_v(:,:,i+1,:) = V(:,:,i,neg_y+1);
        norm_u(:,:,i+1,:) = U(:,:,i,norm_y+1);
        norm_v(:,:,i+1,:) = V(:,:,i,norm_y+1);        
    end
    
    [Lon1,Lat1] = meshgrid(lon,lat);
    Lon1 = Lon1';
    Lat1 = Lat1';
    %%
    
    cmin = -10;
    cmax =  10;
        
    scrsz = get(0,'screensize');
    figure('position',scrsz);
    
    tit = {'(a) nIOD Dec(-1)';'(b) nIOD Feb(0)';'(c) nIOD Apr(0)';'(d) nIOD Jun(0)';'(e) nIOD Aug(0)';'(f) nIOD Oct(0)';...
           '(g) pIOD Dec(-1)';'(h) pIOD Feb(0)';'(i) pIOD Apr(0)';'(j) pIOD Jun(0)';'(k) pIOD Aug(0)';'(l) pIOD Oct(0)'};
    
    for j = 1:2
        if j == 1
            tmp1 = neg_d20;
            tmp3 = neg_u;
            tmp4 = neg_v;
            tmp5 = neg_y;
        else
            tmp1 = pos_d20;
            tmp3 = pos_u;
            tmp4 = pos_v;
            tmp5 = pos_y;            
            
        end
        
        for i = 1:6
           
            x1 = reshape(tmp1(:,:,(i-1)*2+1,:),size(tmp1,1),size(tmp1,2),[]);
            x2 = reshape(norm_d20,size(tmp1,1),size(tmp1,2),[]);
            h1 = reshape(ttest2(x1,x2,'dim',3,'alpha',0.05),size(tmp1,1),size(tmp1,2));
            tmp2 = nanmean(x1,3);
            tmp2(h1==0) = nan; 
            
%             x1 = squeeze(tmp1(:,:,(i-1)*2+1,:));
%             x2 = reshape(norm_d20(:,:,(i-1)*2+1,:),size(x1,1),size(x1,2),[]);
%             n1 = length(tmp5);
%             n2 = size(x2,3);
%             s2 = (sum(bsxfun(@minus,x1,nanmean(x1,3)).^2,3)+sum(bsxfun(@minus,x2,nanmean(x2,3)).^2,3))/(n1+n2-2);
%             test = (nanmean(x1,3)-nanmean(x2,3))./sqrt(s2)/sqrt(1/n1+1/n2);
%             tmp2 = nanmean(x1,3);
%             tmp2(abs(test)<2.0) = nan;                      
            
            subplot(4,3,i+(j-1)*6);
            h = pcolor(Lon,Lat,nanmean(x1,3)); 
            set(h,'linestyle','none');
            caxis([cmin,cmax]);
            colormap(BlRe);
            title(tit{i+(j-1)*6},'fontsize',10);
            axis([30,120,-20,20]);
            hold on;
            scatter(Lon(isnan(tmp2)==0),Lat(isnan(tmp2)==0),0.3,'g.');
            hold on;
            
            x1 = reshape(sqrt(tmp3(:,:,(i-1)*2+1,:).^2+tmp4(:,:,(i-1)*2+1,:).^2),size(tmp3,1),size(tmp3,2),[]);
            x2 = reshape(sqrt(norm_u.^2+norm_v.^2),size(tmp3,1),size(tmp3,2),[]);
            h1 = reshape(ttest2(x1,x2,'dim',3,'alpha',0.05),size(tmp3,1),size(tmp3,2));
            tmp33 = nanmean(squeeze(tmp3(:,:,(i-1)*2+1,:)),3);
            tmp333 = tmp33; 
            tmp333(h1==0) = nan;
            tmp44 = nanmean(squeeze(tmp4(:,:,(i-1)*2+1,:)),3); 
            tmp444 = tmp44;
            tmp444(h1==0) = nan;
            
%             x1 = squeeze(tmp3(:,:,(i-1)*2+1,:));
%             x2 = reshape(norm_u(:,:,(i-1)*2+1,:),size(x1,1),size(x1,2),[]);
%             n1 = length(tmp5);
%             n2 = size(x2,3);
%             s2 = (sum(bsxfun(@minus,x1,nanmean(x1,3)).^2,3)+sum(bsxfun(@minus,x2,nanmean(x2,3)).^2,3))/(n1+n2-2);
%             test = (nanmean(x1,3)-nanmean(x2,3))./sqrt(s2)/sqrt(1/n1+1/n2);
%             tmp22 = nanmean(x1,3);
%             tmp22(abs(test)<2.0) = 0;                
%             
%             x1 = squeeze(tmp4(:,:,(i-1)*2+1,:));
%             x2 = reshape(norm_v(:,:,(i-1)*2+1,:),size(x1,1),size(x1,2),[]);
%             n1 = length(tmp5);
%             n2 = size(x2,3);
%             s2 = (sum(bsxfun(@minus,x1,nanmean(x1,3)).^2,3)+sum(bsxfun(@minus,x2,nanmean(x2,3)).^2,3))/(n1+n2-2);
%             test = (nanmean(x1,3)-nanmean(x2,3))./sqrt(s2)/sqrt(1/n1+1/n2);
%             tmp33 = nanmean(x1,3);
%             tmp33(abs(test)<2.0) = 0;   
            quiver(Lon1(1:5:end,1:5:end),Lat1(1:5:end,1:5:end),tmp33(1:5:end,1:5:end),tmp44(1:5:end,1:5:end),1,'color','k');
            hold on;
            quiver(Lon1(1:5:end,1:5:end),Lat1(1:5:end,1:5:end),tmp333(1:5:end,1:5:end),tmp444(1:5:end,1:5:end),1,'color','m');
            geoshow('landareas.shp','FaceColor', [0.8 0.8 0.8],'linestyle','none');
            set(gca,'ytick',-20:10:20,'yticklabel',{'20°„S';'10°„S';'0';'10°„N';'20°„N'},...
                'xtick',30:30:120,'xticklabel',{'30°„E';'60°„E';'90°„E';'120°„E'},'fontsize',8);
        end
    end
    
    h = colorbar;
    set(h,'position',[0.94,0.2,0.02,0.6]);

    h = axes;
    quiver(65,0,20,0,1,'color','k','linewidth',2);
    set(h,'position',[0.85,0.88,0.2,0.02]);
    axis([30,120,-20,20]);
    axis off;
    text(65,-35,'20m/s','fontsize',10);
    text(72,-140,'m','fontsize',10);
    
    print(gcf,'-dpng', 'FigureS11.png');
    print(gcf,'-depsc', 'FigureS11.eps');
    
end

if figS12 == 1

    load('../../Research/color/BlueWhiteOrangeRed.rgb');
    BlRe = BlueWhiteOrangeRed/255;

    load('min_K_v2_T150.mat');
    min_K(min_K>th_K) = nan;
    T = 1984:1:2020;
    t = 1984:1/12:2020.999;
    
    neg_y = [8,12,14,21,32];
    pos_y = [10,13,28,33,34,35];
    
    neg_m = min_loc([8,12,14,21,32]);
    pos_m = min_loc([10,13,28,33,34,35]);
        
    load('cal07_10Wind.mat','U','lon','lat');
    
    [Lon1,Lat1] = meshgrid(lon,lat);
    Lon1 = Lon1';
    Lat1 = Lat1';
    
    
    cmin = -3;
    cmax =  3;
    
    figure(1)

    h = pcolor(Lon1,Lat1,U(:,:,12,37));
    set(h,'linestyle','none');
    caxis([cmin,cmax]);
    hold on;
    colormap(BlRe);
    axis([30,120,-20,20]);
    hold on;
    geoshow('landareas.shp','FaceColor', 'none');
    set(gca,'ytick',-20:10:20,'yticklabel',{'20°„S';'10°„S';'0';'10°„N';'20°„N'},...
        'xtick',30:10:120,'xticklabel',{'30°„E';'40°„E';'50°„E';'60°„E';'70°„E';'80°„E';'90°„E';'100°„E';'110°„E';'120°„E'},'fontsize',12);
    title('10m U Wind in Dec 2020','fontsize',16);
    
    line([70,70,100,100,70],[-5,5,5,-5,-5],'linewidth',2);
    h = colorbar;
    set(gca,'position',[0.1,0.3,0.7,0.4]);
    set(h,'ytick',-3:1:3,'yticklabel',{'-3';'-2';'-1';' 0';' 1';' 2';' 3'},'fontsize',12);
    text(124,23,'m/s','fontsize',12);
    
    print(gcf,'-dpng', 'FigureS12.png');
    print(gcf,'-depsc','FigureS12.eps');    

    
    
end


