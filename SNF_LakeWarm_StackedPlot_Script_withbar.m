close all
clear 

% Plot parameters
height = 0.7
bottom = 0.15
nrplots = 20 
limits = [13500 15000]
% limits = [5000 5700]
horizontal = 1

datapath = "X:\04_PROJECTS\2022_Lake_response_warming_SNF\Scripts\Finalized_working_scripts\Masterfile_Processed_Soppensee.xlsx"
sheetnames = ["Lowres_Geochemistry","Stats", "XRF", "HSI","Carotenoids","Pollen"]
% sheetnames = ['GERZ',"ALPS",'XRFRED',"Pollen",'RoCs',"HSIRED","Carotenoids_Wet",'Green_pigments_wet','NGRIP','InsD','InsJ']
% sheetnames = ['Lowstats',"Carotenoids_wet"]
for m = 1:length(sheetnames)
    % Data 1 read table Low resolution usually XRF
    data_raw = readcell(datapath,Sheet = sheetnames{m});                       % Loading the dataframe
    datain.headers{m} = data_raw(3:5,:);                                       % Loading the headers of the data (plot titles)
    datain.info{m} = data_raw(6:8,:);                                          % Loading the info of the data (plot type, color, in-or out criterion)
    datain.data{m} = readmatrix(datapath,Sheet = sheetnames{m});               % Loading the data
    datain.data{m} = datain.data{m}(4:size(datain.data{m},1),:);               % Cutting off the front columns containing the ages
    clear data_raw
end

for m = 1:length(datain.data)
    sel.data{m} = datain.data{m}(3:size(datain.data{m},1),cell2mat(datain.info{m}(3,:)) == 1);
    sel.headers{m} = datain.headers{m}(:,cell2mat(datain.info{m}(3,:)) == 1);
    sel.info{m} = datain.info{m}(:,cell2mat(datain.info{m}(3,:)) == 1);
end

% Determining the appropriate placings of the plots
width = (0.8/nrplots)
left = [0.01+width:width:0.81]

% Starting up the loop over the number of figures
figure,
b = 1
% Looping over the different datasets
for k = 1:length(sel.data)

    % Looping over the number of proxies within one dataset
    for i = 1:size(sel.data{k},2)-1

        % Read the data
        x = sel.data{k}(:,1) 
        y = sel.data{k}(:,i+1) 

        % Read headers
        yhead = append(sel.headers{k}(1,i+1)," ",sel.headers{k}(2,i+1)," ",sel.headers{k}(3,i+1))
        xhead = append(sel.headers{k}(2,1)," ",sel.headers{k}(3,1));

        % Read plotting info
        typ = sel.info{k}{1,i+1}
        if ismissing(sel.info{k}{2,i+1})
            col = [0.2+randi(200)/500, 0.2+randi(200)/500, 0.2+randi(200)/500]  % if color column is empty generate random color
        else
            col = str2num(sel.info{k}{2,i+1})
        end

        % Create new figure and reset counter when number of plots are full
        if b == nrplots+1
            figure,
            b = 1
        end
        
        % Plot a new plot with the right settings
        sp(i) = subplot(1,nrplots,b);set(gca,'position',[left(b) bottom width height]);
        
        % Select the plotting type
        if typ == 'b'
            bar(x,y,50,'FaceColor',[col],'FaceAlpha',0.4,'EdgeColor','none',"Horizontal","on",'ShowBaseLine','off') % Edgecolor is set to black but can be changed
%             hold on, plot(y,x,"LineWidth",1,"Color",[col],"Marker",".","MarkerEdgeColor",'k')
%             xst = x+((x(5)-x(4))/2)
%             hold on, stairs(y,xst,"LineWidth",0.5,"LineStyle",'-',"Color",0+(col*0.3))
            hold on, scatter(y,x,60,".",MarkerEdgeColor="k")
        elseif typ == 'a'
            y = fillmissing(y,"nearest")                                    % Filling data
            y2=zeros(length(y),1)+(min(y));                            % create second curve
            X=[x',fliplr(x')];                                              % create continuous x value array for plotting
            Y=[y', fliplr(y2')];                                            % create y values for out and then back
            fill(Y,X,col,'EdgeColor',0+(col*0.3),'FaceAlpha',0.3);                          % plot filled area
        elseif typ == 's'
            scatter(y,x,50,y,'filled',"LineWidth",1,"Marker","|","MarkerEdgeColor",'k')
            hold on, plot(y,x,"LineWidth",1,"Color",[col])
        elseif typ == 't'
            xst = x+((x(5)-x(4))/2)
            stairs(y,xst,"LineWidth",2,"Color",[col]), hold on,
            scatter(y,x,60,".",MarkerEdgeColor="k")
        elseif typ == 'r'
            % Define rectangle center and +/-
            % Define rectangle center and +/-d
            xst = x+((x(5)-x(4))/2)
            x = xst
            nrofswitch = sum(abs(diff(y))>0)
            fY = y(abs(diff(y))>0)
            fY(2:nrofswitch+1) = y(find(abs(diff(y))>0)+1)
            rectXarr(1) = min(x)
            rectXarr(2:nrofswitch+1) = x(abs(diff(y))>0)
            rectXarr(nrofswitch+2) = max(x)
            Clr = copper(max(y))
            for z = 1:length(rectXarr)-1
                rectX = [rectXarr(z),rectXarr(z+1)]
                rectY = xlim([sp(i)]);
                pch(z) = patch(sp(i), rectY([1 1 2 2]), rectX([1,2,2,1]),Clr(fY(z),:),'EdgeColor', 'k', 'FaceAlpha', 0.3); % FaceAlpha controls transparency
                hold on,
            end
            clear rectXarr fY rectY rectX 
        else
            plot(y,x,"LineWidth",2,"Color",[col],"LineStyle","none","Marker",".","MarkerSize",10,"MarkerEdgeColor",'k')

%             stddev = movstd(smooth(x,y,0.05,'rloess'),5)
%             % hold on, plot(sortscores(:,1),smooth(sortscores(:,1),sortscores(:,2),0.6,'rloess')
% %             stddev = smooth(x,stddev,0.9,'rloess')
% 
%             y1 = movmean(y,5)-2*stddev;
%             y2 = movmean(y,5)+2*stddev;
% 
%             hold on, patch([y1; flip(y2)],[x; flip(x)], 'k', 'FaceAlpha',0.25, 'EdgeColor','none')
            hold on, plot(movmean(y,3),x,"LineWidth",2,"Color",[col],"Marker","none","MarkerEdgeColor",'k')
        end
        
        % Adjust general plotting settings
        set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2,"YDir","reverse");
        box off;
        xlabel({sel.headers{k}{1,i+1},sel.headers{k}{2,i+1},sel.headers{k}{3,i+1}}); 
        ylim(limits); % Depths
        
        % Switch axes top bottom and with or without depth scale
        if b == 1;
        set(gca,'YminorTick','on','YGrid','on','Color','none','YAxisLocation','left','XAxisLocation','bottom');ylabel(xhead);
        elseif and(k==length(sel.data),i==size(sel.data{k},2)-1)
        set(gca,'YminorTick','on','YGrid','on','Color','none','XAxisLocation','top','YAxisLocation','right');ylabel(xhead);
        elseif b == nrplots;
        set(gca,'YminorTick','on','YGrid','on','Color','none','XAxisLocation','top','YAxisLocation','right');ylabel(xhead);
        elseif rem(b,2) == 1;
        set(gca,'YminorTick','on','YGrid','on','YColor','none','Color','none','XAxisLocation','bottom');
        else rem(b,2) == 0;
        set(gca,'YminorTick','on','YGrid','on','YColor','none','Color','none','XAxisLocation','top'); 
        end

        b = b+1
    end
end

%% Adding ROCs and PCs -->

% ROC
x = sel.data{5}(:,1)
y0 = sel.data{5}(:,2)
y1 = sel.data{5}(:,3)
y2 = sel.data{5}(:,4)
figure, subplot(1,5,3);plot(y0,x,"Color",[0.4,0.2,0.2],LineWidth=2)
hold on, patch([y1; flip(y2)], [x; flip(x)], 'k', 'FaceAlpha',0.25, 'EdgeColor','none')
hold on, scatter(y0(sel.data{5}(:,5) == 1),x(sel.data{5}(:,5) == 1),20,"green"',"filled")
set(gca, "YDir", "reverse");ylim([limits])
set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2,'YminorTick','on','YGrid','on','YColor','none','Color','none');box off

% PCS
x = datain.data{2}(3:end,4)
y0 = datain.data{2}(3:end,12)
y1 = datain.data{2}(3:end,13)
T = datain.data{2}(3:end,15)

subplot(1,5,1),scatter(y0,x,20,'x');
hold on,
y = T
nrofswitch = sum(abs(diff(y))>0)
fY = y(abs(diff(y))>0)
fY(2:nrofswitch+1) = y(find(abs(diff(y))>0)+1)
rectXarr(1) = min(x)
rectXarr(2:nrofswitch+1) = x(abs(diff(y))>0)
rectXarr(nrofswitch+2) = max(x)
Clr = bone(max(y))
for z = 1:length(rectXarr)-1
    rectX = [rectXarr(z),rectXarr(z+1)]
    rectY = xlim([subplot(1,5,1)]);
    pch(z) = patch(subplot(1,5,1), rectY([1 1 2 2]), rectX([1,2,2,1]),Clr(fY(z),:),'EdgeColor', 'none', 'FaceAlpha', 0.7); % FaceAlpha controls transparency
    hold on,
end
clear rectXarr fY rectY rectX 
hold on, plot(y0,x,'k'); title('PC1');xlabel('Depth (mm)');ylabel('Score on PC');
set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2,'YminorTick','on','YGrid','on','YColor','none','Color','none')
set(gca, "YDir", "reverse");ylim([limits])

subplot(1,5,2),scatter(y1,x,20,'x');
hold on,
y = T
nrofswitch = sum(abs(diff(y))>0)
fY = y(abs(diff(y))>0)
fY(2:nrofswitch+1) = y(find(abs(diff(y))>0)+1)
rectXarr(1) = min(x)
rectXarr(2:nrofswitch+1) = x(abs(diff(y))>0)
rectXarr(nrofswitch+2) = max(x)
Clr = bone(max(y))
for z = 1:length(rectXarr)-1
    rectX = [rectXarr(z),rectXarr(z+1)]
    rectY = xlim([subplot(1,5,2)]);
    pch(z) = patch(subplot(1,5,2), rectY([1 1 2 2]), rectX([1,2,2,1]),Clr(fY(z),:),'EdgeColor', 'none', 'FaceAlpha', 0.7); % FaceAlpha controls transparency
    hold on,
end
clear rectXarr fY rectY rectX 
hold on, plot(y1,x,'k'); title('PC2');xlabel('Depth (mm)');ylabel('Score on PC');
set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2,'YminorTick','on','YGrid','on','YColor','none','Color','none');box off
set(gca, "YDir", "reverse");ylim([limits])

datain.data{2}(3:end,4)

%% Adding stackplots -->


% Adding some stackplots
colorarrayP = flipud([223 191 159;204 153 102;153 102 51;96 64 31]/255)
colorarrayMn = flipud([204 204 255;230 204 255; 181 102 255; 119 0 230;79 0 153]/255)
colorarrayFe = (abs([255 217 179;255 179 102;255 128 0;153 77 0;77 38 0]-255))/255
x = sel.data{1}(:,1)
y = 100*sel.data{1}(:,2:5)./sum(sel.data{1}(:,2:5),2)
% Phosphorous
clear st
figure,
subplot(1,6,1),
b = bar(x,sel.data{1}(:,2:5),1,'stacked','EdgeColor',[0 0 0],"Horizontal","on",'FaceAlpha',0.4) % Edgecolor is set to black but can be changed
hold on, plot(sum(sel.data{1}(:,2:5),2),x,"LineWidth",1,"Color",[col],"Marker",".","MarkerEdgeColor",'k') 
set(gca,'YminorTick','on','YGrid','on','Color','none','XAxisLocation','top','YAxisLocation','right');ylabel(xhead);
set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2,"YDir","reverse");
for h = 1:4
    st(h) = append(string(sel.headers{1}{1,1+h})," ",string(sel.headers{1}{2,1+h})," ",string(sel.headers{1}{3,1+h}))
    b(h).FaceColor = colorarrayP(h,:);
end
legend(st)
xlabel("Phosphorous {\mu}g/g")
ylim(limits); % Depths
box off; hold off


x = sel.data{1}(:,1)
y = 100*sel.data{1}(:,2:5)./sum(sel.data{1}(:,2:5),2)
% Phosphorous
clear st
subplot(1,6,2),
b = bar(x,y,1,'stacked','EdgeColor',[0 0 0],"Horizontal","on",'FaceAlpha',0.4) % Edgecolor is set to black but can be changed
set(gca,'YminorTick','on','YGrid','on','Color','none','XAxisLocation','top','YAxisLocation','right');ylabel(xhead);
set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2,"YDir","reverse");
for h = 1:4
    st(h) = append(string(sel.headers{1}{1,1+h})," ",string(sel.headers{1}{2,1+h})," ",string(sel.headers{1}{3,1+h}))
    b(h).FaceColor = colorarrayP(h,:);
end
legend(st)
xlabel("Phosphorous {\mu}g/g")
ylim(limits); % Depths
box off; hold off


% Manganese of XRF

% subplot(1,6,2),
% plot(sel.data{1,1}(:,6),sel.data{1,2}(:,1),"LineWidth",0.8,"LineStyle",":","Color",'r',"Marker",".","MarkerEdgeColor",'k'); hold on,
% plot(movmean(sel.data{1,2}(:,6),5),sel.data{1,2}(:,1),"LineWidth",2,"Color",'b',"Marker",".","MarkerEdgeColor",'k') 
% set(gca,'YminorTick','on','YGrid','on','Color','none','XAxisLocation','top','YAxisLocation','right');ylabel(xhead);
% set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2,"YDir","reverse");
% xlabel("Phosphorous xrf normalized cps")
% ylim(limits); % Depths
% box off;

x = sel.data{1}(:,1)
y = 100*sel.data{1}(:,6:10)./sum(sel.data{1}(:,6:10),2)
% Manganese
subplot(1,6,3),
b2 = bar(x,sel.data{1}(:,6:10),1,'stacked','EdgeColor',[0 0 0],"Horizontal","on",'FaceAlpha',0.4) % Edgecolor is set to black but can be changed
hold on, plot(sum(sel.data{1}(:,6:10),2),x,"LineWidth",1,"Color",[col],"Marker",".","MarkerEdgeColor",'k') 
set(gca,'YminorTick','on','YGrid','on','Color','none','XAxisLocation','top','YAxisLocation','right');ylabel(xhead);
set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2,"YDir","reverse");
for h = 1:5
    st(h) = append(string(sel.headers{1}{1,5+h})," ",string(sel.headers{1}{2,5+h})," ",string(sel.headers{1}{3,5+h}))
    b2(h).FaceColor = colorarrayMn(h,:);
end
legend(st)
xlabel("Manganese {\mu}g/g")
ylim(limits); % Depths
box off;

x = sel.data{1}(:,1)
y = 100*sel.data{1}(:,6:10)./sum(sel.data{1}(:,6:10),2)
% Manganese
subplot(1,6,4),
b2 = bar(x,y,1,'stacked','EdgeColor',[0 0 0],"Horizontal","on",'FaceAlpha',0.4) % Edgecolor is set to black but can be changed
set(gca,'YminorTick','on','YGrid','on','Color','none','XAxisLocation','top','YAxisLocation','right');ylabel(xhead);
set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2,"YDir","reverse");
for h = 1:4
    st(h) = append(string(sel.headers{1}{1,5+h})," ",string(sel.headers{1}{2,5+h})," ",string(sel.headers{1}{3,5+h}))
    b2(h).FaceColor = colorarrayMn(h,:);
end
legend(st)
xlabel("Manganese {\mu}g/g")
ylim(limits); % Depths
box off;

% Manganese of XRF
% subplot(1,6,4),
% plot(sel.data{1,2}(:,5),sel.data{1,2}(:,1),"LineWidth",0.8,"LineStyle",":","Color",'r',"Marker",".","MarkerEdgeColor",'k'); hold on,
% plot(movmean(sel.data{1,2}(:,5),5),sel.data{1,2}(:,1),"LineWidth",2,"Color",'b',"Marker",".","MarkerEdgeColor",'k') 
% set(gca,'YminorTick','on','YGrid','on','Color','none','XAxisLocation','top','YAxisLocation','right');ylabel(xhead);
% set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2,"YDir","reverse");
% xlabel("Manganese xrf normalized cps")
% ylim(limits); % Depths
% box off;

% Iron
subplot(1,6,5),
b3 = bar(x,sel.data{1}(:,11:15),1,'stacked','EdgeColor',[0 0 0],"Horizontal","on",'FaceAlpha',0.4) % Edgecolor is set to black but can be changed
hold on, plot(sum(sel.data{1}(:,11:15),2),x,"LineWidth",1,"Color",[col],"Marker",".","MarkerEdgeColor",'k') 
set(gca,'YminorTick','on','YGrid','on','Color','none','XAxisLocation','top','YAxisLocation','right');ylabel(xhead);
set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2,"YDir","reverse");
for h = 1:5
    st(h) = append(string(sel.headers{1}{1,10+h})," ",string(sel.headers{1}{2,10+h})," ",string(sel.headers{1}{3,10+h}))
    b3(h).FaceColor = colorarrayFe(h,:);
end
legend(st)
xlabel("Iron {\mu}g/g")
ylim(limits); % Depths
box off;

x = sel.data{1}(:,1)
y = 100*sel.data{1}(:,11:15)./sum(sel.data{1}(:,11:15),2)
% Iron
subplot(1,6,6),
b3 = bar(x,y,1,'stacked','EdgeColor',[0 0 0],"Horizontal","on",'FaceAlpha',0.4) % Edgecolor is set to black but can be changed
set(gca,'YminorTick','on','YGrid','on','Color','none','XAxisLocation','top','YAxisLocation','right');ylabel(xhead);
set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2,"YDir","reverse");
for h = 1:5
    st(h) = append(string(sel.headers{1}{1,10+h})," ",string(sel.headers{1}{2,10+h})," ",string(sel.headers{1}{3,10+h}))
    b3(h).FaceColor = colorarrayFe(h,:);
end
legend(st)
xlabel("Iron {\mu}g/g")
ylim(limits); % Depths
yticks([5000:100:5700])
yticklabels({'11300 +/- 300','11950 +/- 150','13050 +/- 100','13870 +/- 220','14810 +/-300','15780 +/- 800','16530 +/- 1000'})
box off;

% % Iron of xrf
% subplot(1,6,6),
% plot(sel.data{1,2}(:,4),sel.data{1,2}(:,1),"LineStyle","none","Color",'r',"Marker",".","MarkerEdgeColor",'k'); hold on,
% plot(movmean(sel.data{1,2}(:,4),5),sel.data{1,2}(:,1),"LineWidth",2,"Color",'b',"Marker",".","MarkerEdgeColor",'k') 
% set(gca,'YminorTick','on','YGrid','on','Color','none','XAxisLocation','top','YAxisLocation','right');ylabel(xhead);
% set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2,"YDir","reverse");
% xlabel("Iron xrf normalized cps")
% ylim(limits); % Depths
% box off;

