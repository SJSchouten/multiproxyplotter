% Code produced for the SNF lake warming effect
% Paleo-proxy plotter for MATLAB
% Horizontal plotter preferrably used with data on age

close all
clear 

% CTRL+R to comment out sections, CTRL+T to activate commented code

% Plot parameters
width = 0.6
left = 0.3
nrplots = 10 
limits = [11000 15000]
% limits = [5000 5700]

datapath = "X:\04_PROJECTS\2022_Lake_response_warming_SNF\Scripts\Finalized_working_scripts\Masterfile_Processed_Soppensee.xlsx"
sheetnames = ["Lowres_Geochemistry","Stats", "XRF", "HSI","Carotenoids","Pollen","MOS","Gerz","NGRIP","insJ","insD","AMS","Burg","Alps","7H"]
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
height = (0.8/nrplots)
bottom = flip([0.01+height:height:0.81])

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
        sp(i) = subplot(nrplots,1,b); 
        initpos = sp(i).Position; set(gca,'position',initpos+[0.07 0 -0.15 0]) 
        % Select the plotting type
        if typ == 'b'
            barwidth = mean(diff(x))
            bar(x,y,barwidth,'FaceColor',[col],'FaceAlpha',0.4,'EdgeColor','none','ShowBaseLine','off') % Edgecolor is set to black but can be changed
%             hold on, plot(y,x,"LineWidth",1,"Color",[col],"Marker",".","MarkerEdgeColor",'k')
%             xst = x+((x(5)-x(4))/2)
%             hold on, stairs(y,xst,"LineWidth",0.5,"LineStyle",'-',"Color",0+(col*0.3))
            hold on, scatter(x,y,60,".",MarkerEdgeColor="k")
        elseif typ == 'a'
            y = fillmissing(y,"nearest")                                    % Filling data
            y2=zeros(length(y),1)+(min(y));                            % create second curve
            X=[x',fliplr(x')];                                              % create continuous x value array for plotting
            Y=[y', fliplr(y2')];                                            % create y values for out and then back
            fill(X,Y,col,'EdgeColor',0+(col*0.3),'FaceAlpha',0.3);                          % plot filled area
        elseif typ == 's'
            scatter(x,y,50,y,'filled',"LineWidth",1,"Marker","|","MarkerEdgeColor",'k')
            hold on, plot(x,y,"LineWidth",1,"Color",[col])
        elseif typ == 't'
            xst = x+((x(5)-x(4))/2)
            stairs(xst,y,"LineWidth",2,"Color",[col]), hold on,
            scatter(x,y,60,".",MarkerEdgeColor="k")
        elseif typ == 'r'
            % Define rectangle center and +/-
            % Define rectangle center and +/-d
            xst = x+((x(3)-x(2))/2)
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
                pch(z) = patch(sp(i), rectX([1 1 2 2]), rectY([1,2,2,1]),Clr(fY(z),:),'EdgeColor', 'k', 'FaceAlpha', 0.3); % FaceAlpha controls transparency
                hold on,
            end
            clear rectXarr fY rectY rectX 
        else
            plot(x,y,"LineWidth",2,"Color",[col],"LineStyle","--","Marker",".","MarkerSize",10,"MarkerEdgeColor",'k')

%             stddev = movstd(smooth(x,y,0.05,'rloess'),5)
%             % hold on, plot(sortscores(:,1),smooth(sortscores(:,1),sortscores(:,2),0.6,'rloess')
% %             stddev = smooth(x,stddev,0.9,'rloess')
% 
%             y1 = movmean(y,5)-2*stddev;
%             y2 = movmean(y,5)+2*stddev;
% 
%             hold on, patch([y1; flip(y2)],[x; flip(x)], 'k', 'FaceAlpha',0.25, 'EdgeColor','none')
%             hold on, plot(x,movmean(y,3),"LineWidth",2,"Color",[col],"Marker","none","MarkerEdgeColor",'k')
        end
        
        % Adjust general plotting settings
        set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2);
        box off;
        ylabel({sel.headers{k}{1,i+1},sel.headers{k}{2,i+1},sel.headers{k}{3,i+1}}); 
        xlim(limits); % Depths
        
        % Switch axes top bottom and with or without depth scale
        if b == 1;
        set(gca,'XminorTick','on','XGrid','off','Color','none','XAxisLocation','top','YAxisLocation','left');ylabel(xhead);
        elseif and(k==length(sel.data),i==size(sel.data{k},2)-1)
        set(gca,'XminorTick','on','XGrid','off','Color','none','YAxisLocation','left','XAxisLocation','bottom');ylabel(xhead);
        elseif b == nrplots;
        set(gca,'XminorTick','on','XGrid','off','Color','none','YAxisLocation','right','XAxisLocation','bottom');ylabel(xhead);
        elseif rem(b,2) == 1;
        set(gca,'XminorTick','on','XGrid','off','XColor','none','Color','none','YAxisLocation','left');
        else rem(b,2) == 0;
        set(gca,'XminorTick','on','XGrid','off','XColor','none','Color','none','YAxisLocation','right'); 
        end

        b = b+1
    end
end