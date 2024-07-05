%% --- Statistics based on excel sheet of collected data ---

% Script which calculates basic statistical plots for selected proxies
% SNF Lake Warming project - S.J.Schouten - Initiated: 06-09-2023 
% Last updated: 
% 15-09-2023 By: Stan - claycut function
% 08-01-2024 By: Stan - Beautiful clustering PCA
% 08-02-2024 By: Stan - MDS added

% Still to improve:
% - Spacing between header lines, i.e. the titles of the plot
% - Adding a box indicating the clusters to the plot
% - Setting line of threshold for pigments and Fe-Mn in the bar plots

% Default settings:
% height = 0.4
% bottom = 0.3
% nrplots = 20 

close all
clear 

height = 0.7;
bottom = 0.15;
nrplots = 20; 
resolution = 10;    %Resampling resolution
claycut = 0;
% claycut = [80,195]

limits = [4907 5797];
datapath = "X:\04_PROJECTS\2022_Lake_response_warming_SNF\Scripts\Finalized_working_scripts\Masterfile_Processed_Soppensee.xlsx";
cd X:\04_PROJECTS\2022_Lake_response_warming_SNF\Scripts\Finalized_working_scripts;
sheetnames = ["HSI","Lowres_Pigments10"];
depthheaders = ["Depthhsi","Depthlow"];

% Data table building
for m = 1:length(sheetnames)
    data_raw = readcell(datapath,Sheet = sheetnames{m});
    datainP.headers{m} = data_raw(3:5,:);
    datainP.info{m} = data_raw(6:8,:);
    datainP.data{m} = readmatrix(datapath,Sheet = sheetnames{m}); datainP.data{m} = datainP.data{m}(4:size(datainP.data{m},1),:);
    clear data_raw;
end

for m = 1:length(datainP.data)
    sel.data{m} = datainP.data{m}(3:size(datainP.data{m},1),datainP.data{m}(2,:) == 1);
    sel.headers{m} = datainP.headers{m}(:,datainP.data{m}(2,:) == 1);
    sel.info{m} = datainP.info{m}(:,datainP.data{m}(2,:) == 1);
    datain{m} = array2table(sel.data{m},"VariableNames",sel.headers{m}(1,:));
end

subheaders_in = sel.headers{1}(2,2:size(sel.headers{1},2));
units_in = sel.headers{1}(3,2:size(sel.headers{1},2));
types_in = sel.info{1}(1,2:size(sel.info{1},2));   
colors_in = sel.info{1}(2,2:size(sel.info{1},2));

for h = 2:length(sheetnames)
    if length(sheetnames) >= h-1
    subheaders_in = [subheaders_in, sel.headers{h}(2,2:size(sel.headers{h},2))];
    units_in = [units_in, sel.headers{h}(3,2:size(sel.headers{h},2))];
    types_in = [types_in, sel.info{h}(1,2:size(sel.info{h},2))];
    colors_in = [colors_in, sel.info{h}(2,2:size(sel.info{h},2))];
    end
end

[yt,ytab] = ResamplingMULTIPROX(datain,depthheaders,resolution);
% writetable(ytab,"G:\Research Projects\2022_Soppensee\BA_data\Output\Soppensee_data_resampled.csv");

%% Stacked plot of subsampled proxies

% Find where the resmapling is giving NaN (This is an artefact of the
% linear interpolation method
cutoff = [find(sum(isnan(yt),2)==0,1),find(sum(isnan(yt),2)==0,1,'last')];

% Calculate the number of tiles
nroftiles = size(yt,2);
x = yt(cutoff(1):cutoff(2),1);
yt = zscore(yt(cutoff(1):cutoff(2),2:size(yt,2)))+1;

% Optional clay cutoff
if claycut > 1
    yt = yt(claycut(1):claycut(2),:)
    x = x(claycut(1):claycut(2),:)
end
xlab = x
figure,
b = 1;
k = 1;

for i = 1:size(yt,2)

    if b == nrplots+1
        figure,
        b = 1;
        k = k+1;
    end

    y = yt(:,i);
    yhead = append(ytab.Properties.VariableNames(i+1)," ",subheaders_in(i)," ",units_in(i));
    xhead = append(sel.headers{1}(2)," ",sel.headers{1}(3));
    typ = types_in{i};

    if ismissing(colors_in{i})
        col = [0.2+randi(200)/500, 0.2+randi(200)/500, 0.2+randi(200)/500];  % if color column is empty generate random color
    else
        col = str2num(colors_in{i});
    end
    
    width = (0.8/nrplots);
    left = 0.01+width:width:0.81;
    subplot(1,nrplots,b);set(gca,'position',[left(b) bottom width height]);
    
    % Select the plotting type
    if typ == 'b'
        bar(x,y,1,'FaceColor',col,'EdgeColor',[0 0 0],"Horizontal","on"); % Edgecolor is set to black but can be changed
    elseif typ == 'a'
        y = fillmissing(y,"nearest");                                   % Filling data
        y2=zeros(length(y),1)+(min(y)*0.99);                            % create second curve
        X=[x',fliplr(x')];                                              % create continuous x value array for plotting
        Y=[y', fliplr(y2')];                                            % create y values for out and then back
        fill(Y,X,col,'EdgeColor',0+(col*0.3));                          % plot filled area
    elseif typ == 's'
        scatter(y,x,"LineWidth",2,"Color",col,"Marker",".","MarkerEdgeColor",'k');
        hold on, plot(y,x,"LineWidth",1,"Color",col,"Marker","x","MarkerEdgeColor",'k');
    elseif typ == 't'
        stairs(y,x,"LineWidth",2,"Color",col);
    else
        plot(y,x,"LineWidth",2,"Color",col,"Marker",".","MarkerEdgeColor",'k');
    end

    % Adjust general plotting settings
    set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2,"YDir","reverse");
    box off;
    xlabel({ytab.Properties.VariableNames{i+1},subheaders_in{i},units_in{i}}); 
    ylim(limits); % Depths
    
    % Switch axes top bottom and with or without depth scale
    if b == 1
    set(gca,'YminorTick','on','YGrid','on','Color','none','YAxisLocation','left','XAxisLocation','bottom');ylabel(xhead);
    elseif and(k==length(sheetnames),i==size(yt,2))
    set(gca,'YminorTick','on','YGrid','on','Color','none','XAxisLocation','top','YAxisLocation','right');ylabel(xhead);
    elseif b == nrplots
    set(gca,'YminorTick','on','YGrid','on','Color','none','XAxisLocation','top','YAxisLocation','right');ylabel(xhead);
    elseif rem(b,2) == 1
    set(gca,'YminorTick','on','YGrid','on','YColor','none','Color','none','XAxisLocation','bottom');
    else 
    set(gca,'YminorTick','on','YGrid','on','YColor','none','Color','none','XAxisLocation','top'); 
    end

    b = b+1;
end

%% PCA

figure('Units','normalized','Position',[0.3 0.3 0.3 0.5]);

% PCA biplot
[coefs,scores,latents] = pca(yt);

variables = string([ytab.Properties.VariableNames(2:size(ytab,2))]);
ax1 = subplot(2,2,1); % Top subplot
biplot(ax1,coefs(:,1:2),'Scores',scores(:,1:2),'VarLabels',variables);
xlabel('Component 1');
ylabel('Component 2');
title('PCs 1-2');

subplot(2,2,2); % Top subplot
scatter(scores(:,1),scores(:,2),"filled"); 
sortscores = sortrows(scores,1)
% std = movstd(sortscores(:,2),5)
stddev = movstd(smooth(sortscores(:,1),sortscores(:,2),0.6,'rloess'),25)
hold on, plot(sortscores(:,1),smooth(sortscores(:,1),sortscores(:,2),0.6,'rloess'))
stddev = smooth(sortscores(:,1),stddev,0.6,'rloess')

y1 = smooth(sortscores(:,1),sortscores(:,2),0.6,'rloess')-stddev;
y2 = smooth(sortscores(:,1),sortscores(:,2),0.6,'rloess')+stddev;

hold on, patch([sortscores(:,1); flip(sortscores(:,1))], [y1; flip(y2)], 'k', 'FaceAlpha',0.25, 'EdgeColor','none')
xlabel('Component 1');
ylabel('Component 2');
title('PCs 1-2 scores');

ax2 = subplot(2,2,3); % Bottom subplot
biplot(ax2,coefs(:,2:3),'Scores',scores(:,2:3),'VarLabels',variables);
% xlim(ax2,[-1 1])
% ylim(ax2,[-1 1])
xlabel('Component 2');
ylabel('Component 3');
title('PCs 2-3');

% Screeplot
subplot(2,2,4), stem(cumsum(latents/sum(latents)*100));
xlabel('PC number');
ylabel('Variance explained (%)');
title('Scree-plot');

%% Walking trough PCA

% figure,
% h = animatedline("Marker",".");
% for m = 1:length(scores(:,1))
%     addpoints(h,scores(m,2),scores(m,3));
%     title(string(yt(3+m,1)))
%     drawnow;pause(0.1)
% end

%% Correlation matrix

figure,
r = corr(yt);

isupper = logical(triu(ones(size(r)),1));
r(isupper) = NaN;

h = heatmap(r,'MissingDataColor','w','Colormap',parula);
labels = ytab(:,2:size(ytab,2)).Properties.VariableNames;
h.XDisplayLabels = labels;
h.YDisplayLabels = labels;

title('Correlation matrix');xlabel('Proxies'); ylabel('Proxies');

% %% Linear regressions
% 
% ytabred = ytab(cutoff(1):cutoff(2),:);
% fitnames = {'TOC','Crinc.Crcoh / Br'};
% comp1 = zscore(ytabred.("TOC"));
% comp2 = zscore(ytabred.("Crinc.Crcoh / Br")); % Crinc.Crcoh / Br
% mdl = fitlm(comp1, comp2);
% 
% figure, subplot(3,2,1),plot(table2array(ytabred(:,1)),comp1,"Marker",".","MarkerEdgeColor",'k');
% hold on, plot(table2array(ytabred(:,1)),comp2,"Marker",".","MarkerEdgeColor",'k');title("Fitted signals");
% legend(fitnames);
% subplot(3,2,2),plot(mdl); title("Model fit"); xlabel(fitnames{1}); ylabel(fitnames{2});
% subplot(3,2,3),histfit(mdl.Residuals.Standardized); title("Residual distribution");
% subplot(3,2,4),bar(mdl.Diagnostics.CooksDistance); title("Cooks distance");
% subplot(3,2,5),plotResiduals(mdl,'probability'); title("Residual distribution Q-Qplot");

%% Unconstrained clustering (hierarchial) wards method
figure,subplot(2,1,1)
cut = 17;
dendrogram(linkage(yt,'ward'),0,'ColorThreshold',cut);set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2)
T = cluster(linkage(yt,'ward'),'cutoff',cut,'Criterion','distance');
title(append('Dendrogram wards clusters cutoff ',num2str(cut))); xlabel('Point nr (large = downcore - small = upcore)'); ylabel('Eucledian distance');

subplot(2,1,2);
scatter(scores(:,1),scores(:,2),40,T,"filled");
text(scores(:,1),scores(:,2),string([1:1:length(scores(:,1))]))
xlabel('PC 1');
ylabel('PC 2');
title('PCs 1-2 wards clusters');set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2)

figure, subplot(4,1,1),scatter(x,scores(:,1),20,T,"filled");
hold on, plot(x,scores(:,1),'k'); title('PC1');xlabel('Depth (mm)');ylabel('Score on PC');
set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2)
subplot(4,1,2),scatter(x,scores(:,2),20,T,"filled");
hold on, plot(x,scores(:,2),'k'); title('PC2');xlabel('Depth (mm)');ylabel('Score on PC');
set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2)
subplot(4,1,3),scatter(x,scores(:,3),20,T,"filled");
hold on, plot(x,scores(:,3),'k'); title('PC3');xlabel('Depth (mm)');ylabel('Score on PC');
set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2)
subplot(4,1,4),scatter(x,scores(:,4),20,T,"filled");
hold on, plot(x,scores(:,4),'k'); title('PC4');xlabel('Depth (mm)');ylabel('Score on PC');
set(gca,"FontName","Gill Sans MT","FontSize",12,"LineWidth",2)

pcscoretable(:,1) = array2table(x)
pcscoretable(:,2) = array2table(scores(:,1))
pcscoretable(:,3) = array2table(scores(:,2))
pcscoretable(:,4) = array2table(scores(:,3))
pcscoretable(:,5) = array2table(scores(:,4))
pcscoretable(:,6) = array2table(T)



% -- Experimental stuff

% X = yt
% figure;
% plot(X(:,1),X(:,2),'.');
% title 'Data';
% opts = statset('Display','final');
% [idx,C,sumd,D] = kmeans(X,4,'Distance','sqeuclidean',...
%     'Replicates',8,'Options',opts);
% 
% 
% [coefs,scores,latents] = pca(yt);
% 
% figure;
% scatter(scores(:,1),scores(:,2),20,idx,'o','filled')
% hold on
% % plot(scores(idx==2,1),scores(idx==2,2),'.','MarkerSize',12)
% % hold on,
% % plot(scores(idx==3,2),scores(idx==3,3),'.','MarkerSize',12)
% % hold on,
% % plot(scores(idx==4,3),scores(idx==4,4),'.','MarkerSize',12)
% % hold on,
% % plot(scores(idx==5,4),scores(idx==5,5),'.','MarkerSize',12)
% % hold on,
% Cn(:,1) = [mean(scores(idx==1,1)),mean(scores(idx==1,2))]
% Cn(:,2) = [mean(scores(idx==2,1)),mean(scores(idx==2,2))]
% Cn(:,3) = [mean(scores(idx==3,1)),mean(scores(idx==3,2))]
% Cn(:,4) = [mean(scores(idx==4,1)),mean(scores(idx==4,2))]
% Cn(:,5) = [mean(scores(idx==5,1)),mean(scores(idx==5,2))]
% plot(Cn(1,:), Cn(2,:),'kx',...
%      'MarkerSize',15,'LineWidth',3) 
% legend('Cluster 1','Cluster 2','Centroids',...
%        'Location','NW')
% title 'Cluster Assignments and Centroids'
% 
% mu1 = [2 2];          % Mean of the 1st component
% sigma1 = [2 0; 0 1];  % Covariance of the 1st component
% mu2 = [-2 -1];        % Mean of the 2nd component
% sigma2 = [1 0; 0 1];  % Covariance of the 2nd component
% 
% rng('default') % For reproducibility
% r1 = scores(idx==1,1:2)
% r2 = scores(idx==2,1:2)
% X = [scores(:,1),scores(:,2)];
% 
% gm = fitgmdist(X,5);
% 
% figure,
% scatter(X(:,1),X(:,2),20,idx,'o','filled') % Scatter plot with points of size 10
% hold on
% gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(gm,[x0 y0]),x,y);
% fcontour(gmPDF,[-8 8 -4 6])
% hold on
% scatter(coefs(:,1),coefs(:,2))
% 
% close all

% --- Nice Biplots ---

% [coeff,score,latent,tsquared,explained] = pca(X);
% f1=figure(); 
% h=biplot(coeff(:,1:2),'scores',score(:,1:2),'color','k','marker','.','markersize',17,'varlabels',...
%     [ytab.Properties.VariableNames(2:size(ytab,2))],'ObsLabels',x);
% hold on
% %color by site
% hID = get(h,'tag'); %identify handle
% hPt = h(strcmp(hID,'obsmarker')); %isolate handles to scatter points
% [grp, grpID] = findgroups(idx);
% clrMap = winter(5);
% p = 95; % CI level
% for i = 1:max(grpID)
%     set(hPt(grp==i), 'Color', clrMap(i,:), 'DisplayName', sprintf('MSP%d', grpID(i)))
%     
%     % Compute centers (means)
%     allX = arrayfun(@(hh)hh.XData(1), hPt(grp==i)); 
%     allY = arrayfun(@(hh)hh.YData(1), hPt(grp==i)); 
%     centers(1) = mean(allX);  % x mean
%     centers(2) = mean(allY);  % y mean
%     
%     % Plot centers, do they make sense?
%     plot(centers(1), centers(2), 'rp', 'MarkerFaceColor', clrMap(i,:), 'MarkerSize', 20, 'LineWidth', 1)
%     
%     % Compute 95% CI using percentile method
%     CIx = prctile(allX, [(100-p)/2, p+(100-p)/2]); % x CI [left, right]
%     CIy = prctile(allY, [(100-p)/2, p+(100-p)/2]); % y CI [lower, upper]
%     CIrng(1) = CIx(2)-CIx(1); % CI range (x)
%     CIrng(2) = CIy(2)-CIy(1); % CI range (y)
%     
%     % Draw ellipses
%     llc = [CIx(1), CIy(1)]; % (x,y) lower left corners
%     rectangle('Position',[llc,CIrng],'Curvature',[1,1], 'EdgeColor', clrMap(i,:));
%     
% end
% title('Full Deployment (12h)');
% set(gca,'fontsize',18)
% xlabel('Component 1 (55.26%)') % how can I add percent ('%s %', explained(1))')
% ylabel('Component 2 (21.87%)')
% [~, unqIdx] = unique(grp);
% legend(hPt(unqIdx))
% 
% (scores(:,2)-mean(scores(:,2)))/std(scores(:,2))^2
% 
% % Creating data
% %close all

% T = table2array(temp)
clear x y

[coeff,scores,latent,tsquared,explained] = pca(yt);

f1=figure(); 
h=biplot(coeff(:,1:2),'scores',scores(:,1:2),'color',[0.7 0.7 0.7],'marker','.','markersize',17,'varlabels',...
    [ytab.Properties.VariableNames(2:size(ytab,2))],'ObsLabels',num2str(xlab));
title('Gechemical PCA')
hold on,
hID = get(h,'tag'); %identify handle
hVar = h(strcmp(hID,'varmarker'))
set(hVar,'Marker','none','MarkerSize',5)
hPt = h(strcmp(hID,'obsmarker')); %isolate handles to scatter points
[grp, grpID] = findgroups(T);
clrMap = hsv(10);
for i=1:max(T)
    x0=arrayfun(@(hh)hh.XData(1), hPt(grp==i));
    y0=arrayfun(@(hh)hh.YData(1), hPt(grp==i));
    centers(1) = mean(x0);  % x mean
    centers(2) = mean(y0);  % y mean
    % Plot centers, do they make sense?
    plot(centers(1), centers(2),'x','MarkerEdgeColor', clrMap(i,:),'MarkerSize',15,'LineWidth',3)
%     text(centers(1), centers(2),num2str(mean(xlab(grp == i))))
    orientation_rad =0;
    cos_phi = cos( orientation_rad );
    sin_phi = sin( orientation_rad );
    R = [ cos_phi sin_phi; -sin_phi cos_phi ];
    xy=[x0,y0]*R;
    x{i}=xy(:,1);
    y{i}=xy(:,2);
end

arr = [arrayfun(@(hh)hh.XData(1),hPt), arrayfun(@(hh)hh.YData(1),hPt)]
sortscores = sortrows(arr,1)
% std = movstd(sortscores(:,2),5)
stddev = movstd(smooth(sortscores(:,1),sortscores(:,2),0.6,'rloess'),50)
% hold on, plot(sortscores(:,1),smooth(sortscores(:,1),sortscores(:,2),0.6,'rloess')
stddev = smooth(sortscores(:,1),stddev,0.6,'rloess')

y1 = smooth(sortscores(:,1),sortscores(:,2),0.6,'rloess')-stddev;
y2 = smooth(sortscores(:,1),sortscores(:,2),0.6,'rloess')+stddev;

hold on, patch([sortscores(:,1); flip(sortscores(:,1))], [y1; flip(y2)], 'k', 'FaceAlpha',0.25, 'EdgeColor','none')

% Plot the data with ellipse
for i=1:length(x)
set(hPt(grp==i), 'Color', clrMap(i,:), 'DisplayName', append('Cluster ',num2str(grpID(i))))
    hold on
    ym = y{i}
    xm = x{i}
        % Center of the ellipse
        X0=mean(x{i});
        Y0=mean(y{i});
        % Lenght of the axes of the ellipse
        
        X = [x{i}-X0,y{i}-Y0]
        D=bsxfun(@minus, X, mean(X, 1));
        CVD=D'*D/size(X,1);
        [~,I]=sort(eig(CVD),'descend');
        [H,V]=eig(CVD);
        l=H(:,I);
    
        orientation_rad =-atan(l(2,1)/l(1,1));
        cos_phi = cos( orientation_rad );
        sin_phi = sin( orientation_rad );
        R = [cos_phi sin_phi; -sin_phi cos_phi ];
        xy=X*R;
        a=2*std(xy(:,1));
        b=2*std(xy(:,2));
        % Matrix of rotation
        R = [ cos_phi -sin_phi; sin_phi cos_phi ];
        % Drawing the ellipse
        theta_r         = linspace(0,2*pi);
        ellipse_x_r     =  a*cos( theta_r );
        ellipse_y_r     = b*sin( theta_r );
        rotated_ellipse =  [ellipse_x_r;ellipse_y_r]'*R;
        rotated_ellipse(:,1)=rotated_ellipse(:,1) + X0;
        rotated_ellipse(:,2)=rotated_ellipse(:,2) + Y0;
        handle_ellipse =plot( rotated_ellipse(:,1),rotated_ellipse(:,2),'r' );
        hold on,
        handle_elipseface = fill(rotated_ellipse(:,1),rotated_ellipse(:,2),'r','FaceAlpha',0.3)
    handle_ellipse.Color=clrMap(i,:);
    handle_elipseface.FaceColor=clrMap(i,:);
    [~, unqIdx] = unique(grp);
    legend(hPt(unqIdx))
end
% hold on, patch([sortscores(:,1); flip(sortscores(:,1))], [y1; flip(y2)], 'k', 'FaceAlpha',0.25, 'EdgeColor','none')
xlabel(append('Component 1 (~',num2str(round(explained(1)),4),' % of variance)'));ylabel(append('Component 2 (~',num2str(round(explained(2)),4),' % of variance)'))


%% Multidimensional scaling (MDS) Example metric / classical
% Load matlab cities data
ratings = zscore(yt)

% Step 1: Set up our proximity matrix
% First let's create our similarity (proximity) matrix by calculating the
proximities = zeros(size(ratings,1));
for i=1:size(ratings,1)
    for j =1:size(ratings,1)
        proximities(i,j) = pdist2(ratings(i,:),ratings(j,:),'euclidean');
    end
end
[Xt,eigv] = cmdscale(proximities)
figure, plot(cumsum(eigv/sum(eigv)*100),'o')
% We can now look at our cities in 2D:
figure,scatter(Xt(:,1),Xt(:,2),20,T)
title('Example of Classical MDS with M=2 for City Ratings');

% Throw on some labels!
text(Xt(:,1),Xt(:,2), num2str(table2array(ytab(cutoff(1):cutoff(2),1))), 'VerticalAlignment','bottom', ...
                             'HorizontalAlignment','right')

%% Non-metric MDS

figure,

% Most important step choose the right distance measure for your data-type:
dissimilarities = pdist(ratings,"seuclidean")

% Running the ordination
[Y,stress,disparities] = mdscale(dissimilarities,2,'criterion','stress');

% Make a Shepards plot of the results, which displays the goodness of fit
distances = pdist(Y);
[dum,ord] = sortrows([disparities(:) dissimilarities(:)]);
plot(dissimilarities,distances,'bo', dissimilarities(ord),disparities(ord),'r.-', [0 15],[0 15],'k-')
xlabel('Dissimilarities')
ylabel('Distances/Disparities')
legend({'Distances' 'Disparities' '1:1 Line'},'Location','NorthWest');

% plot variables on in the 2 dimensional space
figure, h = scatter(Y(:,1),Y(:,2),30,T,"filled")
% hold on, h = animatedline("Marker",".");
% for m = 1:length(Y(:,1))
%     addpoints(h,Y(m,1),Y(m,2));
%     drawnow;pause(0.1)
% end
text(Y(:,1),Y(:,2), num2str(table2array(ytab(cutoff(1):cutoff(2),1))), 'VerticalAlignment','bottom','HorizontalAlignment','right')

hID = get(h,'tag'); %identify handle
hVar = h(strcmp(hID,'varmarker'))
set(hVar,'Marker','none','MarkerSize',5)
hPt = h(strcmp(hID,'obsmarker')); %isolate handles to scatter points
[grp, grpID] = findgroups(T);
clrMap = copper(12);
for i=1:max(T)
    x0=Y(grp==i,1);
    y0=Y(grp==i,2);
    centers(1) = mean(x0);  % x mean
    centers(2) = mean(y0);  % y mean
    % Plot centers, do they make sense?
    plot(centers(1), centers(2),'x','MarkerEdgeColor', clrMap(i,:),'MarkerSize',15,'LineWidth',3)
%     text(centers(1), centers(2),num2str(mean(xlab(grp == i))))
    orientation_rad =0;
    cos_phi = cos( orientation_rad );
    sin_phi = sin( orientation_rad );
    R = [ cos_phi sin_phi; -sin_phi cos_phi ];
    xy=[x0,y0]*R;
    x{i}=xy(:,1);
    y{i}=xy(:,2);
    hold on,
end

for i=1:length(x)
    hold on
    ym = y{i}
    xm = x{i}
        % Center of the ellipse
        X0=mean(x{i});
        Y0=mean(y{i});
        % Lenght of the axes of the ellipse
        
        X = [x{i}-X0,y{i}-Y0]
        D=bsxfun(@minus, X, mean(X, 1));
        CVD=D'*D/size(X,1);
        [~,I]=sort(eig(CVD),'descend');
        [H,V]=eig(CVD);
        l=H(:,I);
    
        orientation_rad =-atan(l(2,1)/l(1,1));
        cos_phi = cos( orientation_rad );
        sin_phi = sin( orientation_rad );
        R = [cos_phi sin_phi; -sin_phi cos_phi ];
        xy=X*R;
        a=2*std(xy(:,1));
        b=2*std(xy(:,2));
        % Matrix of rotation
        R = [ cos_phi -sin_phi; sin_phi cos_phi ];
        % Drawing the ellipse
        theta_r         = linspace(0,2*pi);
        ellipse_x_r     =  a*cos( theta_r );
        ellipse_y_r     = b*sin( theta_r );
        rotated_ellipse =  [ellipse_x_r;ellipse_y_r]'*R;
        rotated_ellipse(:,1)=rotated_ellipse(:,1) + X0;
        rotated_ellipse(:,2)=rotated_ellipse(:,2) + Y0;
        handle_ellipse =plot( rotated_ellipse(:,1),rotated_ellipse(:,2),'r' );
        hold on,
        handle_elipseface = fill(rotated_ellipse(:,1),rotated_ellipse(:,2),'r','FaceAlpha',0.3)
end
hold on, scatter(Y(:,1),Y(:,2),30,T,"filled")
% Make an iterative scree plot
for k = 1:10;
    opts = statset('Display','final');
    [Y,stress(k)] = mdscale(dissimilarities,k,'criterion','stress','start','random','replicates',5,'Options',opts);
end
figure, plot(stress,'.-k')



