
function [yt,ytab] = ResamplingMULTIPROX(datain,depthheaders,resolution) 

% datain = input table with all the data you want to homogenize separated
% in cells. Depth should be provided in mm
% datain{1} = readtable("G:\Research Projects\2022_Uebeschisee\UBE22-1e_xrf.xlsx")
% measselect = ["position_mm_","Al","Fe","Mn","K","Ti","Ca","Si","S","Br","Sr","CrInc","CrCoh"]
% datain{1} = datain{1}(:,measselect)
% datain{2} = readtable("G:\Research Projects\2022_Uebeschisee\UBE22-1e_hsi.csv")
% datain{2} = datain{2}(:,3:10)
% resolution = subsampling rate on desired mm
% resolution = 1
% headers of the depth collumns are provided as text
% depthheaders = ["position_mm_","Core_Y_mm_"];

clear y yt

for i = 1:length(datain)
 
    headers{i} = datain{i}.Properties.VariableNames;

    % datain{2} = data.hsi{6}
    % head{2} = datain{2}.Properties.VariableNames
    
    depthcol{i} = find(datain{i}.Properties.VariableNames == depthheaders(i));
    datain{i} = table2array(datain{i})
    range(:,i) = [max(datain{i}(:,depthcol{i})),min(datain{i}(:,depthcol{i}))];

end

minval = max(range(2,:)); maxval = min(range(1,:));
gridtime = (minval:resolution:maxval)';

% minval = round(min(datain{1}(:,depthcol{1})),2)
% maxval = round(max(datain{1}(:,depthcol{1})),2)

z = 1;
newhead(1) = depthheaders(1);

for i = 1:length(datain);

    data_cut{i} = datain{i}((logical((minval<datain{i}(:,depthcol{i})).*(datain{i}(:,depthcol{i})<maxval))),:);
    % xrfinit= datain{1}(logical((minval<table2array(data.sel{1,18}(:,"position (mm)"))).*(table2array(data.sel{1,18}(:,"position (mm)"))<maxval)),:)

    for j = 1:size(data_cut{i},2)-1; 
        ts = timeseries(data_cut{i}(:,j+1),data_cut{i}(:,depthcol{i}));
        tsnew = (resample(ts,gridtime)); 
        y(:,j+1) = tsnew.data;
        yt(:,j+z) = tsnew.data;
        newhead(j+z) = headers{i}(j+1);
    end
    y(:,1) = gridtime;
    z = z+j

end
yt(:,1) = gridtime;
ytab = array2table(yt,'VariableNames',[newhead])

end



% writetable(yt,"F://Amsoldingersee/Output/HSIXRFdata_resampled.csv");
