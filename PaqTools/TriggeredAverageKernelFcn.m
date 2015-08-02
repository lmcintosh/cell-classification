function outputfig = TriggeredAverageKernelFcn(trigtraces, rawtraces, zeroedtraces);
%varargin should contain
% - trigtraces
% - rawtraces
% - zeroedtraces

meantrigtraces = mean(trigtraces,1);
meanrawtraces = mean(rawtraces,1);
meanzeroedtraces = mean(zeroedtraces,1);        

data.data{1} = trigtraces;
data.data{2} = meantrigtraces;
data.data{3} = rawtraces;
data.data{4} = meanrawtraces;
data.data{5} = zeroedtraces;
data.data{6} = meanzeroedtraces;
data.names{1} = 'Triggers';
data.names{2} = 'AveragedTriggers';
data.names{3} = 'TriggeredTraces';
data.names{4} = 'AveragedTraces';
data.names{5} = 'ZeroedTraces';
data.names{6} = 'AveragedZeroedTraces';


if isempty(trigtraces);
    errordlg('No Triggers Detected')
    return
end

outputfig = figure('units','normalized','position',[.35 .1 .3 .85]);
TrigAvgExportMenu = uimenu('Label','Export','Parent',outputfig);
uimenu(TrigAvgExportMenu(end),'Label','To Base Workspace','callback',@ExportPeakTrigAvgToBaseWorkspace);
setappdata(outputfig,'data',data)

subplot(6,1,1);
plot(trigtraces');
xlim([0 size(trigtraces,2)])
maxval = max(max(trigtraces));
minval = min(min(trigtraces));
margin = .1*(maxval - minval);
ylim([minval-margin maxval+margin]);
title(['Triggers: ',num2str(size(trigtraces,1)),' total.'])

subplot(6,1,2);
plot(meantrigtraces');
xlim([0 size(meanrawtraces,2)])
maxval = max(max(meantrigtraces));
minval = min(min(meantrigtraces));
margin = .1*(maxval - minval);
ylim([minval-margin maxval+margin]);
title ('Averaged Triggers')

subplot(6,1,3);
plot(rawtraces');
xlim([0 size(rawtraces,2)])
maxval = max(max(rawtraces));
minval = min(min(rawtraces));
margin = .1*(maxval - minval);
ylim([minval-margin maxval+margin]);
title ('Triggered Traces')

subplot(6,1,4);
plot(meanrawtraces');
xlim([0 size(meanrawtraces,2)])
maxval = max(max(meanrawtraces));
minval = min(min(meanrawtraces));
margin = .1*(maxval - minval);
ylim([minval-margin maxval+margin]);
title ('Averaged Triggered Traces')

subplot(6,1,5);
plot(zeroedtraces');
xlim([0 size(zeroedtraces,2)])
maxval = max(max(zeroedtraces));
minval = min(min(zeroedtraces));
margin = .1*(maxval - minval);
ylim([minval-margin maxval+margin]);
title ('Zeroed Triggered Traces.  (Zero set by last point before trigger)')

subplot(6,1,6);
plot(meanzeroedtraces');
xlim([0 size(meanzeroedtraces,2)])
maxval = max(max(meanzeroedtraces));
minval = min(min(meanzeroedtraces));
margin = .1*(maxval - minval);
ylim([minval-margin maxval+margin]);
title ('Averaged Zeroed Traces');