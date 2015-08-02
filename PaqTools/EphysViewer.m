%%
function EphysViewer(varargin)
%this function will read a .daq file and will plot the data in a way
%similar to RecordDaq, but will plot each channel on a separate axes.  It
%will eventually hopefully have some nice data output options. 

warning off

handles.filetypes = {'*.daq;*.bin;*.paq;*.fig'};
handles.position = [0.2266    0.3359    0.5469    0.5469];

if nargin==1;%assume input is the path to the file with the data to view    
    handles.path=varargin{1};%save pathname into handles
elseif nargin==0;%if no inputs
    [filename,pathname,FilterIndex]=uigetfile(handles.filetypes,'Choose a data file');
    if ~FilterIndex
        return
	end
    handles.path=[pathname,filename];
end
ViewerMainFcn(handles);%open figure with all necessary axes etc.

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = ViewerMainFcn(handles)
%function ViewerMainFcn(handles)
%handles must contain .path, containing the path to the file to be opened

warning off
[pathstr, name, ext, versn] = fileparts(handles.path);
handles.filename = name;
if strcmp(ext,'.daq')
    info=bdaqread(handles.path,'info');
elseif strcmp(ext,'.bin');
    info=paq2lab(handles.path,'info');
elseif strcmp(ext,'.paq');
    info=paq2lab(handles.path,'info');
elseif strcmp(ext,'.fig');%for opening Figures representing Triggered Averages from EphysViewer
    ReadFig = hgload(handles.path);
    set(ReadFig,'visible','off');
    try
        ReadFigDat = getappdata(ReadFig,'data');
    catch
        errdlg ('Invalid Figure.  Must have "data" in appdata containing a structure of line data')
        return
    end
    close(ReadFig)
    clear ReadFig
    info.ObjInfo.SampleRate = 1;%assume nothing
    info.HwInfo = 'From Triggered Average Figure';
    info.ObjInfo.SamplesAquired = size(ReadFigDat.data{1},2);
    AveragedChannels = [2 4 6];%these are only channels with single datastrings
    for a = 1:length(AveragedChannels);%for each averaged channel
        thischannel = AveragedChannels(a);
        info.ObjInfo.Channel(a).HwChannel = thischannel;
        info.ObjInfo.Channel(a).Units = 'Units';
        info.ObjInfo.Channel(a).ChannelName = ReadFigDat.names{thischannel};
    end        
end
%Viewer expects the following information in info
% info.ObjInfo.SampleRate = Data Samples (Points) Per Second;
% info.HwInfo = Name of Acquistion system?;
% info.ObjInfo.SamplesAcquired = number of datapoints per channel;
%
% for each channel (channel = a)
% info.ObjInfo.Channel(a).HwChannel = Channel Number in the acquisition hardware;
% info.ObjInfo.Channel(a).Units = ie mV, nA, Volts, etc;
% info.ObjInfo.Channel(a).ChannelName = user-chosen name for each channel;
for a=1:length(info.ObjInfo.Channel)%for each channel
    if strcmp(ext,'.paq');
        filechans{a}=info.ObjInfo.Channel(a).HwChannel;
        nam=info.ObjInfo.Channel(a).ChannelName;
        charcell{a}=strcat(filechans{a},'-',nam);
    else
        filechans(a)=info.ObjInfo.Channel(a).HwChannel;
        nam=info.ObjInfo.Channel(a).ChannelName;
        charcell{a}=[num2str(filechans(a)),' - ',nam];
    end
end
% filechans=sort(filechans);

[handles.includedchannels,ok]=listdlg('ListString',charcell,'PromptString','Plot which channels... All?',...
    'Name','Channel Selection','InitialValue',1:length(charcell));%find out which channels to plot
if isempty(ok) || ~ok;
    return
end

totalaxes=length(handles.includedchannels);%taking data from info
handles.obj=info.ObjInfo;
handles.obj.HwInfo=info.HwInfo;
objuserdata.channelaxes=1:length(handles.includedchannels);

handles.displayfig=figure('Visible','off',...
    'DoubleBuffer','on',...
    'BackingStore','off',...
    'DithermapMode','manual',...
    'RendererMode','manual',...
    'Renderer','OpenGL',...
    'Units','Normalized',...
    'Position',handles.position,...
    'CloseRequestFcn',@CloseViewerFcn,...
    'numbertitle','off',...
    'tag','EphysViewer',...
    'name',handles.filename);%name after the file
%,...%     'menubar','none',...
%     'ResizeFcn',@ViewerResizing);%all may help faster rendering

handles.MainFigDataMenu=uimenu('Label','Data');
uimenu(handles.MainFigDataMenu,'Label','Load Data','Callback',@LocalMenuLoadDataFcn);
uimenu(handles.MainFigDataMenu,'Label','Export Visible Data','Callback',@LocalMenuExportVisibleDataFcn);

handles.MainFigDisplayMenu=uimenu('Label','Display');
uimenu(handles.MainFigDisplayMenu,'Label','Grid','Callback',@LocalMenuGridFcn);
% uimenu(handles.MainFigDisplayMenu,'Label','Background Color','Callback',@LocalMenuBackgroundColorFcn);

handles.MainFigDisplayMenu=uimenu('Label','Analyze');
uimenu(handles.MainFigDisplayMenu,'Label','Peak-Triggered Average','Callback',@PeakTrigAvgFcn);
uimenu(handles.MainFigDisplayMenu,'Label','Dip-Triggered Average','Callback',@DipTrigAvgFcn);
% uimenu(handles.MainFigDisplayMenu,'Label','Background Color','Callback',@LocalMenuBackgroundColorFcn);


handles.HorizScaleLookup=[.5 1 2 5 10];%these are horizontal scales to be chosen from: referring to fractions of a second 
%   (but will be used at different powers of 10)
handles.HorizScaleInButton=uicontrol('Style','PushButton',...
    'parent',handles.displayfig,...
    'Units','Normalized',...
    'Position',[.12 .065 .045 .045],...
    'String','In',...
    'Callback',@HorizScaleInButtonFcn);
handles.HorizScaleOutButton=uicontrol('Style','PushButton',...
    'parent',handles.displayfig,...
    'Units','Normalized',...
    'Position',[.166 .065 .045 .045],...
    'String','Out',...
    'Callback',@HorizScaleOutButtonFcn);
handles.HorizMinEditBox=uicontrol('Style','Edit',...
    'parent',handles.displayfig,...
    'Units','Normalized',...
    'Position',[.12 .0225 .045 .039],...
    'Callback',@HorizMinSetFcn);
handles.HorizMaxEditBox=uicontrol('Style','Edit',...
    'parent',handles.displayfig,...
    'Units','Normalized',...
    'Position',[.166 .0225 .045 .039],...
    'Callback',@HorizMaxSetFcn);
handles.TimeUnitsLabel = uicontrol('style','text',...
    'parent',handles.displayfig,...
    'units','normalized',...
    'position',[.12 .005 .1 .015],...
    'string','Seconds',...
    'backgroundcolor',[.8 .8 .8]);

handles.MoveLeftButton=uicontrol('style','pushbutton',...
    'parent',handles.displayfig,...
    'units','normalized',...
    'position',[.24 .01 .03 .1],...
    'string','<',...
    'callback',@MoveLeftFcn);
handles.MoveRightButton=uicontrol('style','pushbutton',...
    'parent',handles.displayfig,...
    'units','normalized',...
    'position',[.271 .01 .03 .1],...
    'string','>',...
    'callback',@MoveRightFcn);


handles.vertzoomfactor=2;

vertgap=.05;%constants of figure spacing
bottombottom=.11;
toptop=.98;
totalvertspace=toptop-bottombottom;
eachvertspace=totalvertspace/totalaxes;
axesheight=eachvertspace-vertgap;
for a = 1:totalaxes;%create each axes, equally spaced between certain set points

    ainv = totalaxes - a + 1;%last, last-1, last-2... 2, 1... just to make POSITION of axes right
    handles.axeshandles(a)=axes('parent',handles.displayfig,...
        'Units','Normalized',...
        'position',[.12 bottombottom+(eachvertspace*(ainv-1))+vertgap .78 axesheight],...
        'DrawMode','fast',...
        'XLimMode','manual',...%if not, then auto rescales sometimes on one axis
        'YLimMode','manual',...
        'ButtonDownFcn',@ZoomFcn);%immediately activate zooming

    axesposition=get(handles.axeshandles(a),'position');%do relative to these
    middlevert=axesposition(2)+(.5*axesposition(4));    
    widthslides=.047;
    
    try
        handles.VertSlideUp(a)=uicontrol('Style','PushButton',...
            'parent',handles.displayfig,...
            'Units','Normalized',...
            'Position',[axesposition(1)+axesposition(3)+.005 middlevert+.015 widthslides axesposition(4)/2-.015],...
            'String','^',...
            'Callback',@VertSlideUpFcn);%slides up by 1/4 field of view
        handles.VertSlideDown(a)=uicontrol('Style','PushButton',...
            'parent',handles.displayfig,...
            'Units','Normalized',...
            'Position',[axesposition(1)+axesposition(3)+.005 axesposition(2)+.03 widthslides axesposition(4)/2-.015],...
            'String','v',...
            'Callback',@VertSlideDownFcn);%slides down by 1/4 field of view                

        handles.VertScaleIn(a)=uicontrol('Style','PushButton',...
            'parent',handles.displayfig,...
            'Units','Normalized',...
            'Position',[axesposition(1)+axesposition(3)+.005+widthslides middlevert+.015 widthslides axesposition(4)/2-.015],...
            'String','In',...
            'Callback',@VertScaleInFcn);
        handles.VertScaleOut(a)=uicontrol('Style','PushButton',...
            'parent',handles.displayfig,...
            'Units','Normalized',...
            'Position',[axesposition(1)+axesposition(3)+.005+widthslides axesposition(2)+.03 widthslides axesposition(4)/2-.015],...
            'String','Out',...
            'Callback',@VertScaleOutFcn);
    catch
        error ('too many channels selected')
    end
    
    handles.VertMinEditBox(a)=uicontrol('Style','Edit',...
        'parent',handles.displayfig,...
        'Units','Normalized',...
        'Position',[axesposition(1)+axesposition(3)+.005 axesposition(2)-.025 widthslides .05],...
        'Callback',@VertMinSetFcn);
    handles.VertMaxEditBox(a)=uicontrol('Style','Edit',...
        'parent',handles.displayfig,...
        'Units','Normalized',...
        'Position',[axesposition(1)+axesposition(3)+.005+widthslides axesposition(2)-.025 widthslides .05],...
        'Callback',@VertMaxSetFcn);
        
    handles.ChannelsLabelAxes(a)=axes('Units','Normalized',...%an axes for putting names of channels.  Has xylims of 0,1
        'parent',handles.displayfig,...
        'Position',[0 axesposition(2) .08 axesposition(4)],...%to left of the plotting axes, same vert position and height
        'Visible','Off');
    
    handles.ChannelMeasureCheckbox(a) = uicontrol('style','checkbox',...
        'parent',handles.displayfig,...
        'units','normalized',...
        'position',[0 axesposition(2)+(axesposition(4)*0) .08 .02],...
        'backgroundcolor',[.8 .8 .8],...
        'string','Measure',...
        'FontSize',10,...
        'value',1);
end

handles.ResetAxesButton=uicontrol('style','pushbutton',...
    'units','normalized',...
    'position',[.33 .01 .05 .1],...
    'string',{'Reset';'Axes'},...
    'Callback',@ResetAxesButtonFcn);

handles.SingleMeasureButton = uicontrol('style','pushbutton',...
    'parent',handles.displayfig,...
    'units','normalized',...
    'position',[.65 .06 .12 .05],...
    'string','Single Measure',...
    'callback',@SingleMeasureFcn);%see comments below
handles.PairedMeasureButton = uicontrol('style','pushbutton',...
    'parent',handles.displayfig,...
    'units','normalized',...
    'position',[.65 .01 .12 .05],...
    'string','Paired Measure',...
    'callback',@PairedMeasureFcn);%see comments below
handles.SingleMeasureModeButton = uicontrol('style','pushbutton',...
    'parent',handles.displayfig,...
    'units','normalized',...
    'position',[.78 .06 .05 .05],...
    'string','Mode:',...
    'callback',@SingleMeasureModeFcn);
handles.PairedMeasureModeButton = uicontrol('style','pushbutton',...
    'parent',handles.displayfig,...
    'units','normalized',...
    'position',[.78 .01 .05 .05],...
    'string','Mode:',...
    'callback',@PairedMeasureModeFcn);
handles.SingleMeasureMode = 'Click';%default
handles.PairedMeasureMode = 'Click & Click';%default
handles.SingleMeasureModeIndic = uicontrol('style','text',...
    'parent',handles.displayfig,...
    'units','normalized',...
    'position',[.84 .0725 .14 .025],...
    'backgroundcolor',[.9 .9 .9],...
    'string',handles.SingleMeasureMode);
handles.PairedMeasureModeIndic = uicontrol('style','text',...
    'parent',handles.displayfig,...
    'units','normalized',...
    'position',[.84 .0225 .14 .025],...
    'backgroundcolor',[.9 .9 .9],...
    'string',handles.PairedMeasureMode);

handles.SingleMeasuresMatrix.MeasureNumber = [];
handles.SingleMeasuresMatrix.FileName = {};
handles.SingleMeasuresMatrix.ChannelName = {};
handles.SingleMeasuresMatrix.xs = [];
handles.SingleMeasuresMatrix.ys = [];
handles.SingleMeasuresMatrix.MeasureMode = {};
handles.SingleMeasuresMatrix.Comments = {};

handles.PairedMeasuresMatrix.MeasureNumber = [];
handles.PairedMeasuresMatrix.FileName = {};
handles.PairedMeasuresMatrix.ChannelName = {};
handles.PairedMeasuresMatrix.x1s = [];
handles.PairedMeasuresMatrix.y1s = [];
handles.PairedMeasuresMatrix.x2s = [];
handles.PairedMeasuresMatrix.y2s = [];
handles.PairedMeasuresMatrix.xdiffs = [];
handles.PairedMeasuresMatrix.ydiffs = [];
handles.PairedMeasuresMatrix.MeasureMode = {};
handles.PairedMeasuresMatrix.Comments = {};


handles.obj.userdata=objuserdata;

LoadDecimateStoreData(handles.displayfig,handles.path,handles.includedchannels,handles.obj.SampleRate);%load in data when the figure creation is finished

handles.sizedata = getappdata(handles.displayfig,'sizedata');
% handles.displaylength=handles.obj.SamplesAcquired/handles.obj.SampleRate; %total seconds
handles.displaylength=handles.sizedata(1)/handles.obj.SampleRate;%total seconds
set(handles.axeshandles,'XLim',[1/handles.obj.SampleRate handles.displaylength]);%set width to that specified earlier
% handles.HorizScaleTextBox=uicontrol('String',num2str(handles.displaylength));

PlotMinWithinXlims(handles,get(handles.axeshandles(1),'XLim'));

decdata=getappdata(handles.displayfig,'decdata');
for a=1:totalaxes
	channelsthisaxes=find(objuserdata.channelaxes==a);%find which channels are in that axes
    if ~isempty(channelsthisaxes)
        maxy=0;
        miny=0;
        for b=1:length(channelsthisaxes);%for each channel this axes
            maxy=max([maxy,max(decdata.data{1}(:,channelsthisaxes(b)))]);
            miny=min([miny,min(decdata.data{1}(:,channelsthisaxes(b)))]);
        end
        if maxy==0  % Added by Adam 01/18/07 to account for a channel that is ALWAYS zero!
            if miny==0
                maxy=1;
            end
        end
        handles.maxy(a)=maxy+abs(.05*maxy);%go up 5%
        handles.miny(a)=miny-abs(.05*miny);%go up 5%
        
        set(handles.axeshandles(a),'ylim',[handles.miny(a) handles.maxy(a)]);
        set(handles.VertMinEditBox(a),'String',num2str(handles.miny(a),3));
        set(handles.VertMaxEditBox(a),'String',num2str(handles.maxy(a),3));
    end
end
set(handles.HorizMinEditBox,'string','0');
set(handles.HorizMinEditBox,'string',num2str(handles.displaylength));

% colorlist={'blue';'green';'red';'cyan';'magenta';'yellow';'black';'white'};
% colorlookup=[0 0 1; 0 .5 .0; 1 0 0; 0 .75 .75; .75 0 .75; .75 .75 0; 0 0 0; 1 1 1];
for a=1:totalaxes;%for every axes
	channelsthisaxes=find(objuserdata.channelaxes==a);%find which channels are in that axes
    if ~isempty(channelsthisaxes)
%         colororder=[];
        for b=1:length(channelsthisaxes);%for each channel this axes
%             thiscolor=objuserdata.channelcolors{channelsthisaxes(b)};%extract name of chosen color
%             thiscolor=colorlookup(strmatch(thiscolor,colorlist,'exact'),:);%store the rgb value for that color
%             colororder(b,:)=thiscolor;
%             handles.ChannelLabels(a,b)=text(0, 1-((b-1)*.1),handles.obj.Channel(channelsthisaxes(b)).ChannelName,'Parent',handles.ChannelsLabelAxes(a),...
%                 'FontSize',12);%,'color',thiscolor);
            thischannelindex = handles.includedchannels(channelsthisaxes(b));
            handles.ChannelLabels(a,b)=text(0, .9-((b-1)*.1),...
                [handles.obj.Channel(thischannelindex).ChannelName,'(',handles.obj.Channel(thischannelindex).Units,')'],...
                'Parent',handles.ChannelsLabelAxes(a),...
                'FontSize',12);%,'color',thiscolor);
        end
%         set(handles.axeshandles(a),'ColorOrder',colororder);
    end
end
guidata(handles.displayfig,handles);

set(handles.displayfig,'Visible','on');
varargout{1} = 1;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LoadDecimateStoreData(fighandle,filename,includedchannels,samplingrate)
%uses handles to store decimated data into the appdata of the main figure

[pathstr, name, ext, versn] = fileparts(filename);
%data must be output as matrix having dims = pointnumber, channelnumber
if strcmp(ext,'.daq')
    alldata=bdaqread(filename,'Channels',includedchannels);
    %read data back to this workspace from the written file...only read specified channels
elseif strcmp(ext,'.bin')
    alldata=paq2lab(filename,'Channels',includedchannels);
elseif strcmp(ext,'.paq')
    alldata=paq2lab(filename,'Channels',includedchannels);
elseif strcmp(ext,'.fig');%again... assuming the fig is like a triggered average fig.
    ReadFig = hgload(filename);
    set(ReadFig,'visible','off');
    try
        ReadFigDat = getappdata(ReadFig,'data');
    catch
        errdlg ('Invalid Figure.  Must have "data" in appdata containing a structure of line data')
        return
    end
    close(ReadFig)
    clear ReadFig
    averagedchannels = [2 4 6];
    includedchannels = averagedchannels(includedchannels);
    for cidx = 1:length(includedchannels);
        alldata(:,cidx) = ReadFigDat.data{includedchannels(cidx)}';
    end
end

ordmag=floor(log10(length(alldata)));%get number of orders of magnitude occupied by the length 
decdata.data{1}=alldata;%store data for plotting
clear alldata;%save space
% decdata.lengths(1)=length(alldata);
decdata.decfactor(1)=1;
decdata.timepoints{1}=[];%this is just to take up space... don't want another vector as long as alldata
    %instead, for the special case of data{1}, time points will be the raw
    %index numbers of the points.

waithandle = waitbar(0,'Decimating and Storing Data');%for user
for a=1:ordmag-1;%for each order of magnitude
    [decdata.data{a+1},decdata.timepoints{a+1}]=decimatebymaxmin(decdata.data{1},10^a);
    decdata.timepoints{a+1}=decdata.timepoints{a+1}/samplingrate;
%store a version of the data decimated by that order of magnitude
%     decdata.lengths(a+1)=length(decdata.data{a+1});
    decdata.decfactor(a+1)=10^a;%store the decimation factor
    waithandle = waitbar(a/(ordmag-1),waithandle);
end
setappdata(fighandle,'decdata',decdata);
setappdata(fighandle,'sizedata',size(decdata.data{1}));
close(waithandle)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout=decimatebymaxmin(vect,factor)

flipindicator=0;
if length(size(vect))>2;
    error('input must be 2D or less')
elseif size(vect,1)<size(vect,2);%want size to be big x small (data stream by stream number)
    vect=vect';
    flipindicator=1;
end
numchans = size(vect,2);

factor=factor*2;
remain=rem(length(vect),factor);
if remain>0;%if number of points is not divisible by factor
    newvect=vect(1:end-remain,:);
    if remain>2;%if enough points to get a max and a min
        temp=vect(end-remain+1:end,:);
        temp(end+1:factor,:)=repmat(mean(temp,1),[factor-size(temp,1) 1]).*ones(factor-size(temp,1),size(temp,2));
        newvect(end+1:end+factor,:)=temp;
    end
else%if divisible by factor
    newvect=vect;
end
clear vect

newvect=reshape(newvect,[factor length(newvect)/factor  size(newvect,2)]);%set up for factor-fold decimation of data  

[nvmax,maxpt]=max(newvect,[],1);
[nvmin,minpt]=min(newvect,[],1);

newvect=zeros(2*size(nvmax,2),numchans);
newvect(1:2:size(newvect,1)-1,:)=squeeze(nvmax);%odd numbered points are maxs
newvect(2:2:size(newvect,1),:)=squeeze(nvmin);%even points are local mins

if nargout>1;
    temp=((0:2:(size(newvect,1)-1))*(factor/2));
    timepoints=sort(cat(2,temp+round(factor/2/3),temp+round(factor/3)))';
end

if flipindicator==1;
    newvect=newvect';
    timepoints=timepoints';
end

varargout{1}=newvect;
if nargout>1;
    varargout{2}=timepoints;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotMinWithinXlims(handles,xlims)
%this function takes the xlims it's given and takes handles, which has
%decimated data in its appdata (decdata structure) and plots the
%minimum-necessary-resolution data between those xlims.  It also deletes
%previous data and it plots specific lines in specific axes, as indicated
%by channelaxes in the daqobject userdata (handles.obj.userdata).

decdata=getappdata(handles.displayfig,'decdata');
       
set(handles.axeshandles(1),'units','pixels');%set units of an axes to pixels so can measure in pixels
pixelwidth=get(handles.axeshandles(1),'position');%record the pixels occupied by that axes
set(handles.axeshandles(1),'units','normalized');%restore the "units" property to its original value
pixelwidth=pixelwidth(3);%save just the width (in pixels) of the axes (not all about its location)

rationow=(((xlims(2)-xlims(1))*handles.obj.SampleRate)/pixelwidth);
index=max(find(decdata.decfactor<(rationow/2)));%find which is the shortest version of the decimated data that still has 
%at least 2x as many points as the axes has pixels
if isempty(index);%if the needed resolution is greater than exists, just use the highest available
    index=1;
end
objuserdata=handles.obj.userdata;
for a=1:length(handles.axeshandles);%for each axes
    channelsthisaxes=find(objuserdata.channelaxes==a);%find which channels are in that axes
    if ~isempty(channelsthisaxes);
		delete(get(handles.axeshandles(a),'children'));%delete all lines in the window
        temp=decdata.data{index}(:,channelsthisaxes);%get the data for just the channels in this axes
        thisaxes=handles.axeshandles(a);

        if index==1;
            yindices=round(xlims(1)*handles.obj.SampleRate):round(xlims(2)*handles.obj.SampleRate);
            yindices(yindices>length(temp))=[];
            yindices(yindices<1)=[];
            xpoints=yindices/handles.obj.SampleRate;
            line(xpoints,temp(yindices,:),'parent',thisaxes,'HitTest','Off')%hit test off for easier zooming
		    drawnow
        else
            bottime=max(find(decdata.timepoints{index}<xlims(1)));%find last point less than the bottom xlim
            if isempty(bottime);
                bottime=1;
            end
            toptime=min(find(decdata.timepoints{index}>xlims(2)));%find first point greater than the top xlim
            if isempty(toptime);
                toptime=length(decdata.timepoints{index});
            end
            
            line(decdata.timepoints{index}(bottime:toptime),temp(bottime:toptime,:),'parent',thisaxes,'HitTest','Off')%hit test off for easier zooming
		    drawnow
        end
    end    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LocalMenuLoadDataFcn(obj,event)

handles=guidata(obj);%get general info from the figure

[pathstr, name, ext, versn]=fileparts(handles.path);
[filename,pathname,FilterIndex]=uigetfile(handles.filetypes,'Choose a data file',[pathstr,'\']);
if ~FilterIndex
    return
end
%reinitialize handles for next version of Viewer
handles2.filetypes = {'*.daq;*.bin;*.paq;*.fig'};%because ViewerMainFcn expects a handles with very little (ie from top fcn)
handles2.path=[pathname,filename];
handles2.position = get(handles.displayfig,'position');
try
    vieweroutput = ViewerMainFcn(handles2);%vieweroutput it set to 1 at end of function
        % ie if all went well... if not all goes well the expectation of
        % an output creates an error... so go to the catch
catch
    vieweroutput = 0;
end
if vieweroutput
    close(handles.displayfig);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LocalMenuExportVisibleDataFcn(obj,event)
%function LocalMenuExportDataFcn(obj,event);
%This function allows the user to export the data (actual numbers)
%contained in the viewed data (data within the current xlimits).  Data can
%either be exported to the base workspace, to a .mat file or to a binary
%file, which has the same precision as the data was originally encoded in.

handles=guidata(obj);%get general info from the figure
data=getappdata(handles.displayfig,'decdata');%get saved loaded data
if isempty(data)
    helpdlg('Data must be loaded before exporting')
    return
end

charcell={'Base Workspace','Base Workspace (mem efficient)','.Mat File','Binary File'};%set up a list of options for export
[selection,ok]=listdlg('ListString',charcell,'SelectionMode','single','Name','Export List','PromptString','Export to where?');
if isempty(ok) || ~ok;%if user hit cancel or closed the window
    return%exit the export function
end

data=data.data{1};%keep just the full data (not the decimated versions)
xlims=get(handles.axeshandles(1),'xlim');%get the width of the current viewing window

yindices=round(xlims(1)*handles.obj.SampleRate):round(xlims(2)*handles.obj.SampleRate);
if yindices(1)<1;
    yindices(1)=[];
end
if yindices(end)>handles.obj.SamplesAcquired/handles.obj.SampleRate;
    yindices(end)=[];
end
yindices=yindices';
data=data(yindices,:);
% pause(1);
% xpoints=yindices./handles.obj.SampleRate;
yindices=yindices./handles.obj.SampleRate;

switch selection%depending on what the user chose to do
    case 1 %'Workspace'
% 		data=data(yindices,:);
        for a=1:size(data,2);%also corresponds to handles.includedchannels
            eval(['ExpDat.',...
                'chan',num2str(handles.obj.Channel(handles.includedchannels(a)).HwChannel),...
                '_',handles.obj.Channel(handles.includedchannels(a)).ChannelName,...
                '=data(:,a);']);
        end
%         ExpDat.time=xpoints;
        ExpDat.time=yindices;
%         c=clock;
%         c=[num2str(c(4)),num2str(c(5))];
%         name=inputdlg('Enter Name of Output Variable','Enter Name',1,{['ExpDat',c]});
        assignin('base','ExpDat',ExpDat);%export to base workspace
    case 2 %'Base Workspace (mem efficient)'
%         data=data(yindices,:);
        assignin('base','Time',yindices);
        assignin('base','ChannelNames',{handles.obj.Channel(handles.includedchannels).ChannelName});
        assignin('base','DataFromViewer',data);
    case 3 %'.Mat File'
% 		data=data(yindices,:);
        for a=1:size(data,2);%also corresponds to handles.includedchannels
            eval(['ExpDat.',...
                'chan',num2str(handles.obj.Channel(handles.includedchannels(a)).HwChannel),...
                '_',handles.obj.Channel(handles.includedchannels(a)).ChannelName,...
                '=data(:,a);']);
        end
%         ExpDat.time=xpoints;
        ExpDat.time=yindices;
        try
            [pathstr,name,ext,versn] = fileparts(handles.obj.LogFileName);%get name of the .daq file from which data was loaded        
            [FileName,PathName,FilterIndex] = uiputfile('.mat','Save Data to Disk',[pathstr,'\']);%ask user for name... above.mat is default
        catch
            [FileName,PathName,FilterIndex] = uiputfile('.mat','Save Data to Disk');%ask user for name... above.mat is default
        end
        if FilterIndex%if user did not hit cancel
            savefile=[PathName,FileName];
            save(savefile,'ExpDat');%save to mat file
        else
            return
        end
    case 4 %'Binary File'
        msgbox('Not yet implemented','Binary output','Help')

%         if strcmp(class(handles.obj),'analoginput');%depending on whether this is being called from RecordDaq or DaqViewer
%     		datatype=daqhwinfo(handles.obj);%get hwinfo
% 			datatype=datatype.NativeDataType;%ie int16
%         else
%             datatype=handles.obj.HwInfo.NativeDataType;%ie int16,uint16, etc
%         end
%         clear data%dump data... we don't need it any more since we'll read it from the file in native format
%         yindices=[yindices(1) yindices(end)];%this is all we'll need for daqread, this way is better for ram if big stuff later
%         [pathstring,name,ext,versn] = fileparts(handles.path);
%         for a=1:length(handles.obj.Channel);%for every hwchannel
%             hwchans(a)=handles.obj.Channel(a).HwChannel;%record it's number
%         end
%         for a=1:length(handles.includedchannels);%for each channel
%             thishwchannel=hwchans(handles.includedchannels(a));%record the number of the specified hwchannel
%             totalpathstr=[pathstring,'\',name,'_chan',num2str(thishwchannel),'_',...%make a path for a new file to be saved for each channel
%                     handles.obj.Channel(handles.includedchannels(a)).ChannelName,'_',datatype,'.bin'];%including it's HW channel number and name
% 			[FileName,PathName,FilterIndex] = uiputfile('.bin','Save This Channel to Disk?',totalpathstr);
%             if FilterIndex%if user did not hit cancel
%                 data=bdaqread(handles.path,'DataFormat','Native',...%read data only for specified channels for specified samples
%                     'Channels',handles.includedchannels(a),'Samples',yindices);
%                 savefile=[PathName,FileName];
%                 [trash1,trash2,ext,versn] = fileparts(savefile);
%                 if isempty(ext);%if no extension;
%                     savefile=[savefile,'.bin'];%add one
%                 end
%                 fid=fopen(savefile,'w');%create a file with the user-specified name
%                 eval(['fwrite(fid,data,''',datatype,''');']);%save data into it
%                 fclose(fid);
%             else
%                 continue
%             end
%         end
%         totalpathstr=[pathstring,'\',name,'_times_double.bin'];%including it's HW channel number and name
% 		[FileName,PathName,FilterIndex] = uiputfile('.bin','Save This Channel to Disk?',totalpathstr);
%         if FilterIndex%if user did not hit cancel
%             savefile=[PathName,FileName];
%             [trash1,trash2,ext,versn] = fileparts(savefile);
%             if isempty(ext);%if no extension;
%                 savefile=[savefile,'.bin'];%add one
%             end
%             fid=fopen(savefile,'w');%create a file with the user-specified name
%             fwrite(fid,xpoints,'double');%save data into it
%             fclose(fid);
%         end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LocalMenuGridFcn(obj,event)

handles=guidata(obj);%get general info from the figure

if strcmp(get(handles.axeshandles(1),'XGrid'),'off');
    set(handles.axeshandles,'XGrid','on','YGrid','on')
else
    set(handles.axeshandles,'XGrid','off','YGrid','off')
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CloseViewerFcn(obj,event)

delete(obj)%deleting the figure

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function HorizScaleInButtonFcn(obj,ev)

handles=guidata(obj);%get general info from the figure
xlimits=get(handles.axeshandles(1),'XLim');%in seconds
power10=10^floor(log10(handles.displaylength));
m=handles.displaylength/power10;
value=max(find(handles.HorizScaleLookup<m));%find indices of the preset largest value that is less than the current display length
value=handles.HorizScaleLookup(value)*power10;%seconds

if value>handles.sizedata(1)/handles.obj.SampleRate%if wants to zoom out to wider than total num of points
    value = handles.sizedata(1)/handles.obj.SampleRate;
end

if value<5/handles.obj.SampleRate;%don't allow to zoom too much
    value=5/handles.obj.SampleRate;
end
meanpos=mean(xlimits);
newlims=[meanpos-.5*value meanpos+.5*value];
set(handles.HorizMinEditBox,'String',num2str(newlims(1),3));
set(handles.HorizMaxEditBox,'String',num2str(newlims(2),3));
set(handles.axeshandles,'XLim',newlims);

PlotMinWithinXlims(handles,newlims)

handles.displaylength=value;
% set(handles.HorizScaleTextBox,'String',num2str(handles.displaylength));%set value of text box in seconds, and set it late
%     %in the function so it represents the most recent information possible

guidata(handles.displayfig,handles);%pass data back to the figure


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function HorizScaleOutButtonFcn(obj,ev);

handles=guidata(obj);%get general info from the figure

xlimits=get(handles.axeshandles(1),'XLim');
power10=10^floor(log10(handles.displaylength));
m=handles.displaylength/power10;
value=min(find(handles.HorizScaleLookup>m));%find indices of the preset largest value that is less than the current display length
value=handles.HorizScaleLookup(value)*power10;
meanpos=mean(xlimits);
newlims=[meanpos-.5*value meanpos+.5*value];

defaultstart=1/handles.obj.SampleRate;
% defaultend=handles.obj.SamplesAcquired/handles.obj.SampleRate;
defaultend=handles.sizedata(1)/handles.obj.SampleRate;%total seconds
if value>defaultend%both the left and right ends will be over the limits
    newlims=[defaultstart defaultend];%fix that
    value=handles.obj.SamplesAcquired/handles.obj.SampleRate;%save value correctly so later handles.displaylength will be correct
elseif newlims(1)<defaultstart;%if only the left side is too low (but right is fine)
    newlims=[defaultstart value];%plot from zero to the required length
elseif newlims(2)>defaultend;%if only the right side is too high
    newlims=[defaultend-value defaultend];%plot required length back from the end of the data
end
set(handles.HorizMinEditBox,'String',num2str(newlims(1),3));
set(handles.HorizMaxEditBox,'String',num2str(newlims(2),3));
set(handles.axeshandles,'XLim',newlims);

PlotMinWithinXlims(handles,newlims)

handles.displaylength=value;
% set(handles.HorizScaleTextBox,'String',num2str(handles.displaylength));%set value of text box in seconds, and set it late
    %in the function so it represents the most recent information possible

guidata(handles.displayfig,handles);%pass data back to the figure


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function HorizScaleTextBoxFcn(obj,ev);
% 
% handles=guidata(obj);%get general info from the figure
% 
% value=str2double(get(handles.HorizScaleTextBox,'string'));%in seconds
% if ~isempty(value);
%     if value<5/handles.obj.SampleRate;
%         value=5/handles.obj.SampleRate;
%         disp('Entered value is too low')
%     end
% 	xlimits=get(handles.axeshandles(1),'XLim');
%     value=value;%set this value for later
%     meanpos=mean(xlimits);
%     newlims=[meanpos-.5*value meanpos+.5*value];
%     
% 	defaultstart=1/handles.obj.SampleRate;
% 	defaultend=handles.obj.SamplesAcquired/handles.obj.SampleRate;
%     if value>defaultend%both the left and right ends will be over the limits
%         newlims=[defaultstart defaultend];%fix that
%         value=defaultend;%save value correctly so later handles.displaylength will be correct
%         disp('Entered value is too high')
%     elseif newlims(1)<defaultstart;%if only the left side is too low (but right is fine)
%         newlims=[defaultstart value];%plot from zero to the required length
%     elseif newlims(2)>defaultend;%if only the right side is too high
%         newlims=[defaultend-value defaultend];%plot required length back from the end of the data
%     end
% 	set(handles.axeshandles,'XLim',newlims);        
% 	PlotMinWithinXlims(handles,newlims)
%     
%     handles.displaylength=value;
% 	set(handles.HorizScaleTextBox,'String',num2str(handles.displaylength));%set value of text box in seconds, and set it late
%         %in the function so it represents the most recent information possible
% 	guidata(handles.displayfig,handles);%pass data back to the figure
% end
% 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MoveLeftFcn(obj,ev)

handles=guidata(obj);%get general info from the figure
% set(handles.obj,'TimerFcn','');%pause plotting (resume later)
% pause(.1)
xlimits=get(handles.axeshandles(1),'XLim');
dl=xlimits(2)-xlimits(1);

newlims=[xlimits(1)-.5*dl xlimits(2)-.5*dl];
if newlims(1)<0;
    newlims=[1/handles.obj.SampleRate 1/handles.obj.SampleRate+dl];
end
set(handles.axeshandles,'XLim',newlims);
set(handles.HorizMinEditBox,'String',num2str(newlims(1),3));
set(handles.HorizMaxEditBox,'String',num2str(newlims(2),3));
PlotMinWithinXlims(handles,newlims)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MoveRightFcn(obj,ev)

handles=guidata(obj);%get general info from the figure
% set(handles.obj,'TimerFcn','');%pause plotting (resume later)
% pause(.1)
xlimits=get(handles.axeshandles(1),'XLim');
dl=xlimits(2)-xlimits(1);

newlims=[xlimits(1)+.5*dl xlimits(2)+.5*dl];
decdata=getappdata(handles.displayfig,'decdata');
%len=length(decdata.data{1});

% defaultend=handles.obj.SamplesAcquired/handles.obj.SampleRate;
defaultend=handles.sizedata(1)/handles.obj.SampleRate;%total seconds
if newlims(2)>defaultend;%if right side is going beyond the end of the data
    newlims=[defaultend-dl defaultend];
end
set(handles.axeshandles,'XLim',newlims);
set(handles.HorizMinEditBox,'String',num2str(newlims(1),3));
set(handles.HorizMaxEditBox,'String',num2str(newlims(2),3));
PlotMinWithinXlims(handles,newlims)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function VertSlideUpFcn(obj,ev)

handles=guidata(obj);%get general info from the figure

ind=find(handles.VertSlideUp==obj);%find number of axis to be changed

ylims=get(handles.axeshandles(ind),'Ylim');
range=ylims(2)-ylims(1);
newylims=ylims+.3*range;
maxy=handles.maxy(ind);
if max(newylims)>maxy;
    newylims=[maxy-range maxy];
end

set(handles.axeshandles(ind),'Ylim',newylims);%slide up 1/4 field of view
set(handles.VertMinEditBox(ind),'String',num2str(newylims(1),3));
set(handles.VertMaxEditBox(ind),'String',num2str(newylims(2),3));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function VertSlideDownFcn(obj,ev)

handles=guidata(obj);%get general info from the figure

ind=find(handles.VertSlideDown==obj);%find number of axis to be changed
ylims=get(handles.axeshandles(ind),'Ylim');
range=ylims(2)-ylims(1);
newylims=ylims-.3*range;

miny=handles.miny(ind);
if min(newylims)<miny;
    newylims=[miny miny+range];
end

set(handles.axeshandles(ind),'Ylim',newylims);%slide up 1/4 field of view
set(handles.VertMinEditBox(ind),'String',num2str(newylims(1),3));
set(handles.VertMaxEditBox(ind),'String',num2str(newylims(2),3));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function VertScaleInFcn(obj,ev)

handles=guidata(obj);%get general info from the figure
ind=find(handles.VertScaleIn==obj);%find number of axis to be changed

ylims=get(handles.axeshandles(ind),'Ylim');
m=mean(ylims);
siderange=m-ylims(1);
newylims=[m-(1/handles.vertzoomfactor)*siderange m+(1/handles.vertzoomfactor)*siderange];
set(handles.axeshandles(ind),'Ylim',newylims);%slide up 1/4 field of view
set(handles.VertMinEditBox(ind),'String',num2str(newylims(1),3));
set(handles.VertMaxEditBox(ind),'String',num2str(newylims(2),3));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function VertScaleOutFcn(obj,ev)

handles=guidata(obj);%get general info from the figure
if isempty(ev)
    ind=find(handles.VertScaleOut==obj);%find number of axis to be changed
else
    ind=find(handles.axeshandles==obj);%find the axis to be changed
end

ylims=get(handles.axeshandles(ind),'Ylim');
m=mean(ylims);
siderange=m-ylims(1);
newylims=[m-(handles.vertzoomfactor)*siderange m+(handles.vertzoomfactor)*siderange];

miny=handles.miny(ind);%use the index number found above
maxy=handles.maxy(ind);
if siderange*2*handles.vertzoomfactor>maxy-miny;
    newylims=[miny maxy];
elseif newylims(1)<miny && newylims(2)<=maxy;%if top side is greater than allowable, but bottom is ok
    newylims=[miny miny+siderange*2*handles.vertzoomfactor];%keep same range, but set top to maxy
elseif newylims(2)>maxy && newylims(1)>=miny;%if bottom side is less than allowable, but top is ok
    newylims=[maxy-siderange*2*handles.vertzoomfactor maxy];%keep same range, but set bottom to miny
end

set(handles.axeshandles(ind),'Ylim',newylims);%slide up 1/4 field of view
set(handles.VertMinEditBox(ind),'String',num2str(newylims(1),3));
set(handles.VertMaxEditBox(ind),'String',num2str(newylims(2),3));


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function VertMinSetFcn(obj,ev)

handles=guidata(obj);%get general info from the figure
ind=find(handles.VertMinEditBox==obj);%find number of axis to be changed
newmin=str2double(get(handles.VertMinEditBox(ind),'String'));%get user input
if newmin<handles.miny(ind);%if value was beyond possible input values
    newmin=handles.miny(ind);%set to the limit input value
    set(handles.VertMinEditBox(ind),'String',num2str(newmin,3));
end
newylims=get(handles.axeshandles(ind),'Ylim');%get old limits
newylims(1)=newmin;%sub in new value
set(handles.axeshandles(ind),'Ylim',newylims);%set axes


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function VertMaxSetFcn(obj,ev)

handles=guidata(obj);%get general info from the figure
ind=find(handles.VertMaxEditBox==obj);%find number of axis to be changed
newmax=str2double(get(handles.VertMaxEditBox(ind),'String'));%get user input
if newmax>handles.maxy(ind);%if value was beyond possible input values
    newmax=handles.maxy(ind);%set to the limit input value
    set(handles.VertMaxEditBox(ind),'String',num2str(newmax,3));
end
newylims=get(handles.axeshandles(ind),'Ylim');%get old limits
newylims(2)=newmax;%sub in new value
set(handles.axeshandles(ind),'Ylim',newylims);%set axes


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function HorizMinSetFcn(obj,ev)

handles=guidata(obj);%get general info from the figure
xlimits=get(handles.axeshandles(1),'XLim');%in seconds

newlims = xlimits;
newlims(1) = str2double(get(obj,'string'));
if newlims(1)<0;
    newlims(1)=1/handles.obj.SampleRate;
end
set(handles.axeshandles,'XLim',newlims);
PlotMinWithinXlims(handles,newlims)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function HorizMaxSetFcn(obj,ev)

handles=guidata(obj);%get general info from the figure
xlimits=get(handles.axeshandles(1),'XLim');%in seconds

newlims = xlimits;
newlims(2) = str2double(get(obj,'string'));
defaultend=handles.sizedata(1)/handles.obj.SampleRate;%total seconds
if newlims(2)>defaultend;%if right side is going beyond the end of the data
    newlims(2) = defaultend;
end
set(handles.axeshandles,'XLim',newlims);
PlotMinWithinXlims(handles,newlims)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ZoomFcn(obj,ev)

handles=guidata(obj);

point1 = get(handles.displayfig,'CurrentPoint'); % button down detected
if strcmp(get(handles.displayfig,'SelectionType'),'normal');%if a left click
	rect = [point1(1,1) point1(1,2) 0 0];
    r2 = rbbox(rect);%in figure units... need to convert to x and y coords... mult by axes relative position and by xlim and ylim of ax
    if (r2(3)+r2(4))>=.01;%if at least some minimum box was made
		axpos=get(obj,'position');
		oldx=get(obj,'xlim');
		oldy=get(obj,'ylim');
		xlims(1)=(r2(1)-axpos(1))/axpos(3);%convert first point to x start in axes normalized
		xlims(2)=r2(3)/axpos(3)+xlims(1);%convert third point to x stop in axes normalized
		xlims(1)=oldx(1)+xlims(1)*(oldx(2)-oldx(1));
		xlims(2)=oldx(1)+xlims(2)*(oldx(2)-oldx(1));
		ylims(1)=(r2(2)-axpos(2))/axpos(4);
		ylims(2)=r2(4)/axpos(4)+ylims(1);
		ylims(1)=oldy(1)+ylims(1)*(oldy(2)-oldy(1));
		ylims(2)=oldy(1)+ylims(2)*(oldy(2)-oldy(1));
		
		handles.displaylength=xlims(2)-xlims(1);
    	set(handles.axeshandles,'xlim',xlims);%set all xlims the same
		set(obj,'ylim',ylims);
		PlotMinWithinXlims(handles,xlims)
        
    	handles.displaylength=(xlims(2)-xlims(1));%
        
        %set x and y text boxes
        set(handles.HorizMinEditBox,'String',num2str(xlims(1),3));
        set(handles.HorizMaxEditBox,'String',num2str(xlims(2),3));
        ind = find(handles.axeshandles==obj);%find number of axis to be changed
        set(handles.VertMinEditBox(ind),'String',num2str(ylims(1),3));%set text boxes
        set(handles.VertMaxEditBox(ind),'String',num2str(ylims(2),3));

% 		set(handles.HorizScaleTextBox,'String',num2str(handles.displaylength));%set value of text box in seconds, and set it late
            %in the function so it represents the most recent information possible
		guidata(handles.displayfig,handles);%pass data back to the figure

    end
elseif strcmp(get(handles.displayfig,'SelectionType'),'alt');%if a left click
    HorizScaleOutButtonFcn(obj,'notempty');
    VertScaleOutFcn(obj,'notempty');
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ResetAxesButtonFcn(obj,ev)

handles=guidata(obj);
% handles.displaylength=handles.obj.SamplesAcquired/handles.obj.SampleRate;
handles.displaylength=handles.sizedata(1)/handles.obj.SampleRate;%total seconds
for a=1:length(handles.axeshandles)
    set(handles.axeshandles(a),'ylim',[handles.miny(a) handles.maxy(a)])
    set(handles.VertMinEditBox(a),'String',num2str(handles.miny(a),3));
    set(handles.VertMaxEditBox(a),'String',num2str(handles.maxy(a),3));
end
set(handles.HorizMinEditBox,'string',num2str(0,3));
set(handles.HorizMaxEditBox,'string',num2str(handles.displaylength,3));

set(handles.axeshandles,'xlim',[1/handles.obj.SampleRate handles.displaylength]);

PlotMinWithinXlims(handles,[1/handles.obj.SampleRate handles.displaylength])
guidata(handles.displayfig,handles)

%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function ViewerResizing(obj,ev)
% 
% handles=guidata(obj);
% for a=1:size(handles.ChannelLabels,1);
%     set(handles.axeshandles(a),'units','points');
%     pts=get(handles.axeshandles(a),'position');
%     set(handles.axeshandles(a),'units','normalized');
%     pts=pts(4);%get full height of the axes
%     for b=1:size(handles.ChannelLabels,2);
%         ypos=pts-((b-1)*18);
%         set(handles.ChannelLabels(a,b),'units','points','position',[1 ypos]);
%     end
% end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SingleMeasureFcn(obj,ev)

handles = guidata(obj);
measuremode = handles.SingleMeasureMode;
[x,y,channelstomeasure] = MeasurementKernelFcn(handles,measuremode);
if ~isempty(handles.SingleMeasuresMatrix.MeasureNumber)
    measnum = handles.SingleMeasuresMatrix.MeasureNumber(end)+1;
else
    measnum = 1;
end
for chidx = 1:size(y,2);
    handles.SingleMeasuresMatrix.MeasureNumber(end+1) = measnum;
    handles.SingleMeasuresMatrix.FileName{end+1} = handles.filename;
    cnn = handles.includedchannels(channelstomeasure(chidx));
	chname = handles.obj.Channel(cnn).ChannelName;
    handles.SingleMeasuresMatrix.ChannelName{end+1} = chname;
    handles.SingleMeasuresMatrix.xs(end+1) = x;
    handles.SingleMeasuresMatrix.ys(end+1) = y(chidx);
    handles.SingleMeasuresMatrix.MeasureMode{end+1} = measuremode;
    handles.SingleMeasuresMatrix.Comments{end+1} = '';
end

handles = SingleMeasuresBox(handles);

guidata(handles.displayfig,handles);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y,channelstomeasure] = MeasurementKernelFcn(handles,measuremode);
if ~strcmp('Click',measuremode)
    [tempx,tempy] = ginput(1);
    point1 = get(handles.displayfig,'CurrentPoint'); % button down detected
    if strcmp(get(handles.displayfig,'SelectionType'),'normal');%if a left click
        rect = [point1(1,1) point1(1,2) 0 0];
        r2 = rbbox(rect);%in figure units = seconds
        thisaxishand = gca;
        axpos = get(thisaxishand,'position');
        axxlims = get(thisaxishand,'xlim');
        xlims(1)=(r2(1)-axpos(1))/axpos(3);%convert first point to x start in axes normalized
		xlims(2)=r2(3)/axpos(3)+xlims(1);%convert third point to x stop in axes normalized
        

        xstartsec = axxlims(1) + xlims(1)*(axxlims(2)-axxlims(1));
        xstopsec = axxlims(1) + xlims(2)*(axxlims(2)-axxlims(1));
%         if (r2(3)+r2(4))>=.001;%if at least some minimum box was made
%         
            %get the number of the channel in the axes clicked.  Will
            %defnitely get measures from this channel
            thisaxisidx = find(handles.axeshandles == thisaxishand);%grab the number of the axes clicked... 
            channelsthisaxes=find(handles.obj.userdata.channelaxes==thisaxisidx);%get the channel in the axes...
            %will be only one channel in this program
            %find which channels had boxes checked indicating to measure them
            for a = 1:length(handles.ChannelMeasureCheckbox);
                channelschecked(a) = get(handles.ChannelMeasureCheckbox(a),'value');
            end
            channelschecked = find(channelschecked);
            %make final list of all channels to measure from channel
            %clicked on and those with boxes checked
            channelstomeasure = union(channelsthisaxes,channelschecked);%auto sorted in ascending order by union fcn
            
            decdata=getappdata(handles.displayfig,'decdata');
            channeldata=decdata.data{1}(:,channelstomeasure);
            clear decdata;%decdata is big, get rid of it
            
            yindices=round(xstartsec*handles.obj.SampleRate):round(xstopsec*handles.obj.SampleRate);
            yindices(yindices>length(channeldata))=[];
            yindices(yindices<1)=[];
            if isempty(yindices);
                return
            end
            xdata=yindices/handles.obj.SampleRate;
            ydata = channeldata(yindices,:);
            clear channeldata
            %for max and min: needs to be max/min of clicked channel and same time points on
            %other channels.  For mean: take mean across selected region on
            %each channel.
            %Then give just one x output
            if strcmp('Minimum',measuremode)
                selectedchidx = find(channelstomeasure == channelsthisaxes);
                [trash,x] = min(ydata(:,selectedchidx),[],1);%get x&y of min in the CLICKED ON channel
                %use x to get all the y's
                y = ydata(x,:);
                x = xdata(x);
            elseif strcmp('Maximum',measuremode)
                selectedchidx = find(channelstomeasure == channelsthisaxes);
                [trash,x] = max(ydata(:,selectedchidx),[],1);%get x&y of max in the CLICKED ON channel
                %use x to get all the y's
                y = ydata(x,:);
                x = xdata(x);
            elseif strcmp('Mean',measuremode)
                y = mean(ydata,1);
                x = mean(xdata);                
            end            
%         end
    end     
else
    [x,trash] = ginput(1);
    
    thisaxishand = gca;
    thisaxisidx = find(handles.axeshandles == thisaxishand);%grab the number of the axes clicked... 
    channelsthisaxes=find(handles.obj.userdata.channelaxes==thisaxisidx);%get the channel in the axes...
    %will be only one channel in this program
    %find which channels had boxes checked indicating to measure them
    for a = 1:length(handles.ChannelMeasureCheckbox);
        channelschecked(a) = get(handles.ChannelMeasureCheckbox(a),'value');
    end
    channelschecked = find(channelschecked);
    %make final list of all channels to measure from channel
    %clicked on and those with boxes checked
    channelstomeasure = union(channelsthisaxes,channelschecked);%auto sorted in ascending order by union fcn

    decdata=getappdata(handles.displayfig,'decdata');
    channeldata=decdata.data{1}(:,channelstomeasure);
    clear decdata;%decdata is big, get rid of it
    
    yindices=round(x*handles.obj.SampleRate);
    if yindices > length(channeldata) || yindices<1
        return
    end
    x = yindices/handles.obj.SampleRate;
    y = channeldata(yindices,:);
end



%% Creates/handles auxiliary figure to display the Paired Measures Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = SingleMeasuresBox(handles);
%makes an auxiliary figure to display the Single Measures Data
thisfig = findobj('type','figure','name','Single Measures','userdata',handles.displayfig);
if isempty(thisfig);
    handles.SingleMeasuresWindow = figure('units','points',...
    'position',[0 40 325 1],...
    'toolbar','none',...
    'menubar','none',...
    'name','Single Measures',...
    'numbertitle','off',...
    'userdata',handles.displayfig);
else
    handles.SingleMeasuresWindow = thisfig;
end

delete(get(handles.SingleMeasuresWindow,'children'));

handles.SingleMeasuresExportMenu = uimenu('Label','Export','Parent',handles.SingleMeasuresWindow);
uimenu(handles.SingleMeasuresExportMenu,'Label','To Workspace','callback',@ExportSingleMeasuresToWS);
uimenu(handles.SingleMeasuresExportMenu,'Label','To .mat File','callback',@ExportSingleMeasuresToMAT);
uimenu(handles.SingleMeasuresExportMenu,'Label','To .xls File','callback',@ExportSingleMeasuresToXLS);
uimenu(handles.SingleMeasuresExportMenu,'Label','To current Excel Doc','callback',@ExportSingleMeasuresToCurrentExcelBook);
handles.SingleMeasuresDeleteMenu = uimenu('Label','Delete','Parent',handles.SingleMeasuresWindow);
uimenu(handles.SingleMeasuresDeleteMenu,'Label','All','callback',@DeleteSingleMeasuresAll);
uimenu(handles.SingleMeasuresDeleteMenu,'Label','Select','callback',@DeleteSingleMeasuresSelect);

precision = 8;
lineheight = 12;
fontsize = 8;
set(handles.SingleMeasuresWindow,'units','points');
nummeas = length(handles.SingleMeasuresMatrix.xs);
winsz = get(handles.SingleMeasuresWindow,'position');
winsz = [winsz(1) max([0 winsz(2)-lineheight]) winsz(3) (nummeas+2)*lineheight];
set(handles.SingleMeasuresWindow,'position',winsz);

%Column titles
uicontrol('parent',handles.SingleMeasuresWindow,'style','Text','units','points',...
    'position',[10 winsz(4)-12 12 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','#');
uicontrol('parent',handles.SingleMeasuresWindow,'style','Text','units','points',...
    'position',[35 winsz(4)-12 55 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','Channel');
uicontrol('parent',handles.SingleMeasuresWindow,'style','Text','units','points',...
    'position',[120 winsz(4)-12 12 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','X');
uicontrol('parent',handles.SingleMeasuresWindow,'style','Text','units','points',...
    'position',[167 winsz(4)-12 12 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','Y');
uicontrol('parent',handles.SingleMeasuresWindow,'style','Text','units','points',...
    'position',[227 winsz(4)-12 55 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','Comment');

nummeas = length(handles.SingleMeasuresMatrix.xs);
for midx = 1:nummeas;
    uicontrol('parent',handles.SingleMeasuresWindow,'style','frame','units','points',...
        'position',[7  winsz(4)-((lineheight-1)*(midx+1))-1 20 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.SingleMeasuresWindow,'style','text','units','points',...
        'position',[8  winsz(4)-((lineheight-1)*(midx+1)) 18 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',num2str(handles.SingleMeasuresMatrix.MeasureNumber(midx)));

    uicontrol('parent',handles.SingleMeasuresWindow,'style','frame','units','points',...
        'position',[26  winsz(4)-((lineheight-1)*(midx+1))-1 75 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.SingleMeasuresWindow,'style','text','units','points',...
        'position',[27  winsz(4)-((lineheight-1)*(midx+1)) 73 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',handles.SingleMeasuresMatrix.ChannelName{midx});    
    
    uicontrol('parent',handles.SingleMeasuresWindow,'style','frame','units','points',...
        'position',[101  winsz(4)-((lineheight-1)*(midx+1))-1 50 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.SingleMeasuresWindow,'style','text','units','points',...
        'position',[102  winsz(4)-((lineheight-1)*(midx+1)) 48 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',num2str(handles.SingleMeasuresMatrix.xs(midx),precision));

    uicontrol('parent',handles.SingleMeasuresWindow,'style','frame','units','points',...
        'position',[150  winsz(4)-((lineheight-1)*(midx+1))-1 50 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.SingleMeasuresWindow,'style','text','units','points',...
        'position',[151  winsz(4)-((lineheight-1)*(midx+1)) 48 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',num2str(handles.SingleMeasuresMatrix.ys(midx),precision));

    uicontrol('parent',handles.SingleMeasuresWindow,'style','frame','units','points',...
        'position',[199  winsz(4)-((lineheight-1)*(midx+1))-1 119 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.SingleMeasuresWindow,'style','edit','units','points',...
        'position',[201  winsz(4)-((lineheight-1)*(midx+1)) 117 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',handles.SingleMeasuresMatrix.Comments{midx},...
        'tag',num2str(midx),'userdata',handles.displayfig,'callback',@SingleMeasuresCommentFcn);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SingleMeasuresCommentFcn(obj,ev);

EphysFig = get(obj,'userdata');
handles = guidata(EphysFig);
midx = str2double(get(obj,'tag'));
str = get(obj,'string');
handles.SingleMeasuresMatrix.Comments{midx}= str;
guidata(handles.displayfig,handles);

%%
function PairedMeasureFcn(obj,ev)
% take handles.PairedMeasureMode and split with the _&_
handles = guidata(obj);
measuremode = handles.PairedMeasureMode;
andspot = strfind(measuremode,' & ');
measuremode1 = measuremode(1:andspot-1);
measuremode2 = measuremode(andspot+3:end);
[x1,y1,channelstomeasure] = MeasurementKernelFcn(handles, measuremode1);
[x2,y2,channelstomeasure] = MeasurementKernelFcn(handles, measuremode2);

if ~isempty(handles.PairedMeasuresMatrix.MeasureNumber)
    measnum = handles.PairedMeasuresMatrix.MeasureNumber(end)+1;
else
    measnum = 1;
end
for chidx = 1:size(y2,2);
    xdiff = x2-x1;
    ydiff = y2(chidx)-y1(chidx);
    handles.PairedMeasuresMatrix.MeasureNumber(end+1) = measnum;
    handles.PairedMeasuresMatrix.FileName{end+1} = handles.filename;
    cnn = handles.includedchannels(channelstomeasure(chidx));
    chname = handles.obj.Channel(cnn).ChannelName;
    handles.PairedMeasuresMatrix.ChannelName{end+1} = chname;
    handles.PairedMeasuresMatrix.x1s(end+1) = x1;
    handles.PairedMeasuresMatrix.y1s(end+1) = y1(chidx);
    handles.PairedMeasuresMatrix.x2s(end+1) = x2;
    handles.PairedMeasuresMatrix.y2s(end+1) = y2(chidx);
    handles.PairedMeasuresMatrix.xdiffs(end+1) = xdiff;
    handles.PairedMeasuresMatrix.ydiffs(end+1) = ydiff;
    handles.PairedMeasuresMatrix.MeasureMode{end+1} = measuremode;
    handles.PairedMeasuresMatrix.Comments{end+1} = '';
end

handles = PairedMeasuresBox(handles);

guidata(handles.displayfig,handles);

%% Creates/handles auxiliary figure to display the Paired Measures Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = PairedMeasuresBox(handles)
thisfig = findobj('type','figure','name','Paired Measures','userdata',handles.displayfig);
set(0,'units','points');
screenpix = get(0,'ScreenSize');
pixperinch = get(0,'ScreenPixelsPerInch');
inches = screenpix/pixperinch;
screenpts = inches * 72;
boxwidth = 525;

if isempty(thisfig);
    handles.PairedMeasuresWindow = figure('units','points',...
    'position',[screenpts(3)-boxwidth 40 boxwidth 1],...
    'toolbar','none',...
    'menubar','none',...
    'name','Paired Measures',...
    'numbertitle','off',...
    'userdata',handles.displayfig);
else
    handles.PairedMeasuresWindow = thisfig;
end

delete(get(handles.PairedMeasuresWindow,'children'));

precision = 8;
lineheight = 12;
fontsize = 8;
set(handles.PairedMeasuresWindow,'units','points');
nummeas = length(handles.PairedMeasuresMatrix.x1s);
winsz = get(handles.PairedMeasuresWindow,'position');
winsz = [winsz(1) max([0 winsz(2)-lineheight]) winsz(3) (nummeas+2)*lineheight];
set(handles.PairedMeasuresWindow,'position',winsz);

handles.PairedMeasuresExportMenu = uimenu('Label','Export','Parent',handles.PairedMeasuresWindow);
uimenu(handles.PairedMeasuresExportMenu,'Label','To Workspace','callback',@ExportPairedMeasuresToWS);
uimenu(handles.PairedMeasuresExportMenu,'Label','To .mat File','callback',@ExportPairedMeasuresToMAT);
uimenu(handles.PairedMeasuresExportMenu,'Label','To .xls File','callback',@ExportPairedMeasuresToXLS);
uimenu(handles.PairedMeasuresExportMenu,'Label','To current Excel Doc','callback',@ExportPairedMeasuresToCurrentExcelBook);
handles.PairedMeasuresDeleteMenu = uimenu('Label','Delete','Parent',handles.PairedMeasuresWindow);
uimenu(handles.PairedMeasuresDeleteMenu,'Label','All','callback',@DeletePairedMeasuresAll);
uimenu(handles.PairedMeasuresDeleteMenu,'Label','Select','callback',@DeletePairedMeasuresSelect);



%Column titles
uicontrol('parent',handles.PairedMeasuresWindow,'style','Text','units','points',...
    'position',[10 winsz(4)-12 12 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','#');
uicontrol('parent',handles.PairedMeasuresWindow,'style','Text','units','points',...
    'position',[35 winsz(4)-12 55 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','Channel');
uicontrol('parent',handles.PairedMeasuresWindow,'style','Text','units','points',...
    'position',[120 winsz(4)-12 12 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','X1');
uicontrol('parent',handles.PairedMeasuresWindow,'style','Text','units','points',...
    'position',[167 winsz(4)-12 12 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','Y1');
uicontrol('parent',handles.PairedMeasuresWindow,'style','Text','units','points',...
    'position',[220 winsz(4)-12 12 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','X2');
uicontrol('parent',handles.PairedMeasuresWindow,'style','Text','units','points',...
    'position',[267 winsz(4)-12 12 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','Y2');
uicontrol('parent',handles.PairedMeasuresWindow,'style','Text','units','points',...
    'position',[311 winsz(4)-12 30 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','X2-X1');
uicontrol('parent',handles.PairedMeasuresWindow,'style','Text','units','points',...
    'position',[358 winsz(4)-12 30 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','Y2-Y1');
uicontrol('parent',handles.PairedMeasuresWindow,'style','Text','units','points',...
    'position',[430 winsz(4)-12 55 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','Comment');

nummeas = length(handles.PairedMeasuresMatrix.x1s);
for midx = 1:nummeas;
    uicontrol('parent',handles.PairedMeasuresWindow,'style','frame','units','points',...
        'position',[7  winsz(4)-((lineheight-1)*(midx+1))-1 20 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.PairedMeasuresWindow,'style','text','units','points',...
        'position',[8  winsz(4)-((lineheight-1)*(midx+1)) 18 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',num2str(handles.PairedMeasuresMatrix.MeasureNumber(midx)));

    uicontrol('parent',handles.PairedMeasuresWindow,'style','frame','units','points',...
        'position',[26  winsz(4)-((lineheight-1)*(midx+1))-1 75 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.PairedMeasuresWindow,'style','text','units','points',...
        'position',[27  winsz(4)-((lineheight-1)*(midx+1)) 73 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',handles.PairedMeasuresMatrix.ChannelName{midx});    
    
    uicontrol('parent',handles.PairedMeasuresWindow,'style','frame','units','points',...
        'position',[101  winsz(4)-((lineheight-1)*(midx+1))-1 50 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.PairedMeasuresWindow,'style','text','units','points',...
        'position',[102  winsz(4)-((lineheight-1)*(midx+1)) 48 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',num2str(handles.PairedMeasuresMatrix.x1s(midx),precision));
    uicontrol('parent',handles.PairedMeasuresWindow,'style','frame','units','points',...
        'position',[150  winsz(4)-((lineheight-1)*(midx+1))-1 50 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.PairedMeasuresWindow,'style','text','units','points',...
        'position',[151  winsz(4)-((lineheight-1)*(midx+1)) 48 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',num2str(handles.PairedMeasuresMatrix.y1s(midx),precision));

    uicontrol('parent',handles.PairedMeasuresWindow,'style','frame','units','points',...
        'position',[201  winsz(4)-((lineheight-1)*(midx+1))-1 50 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.PairedMeasuresWindow,'style','text','units','points',...
        'position',[202  winsz(4)-((lineheight-1)*(midx+1)) 48 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',num2str(handles.PairedMeasuresMatrix.x2s(midx),precision));
    uicontrol('parent',handles.PairedMeasuresWindow,'style','frame','units','points',...
        'position',[250  winsz(4)-((lineheight-1)*(midx+1))-1 50 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.PairedMeasuresWindow,'style','text','units','points',...
        'position',[251  winsz(4)-((lineheight-1)*(midx+1)) 48 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',num2str(handles.PairedMeasuresMatrix.y2s(midx),precision));    
    
    
    uicontrol('parent',handles.PairedMeasuresWindow,'style','frame','units','points',...
        'position',[301  winsz(4)-((lineheight-1)*(midx+1))-1 50 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.PairedMeasuresWindow,'style','text','units','points',...
        'position',[302  winsz(4)-((lineheight-1)*(midx+1)) 48 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',num2str(handles.PairedMeasuresMatrix.xdiffs(midx),precision));
    uicontrol('parent',handles.PairedMeasuresWindow,'style','frame','units','points',...
        'position',[350  winsz(4)-((lineheight-1)*(midx+1))-1 50 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.PairedMeasuresWindow,'style','text','units','points',...
        'position',[351  winsz(4)-((lineheight-1)*(midx+1)) 48 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',num2str(handles.PairedMeasuresMatrix.ydiffs(midx),precision));    
    
    
    uicontrol('parent',handles.PairedMeasuresWindow,'style','frame','units','points',...
        'position',[401  winsz(4)-((lineheight-1)*(midx+1))-1 119 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.PairedMeasuresWindow,'style','edit','units','points',...
        'position',[403  winsz(4)-((lineheight-1)*(midx+1)) 117 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',handles.PairedMeasuresMatrix.Comments{midx},...
        'tag',num2str(midx),'userdata',handles.displayfig,'callback',@PairedMeasuresCommentFcn);
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SingleMeasureModeFcn(obj,ev)
handles = guidata(obj);
choicelist = {'Click';'Mean';'Maximum';'Minimum'};
initialidx = strmatch(handles.SingleMeasureMode, choicelist,'exact');
[Selection,ok] = listdlg('ListString',choicelist,'SelectionMode','Single',...
    'InitialValue',initialidx,'PromptString','Choose a measurement mode');
if ~ok
    return
end
handles.SingleMeasureMode = choicelist{Selection};
set(handles.SingleMeasureModeIndic,'String',choicelist{Selection});
guidata(handles.displayfig,handles);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PairedMeasureModeFcn(obj,ev)
handles = guidata(obj);
choicelist = {'Click & Click';'Click & Mean';'Click & Maximum';'Click & Minimum';...
    'Mean & Click';'Mean & Mean';'Mean & Maximum';'Mean & Minimum';...
    'Maximum & Click';'Maximum & Mean';'Maximum & Maximum';'Maximum & Minimum';...
    'Minimum & Click';'Minimum & Mean';'Minimum & Maximum';'Minimum & Minimum'};
initialidx = strmatch(handles.PairedMeasureMode, choicelist,'exact');
[Selection,ok] = listdlg('ListString',choicelist,'SelectionMode','Single',...
    'InitialValue',initialidx,'PromptString','Choose a measurement mode');
if ~ok
    return
end
handles.PairedMeasureMode = choicelist{Selection};
set(handles.PairedMeasureModeIndic,'String',choicelist{Selection});
guidata(handles.displayfig,handles);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DeleteSingleMeasuresAll(obj,ev)
thisfig = get(get(obj,'parent'),'parent');%SingleMeasures Display Fig
handles = guidata(get(thisfig,'userdata'));

handles.SingleMeasuresMatrix.MeasureNumber = [];
handles.SingleMeasuresMatrix.FileName = {};
handles.SingleMeasuresMatrix.ChannelName = {};
handles.SingleMeasuresMatrix.xs = [];
handles.SingleMeasuresMatrix.ys = [];
handles.SingleMeasuresMatrix.MeasureMode = {};
handles.SingleMeasuresMatrix.Comments = {};

delete(handles.SingleMeasuresExportMenu);
handles.SingleMeasuresExportMenu = [];
delete(handles.SingleMeasuresDeleteMenu);
handles.SingleMeasuresDeleteMenu = [];

set(thisfig,'userdata',[]);
set(thisfig,'name','Single Measures Deleted')

guidata(handles.displayfig,handles);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DeleteSingleMeasuresSelect(obj,ev)
thisfig = get(get(obj,'parent'),'parent');%SingleMeasures Display Fig
handles = guidata(get(thisfig,'userdata'));

choices = {};
for midx = 1:size(handles.SingleMeasuresMatrix.MeasureNumber,2)
    choices{end+1} = [num2str(handles.SingleMeasuresMatrix.MeasureNumber(midx)),' ',...
        handles.SingleMeasuresMatrix.ChannelName{midx}];
end
[Selection,ok] = listdlg('ListString',choices,'PromptString','Choose measurements to delete');
if ~ok
    return
end

if length(Selection) == length(choices);%if user chose all measures... reinitialize all matrices
    handles.SingleMeasuresMatrix.MeasureNumber = [];
    handles.SingleMeasuresMatrix.FileName = {};
    handles.SingleMeasuresMatrix.ChannelName = {};
    handles.SingleMeasuresMatrix.xs = [];
    handles.SingleMeasuresMatrix.ys = [];
    handles.SingleMeasuresMatrix.MeasureMode = {};
    handles.SingleMeasuresMatrix.Comments = {};

    set(thisfig,'userdata',[]);
    set(thisfig,'name','Single Measures - Deleted')
else%else remove just selected measures
    handles.SingleMeasuresMatrix.MeasureNumber(Selection) = [];
    handles.SingleMeasuresMatrix.FileName(Selection) = [];
    handles.SingleMeasuresMatrix.ChannelName(Selection) = [];
    handles.SingleMeasuresMatrix.xs(Selection) = [];
    handles.SingleMeasuresMatrix.ys(Selection) = [];
    handles.SingleMeasuresMatrix.MeasureMode(Selection) = [];
    handles.SingleMeasuresMatrix.Comments(Selection) = [];

    set(thisfig,'userdata',[]);
    set(thisfig,'name','Single Measures - Partially Deleted')
end

delete(handles.SingleMeasuresExportMenu);
handles.SingleMeasuresExportMenu = [];
delete(handles.SingleMeasuresDeleteMenu);
handles.SingleMeasuresDeleteMenu = [];

guidata(handles.displayfig,handles);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DeletePairedMeasuresAll(obj,ev)
thisfig = get(get(obj,'parent'),'parent');%PairedMeasures Display Fig
handles = guidata(get(thisfig,'userdata'));

handles.PairedMeasuresMatrix.MeasureNumber = [];
handles.PairedMeasuresMatrix.FileName = {};
handles.PairedMeasuresMatrix.ChannelName = {};
handles.PairedMeasuresMatrix.x1s = [];
handles.PairedMeasuresMatrix.y1s = [];
handles.PairedMeasuresMatrix.x2s = [];
handles.PairedMeasuresMatrix.y2s = [];
handles.PairedMeasuresMatrix.xdiffs = [];
handles.PairedMeasuresMatrix.ydiffs = [];
handles.PairedMeasuresMatrix.MeasureMode = {};
handles.PairedMeasuresMatrix.Comments = {};

delete(handles.PairedMeasuresExportMenu);
handles.PairedMeasuresExportMenu = [];
delete(handles.PairedMeasuresDeleteMenu);
handles.PairedMeasuresDeleteMenu = [];

set(thisfig,'userdata',[]);
set(thisfig,'name','Paired Measures Deleted')

guidata(handles.displayfig,handles);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DeletePairedMeasuresSelect(obj,ev)
thisfig = get(get(obj,'parent'),'parent');%PairedMeasures Display Fig
handles = guidata(get(thisfig,'userdata'));

choices = {};
for midx = 1:size(handles.PairedMeasuresMatrix.MeasureNumber,2)
    choices{end+1} = [num2str(handles.PairedMeasuresMatrix.MeasureNumber(midx)),' ',...
        handles.PairedMeasuresMatrix.ChannelName{midx}];
end
[Selection,ok] = listdlg('ListString',choices,'PromptString','Choose measurements to delete');
if ~ok
    return
end

if length(Selection) == length(choices);%if user chose all measures... reinitialize all matrices
    handles.PairedMeasuresMatrix.MeasureNumber = [];
    handles.PairedMeasuresMatrix.FileName = {};
    handles.PairedMeasuresMatrix.ChannelName = {};
    handles.PairedMeasuresMatrix.x1s = [];
    handles.PairedMeasuresMatrix.y1s = [];
    handles.PairedMeasuresMatrix.x2s = [];
    handles.PairedMeasuresMatrix.y2s = [];
    handles.PairedMeasuresMatrix.xdiffs = [];
    handles.PairedMeasuresMatrix.ydiffs = [];
    handles.PairedMeasuresMatrix.MeasureMode = {};
    handles.PairedMeasuresMatrix.Comments = {};

    set(thisfig,'userdata',[]);
    set(thisfig,'name','Paired Measures - Deleted')
else%else remove just selected measures
    handles.PairedMeasuresMatrix.MeasureNumber(Selection) = [];
    handles.PairedMeasuresMatrix.FileName(Selection) = [];
    handles.PairedMeasuresMatrix.ChannelName(Selection) = [];
    handles.PairedMeasuresMatrix.x1s(Selection) = [];
    handles.PairedMeasuresMatrix.y1s(Selection) = [];
    handles.PairedMeasuresMatrix.x2s(Selection) = [];
    handles.PairedMeasuresMatrix.y2s(Selection) = [];
    handles.PairedMeasuresMatrix.xdiffs(Selection) = [];
    handles.PairedMeasuresMatrix.ydiffs(Selection) = [];
    handles.PairedMeasuresMatrix.Comments(Selection) = [];
    handles.PairedMeasuresMatrix.MeasureMode(Selection) = [];
    
    set(thisfig,'userdata',[]);
    set(thisfig,'name','Paired Measures - Partially Deleted')
end

delete(handles.PairedMeasuresExportMenu);
handles.PairedMeasuresExportMenu = [];
delete(handles.PairedMeasuresDeleteMenu);
handles.PairedMeasuresDeleteMenu = [];

guidata(handles.displayfig,handles);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ExportSingleMeasuresToWS(obj,ev)
thisfig = get(get(obj,'parent'),'parent');%SingleMeasures Display Fig
handles = guidata(get(thisfig,'userdata'));
outputmat = getsinglemeasuresoutputmat(handles);

assignin('base','SingleMeasures',outputmat);


%%
function ExportSingleMeasuresToMAT(obj,ev)
thisfig = get(get(obj,'parent'),'parent');%SingleMeasures Display Fig
handles = guidata(get(thisfig,'userdata'));
outputmat = getsinglemeasuresoutputmat(handles);

[FileName,PathName] = uiputfile('.mat','Select a File to Write',['SingleMeasures_',handles.filename]);
if FileName == 0 & PathName == 0
    return
end
save([PathName,FileName],'outputmat');

%%
function ExportSingleMeasuresToXLS(obj,ev)
thisfig = get(get(obj,'parent'),'parent');%SingleMeasures Display Fig
handles = guidata(get(thisfig,'userdata'));
outputmat = getsinglemeasuresoutputmat(handles);

[FileName,PathName] = uiputfile('.xls','Select a File to Write',['SingleMeasures_',handles.filename]);
if FileName == 0 & PathName == 0
    return
end
xlswrite([PathName,FileName],outputmat);

%%
function ExportSingleMeasuresToCurrentExcelBook(obj,ev);
thisfig = get(get(obj,'parent'),'parent');%SingleMeasures Display Fig
handles = guidata(get(thisfig,'userdata'));
outputmat = getsinglemeasuresoutputmat(handles);

answer = inputdlg('Enter sheet cell number to insert at');
if isempty(answer)
    return
end
answer=lower(answer{1});
letter = isstrprop(answer,'alpha');
if sum(letter) > 1
    error('Can only handle single-letter column names');
end
if letter(1) ~= 1;
    error('Please enter a single letter followed by a number');
end
numbers = isstrprop(answer,'digit');
letter = answer(letter);
numbers = answer(numbers);
colnum = letter-'a'+1;
rownum = str2double(numbers);

chan = ddeinit('excel', 'Sheet1');
for midx = 1:size(outputmat,1);
    for iidx = 1:size(outputmat,2);
        rc = ['r',num2str(rownum+midx-1),'c',num2str(colnum+iidx-1)];%start/stop cell number
        rc = [rc,':',rc];
        ddepoke(chan,rc,outputmat{midx,iidx})
    end
end
%go thru each element and ddepoke
ddeterm(chan);

%% subfunction for export functions
function outputmat = getsinglemeasuresoutputmat(handles);

for midx = 1:length(handles.SingleMeasuresMatrix.MeasureNumber);
    outputmat{midx,1} = handles.SingleMeasuresMatrix.FileName{midx};
    outputmat{midx,2} = handles.SingleMeasuresMatrix.MeasureNumber(midx);
    outputmat{midx,3} = handles.SingleMeasuresMatrix.ChannelName{midx};
    outputmat{midx,4} = handles.SingleMeasuresMatrix.xs(midx);
    outputmat{midx,5} = handles.SingleMeasuresMatrix.ys(midx);
    outputmat{midx,6} = handles.SingleMeasuresMatrix.MeasureMode{midx};
    outputmat{midx,7} = handles.SingleMeasuresMatrix.Comments{midx};
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ExportPairedMeasuresToWS(obj,ev)
thisfig = get(get(obj,'parent'),'parent');%PairedMeasures Display Fig
handles = guidata(get(thisfig,'userdata'));
outputmat = getpairedmeasuresoutputmat(handles);

assignin('base','PairedMeasures',outputmat);


%%
function ExportPairedMeasuresToMAT(obj,ev)
thisfig = get(get(obj,'parent'),'parent');%PairedMeasures Display Fig
handles = guidata(get(thisfig,'userdata'));
outputmat = getpairedmeasuresoutputmat(handles);

[FileName,PathName] = uiputfile('.mat','Select a File to Write',['PairedMeasures_',handles.filename]);
if FileName == 0 & PathName == 0
    return
end
save([PathName,FileName],'outputmat');

%%
function ExportPairedMeasuresToXLS(obj,ev)
thisfig = get(get(obj,'parent'),'parent');%PairedMeasures Display Fig
handles = guidata(get(thisfig,'userdata'));
outputmat = getpairedmeasuresoutputmat(handles);

[FileName,PathName] = uiputfile('.xls','Select a File to Write',['PairedMeasures_',handles.filename]);
if FileName == 0 & PathName == 0
    return
end
xlswrite([PathName,FileName],outputmat);

%%
function ExportPairedMeasuresToCurrentExcelBook(obj,ev);
thisfig = get(get(obj,'parent'),'parent');%SingleMeasures Display Fig
handles = guidata(get(thisfig,'userdata'));
outputmat = getpairedmeasuresoutputmat(handles);

answer = inputdlg('Enter sheet cell number to insert at');
if isempty(answer)
    return
end
answer=lower(answer{1});
letter = isstrprop(answer,'alpha');
if sum(letter) > 1
    error('Can only handle single-letter column names');
end
if letter(1) ~= 1;
    error('Please enter a single letter followed by a number');
end
numbers = isstrprop(answer,'digit');
letter = answer(letter);
numbers = answer(numbers);
colnum = letter-'a'+1;
rownum = str2double(numbers);

chan = ddeinit('excel','Sheet1');
for midx = 1:size(outputmat,1);
    for iidx = 1:size(outputmat,2);
        rc = ['r',num2str(rownum+midx-1),'c',num2str(colnum+iidx-1)];%start/stop cell number
        rc = [rc,':',rc];
        ddepoke(chan,rc,outputmat{midx,iidx});
    end
end
%go thru each element and ddepoke
ddeterm(chan);

%% subfunction for export functions
function outputmat = getpairedmeasuresoutputmat(handles);

for midx = 1:length(handles.PairedMeasuresMatrix.MeasureNumber);
    outputmat{midx,1} = handles.PairedMeasuresMatrix.FileName{midx};
    outputmat{midx,2} = handles.PairedMeasuresMatrix.MeasureNumber(midx);
    outputmat{midx,3} = handles.PairedMeasuresMatrix.ChannelName{midx};
    outputmat{midx,4} = handles.PairedMeasuresMatrix.x1s(midx);
    outputmat{midx,5} =  handles.PairedMeasuresMatrix.y1s(midx);
    outputmat{midx,6} =  handles.PairedMeasuresMatrix.x2s(midx);
    outputmat{midx,7} =  handles.PairedMeasuresMatrix.y2s(midx);
    outputmat{midx,8} =  handles.PairedMeasuresMatrix.xdiffs(midx);
    outputmat{midx,9} =  handles.PairedMeasuresMatrix.ydiffs(midx);
    outputmat{midx,10} = handles.PairedMeasuresMatrix.MeasureMode{midx};
    outputmat{midx,11} = handles.PairedMeasuresMatrix.Comments{midx};
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PeakTrigAvgFcn(obj,ev)

answer = inputdlg({'Trigger Channel (Axis Number)',...
    'Trigger Level (Triggers above this level)',...
    'Trigger Minimum Duration (ms)',...
    'Trigger Maximum Duration (ms)',...
    'Channel to Average (Axis Number.  If multiple channels separate by comma.)',...
    'Milliseconds Before Trigger to Grab','Milliseconds After Trigger to Grab',...
    'Milliseconds After Trigger To Ignore New Triggers',...
    'Start of data range to use (seconds)','End of data range to use (seconds)',...
    'Ignore Deviations Less Than Following # of Data Points'},...
    'Peak-Triggered Average',...% dialog title
    1,...%default number of lines per box
    {'1','1','0','Inf','2','10','50','0','0','Inf','0'});%default answers.  Default range is all of data
if length(answer)~=11
    return
end
if isempty(answer{1}) || isempty(answer{2}) ||...
        isempty(answer{3}) || isempty(answer{4}) || isempty(answer{5}) ||...
        isempty(answer{6}) || isempty(answer{7}) || isempty(answer{8}) ||...
        isempty(answer{9}) || isempty(answer{10})|| isempty(answer{11})
    return
end
handles = guidata(obj);
decdata=getappdata(handles.displayfig,'decdata');
allrawdata = decdata.data{1};
clear decdata

trigchan = str2double(answer{1});
abovethresh = str2double(answer{2});
mindur = str2double(answer{3});
maxdur = str2double(answer{4});
beforems = str2double(answer{6});
afterms = str2double(answer{7});
trigignoredelayms = str2double(answer{8});
starttime = str2double(answer{9});
stoptime = str2double(answer{10});
ignoredeviations = str2double(answer{11});

%if multiple channels entered to be averaged
commas = findstr(answer{5},',');
if isempty(commas)
    avgchan = str2double(answer{5});
else%hopefully resistant to comma+space separation
    for cidx = 1:length(commas);
        avgchan(cidx) = str2double(answer{5}(commas(cidx)-1));
    end
    avgchan(end+1) =  str2double(answer{5}(end));
end

beforepts = round(beforems * handles.obj.SampleRate / 1000);
afterpts = round(afterms * handles.obj.SampleRate / 1000);
trigignoredelaypts = round(trigignoredelayms * handles.obj.SampleRate / 1000);
mindur = round(mindur * handles.obj.SampleRate / 1000);
maxdur = round(maxdur * handles.obj.SampleRate / 1000);
if starttime > stoptime
    errordlg 'Start time must be less than stop time'
    return
end
startpt = starttime * handles.obj.SampleRate;
startpt = max([startpt 1]);%startpt is at least point number 1
stoppt = stoptime * handles.obj.SampleRate;
stoppt = min([stoppt size(allrawdata,1)]);%stoppt not greater than samples acquired

%get triggers
trigdata = allrawdata(:,trigchan);
trigtimes = pt_continuousabove(trigdata,zeros(size(trigdata)),abovethresh,mindur,maxdur,ignoredeviations);
trigtimes = trigtimes(:,1);
trigtimes(trigtimes > stoppt) = [];
trigtimes(trigtimes < startpt) = [];
% trigtraces = zeros(length(trigtimes),(beforepts+afterpts+1));
lastsavedtrig = -Inf;
trigtraces = [];
newtrigtimes = [];
for tidx = 1:length(trigtimes)
    if trigtimes(tidx)>(lastsavedtrig+trigignoredelaypts)
        trigtraces(end+1,:) = trigdata(trigtimes(tidx)-beforepts:trigtimes(tidx)+afterpts);
        lastsavedtrig = trigtimes(tidx);
        newtrigtimes(end+1) = lastsavedtrig;
    end
end
trigtimes = newtrigtimes;

for cidx = 1:length(avgchan);
    thischan = avgchan(cidx);
    avgdata = allrawdata(:,thischan);
    rawtraces = zeros(length(trigtimes),(beforepts+afterpts+1));
    for tidx = 1:length(trigtimes)
        rawtraces(tidx,:) = avgdata(trigtimes(tidx)-beforepts:trigtimes(tidx)+afterpts);
    end
    zeroedtraces = rawtraces - repmat(rawtraces(:,beforepts), [1 (beforepts+afterpts+1)]);

    outputfig = TriggeredAverageKernelFcn(trigtraces, rawtraces, zeroedtraces);
    trigname = handles.obj.Channel(handles.includedchannels(trigchan)).ChannelName;
    grabbedname = handles.obj.Channel(handles.includedchannels(avgchan(cidx))).ChannelName;
    figname = [handles.filename,': "',trigname,'" triggering "',grabbedname,'"'];
    set(outputfig,'name',figname);
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DipTrigAvgFcn(obj,ev);

answer = inputdlg({'Trigger Channel (Axis Number)',...
    'Trigger Level (Triggers below this level)',...
    'Trigger Minimum Duration (ms)',...
    'Trigger Maximum Duration (ms)',...
    'Channel to Average (Axis Number.  If multiple channels separate by comma.)',...
    'Milliseconds Before Trigger to Grab','Milliseconds After Trigger to Grab',...
    'Start of data range to use (seconds)','End of data range to use (seconds)',...
    'Ignore Deviations Less Than Following # of Data Points'},...
    'Dip-Triggered Average',...% dialog title
    1,...%default number of lines per box
    {'1','1','0','Inf','2','10','20','0','Inf','0'});%default answers.  Default range is all of data
if length(answer)~=10
    return
end
if isempty(answer{1}) || isempty(answer{2}) ||...
        isempty(answer{3}) || isempty(answer{4}) || isempty(answer{5}) ||...
        isempty(answer{6}) || isempty(answer{7}) || isempty(answer{8}) ||...
        isempty(answer{9}) || isempty(answer{10})
    return
end
handles = guidata(obj);
decdata=getappdata(handles.displayfig,'decdata');
allrawdata = decdata.data{1};
clear decdata

trigchan = str2double(answer{1});
belowthresh = str2double(answer{2});
mindur = str2double(answer{3});
maxdur = str2double(answer{4});
beforems = str2double(answer{6});
afterms = str2double(answer{7});
starttime = str2double(answer{8});
stoptime = str2double(answer{9});
ignoredeviations = str2double(answer{10});

%if multiple channels entered to be averaged
commas = findstr(answer{5},',');
if isempty(commas)
    avgchan = str2double(answer{5});
else
    for cidx = 1:length(commas);
        avgchan(cidx) = str2double(answer{5}(commas(cidx)-1));
    end
    avgchan(end+1) =  str2double(answer{5}(end));
end

beforepts = round(beforems * handles.obj.SampleRate / 1000);
afterpts = round(afterms * handles.obj.SampleRate / 1000);
mindur = round(mindur * handles.obj.SampleRate / 1000);
maxdur = round(maxdur * handles.obj.SampleRate / 1000);
if starttime > stoptime
    errordlg 'Start time must be less than stop time'
    return
end
startpt = starttime * handles.obj.SampleRate;
startpt = max([startpt 1]);%startpt is at least point number 1
stoppt = stoptime * handles.obj.SampleRate;
stoppt = min([stoppt size(allrawdata,1)]);%stoppt not greater than samples acquired

%get triggers
trigdata = allrawdata(:,trigchan);
trigtimes = pt_continuousbelow(trigdata,zeros(size(trigdata)),belowthresh,mindur,maxdur,ignoredeviations);
trigtimes = trigtimes(:,1);
trigtimes(trigtimes > stoppt) = [];
trigtimes(trigtimes < startpt) = [];
trigtraces = zeros(length(trigtimes),(beforepts+afterpts+1));
for tidx = 1:length(trigtimes)
    trigtraces(tidx,:) = trigdata(trigtimes(tidx)-beforepts:trigtimes(tidx)+afterpts);
end


for cidx = 1:length(avgchan);
    thischan = avgchan(cidx);
    avgdata = allrawdata(:,thischan);
    rawtraces = zeros(length(trigtimes),(beforepts+afterpts+1));
    for tidx = 1:length(trigtimes)
        rawtraces(tidx,:) = avgdata(trigtimes(tidx)-beforepts:trigtimes(tidx)+afterpts);
    end
    zeroedtraces = rawtraces - repmat(rawtraces(:,beforepts), [1 (beforepts+afterpts+1)]);

    outputfig = TriggeredAverageKernelFcn(trigtraces, rawtraces, zeroedtraces);
end