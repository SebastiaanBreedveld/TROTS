% TROTSViewPatient is a simple patient visualisation tool. It shows the CT,
%   structures, dose wash, isodose lines and beam configuration in sliced
%   axial, coronal and sagittal mode.
%
%   NOTE!   For certain Matlab versions, you need to left-click 
%           in the grey area of the figure window (next tot the
%           image) first, before you can scroll through the slices.
%
% TROTSViewPatient(patient) only shows the CT, structures and beam
%   configuration.
%
% TROTSViewPatient(patient, Dose) shows the dose overlayed on the CT
%
% TROTSViewPatient(patient, Dose[, Arg1, Val1[, Arg2, Val2...) sets
%   options:
%
%   Figure:         Use this figure number
%
%   Contours:       Is a vector of contour indices to display 
%                   corresponding to patient.StructureNames.
%
%   Isodose:        Is a vector of isodose levels to display.
%
%   Dosewash:       Boolean on whether or not to show the dose wash
%
%   ShowLabels:     Boolean on whether or not to show structure 
%                   names/isodose levels
%
%   Presentation:   Presentation mode, shows maximixed view, white CT
%                   background, thicker lines and lighter colours for
%                   contours/isodose lines, so they light up when printed
%                   in black and white.
%
%   ShowBeams:      Boolean on wheter or not to show the beam setup.
%
%
%   Keys during viewing
%
%   Scrollwheel         Move to next/previous slice
%   up/down, Arrows 
%   up/down, n/p         
%
%   Arrows left/right:  Change view: axial, coronal, sagittal
%
%   b:  Toggle beams
%
%   d:  Toggle dose wash
%
%   i:  Toggle isodose lines
%
%   j:  Set isodose levels
%
%   l:  Toggle labels
%
%   k:  Cycle through CT colourmaps (just for fun)
%
%   f:  Toggle fullscreen
%
%   q/ESC:  exit viewer
%   
%
% Literature:
%
%   Breedveld S. & Heijmen B.
%   Data for TROTS - The Radiotherapy Optimisation Test Set
%   Data in Brief (2017) 12:143-149. doi: 10.1016/j.dib.2017.03.037
%
%   Breedveld S., Van den Berg B. & Heijmen B. 
%   An interior-point implementation developed and tuned for radiation 
%   therapy treatment planning 
%   Comput. Optim. Appl. (2017) 68:209-242. doi: 10.1007/s10589-017-9919-4
%
%   Breedveld S., Craft D., van Haveren R. & Heijmen B.
%   Multi-criteria Optimisation and Decision Making in Radiotherapy
%   Eur. J. Oper. Res. (2018) Accepted. doi: 10.1016/j.ejor.2018.08.019
%
% Website:
%
%   http://www.erasmusmc.nl/radiotherapytrots/
%
%
% (c) 2016/2018 Sebastiaan Breedveld / Erasmus University Medical Center
%   Rotterdam, The Netherlands
%
function TROTSViewPatient(patient, Dose, varargin)

% Check input
ArgumentNames = {'Figure', 'Isodose', 'Dosewash', 'Contours', 'ShowLabels', 'Presentation', 'ShowBeams', 'MaxDose'};


if isempty(patient)
    error('I need at least a Patient as input');
end
Beams = patient.Beams;

% Process other input
show_cont = [];
Figure = [];
Isodose = [];
Dosewash = 0;
Contours = [];
ShowLabels = 1;
Presentation = 0;
ViewMode = 'axial';
ContourIdx = 1;
ShowBeams = 1;
MaxDose = 0;
process_input(ArgumentNames, varargin{:});


% Initialize arrays
axhandle = 0;
icth = 0;
idoseh = 0;
cb = 0;
DoseIso = 0;

% Figure
if isempty(Figure)
    Figure = 3;
end

figure(Figure)
clf(Figure, 'reset');    
shg


% Dose wash and/or isodose
if isempty(Isodose) && isequal(Dosewash, 0)
    Dosewash = 1;
end
if nargin==1 || isempty(Dose)
    Dosewash = 0;
    Isodose = 0;
    Dose = [];
    idoseh = -1;
end

% Contours
if isempty(show_cont) || ~isempty(Contours)
    show_cont = Contours;
end
if isempty(show_cont)
    show_cont = 1:size(patient.Contours{1}, 2);    
end
show_cont = intersect(1:size(patient.Contours{1}, 2), show_cont);


% Reindex CT for visualisation
CTVisualise = patient.CT;
CTVisualise = patient.CT + 1024;
CTVisualise = CTVisualise/50;
CTVisualise(CTVisualise>63) = 63;
CTVisualise(CTVisualise<0) = 0;
if Presentation
    % Make a white background
    CTVisualise(patient.CT==-1024) = 63;
end



% First process dose matrix input and determine maximum dose
if ~isempty(Dose) && isempty(MaxDose) || MaxDose<=0 
    
    MaxDose = double(max(Dose(:)));
    
    if MaxDose > 200
        % Maximum dose is ridiculously high, so probably a brachytherapy
        % dose. Set maximum  to V99.8% of the dose
        Dsrt = sort(Dose(:), 'ascend');
        MaxDose = Dsrt(floor(0.998*length(Dsrt)));
    end
end
DoseInterval = [0 MaxDose];

% Toggle isodose
if isempty(Isodose) || isequal(Isodose, 0)
    ShowIsodose = 0;
    Isodose = 1;
else
    ShowIsodose = 1;
end

axhandle = gca;
icth = image(CTVisualise(:,:,1)');
    
axis image

if ~Presentation
    % Make nice ticks
    ts = 50; % tick 50 mm
    Corner1 = patient.Offset;
    Corner2 = (size(patient.CT)-1).*patient.Resolution + patient.Offset;
    Mid = patient.Isocentre(1,:);

    Start = ceil(Corner1(1)/ts)*ts;
    End = floor(Corner2(1)/ts)*ts;
    RangeX = Start:ts:End;
    % Remove ticks close to Mid to prevent cluttering of ticks
    RangeX(abs(RangeX-Mid(1))<25) = [];
    RangeX = sort([RangeX Mid(1)]);
    RangeXVox = (RangeX-patient.Offset(1))/patient.Resolution(1)+1;

    Start = ceil(Corner1(2)/ts)*ts;
    End = floor(Corner2(2)/ts)*ts;
    RangeY = Start:ts:End;
    % Remove ticks close to Mid to prevent cluttering of ticks
    RangeY(abs(RangeY-Mid(2))<10) = [];
    RangeY = sort([RangeY Mid(2)]);
    RangeYVox = (RangeY-patient.Offset(2))/patient.Resolution(2)+1;

    Start = ceil(Corner1(3)/ts)*ts;
    End = floor(Corner2(3)/ts)*ts;
    RangeZ = Start:ts:End;
    % Remove ticks close to Mid to prevent cluttering of ticks
    RangeZ(abs(RangeZ-Mid(3))<10) = [];
    RangeZ = sort([RangeZ Mid(3)]);
    RangeZVox = (RangeZ-patient.Offset(3))/patient.Resolution(3)+1;                
else
    axis off               
end
hold on

if ~isempty(Dose)         
    DoseOrg = Dose;
    % Make separate isodose array
    if ~isempty(Isodose)
        DoseIso = smooth3(Dose, 'gaussian', [3 3 1]);        
        DoseIso = DoseIso/MaxDose*64+65;
    end

    % Scale dose properly in intervals
    Dose = max(0, ((Dose-DoseInterval(1)))/(DoseInterval(2)-DoseInterval(1))*64)+65;
       
    idoseh = image(0, 'Parent', axhandle);
    alpha(idoseh, 0.4);

    % Show colour bar
    if (DoseInterval(2)-DoseInterval(1))>20
        Indices = sort(unique([DoseInterval(1):10:DoseInterval(2) DoseInterval(2)]));
    else
        Indices = sort(unique(round([DoseInterval(1):2:DoseInterval(2) DoseInterval(2)])));
    end
    caxis(axhandle, [65 129]);
    axes(axhandle);
    cb = colorbar('EastOutside', 'YLim', [64.5 128.5], 'YTick', ((Indices-DoseInterval(1))/(DoseInterval(2)-DoseInterval(1)))*64+64.5, 'YTickLabel', Indices, 'ButtonDownFcn', '', 'Visible', 'off');
    if ShowLabels
        set(cb, 'Visible', 'on')
        drawnow
    else
        set(cb, 'Visible', 'off')
    end
else
    DoseInterval = [0 0];
    Isodose = [];
end

% Set colour map, and plant easter eggs for normal user (don't bother
% clinical mode)
cmapidx = 0;
if (isequal(datestr(now, 'dd/mm'), '24/12') && str2double(datestr(now, 'HH'))>18) || isequal(datestr(now, 'dd/mm'), '25/12') || isequal(datestr(now, 'dd/mm'), '26/12')
    % Christmas
    ChangeColourMap(1);
    uiwait(msgbox('Because it is Christmas, and you are working, the CT is now in Christmas colours! (press k to toggle)', 'Christmas!', 'help'));
elseif (isequal(datestr(now, 'dd/mm'), '31/12') && str2double(datestr(now, 'HH'))>20) || isequal(datestr(now, 'dd/mm'), '01/01')
    % New Year
    ChangeColourMap(2);
    uiwait(msgbox('Because it is New Year''s day, and you are working, the CT is now in Evergreen colours! (press k to toggle)'));
elseif isequal(datestr(now, 'dd/mm'), '02/08')  
    % Sebas
    ChangeColourMap(0);
    uiwait(msgbox('Because it is Sebastiaan''s birthday, let''s go and remember him that today he is one year older', 'Birthday!', 'help'));
elseif (weekday(now)==7 && isequal(datestr(now, 'dd/mm/yy'), '26/04/14')) || isequal(datestr(now, 'dd/mm'), '27/04')
    % King's Day 
    ChangeColourMap(7);
    uiwait(msgbox('Because it is King''s Day in the Netherlands (Long Live the King!), and you are working, the CT is now in Royal colours! (press k to toggle)', 'King''s Day!', 'help'));
else
    ChangeColourMap(0);
end

% Fix a nice figure, but separate for normal and presentation mode
if Presentation
    Positioning = 'Position';
    % Background colour
    set(Figure, 'Color', [1 1 1]);
else
    Positioning = 'OuterPosition';
end

[~, SliceIdxX] = min(abs((1:size(patient.CT, 1))-((patient.Isocentre(1)-patient.Offset(1))/patient.Resolution(1)+1)));
[~, SliceIdxY] = min(abs((1:size(patient.CT, 2))-((patient.Isocentre(2)-patient.Offset(2))/patient.Resolution(2)+1)));
[~, SliceIdxZ] = min(abs((1:size(patient.CT, 3))-((patient.Isocentre(3)-patient.Offset(3))/patient.Resolution(3)+1)));

if ~isempty(Isodose)
    if isequal(Isodose, 1)
        Isodose = round(linspace(1, MaxDose-MaxDose/10, 10));
    elseif length(Isodose)==1
        Isodose = round(linspace(1, MaxDose-MaxDose/Isodose, Isodose));
    end
    Isodose = double(sort(Isodose));
    LevelsInImage = Isodose/MaxDose*64+65;   
end

% Set colours
DoseColours = jet(64);
ContColours = 1-jet(length(show_cont));
if Presentation
    % Turn up the saturation for B/W prints
    DoseColours = rgb2hsv(DoseColours);
    DoseColours(:, 2) = DoseColours(:, 2)/2;
    DoseColours = hsv2rgb(DoseColours);

    ContColours = rgb2hsv(ContColours);
    ContColours(:, 2) = ContColours(:, 2)/2;
    ContColours = hsv2rgb(ContColours);
end

% Label axes
labelhandle = 0;
contnameh = zeros(1, length(show_cont));
isodosenameh = zeros(1, length(Isodose));

labelhandle = axes('Position', get(axhandle, 'Position'), 'ActivePositionProperty', 'Position', 'Visible', 'off');

CurrentAxis = [get(labelhandle, 'XLim') get(labelhandle, 'YLim')];
for w=1:length(show_cont)
    XPos = 0.7*CurrentAxis(2);
    YPos = 1-0.03*CurrentAxis(4)*w;
    contnameh(w) = text(XPos,YPos,num2str(patient.StructureNames{show_cont(w)}), 'Parent', labelhandle, 'Color', ContColours(w,:), 'FontWeight', 'bold', 'Visible', 'off');
end

if idoseh~=-1
    Isodose(Isodose>MaxDose) = [];
    Isodose(Isodose<0) = [];
    for w=1:length(Isodose)
        XPos = 0.2*CurrentAxis(2);
        YPos = 1-0.05*CurrentAxis(4)*(length(LevelsInImage)-w+1);
        zcidx = max(1, round(Isodose(w)/MaxDose*size(DoseColours, 1)));
        isodosenameh(w) = text(XPos,YPos,num2str(Isodose(w)), 'FontSize', 16, 'Parent', labelhandle, 'Color', DoseColours(zcidx,:));
        if ShowIsodose && ShowLabels
            set(isodosenameh(w), 'Visible', 'on');
        else
            set(isodosenameh(w), 'Visible', 'off');
        end
    end
end

% Reserve appropriate size (can grow for extreme number of contours)
Level = zeros(100, 1);
NumberOfPoints = zeros(100, 1);
ContourDose = cell(100, 1);
    
% Init indices
if isequal(ViewMode, 'axial')
    SliceIdx = SliceIdxZ;
elseif isequal(ViewMode, 'coronal')
    SliceIdx = SliceIdxY;
elseif isequal(ViewMode, 'sagittal')
    SliceIdx = SliceIdxX;
end
ChangeView



% So, that's the init, make it interactive
Amount = 0; % Scroll amount

% Draw once and set callbacks
ScrollSlice    
set(Figure, 'WindowButtonDownFcn', @ScrollSlice)
set(Figure, 'WindowScrollWheelFcn', @ScrollSlice)
set(Figure, 'WindowKeyPressFcn', @ScrollSlice)
if ~Presentation
    set(Figure, 'WindowButtonMotionFcn', @UpdateCoordinates)
    set(Figure, 'Pointer', 'cross')
end


% Make figure active to regain focus
% See: http://www.mathworks.com/matlabcentral/newsreader/view_thread/235825
figure(Figure);
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
% Add try catch block beacuse these calls fail if Matlab has no GUI window.
try
    javaFrame = get(Figure,'JavaFrame');
    javaFrame.getAxisComponent.requestFocus;
catch
end
return


% Here are the different embedded subfunctions

    % Set different colour map. This is just for fun, why should the world
    % be grey? ;)
    function ChangeColourMap(ThisOne)
        if nargin==0
            cmapidx = mod(cmapidx+1,8);
        else
            cmapidx = ThisOne;
        end
        
        if cmapidx==0
            % Default
            CTMap = gray(64);
        elseif cmapidx==1
            % Red
            CTMap = zeros(64,3);
            CTMap(1:12, 1) = linspace(0, 0.25, 12)';
            CTMap(13:64, 1) = 1; CTMap(:, 2) = linspace(0, 1, 64)'; CTMap(:, 3) = CTMap(:, 2);
        elseif cmapidx==2
            % Green
            CTMap = zeros(64,3);
            CTMap(1:12, 2) = linspace(0, 0.25, 12)';
            CTMap(13:64, 2) = 1; CTMap(:, 1) = linspace(0, 1, 64)'; CTMap(:, 3) = CTMap(:, 1);
        elseif cmapidx==3
            % Blue
            CTMap = zeros(64,3);
            CTMap(1:12, 3) = linspace(0, 0.25, 12)';
            CTMap(13:64, 3) = 1; CTMap(:, 1) = linspace(0, 1, 64)'; CTMap(:, 2) = CTMap(:, 1);
        elseif cmapidx==4
            % Yellow
            CTMap = zeros(64,3);
            CTMap(1:12, 1) = linspace(0, 0.25, 12)'; CTMap(1:12, 2) = linspace(0, 0.25, 12)';
            CTMap(13:64, 1:2) = 1; CTMap(:, 3) = linspace(0, 1, 64)';
        elseif cmapidx==5
            % Purple
            CTMap = zeros(64,3);
            CTMap(1:12, 1) = linspace(0, 0.25, 12)'; CTMap(1:12, 3) = linspace(0, 0.25, 12)';
            CTMap(13:64, [1 3]) = 1; CTMap(:, 2) = linspace(0, 1, 64)';
        elseif cmapidx==6
            % ?
            CTMap = zeros(64,3);
            CTMap(1:12, 2) = linspace(0, 0.25, 12)'; CTMap(1:12, 3) = linspace(0, 0.25, 12)';
            CTMap(13:64, 2:3) = 1; CTMap(:, 1) = linspace(0, 1, 64)';
        elseif cmapidx==7
            % Orange
            CTMap = zeros(64,3);
            CTMap(1:12, 1) = linspace(0, 0.25, 12)';
            CTMap(13:64, 1) = 1; CTMap(:, 2) = linspace(0, 1, 64)'; CTMap(:, 3) = 0;    
        end
        
        % Set map
        if isempty(Dose)
            % CT only
            colormap(CTMap);
        else
            DoseCMap = jet(64);            
            if Presentation
                colormap([CTMap; 1 1 1; DoseCMap]);
            else
                colormap([CTMap; 0 0 0; DoseCMap]);
            end
        end

        if (DoseInterval(2)-DoseInterval(1))>20
            Indices = sort(unique([DoseInterval(1):10:DoseInterval(2) DoseInterval(2)]));
        else
            Indices = sort(unique(round([DoseInterval(1):2:DoseInterval(2) DoseInterval(2)])));
        end
        if cb~=0
            if ShowLabels
                set(cb, 'YLim', [64.5 128.5], 'YTick', ((Indices-DoseInterval(1))/(DoseInterval(2)-DoseInterval(1)))*64+64.5, 'YTickLabel', Indices, 'ButtonDownFcn', '', 'Visible', 'on');
            else
                set(cb, 'YLim', [64.5 128.5], 'YTick', ((Indices-DoseInterval(1))/(DoseInterval(2)-DoseInterval(1)))*64+64.5, 'YTickLabel', Indices, 'ButtonDownFcn', '', 'Visible', 'off');
            end
        end
    end



    % Changes or set current view mode
    function ChangeView(how)
        if nargin==1
            if isequal(ViewMode, 'axial')
                SliceIdxZ = SliceIdx;
                if how==1
                    ViewMode = 'coronal';
                    SliceIdx = SliceIdxY;
                else
                    ViewMode = 'sagittal';
                    SliceIdx = SliceIdxX;
                end                
            elseif isequal(ViewMode, 'coronal')
                SliceIdxY = SliceIdx;
                if how==1
                    ViewMode = 'sagittal';
                    SliceIdx = SliceIdxX;
                else
                    ViewMode = 'axial';
                    SliceIdx = SliceIdxZ;
                end                
            elseif isequal(ViewMode, 'sagittal')
                SliceIdxX = SliceIdx;
                if how==1
                    ViewMode = 'axial';
                    SliceIdx = SliceIdxZ;
                else
                    ViewMode = 'coronal';
                    SliceIdx = SliceIdxY;
                end                
            end
        end    
            
            
        
        if isequal(ViewMode, 'axial')
            % Contour set to use
            ContourIdx = 1;
            
            % Set proper aspect ratio and axes
            set(axhandle, 'DataAspectRatio', patient.Resolution([2 1 3]), 'XDir', 'normal', 'YDir', 'reverse')
            xlim(axhandle, [0 size(patient.CT, 1)+1])
            ylim(axhandle, [0 size(patient.CT, 2)+1])
            if ~Presentation
                set(axhandle, 'XTick', RangeXVox, 'XTickLabel', RangeX);
                set(axhandle, 'YTick', RangeYVox, 'YTickLabel', RangeY);
                xlabel('right-left (mm)')
                ylabel('dorsal-ventral (mm)')
            end
        elseif isequal(ViewMode, 'coronal')
            % Contour set to use
            ContourIdx = 2;
            
            % Set proper aspect ratio
            set(axhandle, 'DataAspectRatio', patient.Resolution([3 1 2]), 'XDir', 'normal', 'YDir', 'normal')
            xlim(axhandle, [0 size(patient.CT, 1)+1])
            ylim(axhandle, [0 size(patient.CT, 3)+1])
            if ~Presentation
                set(axhandle, 'XTick', RangeXVox, 'XTickLabel', RangeX);
                set(axhandle, 'YTick', RangeZVox, 'YTickLabel', RangeZ);
                xlabel('right-left (mm)')
                ylabel('caudal-cranial (mm)')
            end
        elseif isequal(ViewMode, 'sagittal')
            % Contour set to use
            ContourIdx = 3;
            
            % Set proper aspect ratio
            set(axhandle, 'DataAspectRatio', patient.Resolution([3 2 1]), 'XDir', 'reverse', 'YDir', 'normal');
            xlim(axhandle, [0 size(patient.CT, 2)+1])
            ylim(axhandle, [0 size(patient.CT, 3)+1])
            if ~Presentation
                set(axhandle, 'XTick', RangeYVox, 'XTickLabel', RangeY);
                set(axhandle, 'YTick', RangeZVox, 'YTickLabel', RangeZ);            
                xlabel('dorsal-ventral (mm)')
                ylabel('caudal-cranial (mm)')
            end
        end
        
        % Get ratio of figure
        AR = get(axhandle(1), 'DataAspectRatio'); 
        Ratio = [get(axhandle(1), 'XLim')*AR(2); get(axhandle(1), 'YLim')*AR(1)];
        Ratio = Ratio(:, 2) - Ratio(:, 1);
        Ratio = Ratio(2)/Ratio(1);

        colidx = 0;
        rowidx = 0;

        set(axhandle, 'Units', 'normalized');
        set(axhandle, Positioning, [rowidx 1-(1+colidx) 1 1])


        % Fix figure
        ScreenSize = get(0, 'ScreenSize');
        Pos = get(Figure, 'Position');
        if Presentation
            HeightF = 512;
            WidthF = HeightF/Ratio;
        else
            HeightF = 0;
            WidthF = 0;
        end
        
        if Pos(1)+WidthF>ScreenSize(3)
            Pos(1) = max(0, ScreenSize(3)-WidthF);
        end
        if Pos(2)+HeightF+128>ScreenSize(4)
            Pos(2) = ScreenSize(4)-HeightF-128;
        end

        % Leave figure window when unset
        if WidthF~=0 && HeightF~=0
            set(Figure, 'Position', [Pos(1:2) WidthF HeightF], 'PaperPositionMode', 'auto');    
        end
        set(Figure, 'Units', 'pixels');

        return
    end

    % Display a new slice and update all information
    function ScrollSlice(FigSrc, event)            
        % Get event
        if nargin>0
            if isfield(event, 'VerticalScrollCount') || (isequal(class(event), 'matlab.ui.eventdata.ScrollWheelData') && isequal(event.EventName, 'WindowScrollWheel'))
                % Scroll wheel - move in this direction
                Amount = -event.VerticalScrollCount;
            elseif isfield(event, 'Key') || (isequal(class(event), 'matlab.ui.eventdata.KeyData') && isequal(event.EventName, 'WindowKeyPress'))
                % Key press
                if isequal(event.Key, 'uparrow') || isequal(event.Key, 'n')
                    Amount = 1;
                elseif isequal(event.Key, 'downarrow') || isequal(event.Key, 'p')
                    Amount = -1;
                elseif isequal(event.Key, 'leftarrow')
                    Amount = 0;
                    ChangeView(-1);
                    ScrollSlice;
                    return
                elseif isequal(event.Key, 'rightarrow')                    
                    Amount = 0;
                    ChangeView(1);
                    ScrollSlice;
                    return
                elseif isequal(event.Key, 'escape') || isequal(event.Key, 'q')
                    % If we wait, resume first without closing the window
                    close(FigSrc);                    
                    return
                elseif isequal(event.Key, 'b')
                    ShowBeams = 1-ShowBeams;
                    Amount = 0;
                elseif isequal(event.Key, 'd')
                    Dosewash = 1-Dosewash;
                    Amount = 0;
                    if Dosewash==0 && ~isempty(Dose)
                        set(idoseh, 'CData', 0);
                    end
                elseif isequal(event.Key, 'f')
                    jFig = get(handle(Figure), 'JavaFrame');
                    jFig.setMaximized(1-jFig.isMaximized);                                  
                elseif isequal(event.Key, 'i')
                    if ~isempty(Dose)
                        ShowIsodose = 1-ShowIsodose;
                    else
                        ShowIsodose = 0;
                    end
                    
                    Amount = 0;
                    if ~ShowIsodose || ~ShowLabels
                        for z=1:length(Isodose)
                            set(isodosenameh(z), 'Visible', 'off');
                        end                        
                    elseif ShowIsodose && ShowLabels
                        for z=1:length(Isodose)
                            set(isodosenameh(z), 'Visible', 'on');
                        end
                    end
                elseif isequal(event.Key, 'j')
                    % Set isodose levels
                    NewLevels = inputdlg('Isodose levels in Gy:', 'Isodose levels');
                    if ~isempty(NewLevels)
                        NewLevels = str2num(NewLevels{1});
                        for z=1:length(Isodose)                            
                            delete(isodosenameh(z));
                        end
                        Isodose = sort(NewLevels);
                        Isodose(Isodose>MaxDose) = [];
                        Isodose(Isodose<0) = [];
                        
                        LevelsInImage = Isodose/MaxDose*64+65;
                        CurrentAxis = [get(labelhandle, 'XLim') get(labelhandle, 'YLim')];
                        for z=1:length(Isodose)
                            XPos = 0.2*CurrentAxis(2);
                            YPos = 1-0.05*CurrentAxis(4)*(length(LevelsInImage)-z+1);
                            cidx = round(Isodose(z)/MaxDose*size(DoseColours, 1));
                            isodosenameh(z) = text(XPos,YPos,num2str(Isodose(z)), 'FontSize', 16, 'Parent', labelhandle, 'Color', DoseColours(cidx,:));
                            if ShowIsodose && ShowLabels
                                set(isodosenameh(z), 'Visible', 'on');                                    
                            else
                                set(isodosenameh(z), 'Visible', 'off');
                            end
                        end
                    end                    
                elseif isequal(event.Key, 'k')
                    % Change colour map
                    ChangeColourMap;
                    Amount = 0;
                elseif isequal(event.Key, 'l')
                    ShowLabels = 1-ShowLabels;
                    Amount = 0;
                    if ShowLabels
                        if ShowIsodose
                            for z=1:length(Isodose)
                                set(isodosenameh(z), 'Visible', 'on');
                            end
                        end                        
                        set(cb, 'Visible', 'on');
                    else
                        if ShowIsodose
                            for z=1:length(Isodose)
                                set(isodosenameh(z), 'Visible', 'off');
                            end
                        end
                        set(cb, 'Visible', 'off');
                    end
                else
                    return;
                end                   
            elseif isempty(event)
                if strcmp(get(FigSrc, 'SelectionType'), 'normal')
                    Amount = -1;
                elseif strcmp(get(FigSrc, 'SelectionType'), 'alt')
                    Amount = 1;
                elseif strcmp(get(FigSrc, 'SelectionType'), 'open')
                    % do previous amount
                else
                    return;
                end
            else
                return;
            end
            SliceIdx = mod(SliceIdx-1+Amount, size(patient.CT, 4-ContourIdx))+1;
        end


        % Set Slice
        j = SliceIdx;

        % Remove old contour data
        Childs = get(axhandle, 'Children');
        for z=1:length(Childs)
            if isequal(get(Childs(z), 'Type'), 'line')
                delete(Childs(z));
            end
        end

        if isequal(ViewMode, 'axial')
            if j<=size(patient.CT, 3)
                set(icth, 'CData', CTVisualise(:,:,j)');
            end
        elseif isequal(ViewMode, 'coronal')
            if j<=size(patient.CT, 2)
                set(icth, 'CData', squeeze(CTVisualise(:,j,:))');
            end
        elseif isequal(ViewMode, 'sagittal')
            if j<=size(patient.CT, 1)
                set(icth, 'CData', squeeze(CTVisualise(j,:,:))');
            end
        end                       
            
        if ~isempty(Dose)
            if isequal(ViewMode, 'axial')
                if j<=size(patient.CT, 3)
                    DoseDisplay = Dose(:,:,j)';
                end
            elseif isequal(ViewMode, 'coronal')
                if j<=size(patient.CT, 2)
                    DoseDisplay = squeeze(Dose(:,j,:))';
                end
            elseif isequal(ViewMode, 'sagittal')
                if j<=size(patient.CT, 1)
                    DoseDisplay = squeeze(Dose(j,:,:))';
                end
            end
            if Dosewash
                set(idoseh, 'CData', DoseDisplay);                
            end
            if ShowIsodose
                if any(DoseDisplay(:)>Isodose(1))
                    clear ContData
                    if isequal(ViewMode, 'axial')
                        DoseDisplayIso = DoseIso(:,:,j)';
                    elseif isequal(ViewMode, 'coronal')
                        DoseDisplayIso = squeeze(DoseIso(:,j,:))';
                    elseif isequal(ViewMode, 'sagittal')
                        DoseDisplayIso = squeeze(DoseIso(j,:,:))';
                    end
                    ContData = contourc(double(DoseDisplayIso), double(LevelsInImage));
                    teller=1;
                    NumContours = 0;
                    while teller < size(ContData,2)
                        NumContours=NumContours+1;
                        Level(NumContours) = ContData(1,teller);
                        NumberOfPoints(NumContours) = ContData(2,teller);
                        ContourDose{NumContours} = ContData(:,teller+1:teller+NumberOfPoints(NumContours));
                        teller = teller+NumberOfPoints(NumContours)+1;
                    end                    

                    if NumContours>0
                        for iii=1:NumContours
                            % Get closest (rounding off errors!)
                            LineColour = DoseColours(round(Level(iii)-65),:);
                            plot(ContourDose{iii}(1,:),ContourDose{iii}(2,:),'Color',LineColour,'LineStyle','-', 'LineWidth', 1*(Presentation+1), 'Parent', axhandle)
                        end
                    end
                end
            end
        end
        if ~Presentation                
            title(axhandle, {sprintf('%s (%s), Slice %d at %f', patient.Identifier, patient.PatientPosition, j, patient.Offset(4-ContourIdx)+(j-1)*patient.Resolution(4-ContourIdx)) ''})
        end

        if ~isempty(show_cont)
            teller3=0;

            for p=1:length(show_cont)
                z = show_cont(p);
                teller3=teller3+1;
                if j<=size(patient.Contours{ContourIdx}, 1) && z<=size(patient.Contours{ContourIdx}, 2)
                    ContourData = patient.Contours{ContourIdx}{j,z};                        
                else
                    ContourData = [];
                end

                if ~isempty(ContourData)                       
                    for l = 1:length(ContourData)                        
                        if isequal(ViewMode, 'axial')                                
                            Xnodes = [ContourData{l}(:,1); ContourData{l}(1,1)];
                            Ynodes = [ContourData{l}(:,2); ContourData{l}(1,2)];
                            line((Xnodes-patient.Offset(1))/patient.Resolution(1)+1.0,... 
                                (Ynodes-patient.Offset(2))/patient.Resolution(2)+1.0, 'Color', ContColours(teller3,:), 'LineWidth', 2*(Presentation+1), 'Parent', axhandle);                                
                        elseif isequal(ViewMode, 'coronal')
                            Xnodes = [ContourData{l}(:,1); ContourData{l}(1,1)];
                            Znodes = [ContourData{l}(:,2); ContourData{l}(1,2)];
                            line((Xnodes-patient.Offset(1))/patient.Resolution(1)+1.0,...
                                (Znodes-patient.Offset(3))/patient.Resolution(3)+1.0, 'Color', ContColours(teller3,:), 'LineWidth', 2*(Presentation+1), 'Parent', axhandle);                                                               
                        elseif isequal(ViewMode, 'sagittal')
                            Ynodes = [ContourData{l}(:,1); ContourData{l}(1,1)];
                            Znodes = [ContourData{l}(:,2); ContourData{l}(1,2)];
                            line((Ynodes-patient.Offset(2))/patient.Resolution(2)+1.0,...
                                (Znodes-patient.Offset(3))/patient.Resolution(3)+1.0, 'Color', ContColours(teller3,:), 'LineWidth', 2*(Presentation+1), 'Parent', axhandle);
                        end
                    end                
                    if ShowLabels
                        set(contnameh(p), 'Visible', 'on');
                    else
                        set(contnameh(p), 'Visible', 'off');
                    end
                else
                    set(contnameh(p), 'Visible', 'off');                    
                end
            end
        end

        % Show Isocentre
        if isequal(ViewMode, 'axial')
            if abs(j-(patient.Isocentre(3)-patient.Offset(3))/patient.Resolution(3)-1)<1
                plot((patient.Isocentre(1)-patient.Offset(1))/patient.Resolution(1)+1, (patient.Isocentre(2)-patient.Offset(2))/patient.Resolution(2)+1, 'Marker', 'x', 'MarkerSize', 20, 'Color', [1 1 1], 'LineWidth', 3, 'Parent', axhandle);
            else
                plot((patient.Isocentre(1)-patient.Offset(1))/patient.Resolution(1)+1, (patient.Isocentre(2)-patient.Offset(2))/patient.Resolution(2)+1, 'Marker', 'x', 'MarkerSize', 10, 'Color', [1 1 1], 'Parent', axhandle);
            end
        elseif isequal(ViewMode, 'coronal')
            if abs(j-(patient.Isocentre(2)-patient.Offset(2))/patient.Resolution(2)-1)<1
                plot((patient.Isocentre(1)-patient.Offset(1))/patient.Resolution(1)+1, (patient.Isocentre(3)-patient.Offset(3))/patient.Resolution(3)+1, 'Marker', 'x', 'MarkerSize', 20, 'Color', [1 1 1], 'LineWidth', 3, 'Parent', axhandle);
            else
                plot((patient.Isocentre(1)-patient.Offset(1))/patient.Resolution(1)+1, (patient.Isocentre(3)-patient.Offset(3))/patient.Resolution(3)+1, 'Marker', 'x', 'MarkerSize', 10, 'Color', [1 1 1], 'Parent', axhandle);
            end
        elseif isequal(ViewMode, 'sagittal')
            if abs(j-(patient.Isocentre(1)-patient.Offset(1))/patient.Resolution(1)-1)<1
                plot((patient.Isocentre(2)-patient.Offset(2))/patient.Resolution(2)+1, (patient.Isocentre(3)-patient.Offset(3))/patient.Resolution(3)+1, 'Marker', 'x', 'MarkerSize', 20, 'Color', [1 1 1], 'LineWidth', 3, 'Parent', axhandle);
            else
                plot((patient.Isocentre(2)-patient.Offset(2))/patient.Resolution(2)+1, (patient.Isocentre(3)-patient.Offset(3))/patient.Resolution(3)+1, 'Marker', 'x', 'MarkerSize', 10, 'Color', [1 1 1], 'Parent', axhandle);
            end
        end                

        % Draw beam setup
        if ShowBeams && ~isempty(Beams)                
            if Presentation
                rgb = ones(Beams.Num, 3)*0.8;
            else
                hsv = ones(Beams.Num,3);
                hsv(:,1) = ([1:Beams.Num]'-1)/Beams.Num;
                hsv(:,2) = 0.5;
                hsv(:,3) = 1;
                rgb = hsv2rgb(hsv);
            end

            cnt = 1;
            for beamidx = 1:Beams.Num
                if isfield(Beams.BeamConfig(beamidx), 'Gantry')
                    % For external beam therapy

                    
                    % Compute spatial point
                    BeamPoints = [sind(Beams.BeamConfig(beamidx).Gantry)*cosd(Beams.BeamConfig(beamidx).Couch) cosd(Beams.BeamConfig(beamidx).Gantry) sind(Beams.BeamConfig(beamidx).Gantry)*sind(Beams.BeamConfig(beamidx).Couch)];
                    if isequal(ViewMode, 'axial')
                        line((patient.Isocentre(1)-patient.Offset(1))/patient.Resolution(1) + 1 + BeamPoints(1)/patient.Resolution(1)*[0 250], ...
                             (patient.Isocentre(2)-patient.Offset(2))/patient.Resolution(2) + 1 - BeamPoints(2)/patient.Resolution(2)*[0 250], ...
                             (patient.Isocentre(3)-patient.Offset(3))/patient.Resolution(3) + 1 + BeamPoints(3)/patient.Resolution(3)*[0 250], 'LineWidth',1.5*(Presentation+1),'Color',rgb(cnt,:), 'Parent', axhandle)
                    elseif isequal(ViewMode, 'coronal')
                        line((patient.Isocentre(1)-patient.Offset(1))/patient.Resolution(1) + 1 + BeamPoints(1)/patient.Resolution(1)*[0 250], ...
                             (patient.Isocentre(3)-patient.Offset(3))/patient.Resolution(3) + 1 + BeamPoints(3)/patient.Resolution(3)*[0 250], ...
                             (patient.Isocentre(2)-patient.Offset(2))/patient.Resolution(2) + 1 + BeamPoints(2)/patient.Resolution(2)*[0 250], 'LineWidth',1.5*(Presentation+1),'Color',rgb(cnt,:), 'Parent', axhandle)
                    elseif isequal(ViewMode, 'sagittal')
                        line((patient.Isocentre(2)-patient.Offset(2))/patient.Resolution(2) + 1 - BeamPoints(2)/patient.Resolution(2)*[0 250], ...
                             (patient.Isocentre(3)-patient.Offset(3))/patient.Resolution(3) + 1 + BeamPoints(3)/patient.Resolution(3)*[0 250], ...
                             (patient.Isocentre(1)-patient.Offset(1))/patient.Resolution(1) + 1 + BeamPoints(1)/patient.Resolution(1)*[0 250], 'LineWidth',1.5*(Presentation+1),'Color',rgb(cnt,:), 'Parent', axhandle)
                    end
                    cnt = cnt + 1;
                elseif isfield(Beams.BeamConfig(beamidx), 'DwellPositions')
                    % For brachytherapy
                    
                    for dwidx=1:size(Beams.BeamConfig(beamidx).DwellPositions, 1)                        
                        if isequal(ViewMode,'axial')
                            if abs(j-(Beams.BeamConfig(beamidx).DwellPositions(dwidx,3)-patient.Offset(3))/patient.Resolution(3) - 1)<1
                                plot((Beams.BeamConfig(beamidx).DwellPositions(dwidx,1)-patient.Offset(1))/patient.Resolution(1) + 1,...
                                     (Beams.BeamConfig(beamidx).DwellPositions(dwidx,2)-patient.Offset(2))/patient.Resolution(2) + 1,'Marker','o','MarkerSize',4,'Color',[0 1 0],'Parent', axhandle);
                            end
                        elseif isequal(ViewMode,'coronal')
                            if abs(j-(Beams.BeamConfig(beamidx).DwellPositions(dwidx,2)-patient.Offset(2))/patient.Resolution(2) - 1)<1
                                plot((Beams.BeamConfig(beamidx).DwellPositions(dwidx,1)-patient.Offset(1))/patient.Resolution(1) + 1,...
                                     (Beams.BeamConfig(beamidx).DwellPositions(dwidx,3)-patient.Offset(3))/patient.Resolution(3) + 1,'Marker','o','MarkerSize',4,'Color',[0 1 0],'Parent', axhandle);
                            end
                        elseif isequal(ViewMode,'sagittal')
                            if abs(j-(Beams.BeamConfig(beamidx).DwellPositions(dwidx,1)-patient.Offset(1))/patient.Resolution(1) - 1)<1
                                plot((Beams.BeamConfig(beamidx).DwellPositions(dwidx,2)-patient.Offset(2))/patient.Resolution(2) + 1,...
                                     (Beams.BeamConfig(beamidx).DwellPositions(dwidx,3)-patient.Offset(3))/patient.Resolution(3) + 1,'Marker','o','MarkerSize',4,'Color',[0 1 0],'Parent', axhandle);
                            end
                        end                  
                    end
                end
            end
        end
           
        if ~Presentation
            UpdateCoordinates
        end
        drawnow   
    end % End ScrollSlice function

    function UpdateCoordinates(~, ~)
        % First, get active axes        
        CurFigP = get(Figure, 'CurrentPoint');
        
        AxHandle = gca;
        set(axhandle, 'Units', 'pixels')
        OP = get(axhandle, 'OuterPosition');
        set(axhandle, 'Units', 'normalized')
        if all(CurFigP>=OP(1:2)) && all(CurFigP<=(OP(1:2)+OP(3:4)))
            AxHandle = axhandle;
        end
                
        cp = get(AxHandle, 'CurrentPoint');
        if isequal(ViewMode, 'axial')
            x = cp(1,1);
            y = cp(1,2);
            k = SliceIdx;
        elseif isequal(ViewMode, 'coronal')
            x = cp(1,1);
            y = SliceIdx;
            k = cp(1,2);
        elseif isequal(ViewMode, 'sagittal')
            x = SliceIdx;
            y = cp(1,1);
            k = cp(1,2);
        end
        SpatialPoint = ([x y k] - 1).*patient.Resolution + patient.Offset;
            
        x = round(x);
        y = round(y);
        k = round(k);
        OrgTitle = get(get(axhandle, 'Title'), 'String');            
                
        if ~isempty(Dose) && ~isempty(DoseOrg)                
            CurSize = size(DoseOrg);
        else
            CurSize = size(patient.CT);
        end
            
        if length(CurSize)==2
            CurSize = [CurSize 1];
        end
        if (all([x y k]>=[1 1 1]) && (all([x y k]<=CurSize)))
            if ~isempty(Dose)                    
                title(axhandle, {OrgTitle{1} sprintf('(%f , %f , %f) mm, %d HU, %.2f Gy', SpatialPoint, round(patient.CT(x,y,k)), DoseOrg(x,y,k))});
            else
                title(axhandle, {OrgTitle{1} sprintf('(%f , %f , %f) mm, %d HU', SpatialPoint, round(patient.CT(x,y,k)))});
            end
        else
            title(axhandle, {OrgTitle{1} sprintf('(%f , %f , %f) mm', SpatialPoint)});
        end
    end 


    function j = process_input(accept, varargin)
        for j=1:2:nargin-1
            if ~ischar(varargin{j})
                return
            else
                if j+1<=nargin-1
                    if isempty(accept)
                        assignin('caller', varargin{j}, varargin{j+1});
                    else

                        % Make case-insensitive match
                        idx = strmatch(lower(varargin{j}), lower(accept), 'exact');
                        if ~isempty(idx)
                            if ~isequal(accept{idx}, varargin{j})
                                fprintf('Case insensitive match: got ''%s'' in stead of ''%s''. Please pay attention to the case when supplying options.\n', varargin{j}, accept{idx});
                            end

                            assignin('caller', accept{idx}, varargin{j+1});
                        else
                            fprintf('Ignoring unknown option ''%s'' and returning.\n', varargin{j});
                            return
                        end
                    end
                else
                    return
                end
            end
        end
    end
end % End view_patient function
