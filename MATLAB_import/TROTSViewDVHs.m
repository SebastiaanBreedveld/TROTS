% TROTSViewDVHs displays the dose-volume histogram
%
% DVH = TROTSViewDVHs(problem, data[, IDX][, solutionX[, solutionX2[, ...]]]) 
%   displays the DVH for PROBLEM and DATA for solution vector X. 
%   Optionally, multiple solutions can be given, which will be plotted 
%   simultaneously.
%
%   Optionally, IDX is a vector of pointers to structure PROBLEM for
%   structures to be displayed.
%
%   If IDX is not provided, heuristics are used to plot the relevant data
%   in order of importance. This is:
%   1) All tumours (PTV/CTV/PZ)
%   2) All objectives for real clinical structures (organs), in order of priority
%   3) All else
%
%   Excluded are mean dose structures and smoothing criteria. In step 2,
%   artificial structures like Rings, Shells, etc. are also excluded and
%   saved for step 3.
%
%   If > 5 solutions are provided, the INTERACTIVE view mode is enabled,
%   where hoovering the mouse pointer over the DVHs will elicitate all
%   curves belonging to that solution/treatment plan.
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
function varargout = TROTSViewDVHs(problem, data, varargin)

% First variable argument could be the indices of the problem for which
% DVHs are to be constructed
if length(varargin{1})<50 && all(varargin{1}==round(varargin{1}))
    ShowDVHIdx = varargin{1};
    if size(ShowDVHIdx, 2) == 1
        ShowDVHIdx = ShowDVHIdx';
    end
    VararginOffset = 1;
    NumberOfPlans = length(varargin)-1;
else
    ShowDVHIdx = [];
    VararginOffset = 0;
    NumberOfPlans = length(varargin);
end
    

% Parameter settings
% SampleResolution: sample per 0.1 Gy point
SampleResolution = 0.1;

% Set colourmap, line and marker styles
Colour = [   0         0    1.0000
             0    0.5000         0
        1.0000         0         0
             0    0.7500    0.7500
        0.7500         0    0.7500
        0.7500    0.7500         0
        0.2500    0.2500    0.2500
             0    1.0000         0
             0    1.0000    1.0000
        0.0431    0.5176    0.7804
        0.1686    0.5059    0.3373
        0.4784    0.0627    0.8941
        0.6824    0.4667         0
        0.7020    0.7804    1.0000
        0.8706    0.4902         0
        1.0000         0    1.0000
        1.0000    0.6000    0.7843
        1.0000    0.6941    0.3922
        1.0000    1.0000         0];
LineStyles = {'-','--',':','-.','-'};

% Interactive mode - if there are many plans to visualise, allow to
% highlight plans by hoovering the mouse over them.
if NumberOfPlans>5
    Interactive = 1;
    LineStyles = {'-'};
else
    Interactive = 0;
end

% Automated "guessing" of the order of representation
% Order is: 
% - PTV(s)
% - Objectives in order of priority, excluding shells
% - Remainder exluding smoothing data

DataIdx = [];
VolumeNames = cell(0,0);

if isempty(ShowDVHIdx)
    % Get PTVs
    for k=1:length(problem)
        if problem(k).Active && ~ismember(problem(k).dataID, DataIdx) && data.matrix(problem(k).dataID).Type==0
            Name = problem(k).Name;
            if (strncmp(Name, 'PTV', 3) || strncmp(Name, 'CTV', 3) || strncmp(Name, 'PZ', 2) || problem(k).Type==4 || problem(k).Minimise==0) && isempty(strfind(Name, 'Ring')) && isempty(strfind(Name, 'Shell')) && isempty(strfind(Name, 'Max'))
                if size(data.matrix(problem(k).dataID).A, 1)>9 % We check >9 here to account for robust proton plans - no structure should have only 9 voxels.
                    DataIdx = [DataIdx problem(k).dataID];
                    VolumeNames{length(DataIdx)} = Name;
                else
                    % If it is only a mean, check if there is another matrix
                    % present somewhere
                    for z=1:length(data.matrix)
                        if isequal(data.matrix(z).Name, Name) && size(data.matrix(z).A, 1)>9 && ~ismember(z, DataIdx)
                            DataIdx = [DataIdx z];
                            VolumeNames{length(DataIdx)} = Name;
                            break
                        end
                    end
                end
            end
        end
    end
    % Get Objectives
    MaxPrio = 0;
    for k=1:length(problem)
        if problem(k).Active
            MaxPrio = max(MaxPrio, problem(k).Priority);
        end
    end
    for Priority = 1:MaxPrio
        for k=1:length(problem)
            if problem(k).Active && ~ismember(problem(k).dataID, DataIdx) && data.matrix(problem(k).dataID).Type==0
                Name = problem(k).Name;
                if isequal(problem(k).Priority, Priority) && isempty(strfind(Name, 'Ring')) && isempty(strfind(Name, 'Shell')) && isempty(strfind(Name, 'Max'))
                    if size(data.matrix(problem(k).dataID).A, 1)>9 % We check >9 here to account for robust proton plans - no structure should have only 9 voxels.
                        DataIdx = [DataIdx problem(k).dataID];
                        VolumeNames{length(DataIdx)} = Name;
                    else
                        % If it is only a mean, check if there is another matrix
                        % present somewhere
                        for z=1:length(data.matrix)
                            if isequal(data.matrix(z).Name, Name) && size(data.matrix(z).A, 1)>9 && ~ismember(z, DataIdx)
                                DataIdx = [DataIdx z];
                                VolumeNames{length(DataIdx)} = Name;
                                break
                            end
                        end
                    end
                end
            end
        end
    end
    % Get remainder
    for k=1:length(problem)
        if problem(k).Active && ~ismember(problem(k).dataID, DataIdx) && data.matrix(problem(k).dataID).Type==0
            if size(data.matrix(problem(k).dataID).A, 1)>9 % We check >9 here to account for robust proton plans - no structure should have only 9 voxels.
                DataIdx = [DataIdx problem(k).dataID];
                VolumeNames{length(DataIdx)} = problem(k).Name;
            else
                % If it is only a mean, check if there is another matrix
                % present somewhere
                for z=1:length(data.matrix)
                    if isequal(data.matrix(z).Name, problem(k).Name) && size(data.matrix(z).A, 1)>9 && ~ismember(z, DataIdx)
                        DataIdx = [DataIdx z];
                        VolumeNames{length(DataIdx)} = problem(k).Name;
                        break
                    end
                end
            end
        end
    end
else
    % Ok, show what is requested, but still have a sane interpetation
    for k = ShowDVHIdx
        if ~ismember(problem(k).dataID, DataIdx) && data.matrix(problem(k).dataID).Type==0
            if size(data.matrix(problem(k).dataID).A, 1)>9 % We check >9 here to account for robust proton plans - no structure should have only 9 voxels.
                DataIdx = [DataIdx problem(k).dataID];
                VolumeNames{length(DataIdx)} = problem(k).Name;
            else
                % If it is only a mean, check if there is another matrix
                % present somewhere
                for z=1:length(data.matrix)
                    if isequal(data.matrix(z).Name, problem(k).Name) && size(data.matrix(z).A, 1)>9 && ~ismember(z, DataIdx)
                        DataIdx = [DataIdx z];
                        VolumeNames{length(DataIdx)} = problem(k).Name;
                        break
                    end
                end
            end
        end
    end
end

% Here we compute the DVHs

% Keep track of maximum dose of all structures
MaxDose = 0;

% Compute the dose for each found structure
NumberOfStructs = length(DataIdx);
DVH = cell(NumberOfPlans, NumberOfStructs);
for p = 1:NumberOfPlans
    x = varargin{p+VararginOffset};
    for j = 1:length(DataIdx)
        DoseForThisVolume = data.matrix(DataIdx(j)).A*x + data.matrix(DataIdx(j)).b;

        % Make a DVH
        NrVox = length(DoseForThisVolume);
        DVH{p, j} = zeros(ceil(max(DoseForThisVolume)/SampleResolution+1), 2);
        idx = 0;
        for z=0:SampleResolution:ceil(max(DoseForThisVolume)+1)
            idx = idx + 1;
            DVH{p, j}(idx, 1) = z;
            DVH{p, j}(idx, 2) = nnz(DoseForThisVolume>=z)/NrVox;
        end
        
        maxForThisVolume = max(DoseForThisVolume);
        if maxForThisVolume > 200
            % Maximum dose is ridiculously high, so probably a brachytherapy
            % dose. Set maximum  to V90% of the dose
            Dsrt = sort(DoseForThisVolume, 'ascend');
            maxForThisVolume = Dsrt(floor(0.90*length(Dsrt)));
        end
        MaxDose = max(MaxDose, maxForThisVolume);
    end
end
    

% Plot all DVHs
FigNum = figure;
set(FigNum, 'Color', [1 1 1])
DVHAxes = gca;
xlabel(DVHAxes, 'Dose (Gy)', 'FontSize', 16);
ylabel(DVHAxes, 'Volume (%)', 'FontSize', 16);
box on
MaxDose = ceil(MaxDose/5)*5;
set(DVHAxes, 'XLim', [0 MaxDose], 'YLim', [0 1.005], 'XTick', 0:5:MaxDose, 'XTickLabel', 0:5:MaxDose, 'YTick', 0:.1:1, 'YTickLabel', 0:10:100, 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 16);
hold on
% This a trick to ensure that the most important structures lie on top of
% the qraph, and still create an understandable legend
for Idx = 1:NumberOfStructs
    plot(DVHAxes, -1, -1, 'Visible', 'on', 'Color', Colour(mod(Idx-1, size(Colour, 1))+1, :), 'LineStyle', '-', 'LineWidth', 4, 'Marker', 'none');
end
if NumberOfPlans>1 && ~Interactive
    plot(DVHAxes, -1, -1, 'Visible', 'on', 'Color', [1 1 1], 'LineStyle', '-', 'LineWidth', 0.1, 'Marker', 'none');
    VolumeNames{end+1} = '';
    for LIdx = 1:NumberOfPlans
        StyleIdx = mod(LIdx - 1, length(LineStyles)) + 1;
        plot(DVHAxes, -1, -1, 'Visible', 'on', 'Color', [0 0 0], 'LineStyle', LineStyles{StyleIdx}, 'LineWidth', 4, 'Marker', 'none');
        VolumeNames{end+1} = sprintf('Plan %d', LIdx);
    end
end

% Create DVHHandle stucture as to set different linewidths later on in
% interactive mode.
DVHHandle = zeros(NumberOfPlans, NumberOfStructs);
for LIdx = 1:NumberOfPlans
    StyleIdx = mod(LIdx - 1, length(LineStyles)) + 1;
    for Idx = NumberOfStructs:-1:1
        % New plot DVHs
        DVHHandle(LIdx, Idx) = plot(DVHAxes, DVH{LIdx, Idx}(:, 1), DVH{LIdx, Idx}(:, 2), 'Color', Colour(mod(Idx-1, size(Colour, 1))+1, :), 'LineStyle', LineStyles{StyleIdx}, 'LineWidth', 4);
    end
end
legend(VolumeNames, 'Location', 'EastOutside', 'FontSize', 14)

% Check the numerical output is requested
if nargout
    varargout{1} = DVH;
end

% Interactive session - toggle by mouse click
if Interactive
    show_plan_int([], [], 1);
    title(DVHAxes, 'Use mouse button to toggle interactive session')
    MouseButton = 0;
    set(FigNum, 'WindowButtonDownFcn', @interactive_toggle)
    CurrentPlan = 1;
end

function interactive_toggle(~, ~)
    % Toggle mouse button
    if MouseButton
        MouseButton = 0;
        set(FigNum, 'WindowButtonMotionFcn', '');
    else
        MouseButton = 1;
        set(FigNum, 'WindowButtonMotionFcn', @pick_plan);
    end
end

% Find plan that is closest to pointer
function pick_plan(src, event)

    cp = get(DVHAxes, 'CurrentPoint');
    Point = cp(1, 1:2);

    % Get closest match
    MinDist = zeros(NumberOfPlans, NumberOfStructs);
    for iter=1:NumberOfPlans
        for struct=1:NumberOfStructs
            Dist = sqrt(sum((bsxfun(@minus, DVH{iter, struct}, Point).^2), 2));
            MinDist(iter, struct) = min(Dist);
        end
    end

    % Find minimum distance by plan and structure
    [MinDist2, ClosestPoint] = min(MinDist);
    [~, ClosestStructure] = min(MinDist2);
    ClosestPlan = ClosestPoint(ClosestStructure);

    % Show plan
    if CurrentPlan~=ClosestPlan            
        CurrentPlan = ClosestPlan;        
        show_plan_int(src, event, ClosestPlan);        
    end
end

function show_plan_int(~, ~, DisplayPlan)
    if DisplayPlan<=NumberOfPlans
        title(DVHAxes, sprintf('Showing plan %d', DisplayPlan));

        for iter=1:NumberOfPlans
            if isequal(DisplayPlan, iter)
                LineW = 4;
            else
                LineW = 0.5;
            end
            for struct=1:NumberOfStructs
                 set(DVHHandle(iter, struct), 'LineWidth', LineW);
            end
        end
    end
end

end
