% The code expects patientFolder to contain a TROTS *.mat file
% See e.g. https://zenodo.org/records/2708302/files/Protons.zip?download=1
% To open the result, if you did not close matRad, just open matRadGUI and
% click Refresh. No need to touch nothing else.
% You can uncomment the lines saving mini.mat or full.mat depending on how
% much information you want the MAT file to contain.
% If you already closed MATLAB, open it again, start matRadGUI, and then
% click on the LoadMat button, and select mini.mat or full.mat, not the
% original TROTS *.mat file
% Authored by:
% S. Tattenberg - TRIUMF and NOSM
% F. Hueso-González - IFIC (CSIC/UV)
% Assistance:
% N. Wahl, T. Becher, A. Hammi
% See https://github.com/e0404/matRad/issues/695

clear, clc, close all
patientFolder = '/tmp/'; % with TROTS mat file
TrotsMatFile = patientFolder + "Protons_01.mat";
load(TrotsMatFile);

%TROTS has: data, patient, problem, problem_lex, solutionX
%matRad needs: cst, ct, patientFolder, pln, stf, dij

%% Define CT structure.
ct.resolution.x = patient.Resolution(1);
ct.resolution.y = patient.Resolution(2);
ct.resolution.z = patient.Resolution(3);
nRows =  size(patient.CT, 2); % DICOM Rows, goes with y, is index 2 because of how matrix is stored
nColumns =  size(patient.CT, 1);% DICOM Columns, goes with x, is index 1 because of how matrix is stored
nSlices =  size(patient.CT, 3);% DICOM slices, goes with z
ct.x = linspace(patient.Offset(1), patient.Offset(1)+(nColumns-1)*patient.Resolution(1), nColumns);
ct.y = linspace(patient.Offset(2), patient.Offset(2)+(nRows-1)*patient.Resolution(2), nRows);
ct.z = linspace(patient.Offset(3), patient.Offset(3)+(nSlices-1)*patient.Resolution(3), nSlices);
ct.cubeDim = [nRows nColumns nSlices];
ct.numOfCtScen = 1;
ct.timeStamp = string(datetime("now"));
ct.cubeHU = {double(permute(patient.CT,[2 1 3]))};

%% Create cst structure
% see https://github.com/e0404/matRad/blob/master/dicom/matRad_createCst.m
% and https://github.com/e0404/matRad/blob/master/dicom/matRad_dummyCst.m
% and https://github.com/e0404/matRad/wiki/The-cst-cell
nStructures = length(patient.StructureNames);
defaultColors = colorcube(nStructures);
cst = cell(nStructures*2,6);
[grx,gry] = ndgrid(ct.x,ct.y);
disp('Calculating contours')
for i = 1:nStructures
    cst{i,1} = i-1;
    cst{i,2} = patient.StructureNames{i};
    if contains(cst{i,2}, 'ctv', 'IgnoreCase', true) ...
    || contains(cst{i,2}, 'ptv', 'IgnoreCase', true) 
        cst{i,3}          = 'TARGET';
    else
        cst{i,3}          = 'OAR';
    end
    
    linvoxs = [];
    disp(cst{i,2})
    for s = 1:nSlices
        csXY = patient.Contours{1,1}{s,i};
        for sc = 1:length(csXY)
            contoursXY = csXY{sc};
            in = inpolygon(grx,gry,contoursXY(:,1),contoursXY(:,2)).';
            ind = find(in) + (s-1)*nRows*nColumns;
            linvoxs = cat(1, linvoxs, ind);
        end
        if ~isempty(csXY)
            disp(s)
        end
    end
    cst{i,4}{1}       = linvoxs;
    if strcmp(cst{i,3}, 'OAR')
        cst{i,5}.Priority = 2;
    else
        cst{i,5}.Priority = 1;
    end
    cst{i,5}.alphaX   = 0.1;
    cst{i,5}.betaX    = 0.05;
    cst{i,5}.Visible  = 1;
    cst{i,5}.visibleColor = defaultColors(i,:);
    cst{i,6}          = [];
end
% Define now sparse structures for optimization
for i = 1:nStructures
    cst{i+nStructures,1} = i+nStructures-1;
    cst{i+nStructures,2} = ['sp', patient.StructureNames{i}];
    if contains(cst{i+nStructures,2}, 'ctv', 'IgnoreCase', true) ...
    || contains(cst{i+nStructures,2}, 'ptv', 'IgnoreCase', true) 
        cst{i+nStructures,3}          = 'TARGET';
    else
        cst{i+nStructures,3}          = 'OAR';
    end
    sampVox = cell2mat(patient.SampledVoxels(i));
    linvoxs = sub2ind(ct.cubeDim, sampVox(2,:), sampVox(1,:), sampVox(3,:));
    disp(cst{i+nStructures,2})
    cst{i+nStructures,4}{1}       = linvoxs;
    if strcmp(cst{i+nStructures,3}, 'OAR')
        cst{i+nStructures,5}.Priority = 2;
    else
        cst{i+nStructures,5}.Priority = 1;
    end
    cst{i+nStructures,5}.alphaX   = 0.1;
    cst{i+nStructures,5}.betaX    = 0.05;
    cst{i+nStructures,5}.Visible  = 1;
    cst{i+nStructures,5}.visibleColor = defaultColors(i,:);
    cst{i+nStructures,6}          = [];
end

%% Define constraints and objectives in cst

for i=1:nStructures
    OARIndex = i + nStructures;
    totalIndices = size(problem, 2);
    objectiveIndex = 1;
    for index=1:totalIndices
        if contains(cst(OARIndex,2),problem(index).Name)
            cst{OARIndex,6}{objectiveIndex} = struct();
            if problem(index).IsConstraint == 1
                cst{OARIndex,6}{objectiveIndex}.className = 'DoseConstraints.matRad_MinMaxDose';
                if problem(index).Minimise == 1
                    cst{OARIndex,6}{objectiveIndex}.parameters = cell({0,problem(index).Objective,1}); % 1 is approx, 2 is voxel-wise
                else
                    cst{OARIndex,6}{objectiveIndex}.parameters = cell({problem(index).Objective,100,1});
                end
                cst{OARIndex,6}{objectiveIndex}.epsilon = 1.0000e-03;
            else
                if problem(index).Minimise == 1
                    cst{OARIndex,6}{objectiveIndex}.className = 'DoseObjectives.matRad_SquaredOverdosing';
                else
                   cst{OARIndex,6}{objectiveIndex}.className = 'DoseObjectives.matRad_SquaredUnderdosing';
                end
                cst{OARIndex,6}{objectiveIndex}.parameters = cell({problem(index).Objective});
                cst{OARIndex,6}{objectiveIndex}.penalty = problem(index).Weight;
            end
            objectiveIndex = objectiveIndex + 1;
        end
    end
end

%% clear and save
%clearvars -except cst ct patientFolder
%save([patientFolder 'mini.mat'], '-v7.3')

%% Definition of pln
pln = struct;
pln.propStf = struct;
pln.propStf.bixelWidth = ct.resolution.x;
numBeams = size(patient.Beams.BeamConfig,2);
pln.propStf.gantryAngles = [];
for beamIndex=1:numBeams
    pln.propStf.gantryAngles(beamIndex) = patient.Beams.BeamConfig(beamIndex).Gantry;
end
numCouches= size(patient.Beams.BeamConfig,2);
pln.propStf.couchAngles = [];
for couchIndex=1:numCouches
    pln.propStf.couchAngles(couchIndex) = patient.Beams.BeamConfig(couchIndex).Couch;
end
pln.propStf.numOfBeams = numBeams;
for i=1:numBeams
    img_isox = patient.Isocentre(1) - ct.x(1);% matrad subtracts offset
    img_isoy = patient.Isocentre(2) - ct.y(1);% matrad subtracts offset
    img_isoz = patient.Isocentre(3) - ct.z(1);% matrad subtracts offset
    pln.propStf.isoCenter(i,:) = [img_isox img_isoy img_isoz]; % not sure what this is
end
pln.voxelDimensions = ct.cubeDim;
pln.numOfVoxels = prod(pln.voxelDimensions);
pln.numOfFractions = 1;
if contains(TrotsMatFile, 'Proton', 'IgnoreCase', true)
    pln.radiationMode = 'protons';
    pln.propOpt.bioOptimization = 'const_RBExD';
else
    pln.radiationMode = 'photons';
    pln.propOpt.bioOptimization = 'none';
end
pln.machine = 'Generic';
pln.propOpt.runDAO = 0;
pln.propOpt.runSequencing = 0;

%% Pass patient TROTS Dij to MATLAB

spots = size(solutionX, 1);

dij.doseGrid.resolution.x = ct.resolution.x;
dij.doseGrid.resolution.y = ct.resolution.y;
dij.doseGrid.resolution.z = ct.resolution.z;
dij.doseGrid.x = ct.x;
dij.doseGrid.y = ct.y;
dij.doseGrid.z = ct.z;
dij.doseGrid.dimensions  = ct.cubeDim;
dij.doseGrid.numOfVoxels = prod(dij.doseGrid.dimensions);

dij.ctGrid.resolution.x = ct.resolution.x;
dij.ctGrid.resolution.y = ct.resolution.y;
dij.ctGrid.resolution.z = ct.resolution.z;
dij.ctGrid.x = ct.x;
dij.ctGrid.y = ct.y;
dij.ctGrid.z = ct.z;
dij.ctGrid.dimensions  = ct.cubeDim;
dij.ctGrid.numOfVoxels = prod(dij.ctGrid.dimensions);

dij.numOfBeams         = pln.propStf.numOfBeams;
dij.numOfScenarios     = 1; % for the moment we exclude the 9 scenarios TROTS
dij.numOfRaysPerBeam   = patient.Beams.ElementIndex;
dij.totalNumOfBixels   = spots;% sum([stf(:).totalNumOfBixels]);
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);
if contains(patientFolder,'Proton', 'IgnoreCase', true)
    dij.RBE = 1.1;
end
dij.bixelNum = ones(spots,1);
dij.rayNum   = [];
for r=1:length(dij.numOfRaysPerBeam)
    dij.rayNum = [dij.rayNum; (1:1:dij.numOfRaysPerBeam(r))'];
end
dij.beamNum   = [];
for b=1:dij.numOfBeams
    dij.beamNum = [dij.beamNum; ones(patient.Beams.ElementIndex(b),1)*b];
end
% Allocate space for dij.physicalDose sparse matrix
doseNames = struct2cell(data.matrix);
for i = 1:dij.numOfScenarios
    dij.physicalDose{i} = spalloc(dij.doseGrid.numOfVoxels,spots,1);
    
    for structIdx = 1:nStructures
        disp(patient.StructureNames(structIdx));
        sampVox = cell2mat(patient.SampledVoxels(structIdx));
        svoxels = size(sampVox, 2);
        if length(find(strcmpi(doseNames(1,:),patient.StructureNames(structIdx)))) > 1
            disp('Warning: Multiple Dij for one struct not yet supported, taking first one')
        end
        doseIdx = find(strcmpi(doseNames(1,:),patient.StructureNames(structIdx)), 1);%TODO handle double ones!
        struct_dij = data.matrix(doseIdx).A;
        voxels = size(struct_dij,1);
        if data.matrix(doseIdx).b ~= 0
            disp('Plan recalculation not yet supported')
            return
        end
        if ~isempty(data.matrix(doseIdx).c)
            disp('Quadratic cost functions not yet supported')
            return
        end
        if mod(voxels,9) == 0 && voxels/9 == svoxels
            disp('Warning: Multiple scenarios not yet supported')
            struct_dij = struct_dij(1:9:end,:);
            voxels = size(struct_dij,1);
        end
        if voxels ~= svoxels
            disp('Wrong number of voxels')
            return
        end
        sInd = sub2ind(ct.cubeDim, sampVox(2,:), sampVox(1,:), sampVox(3,:));
        dij.physicalDose{1}(sInd,:) = struct_dij;
    end
end

%% Optimize

%pln.propOpt.optimizer = 'fmincon';
[resultGUI, optimizer] = matRad_fluenceOptimization(dij, cst, pln);
solutionM = resultGUI.wUnsequenced;
dvh = matRad_calcDVH(cst(nStructures+1:end,:), resultGUI.physicalDose);

figure;
for i = 1:length(dvh)
    plot(dvh(i).doseGrid.', dvh(i).volumePoints.', 'color', defaultColors(i,:));
    hold on;
end
hold off;
matRadGUI;
qi  = matRad_calcQualityIndicators(cst, pln, resultGUI.physicalDose);

%clearvars -except cst ct dij pln resultGUI patientFolder
%save([patientFolder 'full.mat'], '-v7.3')
