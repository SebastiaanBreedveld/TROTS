% TROTSComputeDose computes dose based on pencil-beam matrices
%
% Dose = TROTSComputeDose(solutionX, patient, data) computes the dose based 
%   on the pencil-beam weight vector xopt. The patient structure contains
%   information on which voxels are related to the dose resulting from the
%   individual matrices in data. A linear interpolation is done to assigin
%   dose to the entire CT.
%
% WARNING!!! CAUTION!!! WARNING!!! CAUTION!!! WARNING!!! CAUTION!!!
%
%        Since this dose is computed based on undersampled data
%        (especially outside the delineated structures), this dose
%        is ONLY an approximation to get an idea of what the dose
%        looks like.
%
%        There are deviations in target coverage (which is generally
%        lower in this reconstruction than in reality).
%
%        Nevertheless, the correspondence with the real dose is
%        still pretty awesome, just not at the volumes with 
%        high dose gradients.
%
% WARNING!!! CAUTION!!! WARNING!!! CAUTION!!! WARNING!!! CAUTION!!!
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
function Dose = TROTSComputeDose(solutionX, patient, data)

% Ensure people read help
fprintf('Reconstructing dose. (BTW, did you read the warning in: help TROTSComputeDose ?)\n');

% Because girddata only interpolates inside a convex hull, and we want to
% avoid NaNs, add the points given by to DoseBox as 0 Gy dose. This makes
% the dose near the end of the body contour more realistic.
Edges(:, 1) = double(patient.DoseBox(:, 1) - 1);
Edges(:, 2) = double(patient.DoseBox(:, 2) + 1);
p = [Edges(1, 1) Edges(2, 1) Edges(3, 1);
     Edges(1, 2) Edges(2, 1) Edges(3, 1);
     Edges(1, 1) Edges(2, 2) Edges(3, 1);
     Edges(1, 1) Edges(2, 1) Edges(3, 2);
     Edges(1, 2) Edges(2, 2) Edges(3, 1);
     Edges(1, 2) Edges(2, 1) Edges(3, 2);
     Edges(1, 1) Edges(2, 2) Edges(3, 2);
     Edges(1, 2) Edges(2, 2) Edges(3, 2)];
p = bsxfun(@times, p - 1, patient.Resolution);
SampledPoints = bsxfun(@plus, p, patient.Offset);
SampledDose = zeros(size(p, 1), 1);

% This is more elaborate, fills all 6 sides of the dose-cube. Seems not be
% necessary.
%
% p = [];
% for i=round(linspace(Edges(1,1), Edges(1,2), 10))
%     for j=round(linspace(Edges(2,1), Edges(2,2), 10))
%         p = [p;
%              i j Edges(3,1);
%              i j Edges(3,2)];
%     end
% end
% for i=round(linspace(Edges(1,1), Edges(1,2), 10))
%     for j=round(linspace(Edges(3,1), Edges(3,2), 10))
%         p = [p;
%              i Edges(2,1) j;
%              i Edges(2,2) j];
%     end
% end
% for i=round(linspace(Edges(2,1), Edges(2,2), 10))
%     for j=round(linspace(Edges(3,1), Edges(3,2), 10))
%         p = [p;
%              Edges(1,1) i j;
%              Edges(1,2) i j];
%     end
% end


% Loop over the data matrices, and find structure in patient, and see if
% they match
for MatId = 1:length(data.matrix)
    % Ignore mean-only doses
    if size(data.matrix(MatId).A, 1)>1 && ~any(data.matrix(MatId).A(:)<0)
        % Find structure in patient
        StructIdx = 0;
        Robust = 0;
        for StructId = 1:length(patient.StructureNames)
            if strcmp(data.matrix(MatId).Name, patient.StructureNames{StructId})
                if size(data.matrix(MatId).A, 1)==size(patient.SampledVoxels{StructId}, 2)
                    % Match!
                    StructIdx = StructId;
                    break
                elseif mod(size(data.matrix(MatId).A, 1), 9)==0 && size(data.matrix(MatId).A, 1)/9==size(patient.SampledVoxels{StructId}, 2)
                    % This is (probably) a robust matrix
                    Robust = 1;
                    StructIdx = StructId;
                    break
                else
                    fprintf('Hmm, matrix (%d) dimensions of structure %s (%d) do not match between pencil-beam matrix and sampled voxel indices. Something fishy is going on!\n', MatId, patient.StructureNames{StructId}, StructId);
                end
            end
        end
        
        % Recored sampled voxels and compute dose for this structure
        if StructIdx>0
%             fprintf('Matched structure %s!\n', patient.StructureNames{StructIdx});
            
            % Dose in sample points
            if ~Robust
                SampledDosePart = data.matrix(MatId).A*solutionX + data.matrix(MatId).b;
            else
                % Take the first of the nine scenarios. First scenario is
                % the nominal one.
                SampledDosePart = data.matrix(MatId).A(1:size(data.matrix(MatId).A, 1)/9, :)*solutionX + data.matrix(MatId).b;
            end
            
            % Transform sampled points to spatial coordinates
            p = bsxfun(@times, double(patient.SampledVoxels{StructIdx}' - 1), patient.Resolution);
            SampledPointsPart = bsxfun(@plus, p, patient.Offset);
    
            % Concatenate
            SampledPoints = [SampledPoints; SampledPointsPart];
            SampledDose = [SampledDose; SampledDosePart];
        end
    end
end

% Find the indices for the voxels where the dose should be computed
CTMask = zeros(size(patient.CT), 'uint8');
CTMask(patient.DoseBox(1,1):patient.DoseBox(1,2), patient.DoseBox(2,1):patient.DoseBox(2,2), patient.DoseBox(3,1):patient.DoseBox(3,2)) = 1;
% Intersect the dose box with the patient external
% Note that a CT value of -1024 is the minimum of the DICOM standard, and
% it is assumed here that -1024 is outside the patient (which is enforced
% for the TROTS dataset). Normal mimimum CT value = -1000 HU.
CTMask = and(CTMask, patient.CT>-1024);
CTPositions = findn(CTMask);

% Transform voxels to spatial coordinates
p = bsxfun(@times, CTPositions - 1, patient.Resolution);
CTPoints = bsxfun(@plus, p, patient.Offset);
            
% Make Dose matrix here
Dose = zeros(size(patient.CT), 'single');

% Some voxels appear double, we know that, please do not remember us that
% we could have sampled more efficiently. Fact is: volumes often overlap,
% so double sampled points are intentional.
warning('off', 'MATLAB:TriScatteredInterp:DupPtsAvValuesWarnId')
warning('off', 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')
Dose(CTMask==1) = griddata(SampledPoints(:, 1), SampledPoints(:, 2), SampledPoints(:, 3), double(SampledDose), CTPoints(:, 1), CTPoints(:, 2), CTPoints(:, 3), 'linear'); % 'natural' gives much nicer doses for brachytherapy, but takes forever.
warning('on', 'MATLAB:TriScatteredInterp:DupPtsAvValuesWarnId')
warning('on', 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')
return
end

function ind=findn(arr)
    %FINDN   Find indices of nonzero elements.
    %   I = FINDN(X) returns the indices of the vector X that are
    %   non-zero. For example, I = FINDN(A>100), returns the indices
    %   of A where A is greater than 100. See RELOP.
    %  
    %   This is the same as find but works for N-D matrices using 
    %   ind2sub function
    %
    %   It does not return the vectors as the third output arguement 
    %   as in FIND
    %   
    %   The returned I has the indices (in actual dimensions)
    %
    %   x(:,:,1)            x(:,:,2)            x(:,:,3)
    %       = [ 1 2 3           =[11 12 13        =[21 22 23
    %           4 5 6             14 15 16          24 25 26
    %           7 8 9]            17 18 19]         27 28 29]
    %
    %   I=find(x==25) will return 23
    %   but findn(x==25) will return 2,2,3
    %   
    %   Also see find, ind2sub

    %   Loren Shure, Mathworks Inc. improved speed on previous version of findn
    %   by Suresh Joel Mar 3, 2003

    in=find(arr);
    sz=size(arr);
    if isempty(in), ind=[]; return; end;
    [out{1:ndims(arr)}] = ind2sub(sz,in);
    ind = cell2mat(out);
    return
end