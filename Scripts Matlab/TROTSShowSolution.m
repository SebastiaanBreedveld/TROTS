% TROTSShowSolution shows the attained values/constraints for a given
%   solution.
%
%   TROTSShowSolution(solutionX, problem, data)
%
%   shows the values of each constraint and objective side by side.
%
%   If the problem contains DVH cost-functions, the parameters settled by
%   the solver can also be provided as argument:
%
%   TROTSShowSolution(solutionX, problem, data, DVHParameters);
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
function TROTSShowSolution(solutionX, problem, data, DVHParameters)

    % Check if DVHParameters were passed on
    if nargin<=3
        DVHParameters = [];
    end
    
    % Now, loop over all functions
    fprintf('%-30s %-8s %-17s %10s %10s %10s %10s\n', 'Volume Name', 'Index', 'Type', 'Value', 'Constraint', 'MC Goal', 'DVH Approx');
    fprintf(repmat('-', 1, 103));
    
    fprintf('\n');
    for j=1:length(problem)
        ConValue = '';
        if problem(j).Active
            if problem(j).IsConstraint
                Flag = 'C';
                ConValue = problem(j).Objective;
                if problem(j).Type==5
                    % For DVH, make constraint show full % and set
                    % steepness parameter
                    ConValue = ConValue * 100;
                    if ~isempty(DVHParameters)
                        problem(j).Parameters(2) = DVHParameters(1,1);
                        if size(DVHParameters, 1)>=2
                            problem(j).Parameters(3) = DVHParameters(2, 1);
                        end
                        DVHParameters = DVHParameters(:, 2:end);
                    end
                end
                GoalValue = '';
            else
                Flag = 'O';
                ConValue = '';
                GoalValue = problem(j).Objective;
                if problem(j).Type==5
                    % For DVH, make constraint show full % and set
                    % steepness parameter
                    ConValue = ConValue * 100;
                    GoalValue = GoalValue * 100;
                    if ~isempty(DVHParameters)                        
                        problem(j).Parameters(2) = DVHParameters(1, 1);
                        if size(DVHParameters, 1)>=2
                            problem(j).Parameters(3) = DVHParameters(2, 1);
                        end
                        DVHParameters = DVHParameters(2:end);
                    end
                end
            end
        else
            Flag = 'I';           
        end
        
        [fval, Type, DVHApprox] = EvaluateFunction(solutionX, data, problem, j);
        fprintf('%-30s (%s %3d): %-17s %10.2f %10s %10s %10.2f\n', problem(j).Name, Flag, j, Type, fval, sprintf('%10.2f', ConValue), sprintf('%10.2f', GoalValue), DVHApprox);
    end

end


function [fval, Type, DVHApprox] = EvaluateFunction(x, data, problem, j)
    DVHApprox = '';
    MatrixStruct = data.matrix(problem(j).dataID);
    d = MatrixStruct.A*x + MatrixStruct.b;
    Params = problem(j).Parameters;
    switch problem(j).Type
        case 1  % Linear        
            if problem(j).Minimise==1
                fval = max(d);
                Type = 'linear (maximum)';
            else
                fval = min(d);
                Type = 'linear (minimum)';
            end
            if isscalar(d)
                Type = 'linear (mean)';
            end
        case 2  % Quadratic
            Type = 'quadratic';
            % Catch a minor "misunderstanding" between me and Matlab: a
            % scalar 0 implies a full 0-vector
            if length(MatrixStruct.b)==1 && MatrixStruct.b==0
                fval = 0.5*(x'*MatrixStruct.A*x) + MatrixStruct.c;
            else
                fval = 0.5*(x'*MatrixStruct.A*x) + MatrixStruct.b'*x + MatrixStruct.c;
            end
        case 3  % EUD - generalised Equivalent Uniform Dose or generalised mean
            Type = 'EUD';
            fval = (sum(d.^Params(1))/length(d)).^(1/Params(1));
        case 4  % LTCP - Logarithmic Tumour Control Probability
            Type = 'LTCP';
            fval = sum(exp(-Params(2)*(d - Params(1)))/length(d));
        case 5  % Dose-Volume - Partial Volume
            Type = 'DVH';
            fval = sum(d >= Params(1))/length(d)*100;
            if length(Params)==2
                K = double(d)/Params(1);
                DVHApprox = sum(K.^Params(2)./(1+K.^Params(2)))/length(d)*100;
            elseif length(Params)==3
                K = double(d)/(Params(1) + Params(3));
                DVHApprox = sum(K.^Params(2)./(1+K.^Params(2)))/length(d)*100;
            else
                DVHApprox = [];
            end
        case 6 % Chain constraint
            Type = 'Chain';
            fval = 0;
            for k=1:size(problem(j).Chain, 1)
                fval = fval + problem(j).Chain(k, 2)*EvaluateFunction(x, data, problem, problem(j).Chain(k, 1)) - x(problem(j).Chain(k, 3));
            end
        otherwise
            Type = 'UNKNOWN';
            fval = NaN;
    end
end
