function [constraints, dataopt] = TROTStoYARTOS(problem, dataopt)

for j=1:length(problem)
    constraints(j).datID = problem(j).dataID;
    constraints(j).VolName = problem(j).Name;
    constraints(j).Minimize = problem(j).Minimise;
    constraints(j).Type = problem(j).Type;
    constraints(j).Objective = problem(j).Objective;
    constraints(j).ObjectiveMax = [];
    constraints(j).ObjectiveSufficient = problem(j).Sufficient;
    constraints(j).Priority = problem(j).Priority;
    constraints(j).Weight = problem(j).Weight;    
    constraints(j).Parameters = problem(j).Parameters;
    constraints(j).Active = problem(j).Active;
    constraints(j).Bounds = problem(j).IsConstraint;
    if constraints(j).Type>1
        constraints(j).numcons = 1;
    else
        constraints(j).numcons = size(dataopt.matrix(constraints(j).datID).A, 1);
    end
end

for j=1:length(dataopt.matrix)
    if ~issparse(dataopt.matrix(j).A)
        dataopt.matrix(j).A = double(dataopt.matrix(j).A);
    end
    dataopt.matrix(j).b = double(dataopt.matrix(j).b);
    dataopt.matrix(j).c = double(dataopt.matrix(j).c);
    dataopt.matrix(j).numvox = size(dataopt.matrix(j).A, 1);
end
dataopt.misc.real = 1:dataopt.misc.real;

dataopt.misc.metadata.Patient.ActiveTargets = 1:length(dataopt.misc.InitialiseMatrixID);
for j=1:length(dataopt.misc.InitialiseMatrixID)
    dataopt.misc.metadata.Dose.Config(j).DataMatrix = dataopt.misc.InitialiseMatrixID(j);
    dataopt.misc.metadata.Dose.Config(j).DataMatrixMean = [];
end
dataopt.misc.metadata.Dose.Reg2 = dataopt.misc.InitialiseRegularisationMatrixID;

dataopt.misc.metadata.Dose.Reg1 = [];
dataopt.misc.metadata.Dose.Reg1Cons = [];
for j=1:length(constraints)
    if isequal(constraints(j).VolName, 'Smoothing Linear')
        dataopt.misc.metadata.Dose.Reg1 = unique([dataopt.misc.metadata.Dose.Reg1 constraints(j).datID]);
        dataopt.misc.metadata.Dose.Reg1Cons = [dataopt.misc.metadata.Dose.Reg1Cons j];
    end
end

