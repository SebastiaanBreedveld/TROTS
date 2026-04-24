% TROTSReadOutput reads numerical output from solver from disk
%
% solutionX = TROTSReadOutput(FILE)
%
%   reads SolutionX from text file FILE and returns the vector. Optionally,
%   for problems containing DVH cost-functions, it will also return the
%   final parameters:
%
% [solutionX, DVHParameters] = TROTSReadOutput(FILE);
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
function [solutionX, DVHParameters] = TROTSReadOutput(OutputFile)

fid = fopen(OutputFile, 'r');
if fid==-1
    error('Error opening file: %s', OutputFile);
end

tline = fgetl(fid);
HeaderLines = 1;
DVHParameters = [];
if ~strncmpi(tline, 'Status', 6)
    fclose(fid);
    error('File %s is probably not an optimisation output!', OutputFile);
else
    % do stuff
    [v, pos] = textscan(tline, '%25c %f');
    if v{2}==0
        fprintf('Problem %s converged! ', OutputFile);
    else
        fprintf('Problem %s FAILED! ', OutputFile);
    end
    while (~strncmpi(tline, 'Optimised result', 16))
        if strncmpi(tline,  'DVH Parameters Offset', 21)
            DVHParameters(2, :) = sscanf(tline(23:end), '%f');
        elseif strncmpi(tline,  'DVH Parameters', 14)
            DVHParameters(1, :) = sscanf(tline(23:end), '%f');
        end
        tline = fgetl(fid);
        HeaderLines = HeaderLines + 1;
        if ~ischar(tline)
            fclose(fid);
            error('File ended prematurely');
        end
    end
    fclose(fid);    
    
    
    % This line holds the number of elements - since Matlab2015a, and extra
    % element is read
    numberOfElements = str2double(tline(regexp(tline, '[\d]')))

    % Load result vector and original data
    solutionX = dlmread(OutputFile, '\n', HeaderLines, 0);    
    
    if ~isnan(numberOfElements)
        solutionX = solutionX(1:numberOfElements);
    elseif ~verLessThan('matlab', '8.6') % R2015b
        % This is just a guess - read old style output with Matlab 2015a+
        solutionX = solutionX(1:end-1);
    end
    solutionX = solutionX(:);

    fprintf('Successfully read output.\n');
end
