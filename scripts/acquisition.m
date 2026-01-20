function [u, y, sel, sig, realizations, power_levels] = acquisition(filename, prefix)
% [u, y, sel, sig, realizations, power_levels] = acquisition(filename, prefix)
%
% Load the acquired input and output signals from .mat files
%
% Inputs:
%   filename - Folder name inside results/ (e.g., 'results_3_power/')
%   prefix - Optional file prefix (default: '' for no prefix)
% Outputs:
%   u - Input signal tensor (samples x realizations x power levels)
%   y - Output signal tensor (samples x realizations x power levels)
%   sel - excited bins
%   sig - input signal
%   realizations - Number of realizations, determined using available files
%   power_levels - Number of power levels, determined using available files
%
% It assumes the measurements name format is:
% '{prefix}ACQ_R{realization}_P{power level}_E0_M0_F0.mat'

    if nargin < 2
        prefix = '';
    end

    currentScriptPath = fileparts(mfilename('fullpath'));

    resultsPath = fullfile(currentScriptPath, '../results/');
    [realizations, power_levels] = detect_counts(resultsPath, filename, prefix);

    % load a single file to get N
    file = strcat(resultsPath, filename, prefix, "ACQ_R0_P0_E0_M0_F0.mat");
    load(file);
    N = length(YR0);

    % u and y are tensor with following dimensions:
    % (sample, realization, power level)
    u = zeros(N, realizations, power_levels);
    y = zeros(N, realizations, power_levels);

    % load measured signals
    for r = 0:realizations-1
        for p = 0:power_levels-1
            file = strcat(resultsPath, filename, prefix, "ACQ_R", num2str(r), "_P", num2str(p), "_E0_M0_F0.mat");
            load(file);
            u(:, r+1, p+1) = YR0(:);
            y(:, r+1, p+1) = YR1(:);
        end
    end

    % load sig
    sigFile = strcat(resultsPath, filename, prefix, "AWG_R0_E0.mat");
    load(sigFile);
    sig = SigR0;

    % load sel
    selFile = strcat(resultsPath, filename, prefix, "AWGEXC_R0_E0.mat");
    load(selFile);
    sel = SelExc0;

end

function [realizations, power_levels] = detect_counts(resultsPath, filename, prefix)
% Detect number of realizations and power levels by probing for existing files.
% Returns counts (non-negative integers). Filenames checked:
% [resultsPath]/[filename][prefix]ACQ_R{r}_P{p}_E0_M0_F0.mat
%
% Example: [realizations, power_levels] = detect_counts('../results/', 'out2k', '');

    % detect realizations by probing P0 files starting at R0
    r = 0;
    while true
        fname = fullfile(resultsPath, strcat(filename, prefix, "ACQ_R", num2str(r), "_P0_E0_M0_F0.mat"));
        if isfile(fname)
            r = r + 1;
        else
            break;
        end
    end
    realizations = r; % may be 0 if no files found

    % detect power levels: find minimum across all realizations
    power_levels = 0;
    if realizations > 0
        min_p = inf;
        for ri = 0:realizations-1
            p = 0;
            while true
                fname = fullfile(resultsPath, strcat(filename, prefix, "ACQ_R", num2str(ri), "_P", num2str(p), "_E0_M0_F0.mat"));
                if isfile(fname)
                    p = p + 1;
                else
                    break;
                end
            end
            min_p = min(min_p, p);
        end
        power_levels = min_p;
    end
end

