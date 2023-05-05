% INPUTS:      
%   - t_nn: time vector of RR interval (in s)
%   - nn: normal to normal interval (in s)
%   - opt: 'normal', 'sqi', 'af', 'mse', 'dfa', 'hrt'
% 
% OUTPUT:     
%   - w : start time (in s) of each window to be analyzed
% 
% Cedric Cannard, 2022

function w = create_window(t_nn, nn, opt)

% win_tol = 0.150; % tolerance

switch opt
    case 'normal'
        increment = 120;
        winLength = 120;

    case 'af'
        increment = 30;
        winLength = 30;
%     case 'mse'
%         increment = HRVparams.MSE.increment * 3600;
%         windowlength = HRVparams.MSE.windowlength * 3600;
%         if isempty(increment)
%             windowRRintervals = 0;
%             return % no need to crate windows , use entair signal
%         end 
%     case 'dfa'
%         increment = HRVparams.DFA.increment * 3600;
%         windowlength = HRVparams.DFA.windowlength * 3600;
%         if isempty(increment)
%             windowRRintervals = 0;
%             return  % no need to crate windows , use entair signal
%         end       
    case 'sqi'
        increment = 1;
        winLength = 10;
    case 'hrt'
        increment = 24;
        winLength = 24*3600;
end

if winLength > t_nn(end)
    w = 0; % use whole file 
else
    nx = floor(t_nn(end));           % length of sequence
    overlap = winLength-increment;   % number of overlapping elements

    % starting index of each windows
    Nwinds = fix((nx-overlap)/(winLength-overlap));    % number of sliding windows
    w = (0:Nwinds-1)*(winLength-overlap);  

%     tstart = 0;  % Window Start Time
%     i = 1;       % counter
% 
%     if ~strcmp(opt,'af') && ~strcmp(opt,'sqi')
%         
%         while tstart <= t_nn(end) - winLength + increment
%     
%             % Find indices of time values in this segment
% %             t_win = t_nn(t_nn >= tstart & t_nn < tstart + winLength);
%     
%             % if nn intervals are supplied, assign them to the current window
%             % if not, put in place a vector of ones as a place holder
% %             if ~isempty(nn)
%             nn_win = nn(t_nn >= tstart & t_nn < tstart + winLength);
% %             else
% %                 nn_win = (winLength/length(t_win))* ones(length(t_win),1);
% %             end
%             
%             % Store the begin time of window
%             w(i) = tstart;
%     
%             % Increment time by sliding segment length (sec)
%             tstart = tstart + increment;
%     
%             % Check Actual Window Length and mark windows that do not meet the
%             % crieria for Adequate Window Length
%             % First remove unphysiologic beats from candidates for this
%             % measurement:
%             idxhi = find(nn_win > HRVparams.preprocess.upperphysiolim);
%             idxlo = find(nn_win < HRVparams.preprocess.lowerphysiolim);
%             comb = [idxhi(:) ; idxlo(:)];
%             nn_win(comb) = [];
%             % Now query the true length of the window by adding up all of the nn
%             % intervals
%             truelength = sum(nn_win(:));
%             if truelength < (winLength * (1-win_tol))
%                 w(i) = NaN; 
%             end
%     
%             % Increment loop index
%             i = i + 1;
%         end
%     end
end



