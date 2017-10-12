function is_line_skipped = step_notification(step_i, Nstep_notif, Nstep_skipline)
%STEP_NOTIFICATION
% step_i : iteration index
% Convenient to display some iterations
% Displays only 1 iterations ever Nstep_notif
% New line every Nstep_skipline iteration

% Author: Christophe Gardella
% Tested on Matlab 2014b
% History:
%   Original: 10/10/2017

if ~mod(step_i, Nstep_notif)
    fprintf([int2str(step_i) ' - ']);
end
if ~mod(step_i, Nstep_skipline)
    fprintf('\n');
    is_line_skipped = 1==1;
else
    is_line_skipped = 0==1;
end

end

