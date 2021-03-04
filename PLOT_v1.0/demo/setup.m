% set up matlab path for PLOT

% look for directory where this file setup.m is installed
if ~exist('PLOTdir', 'var')
	PLOTdir = which('setup');
	PLOTdir = fileparts(PLOTdir);
end

if PLOTdir(end) ~= filesep % make sure there is a '/' at end of directory
	PLOTdir = [PLOTdir filesep];
end

path([PLOTdir 'demo'], path);       % demo
path([PLOTdir 'projector'], path);  % projector
path([PLOTdir 'penalty'], path);    % regularization
path([PLOTdir 'eml'], path);        % emission image reconstruction
path([PLOTdir 'utils'], path);      % utilities
path([PLOTdir 'work'], path);       % work: my own work folder

