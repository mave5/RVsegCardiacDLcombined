% savefig(fig_handle, filename, format_list)
%
% Convenience function for rapid saving of figures into a central location,
% with files grouped into subdirectories by date.
%
%
% Inputs:
%
% curplot  -  The figure number of the plot to be saved. To save the
%             currently selected plot, 'gcf' can be used here.
%
% filename -  Name of the file. No filetype suffix necessary, this is added
%             automatically.
%
% formats  -  (Optional - default value {'png', 'fig'} )
%             String controlling which format(s) to save the figure in.
%             This is done as a cell array of format codes which are the
%             same as used for the saveas command. Type '>> doc saveas'
%             to see a list of these.
%
% Example:    To save a figure with filename 'my figure' in .fig and .eps
%             formats, type:
%             >> savefig(gcf, 'my figure', {'fig', 'eps'})
%
% Matthew Warden
% mattwarden@gmail.com
% 15 June 2010
%

function savefig(fig_handle, filename, format_list)

%Set default file formats to save as
if nargin < 3
    format_list = {'png', 'fig'};
end

%Location to save figures to
%basedir = 'C:/Matlab Figures/';

%Files are saved in subdirectories according to the date.
%yeardir = datestr(now, 'yyyy/');
%subdir = datestr(now, 'yyyy_mm_dd/');

%Construct a string with the entire path to the save directory
%savedir = [basedir yeardir subdir];

%Make this directory (if it's already there, this line does nothing)
%[s, mess] = mkdir(savedir);
%if s~=1
%    error(['savefig failed when trying to create directory with the ' ...
%        'following message: ' mess]);
%end

for format = format_list
    saveas(fig_handle, [savedir filename], format{1});
end
