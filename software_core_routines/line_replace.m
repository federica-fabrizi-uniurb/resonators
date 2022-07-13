
function [] = line_replace(InputFile, OutputFile, line_index, ReplaceString)

    % Replaces the line of InputFile indicated by line_index with ReplaceString,
    % generating OutputFile as a result.
    % E.g. if line_index = 1, then the first line gets replaced.    
    % If the same file is specified as both the input and output file, the file is overwritten.
    %
    % InputFile - string.
    % OutputFile - string.
    % line_index - integer.
    % ReplaceString - string.
    
    % read whole file into cell array
    fid = fopen(InputFile);
    data = textscan(fid, '%s', 'Delimiter', '\n', 'whitespace', '', 'CollectOutput', true);   % 
    fclose(fid);
    
    % modify the cell array
    data{1}{line_index} = ReplaceString; % replace the whole line with ReplaceString

    % write the modified cell array into the text file
    fid = fopen(OutputFile, 'w');
    for ii = 1:length(data{1})
        fprintf(fid, '%s\n', char(data{1}{ii}));
    end
    fclose(fid);
    
end




% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

% F. Fabrizi (2022); federica.fabrizi@uniurb.it

    