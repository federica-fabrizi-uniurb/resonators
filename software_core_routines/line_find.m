
function line_text = line_find(InputFile, line_index)

    % Finds the line of InputFile indicated by line_index,
    % and returns the text included in that line as a string.
    % E.g. if line_index = 1, then the first line gets returned.
    % 
    % InputFile - string.
    % line_index - integer.
    % line_text - string.
    
    % read whole file into cell array
    fid = fopen(InputFile);
    data = textscan(fid, '%s', 'Delimiter', '\n', 'whitespace', '', 'CollectOutput', true);   
    fclose(fid);
    
    line_text = data{1}{line_index};
    
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

    
    

    