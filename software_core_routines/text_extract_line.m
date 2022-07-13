
function line_text = text_extract_line(InputFile, SearchString)

    % Scans for the specified line "SearchString" in the InputFile;
    % returns the text of the whole line in which SearchString is found.
    % It is case sensitive.
    % If the string is found multiple times in InputFile, each of the lines
    % is saved as one cell in the ourput cell array.
    %
    % InputFile - string.
    % SearchString - string.
    % line_text - cell array of strings.
    
    % read whole file into cell array
    fid = fopen(InputFile);
    data = textscan(fid, '%s', 'Delimiter', '\n', 'whitespace', '', 'CollectOutput', true);   
    fclose(fid);
    
    % find
    line_text = {};
    for ii = 1:length(data{1})
        tf = strfind(data{1}{ii}, SearchString); % search for this line in the array
        if ~isempty(tf)
            line_text{end+1} = data{1}{ii};
        end
    end
    
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

    

    