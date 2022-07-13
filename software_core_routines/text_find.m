
function line_indexes = text_find(InputFile, SearchString)

    % Scans for the specified line "SearchString" in the InputFile;
    % returns the list of all line numbers where it was found; returns an empty list if not found at any line.
    % It is case sensitive.
    %
    % InputFile - string.
    % SearchString - string.
    % line_indexes - array of integers.
    
    % read whole file into cell array
    fid = fopen(InputFile);
    data = textscan(fid, '%s', 'Delimiter', '\n', 'whitespace', '', 'CollectOutput', true);   
    fclose(fid);
    
    % find
    line_indexes = [];
    for ii = 1:length(data{1})
        tf = strfind(data{1}{ii}, SearchString); % search for this line in the array
        if ~isempty(tf)
            line_indexes(end+1) = ii;
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

        


    