
function [] = text_replace(InputFile, OutputFile, SearchStringList, ReplaceStringList)

    % Scans for the specified "SearchString" in every line of InputFile; 
    % if found, it replaces the WHOLE line in which it was found with "ReplaceString" in the OutputFile. 
    % If not found, it does nothing (and, if the search and replace lists are cell arrays, proceeds to the next pair of serach/replace strings).
    % If multiple lines containing one SearchString are found, ALL of them are replaced with ReplaceString. 
    % It is NOT case sensitive.
    % If the same file is specified as both the input and output file, the file is overwritten.
    %
    % InputFile - string.
    % OutputFile - string.
    % SearchStringList - string or cell array of strings.
    % ReplaceStringList - string or cell array of strings (must be of same size as SearchString).
    
    % read whole file into cell array
    fid = fopen(InputFile);
    data = textscan(fid, '%s', 'Delimiter', '\n', 'whitespace', '', 'CollectOutput', true);   
    fclose(fid);
    
    if iscell(SearchStringList) && iscell(ReplaceStringList)
        for jj = [1:size(SearchStringList,2)]
            SearchString = SearchStringList{jj};
            ReplaceString = ReplaceStringList{jj};
            % find the position where changes need to be applied and insert new data
            for ii = 1:length(data{1})
                tf = strfind(lower(data{1}{ii}), lower(SearchString)); % search in every line of the array
                if ~isempty(tf)
                    data{1}{ii} = ReplaceString; % replace the whole line with ReplaceString
                end
            end
        end
    else
        % find the position where changes need to be applied and insert new data
        for ii = 1:length(data{1})
            tf = strfind(data{1}{ii}, SearchStringList); % search in every line of the array
            if ~isempty(tf)
                data{1}{ii} = ReplaceStringList; % replace the whole line with ReplaceString
            end
        end
    end
    
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



    