
function [] = text_add(InputFile, OutputFile, SearchString, InsertFile)

    % Searches for SearchString in InputFile, then 
    % inserts the text contained in InsertFile,
    % starting at the line following the one where SearchString is found; 
    % generates OutputFile as a result.
    % If the same file is specified as both the input and output file, the file is overwritten.    
    % The string search is case sensitive.
    % If the string is found multiple times in InputFile, the text gets
    % inserted after the first occurrence; this can be changed by editing the
    % code as indicated by comment, below.
    %
    % InputFile - string.
    % OutputFile - string.
    % SearchString - string.
    % ReplaceString - string.
    
    % read whole file into cell array
    fid = fopen(InputFile);
    data = textscan(fid, '%s', 'Delimiter', '\n', 'whitespace', '', 'CollectOutput', true);    
    fclose(fid);
    
    fin = fopen(InsertFile);
    datain = textscan(fin, '%s', 'Delimiter', '\n', 'whitespace', '', 'CollectOutput', true);    
    fclose(fin);
    
    % find the position to insert new data
    insert_place = [];
    for ii = 1:length(data{1})
        tf = strfind(data{1}{ii}, SearchString); % search in every line of the array
        if ~isempty(tf)
            insert_place(end+1) = ii; 
        end
    end

    % write the cell array + added text into the text file
    insert_here = insert_place(1);  % choose insert_place(end) for last occurrence
    fid = fopen(OutputFile, 'w');
    for ii = 1:insert_here
        fprintf(fid, '%s\n', char(data{1}{ii}));
    end
    for ii = 1:length(datain{1})
        fprintf(fid, '%s\n', char(datain{1}{ii}));
    end
    for ii = insert_here+1:length(data{1})
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

    

    
    