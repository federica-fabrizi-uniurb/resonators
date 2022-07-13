
function obj = load_element_numbers(obj)

    % Loads into calculation_results three integers containing the number
    % of elements used by Ansys for the analysis of the substrate, the 1st coaring, and the 2nd coating, respectively.
    % The number for the two coatings is expected to be the same.
    % This function is used as of now for resonators of type DoublyCoatedDisc;
    % it may be useful to extend it in the future also to simpler resonator types.
    %
    % obj - instance of class DoublyCoatedDisc.

    OutputFolder = obj.calculation_settings.apdl_output_folder;
    OutputFile = fullfile(OutputFolder,'output_file.out');
    
    SearchString = 'Number of total elements subs';
    line_text = text_extract_line(OutputFile, SearchString);
    element_number_substrate = str2num(regexprep(line_text{end},{'\D*([\d\.]+\d)[^\d]*', '[^\d\.]*'}, {'$1 ', ' '}));
    
    SearchString = 'Number of total elements coat_one';
    line_text = text_extract_line(OutputFile, SearchString);
    element_number_coating1 = str2num(regexprep(line_text{end},{'\D*([\d\.]+\d)[^\d]*', '[^\d\.]*'}, {'$1 ', ' '}));   

    SearchString = 'Number of total elements coat_two';
    line_text = text_extract_line(OutputFile, SearchString);
    element_number_coating2 = str2num(regexprep(line_text{end},{'\D*([\d\.]+\d)[^\d]*', '[^\d\.]*'}, {'$1 ', ' '}));   
    
    obj.calculation_results.element_number_substrate = element_number_substrate;
    obj.calculation_results.element_number_coating1 = element_number_coating1;
    obj.calculation_results.element_number_coating2 = element_number_coating2;
    
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




