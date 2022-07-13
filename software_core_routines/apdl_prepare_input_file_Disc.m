
function apdl_prepare_input_file_Disc(discobj)

    % Edits the mechanical apdl input file "input_file.dat" (within input_folder)
    % and saves it as "input_file_edited.dat" in the same folder,
    % by setting the material properties, specifying the solver settings, and appending the post-processing commands to save the relevant data.
    %
    % discobj - instance of class Disc or CantileverFiber.

    InputFolder = discobj.calculation_settings.apdl_input_folder;
    InputFileName = discobj.calculation_settings.apdl_input_file_template;
    InputFileEditedName = discobj.calculation_settings.apdl_input_file_edited;
    OutputFolder = discobj.calculation_settings.apdl_output_folder;

    InputFile = fullfile(InputFolder,InputFileName);
    InputFileEdited = fullfile(InputFolder,InputFileEditedName);
    copyfile(InputFile,InputFileEdited);

    % CHECK UNITS

    line_indexes = text_find(InputFileEdited, '/units,MKS');
    if isempty(line_indexes) == false
        disp('Units MKS');
    else
        disp('Warning: Units not MKS');
        return
    end

    % EDIT DENSITY, YOUNG, POISSON

    SearchStringList = {};
    ReplaceStringList = {};

    jj_body = 1;
    density = discobj.substrate.material.density;
    SearchStringList{end+1} = ['MP,DENS,',num2str(jj_body),','];
    ReplaceStringList{end+1} = ['MP,DENS,',num2str(jj_body),',',num2str(density),', ! kg m^-3   - EDITED BY MATLAB'];
    % if material is amorphous (isotropic)
    if strcmp(discobj.substrate.material.flag_anisotropic,'amorph') 
        if ~isempty(discobj.substrate.material.young_amorph)
            young = discobj.substrate.material.young_amorph;
        else
            young = discobj.substrate.material.young_poly;
        end
        SearchStringList{end+1} = ['MP,EX,',num2str(jj_body),','];
        ReplaceStringList{end+1} = ['MP,EX,',num2str(jj_body),',',num2str(young),', ! Pa   - EDITED BY MATLAB'];
        if ~isempty(discobj.substrate.material.poisson_amorph)
            poisson = discobj.substrate.material.poisson_amorph;
        else
            poisson = discobj.substrate.material.poisson_poly;
        end    
        SearchStringList{end+1} = ['MP,NUXY,',num2str(jj_body),','];
        ReplaceStringList{end+1} = ['MP,NUXY,',num2str(jj_body),',',num2str(poisson),', ! - EDITED BY MATLAB'];    
    else
        % if material is polycrystalline (isotropic)    
        if strcmp(discobj.substrate.material.flag_anisotropic,'poly') 
            if ~isempty(discobj.substrate.material.young_poly) 
                young = discobj.substrate.material.young_poly;
            else
                young = discobj.substrate.material.young_poly_from_stiffness_matrix();
            end
            SearchStringList{end+1} = ['MP,EX,',num2str(jj_body),','];
            ReplaceStringList{end+1} = ['MP,EX,',num2str(jj_body),',',num2str(young),', ! Pa   - EDITED BY MATLAB'];       
            if ~isempty(discobj.substrate.material.poisson_poly) 
                poisson = discobj.substrate.material.poisson_poly;
            else
                poisson = discobj.substrate.material.poisson_poly_from_stiffness_matrix();
            end        
            SearchStringList{end+1} = ['MP,NUXY,',num2str(jj_body),','];
            ReplaceStringList{end+1} = ['MP,NUXY,',num2str(jj_body),',',num2str(poisson),', ! - EDITED BY MATLAB'];  
        else
            % if material is single crystal (anisotropic)    
            if strcmp(discobj.substrate.material.flag_anisotropic,'single_crystal')  
                if ~isempty(discobj.substrate.material.stiffness_matrix)
                    stiffness_matrix = discobj.substrate.material.stiffness_matrix;
                    SearchStringList{end+1} = ['MP,EX,',num2str(jj_body),','];
                    SearchStringList{end+1} = ['MP,NUXY,',num2str(jj_body),','];
                    ReplaceStringList{end+1} = '';
                    % Ansys uses an unusual order in Voigt's vector components:
                    % (xx, yy, zz, xy, yz, xz) instead of usual (xx, yy, zz, yz, xz, xy)
                    ReplaceStringList{end+1} = [newline 'TB,ELASTIC,' num2str(jj_body) ',,21,AELS' newline ...
                    'TBDATA,1,' num2str(stiffness_matrix(1,1)) ',' num2str(stiffness_matrix(2,1)) ',' num2str(stiffness_matrix(3,1)) ',' num2str(stiffness_matrix(6,1)) ',' num2str(stiffness_matrix(4,1))  ',' num2str(stiffness_matrix(5,1)) newline ...
                    'TBDATA,7,' num2str(stiffness_matrix(2,2)) ',' num2str(stiffness_matrix(3,2)) ',' num2str(stiffness_matrix(6,2)) ',' num2str(stiffness_matrix(4,2)) ',' num2str(stiffness_matrix(5,2))  ',' num2str(stiffness_matrix(3,3)) newline ...
                    'TBDATA,13,' num2str(stiffness_matrix(6,3)) ',' num2str(stiffness_matrix(4,3)) ',' num2str(stiffness_matrix(5,3)) ',' num2str(stiffness_matrix(6,6)) ',' num2str(stiffness_matrix(6,4))  ',' num2str(stiffness_matrix(6,5)) newline ...
                    'TBDATA,19,' num2str(stiffness_matrix(4,4)) ',' num2str(stiffness_matrix(5,4)) ',' num2str(stiffness_matrix(5,5)) newline];
                else
                    young = mean(cell2mat(struct2cell(discobj.substrate.material.young_directional)));
                    SearchStringList{end+1} = ['MP,EX,',num2str(jj_body),','];
                    ReplaceStringList{end+1} = ['MP,EX,',num2str(jj_body),',',num2str(young),', ! Pa   - EDITED BY MATLAB'];                
                    poisson = mean(cell2mat(struct2cell(discobj.substrate.material.poisson_directional)));
                    SearchStringList{end+1} = ['MP,NUXY,',num2str(jj_body),','];
                    ReplaceStringList{end+1} = ['MP,NUXY,',num2str(jj_body),',',num2str(poisson),', ! - EDITED BY MATLAB']; 
                end
            end
        end
    end

    text_replace(InputFileEdited, InputFileEdited, SearchStringList, ReplaceStringList);

    % REMOVE UNNECESSARY COMMANDS

    % They are not used in this kind of analysis anyway, 
    % but they can create confusion as they do not correspond to the material and temperature at hand

    SearchStringList = {'toffst','tref','''TEMP''','MP,ALPX','MP,C','MP,KXX','MP,RSVX','MP,MURX'};
    ReplaceStringList = {'','','','','','','',''};

    text_replace(InputFileEdited, InputFileEdited, SearchStringList, ReplaceStringList);

    % EDIT ELEMENT OPTIONS (NOT RECOMMENDED)

    % Uncomment the following part to
    % override the solver's recommendations and keep Workbench's
    % (necessary to keep full integration in solid186)
    % (IT IS NOT RECOMMENDED TO DO THIS)

    % SearchStringList = {'etcon,set'};
    % ReplaceStringList = {'!etcon,set'};
    % 
    % text_replace(InputFileEdited, InputFileEdited, SearchStringList, ReplaceStringList);
    %
    % % check element type
    % 
    % line_index = text_find(InputFileEdited, 'et,1,186');
    % if isempty(line_index)
    %     disp('Warning: element type does not appear to be SOLID186 - check input_file_edited.dat')
    % end
    % line_index = text_find(InputFileEdited, 'keyo,1,2,1');
    % if isempty(line_index)
    %     disp('Warning: full integration in SOLID186 does not appear to be set - check input_file_edited.dat')
    % end

    % EDIT THE SOLVER SETTINGS

    min_frequency = discobj.calculation_settings.min_frequency;
    max_frequency = discobj.calculation_settings.max_frequency;
    max_number_of_modes = discobj.calculation_settings.max_number_of_modes;

    SearchString = 'outres';
    line_indexes = text_find(InputFileEdited, SearchString);
    for jj_line = line_indexes
        line_replace(InputFileEdited, InputFileEdited, jj_line, '');
    end
    SearchString = 'mxpand';
    line_indexes = text_find(InputFileEdited, SearchString);
    for jj_line = line_indexes
        line_replace(InputFileEdited, InputFileEdited, jj_line, '');
    end
    SearchString = 'dmpopt';
    line_indexes = text_find(InputFileEdited, SearchString);
    for jj_line = line_indexes
        line_replace(InputFileEdited, InputFileEdited, jj_line, '');
    end

    SearchString = 'modopt';

    ReplaceString = ['modopt,lanb,',num2str(max_number_of_modes),',',num2str(min_frequency),',',num2str(max_frequency),newline,...
        'outres,erase',newline,...
        'outres,all,none',newline,...
        'outres,nsol,all',newline,...
        'outres,veng,all',newline,...
        'outres,etmp,all',newline,...
        'outres,strs,all',newline,...
        'outres,epel,all',newline,...
        'outres,eppl,all',newline,...
        'mxpand,,,,yes,,yes             ! expand requested element results; write them to file.mode',newline,...
        'dmpopt,esav,no',newline,...
        'dmpopt,emat,no',newline,...
        'dmpopt,full,no'];

    text_replace(InputFileEdited, InputFileEdited, SearchString, ReplaceString);

    % INSERT TEXT FROM SCRIPT, FOR POST-PROCESSING:

    % edit post-processing script
    ScriptFile = fullfile(discobj.calculation_settings.templates_folder,'postprocessing_script.txt');
    ScriptFileEditedTmp = fullfile(discobj.calculation_settings.support_folder,'postprocessing_script_edited_tmp.txt');
    ScriptFileEdited = fullfile(discobj.calculation_settings.support_folder,'postprocessing_script_edited.txt');
    DataFileName_list = ["data_positionsX.txt","data_positionsY.txt","data_positionsZ.txt",...
            "data_energy_sene.txt",...
            "data_vol.txt",...
            "data_strainXX.txt","data_strainYY.txt","data_strainZZ.txt","data_strainXY.txt","data_strainYZ.txt","data_strainXZ.txt",...
            "data_stressXX.txt","data_stressYY.txt","data_stressZZ.txt","data_stressXY.txt","data_stressYZ.txt","data_stressXZ.txt",...
            "data_frequencies.txt"];

    mode_list = discobj.calculation_settings.mode_list;    

    fclose(fopen(ScriptFileEdited,'w'));
    for mode_number = mode_list
        copyfile(ScriptFile,ScriptFileEditedTmp);
        % mode number
        ReplaceString = ['set,,,,,,,',num2str(mode_number)];
        text_replace(ScriptFile, ScriptFileEditedTmp, 'set,,,,,,,1', ReplaceString);
        % output folder
        for DataFileName1 = DataFileName_list
            DataFileName = char(DataFileName1);
            SearchString = ['*cfopen,insert-path-here--',DataFileName(1:end-4),',txt'];
            ReplaceString1 = fullfile(OutputFolder,num2str(mode_number),DataFileName(1:end-4));
            ReplaceString = ['*cfopen,''',ReplaceString1,''',txt'];
            text_replace(ScriptFileEditedTmp, ScriptFileEditedTmp, SearchString, ReplaceString);
        end
        % append
        text_append(ScriptFileEdited, ScriptFileEdited, ScriptFileEditedTmp);
    end

    % insert script in input file
    text_add(InputFileEdited, InputFileEdited, '/wb,file,end', ScriptFileEdited);
    
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

