
function apdl_out = apdl_calculate(discobj)

    % Executes with the Ansys APDL solver the calculation specified in "input_file_edited.dat" 
    % (file located within input_folder).
    % Stores the results within output_folder, deleting old results if present.
    %
    % discobj - instance of class Disc, DoublyCoatedDisc, CurvedDisc, or CantileverFiber (or any other subclass of Resonator).
    % apdl_out - integer (= 0 for a successful calculation).

    InputFolder = discobj.calculation_settings.apdl_input_folder;
    InputFileEditedName = discobj.calculation_settings.apdl_input_file_edited;
    InputFileEdited = fullfile(InputFolder,InputFileEditedName);

    OutputFolder = discobj.calculation_settings.apdl_output_folder;
    OutputFile = fullfile(OutputFolder,'output_file.out');

    % DELETE ALL CONTENTS IN THE OUTPUT_FILES FOLDER, TO AVOID LEFTOVER DATA FROM PREVIOUS CALCULATIONS
    % Comment out the following lines if you wish to keep them
    rmdir(OutputFolder,'s');
    mkdir(OutputFolder);

    mode_list = discobj.calculation_settings.mode_list;    
    for mode_number = mode_list
        if ~exist(fullfile(OutputFolder,num2str(mode_number)),'dir')
            mkdir(fullfile(OutputFolder,num2str(mode_number)));
        end       
    end

    % ANSYS EXECUTE

    old_wd = pwd;
    cd(discobj.calculation_settings.ansys_folder);

    c = clock();
    disp([num2str(c(4),'%02d') ':' num2str(c(5),'%02d') ':' num2str(round(c(6)),'%02d') '  ' num2str(c(3),'%02d') '/' num2str(c(2),'%02d') '/' num2str(c(1),'%04d')]);

    cmd = ['"',discobj.path_ansys_exe,'" -b -i "',InputFileEdited,'" -o "',OutputFile,'"'];
    disp(cmd);
    apdl_out = system(cmd);

    c = clock();
    disp([char(9) num2str(apdl_out) newline num2str(c(4),'%02d') ':' num2str(c(5),'%02d') ':' num2str(round(c(6)),'%02d') '  ' num2str(c(3),'%02d') '/' num2str(c(2),'%02d') '/' num2str(c(1),'%04d') newline]);

    cd(old_wd);

    if apdl_out ~= 0
        disp('Ansys calculation failed. Consider deleting file.lock in subfolder ansys_files.');
        return;
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









