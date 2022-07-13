
function workbench_out = workbench_calculate(discobj)

    % Uses the support files stored in folder "support_files" to
    % generate the Ansys MAPDL input file "input_file.dat",
    % containing information on the geometry, the mesh, and the crystal orientation.
    % (The file will be further edited by the apdl_prepare_input_file functions,
    % with information on the material properties, the solver settings, and the post-processing commands).
    %
    % discobj - instance of class Disc, DoublyCoatedDisc, CurvedDisc, or CantileverFiber (or any other subclass of Resonator).
    % workbench_out - integer (= 0 for a successful calculation).

    SupportFolder = discobj.calculation_settings.support_folder;
    JournalFileEditedName = discobj.calculation_settings.workbench_input_file_journal_edited;
    JournalFile = fullfile(SupportFolder,JournalFileEditedName);

    c = clock();
    disp([newline num2str(c(4),'%02d') ':' num2str(c(5),'%02d') ':' num2str(round(c(6)),'%02d') '  ' num2str(c(3),'%02d') '/' num2str(c(2),'%02d') '/' num2str(c(1),'%04d')]);

    cmd = ['"',discobj.path_workbench_exe,'" -b -r "',JournalFile,'"'];
    disp(cmd);
    workbench_out = system(cmd);

    c = clock();
    disp([char(9) num2str(workbench_out) newline num2str(c(4),'%02d') ':' num2str(c(5),'%02d') ':' num2str(round(c(6)),'%02d') '  ' num2str(c(3),'%02d') '/' num2str(c(2),'%02d') '/' num2str(c(1),'%04d') newline]);

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

