
function workbench_prepare_input_file_DoublyCoatedDisc(discobj)

    % Edits the support files, taken from the templates folder for a resonator of type DoublyCoatedDisc, 
    % and saves them in the folder "support_files";
    % they will be used to generate the Ansys input file "input_file.dat".
    %
    % discobj - instance of class DoublyCoatedDisc.

    TemplatesFolder = discobj.calculation_settings.templates_folder;
    SupportFolder = discobj.calculation_settings.support_folder;

    JournalFileName = discobj.calculation_settings.workbench_input_file_journal_template;
    JournalFileEditedName = discobj.calculation_settings.workbench_input_file_journal_edited;
    GeometryFileName = discobj.calculation_settings.workbench_input_file_geometry_template;
    GeometryFileEditedName = discobj.calculation_settings.workbench_input_file_geometry_edited;
    MechanicalFileName = discobj.calculation_settings.workbench_input_file_mechanical_template;
    MechanicalFileEditedName = discobj.calculation_settings.workbench_input_file_mechanical_edited;

    JournalFile = fullfile(TemplatesFolder,JournalFileName);
    JournalFileEdited = fullfile(SupportFolder,JournalFileEditedName);
    copyfile(JournalFile,JournalFileEdited);

    GeometryFile = fullfile(TemplatesFolder,GeometryFileName);
    GeometryFileEdited = fullfile(SupportFolder,GeometryFileEditedName);
    copyfile(GeometryFile,GeometryFileEdited);

    MechanicalFile = fullfile(TemplatesFolder,MechanicalFileName);
    MechanicalFileEdited = fullfile(SupportFolder,MechanicalFileEditedName);
    copyfile(MechanicalFile,MechanicalFileEdited);

    % Edit Journal file

    SearchStringList = {'script = open(''insert--here--path--to--geom.js'', ''r'')',...
        'scriptm = open(''insert--here--path--to--mech.py'', ''r'')'};

    ReplaceStringList = {['script = open(''',strrep(GeometryFileEdited,'\', '\\'),''', ''r'')'],...
        ['scriptm = open(''',strrep(MechanicalFileEdited,'\', '\\'),''', ''r'')']};

    text_replace(JournalFileEdited, JournalFileEdited, SearchStringList, ReplaceStringList);

    % Edit Geometry file

    radius = discobj.substrate.diameter/2;
    substrate_thickness = discobj.substrate.thickness;

    SearchStringList = {'  p.Cr7 = Circle(0.00000000, 0.00000000, 38.10000000);',...
        '  p.Cr10 = Circle(0.00000000, 0.00000000, 38.10000000);',...
        'var ext1 = agb.Extrude(agc.Add, ps1.Sk1, agc.DirNormal, agc.ExtentFixed, 0.2,',...
        ' pl4.AddTransform(agc.XformZOffset, 0.2);'};

    ReplaceStringList = {['  p.Cr7 = Circle(0.00000000, 0.00000000, ',num2str(radius),');'],...
        ['  p.Cr10 = Circle(0.00000000, 0.00000000, ',num2str(radius),');'],...
        ['var ext1 = agb.Extrude(agc.Add, ps1.Sk1, agc.DirNormal, agc.ExtentFixed, ',num2str(substrate_thickness),','],...
        [' pl4.AddTransform(agc.XformZOffset, ',num2str(substrate_thickness),');']};

    text_replace(GeometryFileEdited, GeometryFileEdited, SearchStringList, ReplaceStringList);

    % Edit Mechanical file

    coating_thickness = discobj.coating.thickness;
    mesh_method = discobj.mesh_settings.method;
    mesh_divisions = discobj.mesh_settings.divisions;
    mesh_fineness = discobj.mesh_settings.fineness;

    InputFolder = discobj.calculation_settings.apdl_input_folder;
    InputFileName = discobj.calculation_settings.apdl_input_file_template;
    InputFile = fullfile(InputFolder,InputFileName);

    SearchStringList = {'coat1.Thickness = Quantity("0.1 [m]")',...
        'coat2.Thickness = Quantity("0.1 [m]")'};

    ReplaceStringList = {['coat1.Thickness = Quantity("',num2str(coating_thickness),' [m]")'],...
        ['coat2.Thickness = Quantity("',num2str(coating_thickness),' [m]")']};

    if strcmp(discobj.substrate.material.flag_anisotropic,'single_crystal') && ~isempty(discobj.substrate.material.stiffness_matrix)

        % cylindrical symmetry: only one vector needed for orientation,
        % oriented along the symmetry axis (= z-axis in Ansys's global coordinate system)

        cylindrical_symmetry_axis = discobj.substrate.material.orientation_axis_1;
        [v1, v2] = find_orthogonal_vectors(cylindrical_symmetry_axis);

        ansys_euler_angles = find_ansys_euler_angles(v1,v2,cylindrical_symmetry_axis);
        eulerZ = ansys_euler_angles(1);
        eulerX = ansys_euler_angles(2);
        eulerY = ansys_euler_angles(3);

        SearchStringList{end+1} = '# Insert--here--substrate--coordinate--system--for--anisotropic--material';

        ReplaceStringList{end+1} = ['coordinate_system_subst = Model.CoordinateSystems.AddCoordinateSystem()' newline ...
            'coordinate_system_subst.OriginDefineBy =  CoordinateSystemAlignmentType.Fixed' newline ...
            'coordinate_system_subst.AddTransformation(TransformationType.Rotation,CoordinateSystemAxisType.PositiveZAxis)' newline ...
            ['coordinate_system_subst.SetTransformationValue(1,' num2str(rad2deg(eulerZ)) ')'] newline ...
            'coordinate_system_subst.AddTransformation(TransformationType.Rotation,CoordinateSystemAxisType.PositiveXAxis)' newline ...
            ['coordinate_system_subst.SetTransformationValue(2,' num2str(rad2deg(eulerX)) ')'] newline ...
            'coordinate_system_subst.AddTransformation(TransformationType.Rotation,CoordinateSystemAxisType.PositiveYAxis)' newline ...
            ['coordinate_system_subst.SetTransformationValue(3,' num2str(rad2deg(eulerY)) ')'] newline ...
            'subst.CoordinateSystem = coordinate_system_subst'];

    end

    if strcmp(discobj.coating.material.flag_anisotropic,'single_crystal') && ~isempty(discobj.coating.material.stiffness_matrix)

        % cylindrical symmetry: only one vector needed for orientation,
        % oriented along the symmetry axis (= z-axis in Ansys's global coordinate system)
        % in the direction of growth of the coating 
        % (opposite for the two coatings)

        cylindrical_symmetry_axis = discobj.coating.material.orientation_axis_1;
        [v1, v2] = find_orthogonal_vectors(cylindrical_symmetry_axis);

        cylindrical_symmetry_axis_rev = -cylindrical_symmetry_axis;
        [v1_rev, v2_rev] = find_orthogonal_vectors(cylindrical_symmetry_axis_rev);

        ansys_euler_angles = find_ansys_euler_angles(v1,v2,cylindrical_symmetry_axis);
        eulerZ = ansys_euler_angles(1);
        eulerX = ansys_euler_angles(2);
        eulerY = ansys_euler_angles(3);
        ansys_euler_angles_rev = find_ansys_euler_angles(v1_rev,v2_rev,cylindrical_symmetry_axis_rev);
        eulerZ_rev = ansys_euler_angles_rev(1);
        eulerX_rev = ansys_euler_angles_rev(2);
        eulerY_rev = ansys_euler_angles_rev(3);    

        SearchStringList{end+1} = '# Insert--here--coating--coordinate--systems--for--anisotropic--material';

        ReplaceStringList{end+1} = ['coordinate_system_coat1 = Model.CoordinateSystems.AddCoordinateSystem()' newline ...
            'coordinate_system_coat1.OriginDefineBy =  CoordinateSystemAlignmentType.Fixed' newline ...
            'coordinate_system_coat1.AddTransformation(TransformationType.Rotation,CoordinateSystemAxisType.PositiveZAxis)' newline ...
            ['coordinate_system_coat1.SetTransformationValue(1,' num2str(rad2deg(eulerZ)) ')'] newline ...
            'coordinate_system_coat1.AddTransformation(TransformationType.Rotation,CoordinateSystemAxisType.PositiveXAxis)' newline ...
            ['coordinate_system_coat1.SetTransformationValue(2,' num2str(rad2deg(eulerX)) ')'] newline ...
            'coordinate_system_coat1.AddTransformation(TransformationType.Rotation,CoordinateSystemAxisType.PositiveYAxis)' newline ...
            ['coordinate_system_coat1.SetTransformationValue(3,' num2str(rad2deg(eulerY)) ')'] newline ...
            'coat1.CoordinateSystem = coordinate_system_coat1' newline ...
            '' newline ...
            'coordinate_system_coat2 = Model.CoordinateSystems.AddCoordinateSystem()' newline ...
            'coordinate_system_coat2.OriginDefineBy =  CoordinateSystemAlignmentType.Fixed' newline ...
            'coordinate_system_coat2.AddTransformation(TransformationType.Rotation,CoordinateSystemAxisType.PositiveZAxis)' newline ...
            ['coordinate_system_coat2.SetTransformationValue(1,' num2str(rad2deg(eulerZ_rev)) ')'] newline ...
            'coordinate_system_coat2.AddTransformation(TransformationType.Rotation,CoordinateSystemAxisType.PositiveXAxis)' newline ...
            ['coordinate_system_coat2.SetTransformationValue(2,' num2str(rad2deg(eulerX_rev)) ')'] newline ...
            'coordinate_system_coat2.AddTransformation(TransformationType.Rotation,CoordinateSystemAxisType.PositiveYAxis)' newline ...
            ['coordinate_system_coat2.SetTransformationValue(3,' num2str(rad2deg(eulerY_rev)) ')'] newline ...
            'coat2.CoordinateSystem = coordinate_system_coat2'];    

    end

    switch mesh_method

        case 'sweep'

            SearchStringList{end+1} = 'mesh_method.SweepNumberDivisions = 3';
            SearchStringList{end+1} = 'mesh.ElementSize = Quantity(''0 [m]'')';

            ReplaceStringList{end+1} = ['mesh_method.SweepNumberDivisions = ',num2str(mesh_divisions)];
            ReplaceStringList{end+1} = ['mesh.ElementSize = Quantity(''',num2str(mesh_fineness),' [m]'')'];

        case 'non-sweep'

            SearchStringList{end+1} = 'mesh_method = mesh.AddAutomaticMethod()';
            SearchStringList{end+1} = 'mesh_method.Location = my_selection';
            SearchStringList{end+1} = 'mesh_method.Method = MethodType.Sweep';
            SearchStringList{end+1} = 'mesh_method.SweepNumberDivisions = 3';
            SearchStringList{end+1} = 'mesh.ElementSize = Quantity(''0 [m]'')';

            ReplaceStringList{end+1} = ''; 
            ReplaceStringList{end+1} = ''; 
            ReplaceStringList{end+1} = ''; 
            ReplaceStringList{end+1} = ''; 
            ReplaceStringList{end+1} = ['mesh.ElementSize = Quantity(''',num2str(mesh_fineness),' [m]'')'];     

    end

    SearchStringList{end+1} = 'an.WriteInputFile(''insert--here--path--to--inputfile.dat'')';

    ReplaceStringList{end+1} = ['an.WriteInputFile(''',strrep(InputFile,'\', '\\'),''')'];

    text_replace(MechanicalFileEdited, MechanicalFileEdited, SearchStringList, ReplaceStringList);

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





