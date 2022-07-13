
function workbench_prepare_input_file_CurvedDisc(discobj)

    % Edits the support files, taken from the templates folder for a resonator of type CurvedDisc, 
    % and saves them in the folder "support_files";
    % they will be used to generate the Ansys input file "input_file.dat".
    %
    % discobj - instance of class CurvedDisc.

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

    curvature_radius = discobj.curvature_radius;
    thickness = discobj.substrate.thickness;
    diameter = discobj.substrate.diameter;

    ext_radius = curvature_radius + thickness/2;
    int_radius = curvature_radius - thickness/2;

    Omega = discobj.calculation_results.curvature_solid_angle;

    Point1 = [0,0];
    Point2 = [0,-(ext_radius-int_radius)];
    Point3 = - ext_radius/(2*pi) * [sqrt(4*pi*Omega-Omega^2),Omega];
    Point4 = - int_radius/(2*pi) * [sqrt(4*pi*Omega-Omega^2),Omega] - [0,ext_radius-int_radius];
    Point5 = [0,-ext_radius];

    SearchString1 = 'with (p.Sk1)';

    ReplaceString1 = ['with (p.Sk1)',newline,...
        '{',newline,...
        '  p.Cr12 = ArcCtrEdge(',newline,...
        '              0.00000000, ',num2str(Point5(2)),',',newline,...
        '              0.00000000, 0.00000000,',newline,...
        '              ',num2str(Point3(1)),', ',num2str(Point3(2)),');',newline,...
        '  p.Cr13 = ArcCtrEdge(',newline,...
        '              0.00000000, ',num2str(Point5(2)),',',newline,...
        '              0.00000000, ',num2str(Point2(2)),',',newline,...
        '              ',num2str(Point4(1)),', ',num2str(Point4(2)),');',newline,...
        '  p.Ln14 = Line(',num2str(Point4(1)),', ',num2str(Point4(2)),', ',num2str(Point3(1)),', ',num2str(Point3(2)),');',newline,...
        '  p.Ln15 = Line(0.00000000, ',num2str(Point2(2)),', 0.00000000, 0.00000000);',newline,...
        '}'];

    SearchString2 = '//Constraints';

    ReplaceString2 = ['  //Constraints',newline,...
      '  VerticalCon(p.Ln15);',newline,...
      '  CoincidentCon(p.Cr12.Center, 0.00000000, ',num2str(Point5(2)),',',newline,...
      '              p.YAxis, 0.00000000, ',num2str(Point5(2)),');',newline,...
      '  CoincidentCon(p.Cr12.Base, 0.00000000, 0.00000000,',newline,... 
      '              p.Origin, 0.00000000, 0.00000000);',newline,...
      '  CoincidentCon(p.Cr13.Center, 0.00000000, ',num2str(Point5(2)),',',newline,... 
      '                p.Cr12.Center, 0.00000000, ',num2str(Point5(2)),');',newline,...
      '  CoincidentCon(p.Cr13.Base, 0.00000000, ',num2str(Point2(2)),',',newline,... 
      '                p.YAxis, 0.00000000, ',num2str(Point2(2)),');',newline,...
      '  CoincidentCon(p.Ln14.Base, ',num2str(Point4(1)),', ',num2str(Point4(2)),',',newline,... 
      '                p.Cr13.End, ',num2str(Point4(1)),', ',num2str(Point4(2)),');',newline,...
      '  CoincidentCon(p.Ln14.End, ',num2str(Point3(1)),', ',num2str(Point3(2)),',',newline,... 
      '                p.Cr12.End, ',num2str(Point3(1)),', ',num2str(Point3(2)),');',newline,...
      '  CoincidentCon(p.Ln15.Base, 0.00000000, ',num2str(Point2(2)),',',newline,... 
      '                p.Cr13.Base, 0.00000000, ',num2str(Point2(2)),');',newline,...
      '  CoincidentCon(p.Ln15.End, 0.00000000, 0.00000000,',newline,... 
      '                p.Cr12.Base, 0.00000000, 0.00000000);',newline,...
      '}'];

    SearchStringList = {SearchString1,SearchString2};
    ReplaceStringList = {ReplaceString1,ReplaceString2};

    text_replace(GeometryFileEdited, GeometryFileEdited, SearchStringList, ReplaceStringList);

    % Edit Mechanical file

    mesh_method = discobj.mesh_settings.method;
    mesh_divisions = discobj.mesh_settings.divisions;
    mesh_fineness = discobj.mesh_settings.fineness;

    InputFolder = discobj.calculation_settings.apdl_input_folder;
    InputFileName = discobj.calculation_settings.apdl_input_file_template;
    InputFile = fullfile(InputFolder,InputFileName);

    SearchStringList = {};
    ReplaceStringList = {};

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


