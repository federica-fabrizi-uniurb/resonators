
classdef Material 
    
    % Class containing information on the material properties. 
    % Note: detailed explanations on the conventions used in this class for
    % tensorial properties and rotations can be found at the end of this file.
    
    properties
       
        keyword = [];               % string
        density                     % [kg m^-3]
        specific_heat_capacity      % [J K^-1 kg^-1] 
        
        flag_anisotropic (:,:) char {mustBeMember(flag_anisotropic,{'amorph','poly','single_crystal'})} = 'amorph'      % string
        
        % isotropic properties
        young_amorph                                % [Pa]
        young_poly                                  % [Pa]
        poisson_amorph                              % [dimensionless]
        poisson_poly                                % [dimensionless]
        bulk_amorph                                 % [Pa]
        bulk_poly                                   % [Pa]
        thermal_linear_expansion_amorph             % [K^-1]
        thermal_linear_expansion_poly               % [K^-1]
        thermal_conductivity_amorph                 % [W m^-1 K^-1]
        thermal_conductivity_poly                   % [W m^-1 K^-1]
        
        % directional properties
        young_directional                           % dictionary [Pa] 
        poisson_directional                         % dictionary [dimensionless] 
        thermal_linear_expansion_directional        % dictionary [K^-1] 
        thermal_conductivity_directional            % dictionary [W m^-1 K^-1] 
        
        % anisotropic properties
        orientation_axis_1                          % 1x3 or 3x1 [dimensionless] direction of highest symmetry 
        orientation_axis_2                          % 1x3 or 3x1 [dimensionless] direction of second highest symmetry
        stiffness_matrix                            % 6x6 [Pa]
        thermal_linear_expansion_matrix             % 3x3 [K^-1]
        thermal_conductivity_matrix                 % 3x3 [W m^-1 K^-1]

        material_folder

    end
    
    properties (SetAccess = private)
        
        temperature
        % temperature dependence, for the material identified by keyword
        density_data                     % [kg m^-3]
        specific_heat_capacity_data      % [J K^-1 kg^-1] 
        % isotropic properties
        young_amorph_data                % [Pa]
        young_poly_data                  % [Pa] 
        poisson_amorph_data              % [dimensionless]
        poisson_poly_data                % [dimensionless]
        bulk_data                        % [Pa]
        thermal_linear_expansion_amorph_data            % [K^-1]
        thermal_linear_expansion_poly_data              % [K^-1]
        thermal_conductivity_amorph_data                % [W m^-1 K^-1]
        thermal_conductivity_poly_data                  % [W m^-1 K^-1]
        % directional properties
        young_directional_data                          % dictionary structure [Pa] 
        poisson_directional_data                        % dictionary structure [dimensionless] 
        thermal_linear_expansion_directional_data       % dictionary structure [K^-1] 
        thermal_conductivity_directional_data           % dictionary structure [W m^-1 K^-1] 
        % anisotropic properties         
        stiffness_matrix_data                           % 6x6 [Pa]
        thermal_linear_expansion_matrix_data            % 3x3 [K^-1]
        thermal_conductivity_matrix_data                % 3x3 [W m^-1 K^-1]          
        
    end
    
    methods
        
        function obj = Material(keyword)
            
            % constructor function, automatically called at instantiation
            
            % you can call this constructor without any input argument (no keyword);
            % if you define a keyword, the object ges re-instantiated from
            % scratch, to avoid leftover data from previous assignments
            
            parentdir = fileparts(fileparts(mfilename('fullpath')));
            obj.material_folder = fullfile(parentdir,'temperature_dependent_material_properties');
            % Initialise all directional properties as structures
            property_keys = properties(obj);
            for jj = 1:length(property_keys) 
               property = property_keys{jj};
               if contains(property,'directional') 
                   obj.(property) = struct();
               end
            end  
            % Load database if a keyword is specified in input
            if nargin == 1
                obj.keyword = keyword;
                obj = local_set_from_database(obj,keyword);
            end
            
        end
        
        function obj = extend_to_zero_temperature(obj)
            
            obj = local_extend_to_zero_temperature(obj);
            
        end
        
        function obj = set_temperature(obj,temperature)
            
            obj.temperature = temperature;
            if ~isempty(obj.keyword) 
                obj = local_set_temperature(obj,temperature);
            end
            
        end
        
        % Note: in the following, in some variable names and function names, 
        % 'alpha' replaces 'thermal_linear_expansion' and
        % 'k' replaces 'thermal_conductivity' and
        % 'B' replaces 'bulk'
        % for brevity.
        
        function stiffness_matrix_rotated = rotate_stiffness_matrix(obj,eulerZ,eulerY,eulerX)
            
            stiffness_matrix_rotated = local_rotate_stiffness_matrix(obj,eulerZ,eulerY,eulerX);
            %stiffness_matrix_rotated = local_rotate_stiffness_matrix_v2(obj,eulerZ,eulerY,eulerX);  % if a check is needed
            
        end
        
        function out = young_directional_from_stiffness_matrix(obj,direction_vector)
            
            out = local_young_directional_from_stiffness_matrix(obj,direction_vector);
            
        end
        
        function k_matrix_rotated = rotate_k_matrix(obj,eulerZ,eulerY,eulerX)
            
            k_matrix_rotated = local_rotate_k_matrix(obj,eulerZ,eulerY,eulerX);
            
        end
        
        function out = k_directional_from_k_matrix(obj,direction_vector)
            
            out = local_k_directional_from_k_matrix(obj,direction_vector);
            
        end
      
        function alpha_matrix_rotated = rotate_alpha_matrix(obj,eulerZ,eulerY,eulerX)
            
            alpha_matrix_rotated = local_rotate_alpha_matrix(obj,eulerZ,eulerY,eulerX);
            
        end
        
        function out = alpha_directional_from_alpha_matrix(obj,direction_vector)
            
            out = local_alpha_directional_from_alpha_matrix(obj,direction_vector);
            
        end        
        
        function young_poly = young_poly_from_stiffness_matrix(obj)
            
            young_poly = local_young_poly_from_stiffness_matrix(obj);
            
        end
        
        function poisson_poly = poisson_poly_from_stiffness_matrix(obj)
            
            poisson_poly = local_poisson_poly_from_stiffness_matrix(obj);
            
        end
        
        function bulk_poly = bulk_poly_from_stiffness_matrix(obj)
            
            bulk_poly = local_bulk_poly_from_stiffness_matrix(obj);
            
        end      
        
        function bulk_sc = bulk_single_crystal_from_stiffness_matrix(obj)
            
            bulk_sc = local_bulk_single_crystal_from_stiffness_matrix(obj);
            
        end 
        
        function B_amorph = bulk_amorph_from_young_amorph_and_poisson_amorph(obj)
            
            B_amorph = local_bulk_amorph_from_young_amorph_and_poisson_amorph(obj);
            
        end
        
        function B_poly = bulk_poly_from_young_poly_and_poisson_poly(obj)
            
            B_poly = local_bulk_poly_from_young_poly_and_poisson_poly(obj);
            
        end
        
        function alpha_poly = alpha_poly_from_alpha_matrix_and_stiffness_matrix(obj)
            
            alpha_poly = local_alpha_poly_from_alpha_matrix_and_stiffness_matrix(obj);
            
        end
        
    end
    
    methods(Static)
        
        % Static methods are associated with a class, but not with specific instances of that class. 
        % These methods do not require an object of the class as an input argument. 
        
        function compliance_voigt = compliance_tensor_to_voigt(compliance_tensor)
            
           compliance_voigt = local_compliance_tensor_to_voigt(compliance_tensor);
            
        end
        
        function stiffness_voigt = stiffness_tensor_to_voigt(stiffness_tensor)
            
           stiffness_voigt = local_stiffness_tensor_to_voigt(stiffness_tensor);
            
        end
        
        function compliance_tensor = compliance_voigt_to_tensor(compliance_voigt)
            
            compliance_tensor = local_compliance_voigt_to_tensor(compliance_voigt);
            
        end
        
        function stiffness_tensor = stiffness_voigt_to_tensor(stiffness_voigt)
            
            stiffness_tensor = local_stiffness_voigt_to_tensor(stiffness_voigt);
            
        end
        
    end
    
end


% Local functions (only callable from within this same .m file)

function obj = local_set_from_database(obj,keyword)

    DataFolder = fullfile(obj.material_folder,keyword);

    file_info = dir(DataFolder);
    file_keys = {file_info.name};
    for jj = 1:length(file_keys) 
        file = file_keys{jj}; 
        if contains(file,'txt') && ~contains(file,'readme')
            if ~contains(file,'directional')
                if contains(file,'matrix')
                    % anisotropic properties (matrix properties)
                    try
                        opts = detectImportOptions(fullfile(DataFolder,file),'Whitespace',';,[]\t','ConsecutiveDelimitersRule','join','ExtraColumnsRule','ignore');
                        data = readtable(fullfile(DataFolder,file),opts);
                        data = table2array(data);
                        % if data is a scalar, re-construct full matrix for thermal_conductivity and thermal_linear_expansion
                        % (still a single crystal tensor but with isotropic symmetry)
                        if size(data,2)==2 && contains(file,["thermal_conductivity","thermal_linear_expansion"])
                            zz = zeros(size(data,1),1);
                            data(:,[2:10])=[data(:,2),zz,zz,zz,data(:,2),zz,zz,zz,data(:,2)];
                        end
                        obj.([file(1:end-4) '_data']) = data;
                    catch
                        disp(['File ' file ' cannot be matched to any material property']);
                    end
                else
                    % isotropic properties (non-directional, non-matrix properties)
                    try
                        obj.([file(1:end-4) '_data']) = dlmread(fullfile(DataFolder,file)); 
                    catch
                        disp(['File ' file ' cannot be matched to any material property']);
                    end                    
                end
            else
                % directional properties
                try
                    file_parsed = split(file,'_directional_');
                    direction = split(file_parsed{2},'.');
                    direction = direction{1};
                    property = file_parsed{1};
                    [obj.([property '_directional_data']).(['direction_' strrep(direction,'-','m')])] = dlmread(fullfile(DataFolder,file));   
                catch
                    disp(['File ' file ' cannot be matched to any material property']);
                end
            end
        end
    end
    
end

function obj = local_extend_to_zero_temperature(obj)

    property_keys = properties(obj);
    for jj = 1:length(property_keys) 
       property = property_keys{jj};
       if contains(property,'data') && ~isempty(obj.(property)) 
           if ~contains(property,'directional')
               % non-directional properties (isotropic and anisotropic properties)
               obj.(property)(1,1) = 0; % extend Temperature range as (almost) constant
           else
               % directional properties
               directions = fieldnames(obj.(property));
               for kk = 1:length(directions) 
                   direction = directions{kk};
                   obj.(property).(direction)(1,1) = 0; % extend Temperature range as (almost) constant
               end
           end
       end
    end  
    
end

function obj = local_set_temperature(obj,temperature)

    property_keys = properties(obj);
    for jj = 1:length(property_keys) 
       property = property_keys{jj};
       if contains(property,'data') && ~isempty(obj.(property)) 
           property_current_temperature = erase(property,'_data');
           if ~contains(property,'directional')
               if contains(property,'matrix')
                   % anisotropic properties (matrix properties)
                   data = zeros(1,size(obj.(property),2)-1);
                   for ll = 1:(size(obj.(property),2)-1)
                       data(ll) = interp1(obj.(property)(:,1),obj.(property)(:,1+ll),temperature);
                   end
                   data = reshape(data,sqrt(size(obj.(property),2)-1),sqrt(size(obj.(property),2)-1))';
                   obj.(property_current_temperature) = data;
               else
                   % isotropic properties (non-directional, non-matrix properties)
                   obj.(property_current_temperature) = interp1(obj.(property)(:,1),obj.(property)(:,2),temperature);
               end
           else
               % directional properties
               directions = fieldnames(obj.(property));
               for kk = 1:length(directions) 
                   direction = directions{kk};
                   obj.(property_current_temperature).(direction) = interp1(obj.(property).(direction)(:,1),obj.(property).(direction)(:,2),temperature); 
               end
           end
       end
    end     

end

% Functions for anisotropic materials

function stiffness_matrix_rotated = local_rotate_stiffness_matrix(obj,eulerZ,eulerY,eulerX)

    % References:
    % Healy et al., Solid Earth, 11, 259–286
    % Zhang et al., Synchrotron Rad. (2014) 21, 507–517
 
    % These are supposed to be Active Euler angles (in degrees)
    
    % Convert in radians
    eulerZ = deg2rad(eulerZ);
    eulerY = deg2rad(eulerY);
    eulerX = deg2rad(eulerX);
    
    % Active rotation matrix 
    R = eul2rotm([eulerZ eulerY eulerX],'ZYX');
    
    % Note: in Healy et al., Eq. (6),
    % the rotation matrix found ("a") is referred to:
    % - rotations about: 1st) +Z; 2nd) -Y; 3rd) +X  (note the different sign for Y)
    % - intrinsic (from Fig. 1)
    % - angles are defined for a passive transformation (from Fig. 1), but
    % - the matrix looked for is the one of an active transformation (it is then used to change the coordinates of the tensor)
    
    % Given this, "a" is consistent with:
    % matrix Z_1 Y_2 X_3 found at https://en.wikipedia.org/wiki/Euler_angles
    % or
    % a = eul2rotm([-phi_z phi_y -phi_x],'ZYX')
    % where
    % all angles have their sign changed to go from passive to active &
    % phi_y has its sign changed again, to account for rotation by -Y

    % The correctedness of our conventions is confirmed by reproducing the results of
    % Zhang et al., Synchrotron Rad. (2014). 21, 507–517
    % for silicon Young's modulus: E100 = 130, E110 = 169, E111 = 188, E311 = 152 GPa (p. 509).

    % construct full tensor and rotate
    
    C = local_stiffness_voigt_to_tensor(obj.stiffness_matrix);
    
    C_rotated = zeros(3,3,3,3);
    for ii =1:3
        for jj = 1:3
            for kk = 1:3
                for ll = 1:3
                    
                    for mm = 1:3
                        for nn = 1:3
                            for oo = 1:3
                                for pp = 1:3
                                    
                                    C_rotated(ii,jj,kk,ll) = C_rotated(ii,jj,kk,ll) + ...
                                        R(ii,mm)*R(jj,nn)*R(kk,oo)*R(ll,pp)*C(mm,nn,oo,pp); 
                                    
%                                      disp([ii jj kk ll]); 
%                                      disp([mm nn oo pp]); 
%                                      disp(r(ii,mm)*r(jj,nn)*r(kk,oo)*r(ll,pp)*C(mm,nn,oo,pp));
%                                      disp('');
                                    
                                end
                            end
                        end
                    end
                    
                end
            end
        end
    end
    
    % reconvert the full tensor to Voigt notation
    
    stiffness_matrix_rotated = local_stiffness_tensor_to_voigt(C_rotated);

end

function stiffness_matrix_rotated = local_rotate_stiffness_matrix_v2(obj,eulerZ,eulerY,eulerX)

    % This function is an alternative way of performing the same calculation 
    % performed by function "local_aniso_rotate_stiffness_matrix".
    % Here the rotation is done directly on the Voigt coefficients,
    % without transforming to full 4th-rank tensor and back.
    % The purpose is to check the correctedness of our calculations and 
    % the consistency of the conventions used.
    
    % References:
    % Minhang Bao, Analysis and Design Principles of MEMS Devices (2005) - Section 6.2    
    
    % These are supposed to be Active Euler angles (in degrres)
    
    % Convert in radians
    eulerZ = deg2rad(eulerZ);
    eulerY = deg2rad(eulerY);
    eulerX = deg2rad(eulerX);
    
    % Active rotation matrix 
    R = eul2rotm([eulerZ eulerY eulerX],'ZYX');
    
    % Note: 
    % the compliance_matrix(I,J) and stiffness_matrix(I,J) in Voigt form (I,J = 1 ... 6)
    % do NOT transform in the same way under rotation.
    % (This is because the tensorial basis of this representation is not normalised, 
    % and therefore the distinction between contra-variant and co-variant components becomes important,
    %  see R.M. Brannon, "Rotation, Reflection, and Frame Changes", IOP Publishing - Chapter 26).
    % Instead, the following transformations can be used (symbol ' denotes the transpose):
    
    % stress_rotated = alpha_matrix * stress
    
    % strain_rotated = beta_matrix  * strain
    
    % stiffness_matrix_rotated = alpha_matrix * stiffness_matrix * alpha_matrix'    =
    %                            alpha_matrix * stiffness_matrix * inv(beta_matrix)
    
    % compliance_matrix_rotated = beta_matrix * compliance_matrix * beta_matrix'    =
    %                             beta_matrix * compliance_matrix * inv(alpha_matrix)
    
    % where: the alpha and beta matrix are defined from the rotation matrix as shown below, 
    % and have the properties:
    % inv(alpha_matrix) =  beta_matrix'
    % inv(beta_matrix)  =  alpha_matrix'
    
    l1 = R(1,1);  m1 = R(1,2);  n1 = R(1,3);
    l2 = R(2,1);  m2 = R(2,2);  n2 = R(2,3);
    l3 = R(3,1);  m3 = R(3,2);  n3 = R(3,3);

    % Voigt notation: for rotation of stress (as voigt 6x1 vector) and stiffness

    alpha_matrix = [ ...
        [ l1^2   m1^2    n1^2    2*m1*n1       2*n1*l1       2*l1*m1     ]; ...
        [ l2^2   m2^2    n2^2    2*m2*n2       2*n2*l2       2*l2*m2     ]; ...
        [ l3^2   m3^2    n3^2    2*m3*n3       2*n3*l3       2*l3*m3     ]; ...
        [ l2*l3  m2*m3   n2*n3   m2*n3+m3*n2   n2*l3+n3*l2   m2*l3+m3*l2 ]; ...
        [ l3*l1  m3*m1   n3*n1   m3*n1+m1*n3   n3*l1+n1*l3   m3*l1+m1*l3 ]; ...
        [ l1*l2  m1*m2   n1*n2   m1*n2+m2*n1   n1*l2+n2*l1   m1*l2+m2*l1 ]; ...
    ];

    alpha_matrix_inv = [ ...
        [ l1^2   l2^2    l3^2    2*l2*l3       2*l3*l1       2*l1*l2     ]; ...
        [ m1^2   m2^2    m3^2    2*m2*m3       2*m3*m1       2*m1*m2     ]; ...
        [ n1^2   n2^2    n3^2    2*n2*n3       2*n3*n1       2*n1*n2     ]; ...
        [ m1*n1  m2*n2   m3*n3   m2*n3+m3*n2   m3*n1+m1*n3   m1*n2+m2*n1 ]; ...
        [ n1*l1  n2*l2   n3*l3   n2*l3+n3*l2   n3*l1+n1*l3   n1*l2+n2*l1 ]; ...
        [ l1*m1  l2*m2   l3*m3   m2*l3+m3*l2   m3*l1+m1*l3   m1*l2+m2*l1 ]; ...
    ];

    % Voigt notation: for rotation of strain (as voigt 6x1 vector) and compliance

    beta_matrix = [ ...
        [ l1^2     m1^2      n1^2      m1*n1         n1*l1         l1*m1     ]; ...
        [ l2^2     m2^2      n2^2      m2*n2         n2*l2         l2*m2     ]; ...
        [ l3^2     m3^2      n3^2      m3*n3         n3*l3         l3*m3     ]; ...
        [ 2*l2*l3  2*m2*m3   2*n2*n3   m2*n3+m3*n2   n2*l3+n3*l2   m2*l3+m3*l2 ]; ...
        [ 2*l3*l1  2*m3*m1   2*n3*n1   m3*n1+m1*n3   n3*l1+n1*l3   m3*l1+m1*l3 ]; ...
        [ 2*l1*l2  2*m1*m2   2*n1*n2   m1*n2+m2*n1   n1*l2+n2*l1   m1*l2+m2*l1 ]; ...
    ];

    beta_matrix_inv = [ ...
        [ l1^2     l2^2      l3^2       l2*l3          l3*l1          l1*l2     ]; ...
        [ m1^2     m2^2      m3^2       m2*m3          m3*m1          m1*m2     ]; ...
        [ n1^2     n2^2      n3^2       n2*n3          n3*n1          n1*n2     ]; ...
        [ 2*m1*n1  2*m2*n2   2*m3*n3   m2*n3+m3*n2   m3*n1+m1*n3   m1*n2+m2*n1 ]; ...
        [ 2*n1*l1  2*n2*l2   2*n3*l3   n2*l3+n3*l2   n3*l1+n1*l3   n1*l2+n2*l1 ]; ...
        [ 2*l1*m1  2*l2*m2   2*l3*m3   m2*l3+m3*l2   m3*l1+m1*l3   m1*l2+m2*l1 ]; ...
    ];


    stiffness_matrix = obj.stiffness_matrix;
    
    stiffness_matrix_rotated = alpha_matrix * stiffness_matrix * inv(beta_matrix);

end

function compliance_voigt = local_compliance_tensor_to_voigt(compliance_tensor)

    S = compliance_tensor;
    compliance_voigt = [ [   S(1,1,1,1)    S(1,1,2,2)    S(1,1,3,3)  2*S(1,1,2,3)  2*S(1,1,1,3)  2*S(1,1,1,2) ]; ...
                         [   S(2,2,1,1)    S(2,2,2,2)    S(2,2,3,3)  2*S(2,2,2,3)  2*S(2,2,1,3)  2*S(2,2,1,2) ]; ...
                         [   S(3,3,1,1)    S(3,3,2,2)    S(3,3,3,3)  2*S(3,3,2,3)  2*S(3,3,1,3)  2*S(3,3,1,2) ]; ...
                         [ 2*S(2,3,1,1)  2*S(2,3,2,2)  2*S(2,3,3,3)  4*S(2,3,2,3)  4*S(2,3,1,3)  4*S(2,3,1,2) ]; ...
                         [ 2*S(1,3,1,1)  2*S(1,3,2,2)  2*S(1,3,3,3)  4*S(1,3,2,3)  4*S(1,3,1,3)  4*S(1,3,1,2) ]; ...
                         [ 2*S(1,2,1,1)  2*S(1,2,2,2)  2*S(1,2,3,3)  4*S(1,2,2,3)  4*S(1,2,1,3)  4*S(1,2,1,2) ]; ...
    ];
    
end

function stiffness_voigt = local_stiffness_tensor_to_voigt(stiffness_tensor)

    C = stiffness_tensor;
    stiffness_voigt = [  [   C(1,1,1,1)    C(1,1,2,2)    C(1,1,3,3)    C(1,1,2,3)    C(1,1,1,3)    C(1,1,1,2) ]; ...
                         [   C(2,2,1,1)    C(2,2,2,2)    C(2,2,3,3)    C(2,2,2,3)    C(2,2,1,3)    C(2,2,1,2) ]; ...
                         [   C(3,3,1,1)    C(3,3,2,2)    C(3,3,3,3)    C(3,3,2,3)    C(3,3,1,3)    C(3,3,1,2) ]; ...
                         [   C(2,3,1,1)    C(2,3,2,2)    C(2,3,3,3)    C(2,3,2,3)    C(2,3,1,3)    C(2,3,1,2) ]; ...
                         [   C(1,3,1,1)    C(1,3,2,2)    C(1,3,3,3)    C(1,3,2,3)    C(1,3,1,3)    C(1,3,1,2) ]; ...
                         [   C(1,2,1,1)    C(1,2,2,2)    C(1,2,3,3)    C(1,2,2,3)    C(1,2,1,3)    C(1,2,1,2) ]; ...
    ];
    
end

function compliance_tensor = local_compliance_voigt_to_tensor(compliance_voigt)

    compliance_tensor = zeros(3,3,3,3);

    for ii =1:3
        for jj = 1:3
            I = local_condense_indexes(ii,jj);
            for kk = 1:3
                for ll = 1:3
                    J = local_condense_indexes(kk,ll);
                    if I <= 3 && J <= 3
                        factor = 1;
                    else
                       if I >= 4 && J >= 4
                           factor = 1/4;
                       else
                           factor = 1/2;
                       end
                    end
                    compliance_tensor(ii,jj,kk,ll) = factor*compliance_voigt(I,J); 
                end
            end
        end
    end
    
end

function stiffness_tensor = local_stiffness_voigt_to_tensor(stiffness_voigt)

    stiffness_tensor = zeros(3,3,3,3);

    for ii =1:3
        for jj = 1:3
            I = local_condense_indexes(ii,jj);
            for kk = 1:3
                for ll = 1:3
                    J = local_condense_indexes(kk,ll);
                    stiffness_tensor(ii,jj,kk,ll) = stiffness_voigt(I,J); 
                end
            end
        end
    end
    
end

function INDEX = local_condense_indexes(index_aa,index_bb)

    % 1 = (1, 1), 2 = (2,2), 3 = (3, 3), 4 = (2, 3) = (3, 2), 5 = (1, 3) = (3, 1), and 6 = (1, 2) = (2, 1)
    
    switch true
        case index_aa==1 && index_bb==1
            INDEX = 1;
        case index_aa==2 && index_bb==2
            INDEX = 2;
        case index_aa==3 && index_bb==3
            INDEX = 3;
        case (index_aa==2 && index_bb==3) || (index_aa==3 && index_bb==2)
            INDEX = 4;
        case (index_aa==1 && index_bb==3) || (index_aa==3 && index_bb==1) 
            INDEX = 5;
        case (index_aa==1 && index_bb==2) || (index_aa==2 && index_bb==1) 
            INDEX = 6;
    end

end

function [index_aa index_bb] = local_expand_indexes(INDEX)

    % 1 = (1, 1), 2 = (2,2), 3 = (3, 3), 4 = (2, 3) = (3, 2), 5 = (1, 3) = (3, 1), and 6 = (1, 2) = (2, 1)
    
    switch INDEX
        case 1
            index_aa = 1;
            index_bb = 1;
        case 2
            index_aa = 2;
            index_bb = 2;
        case 3
            index_aa = 3;
            index_bb = 3;           
        case 4
            index_aa = 2;
            index_bb = 3;        
        case 5
            index_aa = 1;
            index_bb = 3;
        case 6
            index_aa = 1;
            index_bb = 2;            
    end

end

function young_directional = local_young_directional_from_stiffness_matrix(obj,direction_vector)

    % Express direction_vector in coordinates of the Material Reference Frame.
    % It does not need to be normalised.

    % Method: rotate the 6x6 tensor (Active rotation) to make
    % direction_vector coincide with [0 0 1] of the material reference frame;
    % then get component (3,3) of the 6x6 rotated tensor.       
    % (It is easier to extract the Young's modulus along one of the frame axes, from the compliance matrix).
    v3 = direction_vector/norm(direction_vector);
    [v1, v2] = find_orthogonal_vectors(v3);
    % Active rotation that brings vector/tensor component along v3 in [0 0 1]
    R = inv([v1(1) v2(1) v3(1); v1(2) v2(2) v3(2); v1(3) v2(3) v3(3)]);
    
    % Pass the Euler angles of the transformation instead of the matrix,
    % because:
    % - easier to visualise/check;
    % - in accordance with most of the literature.
    
    % Active rotation matrix (angles in radians)
	euler_angles = rotm2eul(R,'ZYX');
    eulerZ = euler_angles(1);
    eulerY = euler_angles(2);
    eulerX = euler_angles(3);

    % Active rotation (pass angles in degrees)
    stiffness_matrix_rotated = obj.rotate_stiffness_matrix(rad2deg(eulerZ),rad2deg(eulerY),rad2deg(eulerX));
    compliance_matrix_rotated = inv(stiffness_matrix_rotated);
    
    young_001 = 1/compliance_matrix_rotated(3,3);
    
    young_directional = young_001;

end

function k_matrix_rotated = local_rotate_k_matrix(obj,eulerZ,eulerY,eulerX)

    % These are supposed to be Active Euler angles (in degrees)
    
    % Convert in radians
    eulerZ = deg2rad(eulerZ);
    eulerY = deg2rad(eulerY);
    eulerX = deg2rad(eulerX);
    
    % Active rotation matrix 
    R = eul2rotm([eulerZ eulerY eulerX],'ZYX');
    
    % Rotated matrix (2nd rank tensor)
    k_matrix_rotated = R*obj.thermal_conductivity_matrix*inv(R);
    
end

function k_directional = local_k_directional_from_k_matrix(obj,direction_vector)

    % Express direction_vector in coordinates of the Material Reference Frame.
    % It does not need to be normalised.
    
    % Method 1: 
    vn = direction_vector/norm(direction_vector);
    
    k_directional = 0;
    k_matrix = obj.thermal_conductivity_matrix;
    for ii =1:3
        for jj = 1:3
            k_directional = k_directional + k_matrix(ii,jj)*vn(ii)*vn(jj);    
        end
    end
     
    %%% as a check (must be k_directional = k_directional2)
    
    % Method 2: rotate the tensor (Active rotation) to make
    % direction_vector coincide with [0 0 1] of the material reference frame;
    % then get component (3,3) of the rotated tensor.   
    v3 = direction_vector/norm(direction_vector);
    [v1, v2] = find_orthogonal_vectors(v3);
    % Active rotation that brings vector/tensor component along v3 in [0 0 1]
    R = inv([v1(1) v2(1) v3(1); v1(2) v2(2) v3(2); v1(3) v2(3) v3(3)]);
    
    % Pass the Euler angles of the transformation instead of the matrix,
    % because:
    % - easier to visualise/check;
    % - in accordance with most of the literature.
    
    % Active Euler angles (in radians)
	euler_angles = rotm2eul(R,'ZYX');
    eulerZ = euler_angles(1);
    eulerY = euler_angles(2);
    eulerX = euler_angles(3);

    % Active rotation (pass angles in degrees)
    thermal_conductivity_matrix_rotated = obj.rotate_k_matrix(rad2deg(eulerZ),rad2deg(eulerY),rad2deg(eulerX));
    
    thermal_conductivity_001 = thermal_conductivity_matrix_rotated(3,3);
    
    k_directional2 = thermal_conductivity_001;
    
end

function alpha_matrix_rotated = local_rotate_alpha_matrix(obj,eulerZ,eulerY,eulerX)

    % These are supposed to be Active Euler angles (in degrees)
    
    % Convert in radians
    eulerZ = deg2rad(eulerZ);
    eulerY = deg2rad(eulerY);
    eulerX = deg2rad(eulerX);
    
    % Active rotation matrix
    R = eul2rotm([eulerZ eulerY eulerX],'ZYX');
    
    % Rotated matrix (2nd rank tensor)
    alpha_matrix_rotated = R*obj.thermal_linear_expansion_matrix*inv(R);
    
end

function alpha_directional = local_alpha_directional_from_alpha_matrix(obj,direction_vector)

    % Express direction_vector in coordinates of the Material Reference Frame.
    % It does not need to be normalised.
    
    % Method 1: 
    vn = direction_vector/norm(direction_vector);
    
    alpha_directional = 0;
    alpha_matrix = obj.thermal_linear_expansion_matrix;
    for ii =1:3
        for jj = 1:3
            alpha_directional = alpha_directional + alpha_matrix(ii,jj)*vn(ii)*vn(jj);    
        end
    end
    
    %%% as a check (must be alpha_directional = alpha_directional2)
    
    % Method 2: rotate the tensor (Active rotation) to make
    % direction_vector coincide with [0 0 1] of the material reference frame;
    % then get component (3,3) of the rotated tensor.
    v3 = direction_vector/norm(direction_vector);
    [v1, v2] = find_orthogonal_vectors(v3);
    % Active rotation that brings vector/tensor component along v3 in [0 0 1]
    R = inv([v1(1) v2(1) v3(1); v1(2) v2(2) v3(2); v1(3) v2(3) v3(3)]);
    
    % Pass the Euler angles of the transformation instead of the matrix,
    % because:
    % - easier to visualise/check;
    % - in accordance with most of the literature.

    % Active Euler angles (in radians)
	euler_angles = rotm2eul(R,'ZYX');
    eulerZ = euler_angles(1);
    eulerY = euler_angles(2);
    eulerX = euler_angles(3);

    % Active rotation (pass angles in degrees)
    thermal_linear_expansion_matrix_rotated = obj.rotate_alpha_matrix(rad2deg(eulerZ),rad2deg(eulerY),rad2deg(eulerX));
      
    thermal_linear_expansion_001 = thermal_linear_expansion_matrix_rotated(3,3);
    
    alpha_directional2 = thermal_linear_expansion_001;    

end

function young_poly = local_young_poly_from_stiffness_matrix(obj)
    
    % References:
    % R Hill 1952 Proc. Phys. Soc. A 65 349
    % J. Am. Ceram. Soc., 96 [12] 3891–3900 (2013); DOI: 10.1111/jace.12618
    % Note: Eq. (5) in Reference 2 (for G_reuss) is missing a factor 15
    % as seen from Reference 1 (Hill) and from numerical checks
    
    stiffness_matrix = obj.stiffness_matrix;    
    compliance_matrix = inv(obj.stiffness_matrix);
    
    B_voigt = 1/9*(stiffness_matrix(1,1) + stiffness_matrix(2,2) + stiffness_matrix(3,3)) + ...
              2/9*(stiffness_matrix(1,2) + stiffness_matrix(1,3) + stiffness_matrix(2,3));
          
    G_voigt = 1/15*(stiffness_matrix(1,1) + stiffness_matrix(2,2) + stiffness_matrix(3,3) - stiffness_matrix(1,2) - stiffness_matrix(1,3) - stiffness_matrix(2,3)) + ...
              1/5 *(stiffness_matrix(4,4) + stiffness_matrix(5,5) + stiffness_matrix(6,6));
          
    compressibility_reuss =    compliance_matrix(1,1) + compliance_matrix(2,2) + compliance_matrix(3,3) + ...
                            2*(compliance_matrix(1,2) + compliance_matrix(1,3) + compliance_matrix(2,3));
          
    B_reuss = 1/compressibility_reuss;
    
    % Notice factor 15, corrected from Reference 2 
    G_reuss = 15/( 4*(compliance_matrix(1,1) + compliance_matrix(2,2) + compliance_matrix(3,3)) + ...
                  -4*(compliance_matrix(1,2) + compliance_matrix(1,3) + compliance_matrix(2,3)) + ...
                   3*(compliance_matrix(4,4) + compliance_matrix(5,5) + compliance_matrix(6,6)) );

    B_hill = (B_voigt + B_reuss)/2;
    
    G_hill = (G_voigt + G_reuss)/2;
    
    young_poly = 9*B_hill*G_hill/(3*B_hill + G_hill);    
    

end

function poisson_poly = local_poisson_poly_from_stiffness_matrix(obj)
    
    % References:
    % R Hill 1952 Proc. Phys. Soc. A 65 349
    % J. Am. Ceram. Soc., 96 [12] 3891–3900 (2013); DOI: 10.1111/jace.12618
    
    stiffness_matrix = obj.stiffness_matrix;    
    compliance_matrix = inv(obj.stiffness_matrix);
    
    B_voigt = 1/9*(stiffness_matrix(1,1) + stiffness_matrix(2,2) + stiffness_matrix(3,3)) + ...
              2/9*(stiffness_matrix(1,2) + stiffness_matrix(1,3) + stiffness_matrix(2,3));
          
    G_voigt = 1/15*(stiffness_matrix(1,1) + stiffness_matrix(2,2) + stiffness_matrix(3,3) - stiffness_matrix(1,2) - stiffness_matrix(1,3) - stiffness_matrix(2,3)) + ...
              1/5 *(stiffness_matrix(4,4) + stiffness_matrix(5,5) + stiffness_matrix(6,6));
          
    compressibility_reuss =    compliance_matrix(1,1) + compliance_matrix(2,2) + compliance_matrix(3,3) + ...
                            2*(compliance_matrix(1,2) + compliance_matrix(1,3) + compliance_matrix(2,3));
          
    B_reuss = 1/compressibility_reuss;
    
    G_reuss = 1/( 4*(compliance_matrix(1,1) + compliance_matrix(2,2) + compliance_matrix(3,3)) + ...
                 -4*(compliance_matrix(1,2) + compliance_matrix(1,3) + compliance_matrix(2,3)) + ...
                  3*(compliance_matrix(4,4) + compliance_matrix(5,5) + compliance_matrix(6,6)) );

    B_hill = (B_voigt + B_reuss)/2;
    
    G_hill = (G_voigt + G_reuss)/2;
    
    poisson_poly = (3*B_hill - 2*G_hill)/(2*(3*B_hill + G_hill));    
    

end

function B_poly = local_bulk_poly_from_stiffness_matrix(obj)

    % References:
    % R Hill 1952 Proc. Phys. Soc. A 65 349
    % J. Am. Ceram. Soc., 96 [12] 3891–3900 (2013); DOI: 10.1111/jace.12618

    stiffness_matrix = obj.stiffness_matrix;    
    compliance_matrix = inv(obj.stiffness_matrix);
    
    B_voigt = 1/9*(stiffness_matrix(1,1) + stiffness_matrix(2,2) + stiffness_matrix(3,3)) + ...
              2/9*(stiffness_matrix(1,2) + stiffness_matrix(1,3) + stiffness_matrix(2,3));
          
    compressibility_reuss =    compliance_matrix(1,1) + compliance_matrix(2,2) + compliance_matrix(3,3) + ...
                            2*(compliance_matrix(1,2) + compliance_matrix(1,3) + compliance_matrix(2,3));
          
    B_reuss = 1/compressibility_reuss;
    
    B_hill = (B_voigt + B_reuss)/2;
    
    B_poly = B_hill;

end

function B_sc = local_bulk_single_crystal_from_stiffness_matrix(obj)

    compliance_matrix = inv(obj.stiffness_matrix);
    S = obj.compliance_voigt_to_tensor(compliance_matrix);
    
    compressibility = S(1,1,1,1) + S(1,1,2,2) + S(1,1,3,3) + ...
                      S(2,2,1,1) + S(2,2,2,2) + S(2,2,3,3) + ...
                      S(3,3,1,1) + S(3,3,2,2) + S(3,3,3,3);
    
    compressibility = compliance_matrix(1,1) + compliance_matrix(1,2) + compliance_matrix(1,3) + ...
                      compliance_matrix(2,1) + compliance_matrix(2,2) + compliance_matrix(2,3) + ...
                      compliance_matrix(3,1) + compliance_matrix(3,2) + compliance_matrix(3,3);
                  
    B_sc = 1/compressibility;

end

function B_amorph = local_bulk_amorph_from_young_amorph_and_poisson_amorph(obj)

    young = obj.young_amorph;

    poisson = obj.poisson_amorph;

    B_amorph = young/(3-6*poisson);

end

function B_poly = local_bulk_poly_from_young_poly_and_poisson_poly(obj)

    young = obj.young_poly;

    poisson = obj.poisson_poly;

    B_poly = young/(3-6*poisson);

end

function alpha_poly = local_alpha_poly_from_alpha_matrix_and_stiffness_matrix(obj)

    % Reference:
    % DeWit, Journal OF Mechanics of Materials and Structures, Vol. 3, No. 2, 2008
    % https://msp.org/jomms/2008/3-2/jomms-v3-n2-p01-s.pdf
    
    stiffness_matrix = obj.stiffness_matrix;    
    
    alpha_matrix = obj.thermal_linear_expansion_matrix;
    
    B_voigt = 1/9*(stiffness_matrix(1,1) + stiffness_matrix(2,2) + stiffness_matrix(3,3)) + ...
              2/9*(stiffness_matrix(1,2) + stiffness_matrix(1,3) + stiffness_matrix(2,3));
          
    % Note that the effective thermal expansion coefficient of a polycrystal is coupled to the elastic constants in the Voigt model
    % cf. Eq. (28) in Reference 
    % (note that in this code stiffness_matrix uses Voigt reduced indexes, and alpha_matrix uses full tensor cartesian indexes)
    
    alpha_voigt = 1/(9*B_voigt) * (   (stiffness_matrix(1,1) + stiffness_matrix(1,2) + stiffness_matrix(1,3))*alpha_matrix(1,1) + ...
                                      (stiffness_matrix(1,2) + stiffness_matrix(2,2) + stiffness_matrix(2,3))*alpha_matrix(2,2) + ...
                                      (stiffness_matrix(1,3) + stiffness_matrix(2,3) + stiffness_matrix(3,3))*alpha_matrix(3,3) + ...
                                    2*(stiffness_matrix(1,4) + stiffness_matrix(2,4) + stiffness_matrix(3,4))*alpha_matrix(1,2) + ...
                                    2*(stiffness_matrix(1,5) + stiffness_matrix(2,5) + stiffness_matrix(3,5))*alpha_matrix(1,3) + ...
                                    2*(stiffness_matrix(1,6) + stiffness_matrix(2,6) + stiffness_matrix(3,6))*alpha_matrix(2,3) );
    
    % By contrast, in the Reuss model the overall thermal properties are independent of the elastic properties
    % cf. Eq. (31) in Reference 
                                
    alpha_reuss = 1/3*( alpha_matrix(1,1) + alpha_matrix(2,2) + alpha_matrix(3,3) );
    
    % In addition to the reference, we average between Voigt and Reuss
    
    alpha_poly = ( alpha_voigt + alpha_reuss )/2;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Conventions used for Rotations %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% General Definitions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Some general notes on rotations:
    
    % Assume all rotations about an axis are right-handed if the angle of rotation is positive.
    % Assume all matrices are meant to pre-multiply a colomn vector.
    % Assume all refernce systems are Cartesian (= orthonormal), therefore
    % there is no need to distinguish between contravariant and covariant physical vectors.

    
    % ACTIVE VS PASSIVE
    
    % Active transformation:
    % vectors (= physical objects) are rotated by a transformation R
    % In coordinates:
    % [v'_1 v'_2 v'_3]' = R * [v_1 v_2 v_3]'
    
    % Passive transformation:
    % basis vectors are rotated by the INVERSE transformation R^-1
    % In coordinates of the OLD system:
    % [e1'_1 e1'_2 e1'_3]' = inv(R) * [e1_1 e1_2 e1_3]' = inv(R) * [1 0 0]'
    % [e2'_1 e2'_2 e2'_3]' = inv(R) * [e2_1 e2_2 e2_3]' = inv(R) * [0 1 0]'
    % [e3'_1 e3'_2 e3'_3]' = inv(R) * [e3_1 e3_2 e3_3]' = inv(R) * [0 0 1]'
    
    % We will always call: 
    % R the active rotation, and 
    % r the passive rotation: r = R^-1 = inv(R).    
    
    % In BOTH active transformation by R and passive transformation by r = R^-1:
    % for any arbitrary vector,
    % described by coordinates "coord_vector":
    % new_coord_vector = R * coord_vector
    
    % Finding the rotation matrix from the active transformation picture:
    % Consider the vector [1 0 0]' (as an object to be rotated, NOT as a basis vector)
    % and find that it gets rotated into [v1'_1 v1'_2 v1'_3]'
    % A matrix applied to [1 0 0]' selects the first column, so v1'_1 v1'_2 v1'_3 is 
    % the first column of the matrix
    % Consider the vectors [0 1 0]' and [0 0 1]' and find the second and
    % third column of the matrix.
    % Therefore:
    % R = [ [ v1'_1   v2'_1   v3'_1 ]; ...
    %       [ v1'_2   v2'_2   v3'_2 ]; ...
    %       [ v1'_3   v2'_3   v3'_3 ]  ]
    
    % Finding the rotation matrix from the passive transformation picture:
    % Consider the vector e1 = [1 0 0]' (as a basis vector), rotate THE BASIS,
    % and find that it gets rotated into [e1'_1 e1'_2 e1'_3]' in the OLD coordinate system.
    % Do the same for e2 = [0 1 0]' and e3 = [0 0 1]'
    % Therefore:
    % R_inv = [ [ e1'_1   e2'_1   e3'_1 ]; ...
    %           [ e1'_2   e2'_2   e3'_2 ]; ...
    %           [ e1'_3   e2'_3   e3'_3 ]  ]   
    % R = inv(R_inv);
    
    % Finding the rotation matrix from the active transformation picture (alternative method):
    % We want the coordinates of our object (a vector) to become [0 0 1]' after the transformation.
    % The transformation of interest is therefore:
    % new_coord_vector = R * coord_vector  -->  in this case:
    % [0 0 1]' = R * [coord_vector(1)  coord_vector(2)  coord_vector(3)]'
    % (where coord_vector has been normalised)
    % In other words, keeping the reference frame fixed, imagine rotating the vector so that
    % it becomes "vertical" (aligned along the axis e3).
    % Select two direction vectors v1 and v2, that after R will have coordinates
    % [1 0 0]' and [0 1 0]', respectively.
    % These vectors are arbitrary, as long as they form an orthonormal
    % right-handed set, with coord_vector as v3.
    % To find the matrix associated with R: we have the opposite situation
    % to the one outlined above (in Finding the rotation matrix from the active transformation picture). 
    % There, the old coordinates were [1 0 0]', [0 1 0]' and [0 0 1]'.
    % Here, the new coordinates are [1 0 0]', [0 1 0]' and [0 0 1]'.
    % Therefore we invert the matrix given by the v1, v2 and v3 coordinates as columns.

    % Rotation of higher-order tensors:
    
    % 2nd rank tensor, active rotation:
    % T_new = R * T_old * R^(-1)
    % where R^(-1) = R translated = r
    
    % Nth rank tensor, active rotation:
    % T_new_ijkl.. = sum_(mnop...) R_im R_jn R_ko R_lp ... T_old_mnop...
    % The transformation of a vector and of a 2nd rank tensor are, of course,
    % special cases of this general rule.

    
    % EULER ROTATIONS    
    
    % Specifically, we use rotations about Tait–Bryan angles.
    
    % Note that sign of the rotation angles is defined, for each elemental rotation:
    
    % For an ACTIVE rotation, as positive when the BODY rotates about the
    % rotation axis in a right-handed way;
    
    % For a PASSIVE rotation, as positive when the AXES rotate about the
    % rotation axis in a right-handed way.
    
    
    % EULER ROTATIONS: INTRINSIC VS EXTRINSIC
    
    % In a composition of successive rotations, such as is the case for Euler rotations,
    % the three elemental rotations may be: 
    
    % INTRINSIC rotations = the elemental rotations occur about the axes of a coordinate system that is attached to the moving body. 
    % Therefore, the axes about which the body rotates change their orientation after each elemental rotation. 
    
    % EXTRINSIC rotations = the elemental rotations occur about the axes of a fixed coordinate system 
    % (the original coordinate system, which is assumed to remain motionless). 
    
    
    % EULER ROTATIONS: MATRIX REPRESENTATION
    
    % The matrix composition:
    
    % M_z(phi_z) * M_y(phi_y) * M_x(phi_x)
    
    % describes, equivalently:
    
    % 1) A sequence of 3 extrinsic rorations, in order:
    %    1st abiut x by phi_x
    %    2nd about y by phi_y
    %    3rd about z by phi_z

    % or
    
    % 2) A sequence of 3 intrinsic rorations, in order:
    %    1st abiut z   by phi_z
    %    2nd about y'  by phi_y
    %    3rd about x'' by phi_x
    
    % The elemental matrices are defined as:
    
    % M_x = [[1   0             0          ];
    %        [0   +cos(phi_x)   -sin(phi_x)];
    %        [0   +sin(phi_x)   +cos(phi_x)]]
    
    % M_y = [[+cos(phi_y)   0   +sin(phi_y)];
    %        [0             1   0          ];
    %        [-sin(phi_y)   0   +cos(phi_y)]]    
    
    % M_z = [[+cos(phi_z)   -sin(phi_z)   0];       
    %        [+sin(phi_z)   +cos(phi_z)   0];
    %        [0             0             1]]    
    
    
    % In both 1) and 2), the rotations can be considered as
    % Active rotations: in this case, the angles are defined as
    % positive for right-handed rotation of the BODY with respect to the rotation axis.
    % Let us call such angles THETA_z, THETA_y and THETA_x.
    % R = M_z(THETA_z) * M_y(THETA_y) * M_x(THETA_x)
    % or 
    % Passive rotations: in this case, the angles are defined as
    % positive for right-handed rotation of the AXES with respect to the rotation axis.
    % Let us call such angles theta_z, theta_y and theta_x.   
    % r = M_z(theta_z) * M_y(theta_y) * M_x(theta_x) =
    %   = M_z(-THETA_z) * M_y(-THETA_y) * M_x(-THETA_x) = inv(R)
    
    
    % This matrix corresponds to the matrix Z_1 Y_2 X_3 found at
    % https://en.wikipedia.org/wiki/Euler_angles
    % and to the matrix given by MatLab in response to the command:
    % eul2rotm([phi_z phi_y phi_x],'ZYX')
    
    
%%%%%%%%%%%%%%% Conventions used by MatLab %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % FUNCTIONS EUL2ROTM AND ROTM2EUL
    
    % eul2rotm([phi_z phi_y phi_x],'ZYX')
    
    % extrinsic in order 1st) x by phi_x 2nd) y by phi_y 3rd) z by phi_z
    % OR 
    % intrinsic in order 1st) z by phi_z 2nd) y' by phi_y 3rd) x'' by phi_x
    
    % Can be used to find active R or passive r, depending on the definition used
    % for the positive sign of the angles.
    
    
%%%%%%%%%%%%%%% Conventions used by Ansys %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % In the MAPDL input file, the values for the three angles are referred to
    % (as can be checked from the GUI)

    % A sequence of 3 intrinsic rorations, in order:
    %    1st abiut z   by phi_z
    %    2nd about x'  by phi_x
    %    3rd about y'' by phi_y
    % Active (if we consider the material reference system as an object to be rotated)

    
%%%%%%%%%%%%%%% Conventions used in this software %%%%%%%%%%%%%%%%%%%%%%%%%

    % Assume all rotations about an axis are right-handed if the angle of rotation is positive.
    % Assume all matrices are meant to pre-multiply a colomn vector.
    % Assume all refernce systems are Cartesian (= orthonormal), therefore
    % there is no need to distinguish between contravariant and covariant physical vectors.
    
    % We always call: 
    % R the active rotation, and 
    % r the passive rotation: r = R^-1 = inv(R).
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




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









