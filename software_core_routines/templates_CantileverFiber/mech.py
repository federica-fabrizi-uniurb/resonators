mymodel = ExtAPI.DataModel.Project.Model 
geo = mymodel.Geometry 
part1 = geo.Children[0]   

fiber = part1.Children[0]

# Create a named selection ("NS1") that selects all faces at the zero location of the XY and XZ planes.
# From ANSYS_Scripting_in_Mechanical_Guide.pdf

NS1 = DataModel.Project.Model.AddNamedSelection()
NS1.ScopingMethod = GeometryDefineByType.Worksheet
GenerationCriteria = NS1.GenerationCriteria
Criterion1 = Ansys.ACT.Automation.Mechanical.NamedSelectionCriterion()
Criterion1.Action = SelectionActionType.Add
Criterion1.EntityType = SelectionType.GeoFace
Criterion1.Criterion = SelectionCriterionType.LocationY
Criterion1.Operator = SelectionOperatorType.Equal
Criterion1.Value = Quantity("0 [m]")
GenerationCriteria.Add(Criterion1)
Criterion2 = Ansys.ACT.Automation.Mechanical.NamedSelectionCriterion()
Criterion2.Action = SelectionActionType.Add
Criterion2.EntityType = SelectionType.GeoFace
Criterion2.Criterion = SelectionCriterionType.LocationZ
Criterion2.Operator = SelectionOperatorType.Equal
Criterion2.Value = Quantity("0 [m]")
GenerationCriteria.Add(Criterion2)
### The following three lines (imposing zero tolerance) are necessary for the case of a fibre, else the whole body gets accepted
NS1.ToleranceType = ToleranceType.Manual 
NS1.ZeroTolerance = 0
NS1.RelativeTolerance = 0
###
NS1.Generate()

# Make it active

NS1.Activate()

# Add fixed support to the active object / named selection (try both ways)

analysis1 = Model.Analyses[0]
support = analysis1.AddFixedSupport()
support.Location = NS1

# Insert--here--substrate--coordinate--system--for--anisotropic--material


# Mesh (you must scope to the body not the face; it will automatically find the starting face)

mesh = mymodel.Mesh
mesh_method = mesh.AddAutomaticMethod()

my_selection = ExtAPI.SelectionManager.CreateSelectionInfo(SelectionTypeEnum.GeometryEntities)
my_selection.Ids = [6] 
mesh_method.Location = my_selection
mesh_method.Method = MethodType.Sweep
mesh_method.SweepNumberDivisions = 3

mesh.ElementSize = Quantity('0 [m]')
mesh.GenerateMesh()

an = Model.Analyses[0]
an.WriteInputFile('insert--here--path--to--inputfile.dat')


