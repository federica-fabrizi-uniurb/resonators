mymodel = ExtAPI.DataModel.Project.Model 
geo = mymodel.Geometry 
part1 = geo.Children[0]   

subst = part1.Children[0]

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

