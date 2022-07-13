mymodel = ExtAPI.DataModel.Project.Model 
geo = mymodel.Geometry 
part1 = geo.Children[0]   

subst = part1.Children[0]

coat1 = part1.Children[1]
coat1.Thickness = Quantity("0.1 [m]")

coat2 = part1.Children[2]
coat2.Thickness = Quantity("0.1 [m]")

# Insert--here--substrate--coordinate--system--for--anisotropic--material

# Insert--here--coating--coordinate--systems--for--anisotropic--material


# Mesh (you must scope to the body=substrate not the face; it will automatically find the starting face)

mesh = mymodel.Mesh
mesh.ElementOrder = ElementOrder.Quadratic	
# The mixed order mesh elements warning occurs when you have solid bodies and surface bodies in a single multibody part.
# With the default settings of 'Program Controlled' for Element Order under meshing, 
# it creates higher order elements for solid bodies and lower order elements for surface bodies. 
# As both the bodies share topology, mixed order elements are created.
# To overcome this issue, it is good to set Element Order to either Quadratic or Linear for all bodies. 

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

