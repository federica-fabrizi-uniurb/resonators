function planeSketch1(p)
{

//Plane
p.Plane  = agb.GetActivePlane();
p.Origin = p.Plane.GetOrigin();
p.XAxis  = p.Plane.GetXAxis();
p.YAxis  = p.Plane.GetYAxis();

//Sketch
p.Sk1 = p.Plane.NewSketch();
p.Sk1.Name = "Sketch1";

//Edges
with (p.Sk1)
{
  p.Cr7 = Circle(0.00000000, 0.00000000, 38.10000000);
  // centre coords in plane + radius (e.g. diameter = 38.1*2 = 76.2)
}

//Dimensions and/or constraints
with (p.Plane)
{
  //Dimensions
  var dim;
  dim = DiameterDim(p.Cr7, 21.24990342, -12.83848360, 0);
  // the 2 coords are the approximate location of the centre of the text and do not matter
  // the 0 is the appropriate type for circle (vs ellipse)
  if(dim) dim.Name = "D1";

  //Constraints
  CoincidentCon(p.Cr7.Center, 0.00000000, 0.00000000, 
                p.Origin, 0.00000000, 0.00000000);
}

p.Plane.EvalDimCons(); //Final evaluate of all dimensions and constraints in plane

return p;
} //End Plane JScript function: planeSketch1


//////////////////////////////////////////////////////////////////////////////

var Yes = agc.Yes;
var No = agc.No;

agb.AutoConstraintGlobal(agc.Yes); 

var plxy = agb.GetXYPlane();

agb.SetActivePlane(plxy);

//Call Plane JScript function
var ps1 = planeSketch1 (new Object());
agb.Regen(); //To insure model validity

//Create Extrude of Sketch1, 0.2 units in the +Z direction
var ext1 = agb.Extrude(agc.Add, ps1.Sk1, agc.DirNormal, agc.ExtentFixed, 0.2,
agc.ExtentFixed, 0.0, agc.No, 0.0, 0.0);
agb.Regen(); //To insure model validity



var mypart = agb.FormNewPartFromAllBodies();
mypart.Name = "mypart";

