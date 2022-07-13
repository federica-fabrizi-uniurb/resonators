//DesignModeler JScript, version: ANSYS DesignModeler 2020 R2 (May 29 2020, 11:21:40; 20,2020,150,1) SV4
//Created via: "Write Script: Sketch(es) of Active Plane"
// Written to: C:\Users\piergiovanni\Desktop\Test_Ansys\data analysis\curvature\spherical_cap_template.js
//         On: 04/12/21, 15:37:53
//Using:
//  agb ... pointer to batch interface


//Note:
// You may be able to re-use below JScript function via cut-and-paste;
// however, you may have to re-name the function identifier.
//

function planeSketchesOnly(p)
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
//{
//  p.Cr12 = ArcCtrEdge(
//              0.00000000, 0.00000000,
//              0.00000000, 17.56586627,
//              -12.79050578, 12.04004235);
//  p.Cr13 = ArcCtrEdge(
//              0.00000000, 0.00000000,
//              0.00000000, 16.78331299,
//              -12.26844432, 11.45272321);
//  p.Ln14 = Line(-12.26844432, 11.45272321, -12.79050578, 12.04004235);
//  p.Ln15 = Line(0.00000000, 16.78331299, 0.00000000, 17.56586627);
//}

//Dimensions and/or constraints
with (p.Plane)
{
    //Constraints
//  VerticalCon(p.Ln15);
//  CoincidentCon(p.Cr12.Center, 0.00000000, 0.00000000, 
//                p.Origin, 0.00000000, 0.00000000);
//  CoincidentCon(p.Cr12.Base, 0.00000000, 17.56586627, 
//                p.YAxis, 0.00000000, 17.32591460);
//  CoincidentCon(p.Cr13.Center, 0.00000000, 0.00000000, 
//                p.Cr12.Center, 0.00000000, 0.00000000);
//  CoincidentCon(p.Cr13.Base, 0.00000000, 16.78331299, 
//                p.YAxis, 0.00000000, 16.80385315);
//  CoincidentCon(p.Ln14.Base, -12.26844432, 11.45272321, 
//                p.Cr13.End, -12.26844432, 11.45272321);
//  CoincidentCon(p.Ln14.End, -12.79050578, 12.04004235, 
//                p.Cr12.End, -12.79050578, 12.04004235);
//  CoincidentCon(p.Ln15.Base, 0.00000000, 16.78331299, 
//                p.Cr13.Base, 0.00000000, 16.78331299);
//  CoincidentCon(p.Ln15.End, 0.00000000, 17.56586627, 
//                p.Cr12.Base, 0.00000000, 17.56586627);
//}

p.Plane.EvalDimCons(); //Final evaluate of all dimensions and constraints in plane

return p;
} //End Plane JScript function: planeSketchesOnly

//Call Plane JScript function
var ps1 = planeSketchesOnly (new Object());
agb.Regen();

var Rev1 = agb.Revolve(agc.Add, ps1.Sk1, ps1.YAxis, agc.DirNormal,360.0, 0.0, agc.No, 0.0, 0.0);

//Finish
agb.Regen(); //To ensure model validity
//End DM JScript
