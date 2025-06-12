// Gmsh project created on Wed Jun 11 14:35:33 2025
SetFactory("OpenCASCADE");

//Import Rail_UIC60
v() = ShapeFromFile("Rail_UIC60.stp");
CreateGeometry;
//Printf(v);
//Curve Loop(1) = v;
//Plane Surface ( 1 ) = {1};




//Physical Curve(1) -= {17, 16, 5, 4, 6, 7, 3, 8, 2, 1, 9, 36, 10, 11, 12, 13, 35, 14, 34, 15, 33, 32, 27, 31, 28, 29, 30, 26, 18, 19, 25, 20, 21, 22, 24, 23};
