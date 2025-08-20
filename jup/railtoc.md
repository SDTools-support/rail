(railtoc)=
# Rail : interface between SDT and other components

This documents meshing scripts related to rail applications.

This is an Open Source set of MATLAB/SDT functions needed to generate meshes 


````matlab
% The following commands are used by sdtm.jup for building
toc={...
'railtoc'          1
'railmesh'          2
'mapRailComp'             3
'urnRailMesh'           3 % 
'd_rail' 2
'railjup_Install' 2
'cntc'          2
};
RO.helpDir=sdtu.f.cffile('@sdt/helpj');
RO.ImgSrc= sdtu.f.cffile('@d_rail.m/../tex/plots');
RO.BuildCmd='sdtm.jup(''Buildrail'')';
````