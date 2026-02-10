function [out,out1]=d_cntc(varargin)

% d_cntc: transcription of CONTACT demonstration scripts
%
% <a href="matlab:d_rail('nmap')">List configurations</a>

%  Gaetan Guillet, Etienne Balmes
% Copyright SDTools & SNCF I&R, 2019-2024
% CONTACT wrapper functions
%  Copyright 2008-2023 by Vtech CMCC.
%  Licensed under Apache License v2.0.  See the file "LICENSE_CONTACT.txt" for more information.


%#ok<*ASGLU,*NASGU,*NOSEM,*NBRAK>
obj=[];evt=[];
if nargin==0; help d_cntc
    return; 
else; [CAM,Cam]=comstr(varargin{1},1);carg=2;
end

if comstr(Cam,'script');[CAM,Cam]=comstr(CAM,7);

if comstr(Cam,'v1')
%% #scriptv1

sdtweb t_gae24  CNTv1
'xxx move here'
elseif comstr(Cam,'profileview')
%% ProfileView

LI=cntc.call; if ~isfield(LI,'prw');d_cntc('LoadWheelFlat');end
colM=d_cntc('nmap.Map:Color');

%Show wheel on rail with interpolated surfaces
RP=struct('list',{{
   struct('from','prr','x',-100:5:100,'is',1:50:801,'bas','r','prop','WireICol')
   struct('from','prw','t',-4:2:4,'is',1:50:300,'bas','w','prop','WireICol')
   struct('from','prw','t',-180:1:180, ...
   'is',sdtm.indNearest(LI.prw.ysurf(1,:),0),'bas','w','prop','LineK2')
   struct('from','prw','t',-4:2:4,'is',1:50:300,'bas','wc','prop','WireICol')
   struct('from','prr','x',-500:5:500,'is',sdtm.indNearest(LI.prr.ProfileS,0),'bas','r','prop','LineK2')
   struct('from','bas','Orig',[])}},'j1',1,'bas','tr','DefLen',40,'gf',11);

   cntc.plot('list',RP);% 50 is the timestep

RP=struct('list',{{
    struct('from','prr','x',-100:5:100,'is',1:50:801,'bas','r','prop','RailPro')
    struct('from','prw','t',-10:2:10,'is',1:50:300,'bas','w','prop','WheelPro')
    struct('from','bas','bas',{{'w','r','cp'}})}},'j1',1,'bas','tr','DefLen',10,'gf',16);
 cntc.plot('list',RP);% sdtweb cntc plot.both


else;error('Script%s',CAM)
end
sdtweb('_link','sdtweb(''t_exp19'',''ref'')');

elseif comstr(Cam,'solve');[CAM,Cam]=comstr(CAM,6);


elseif comstr(Cam,'load');[CAM,Cam]=comstr(CAM,5);
%% #Load 

if comstr(Cam,'wheelflat')
%% #LoadWheelFlat : result of WheelFlat example 
  f1=sdtu.f.safe('@cntc.html/CNTCWheelflat.mat');
  if exist(f1,'file')
   cntc.load(f1);
  elseif exist('t_gae24','file'); t_gae24('cntcWheelFlat')
  end
  evalin('caller','LI=cntc.call;');
else; error('Load%s',CAM)
end

elseif comstr(Cam,'nmap');
%% #nmap list of named experiments used for demos and non-regression tests ---

[key,nmap,uo,carg]=sdtm.stdNmapArgs(varargin,CAM,carg);RO.inKeys=nmap.keys;
projM=nmap;

%% #Cin: #Map:Cin parameter formating and tooltip  ----2

%cinM=nmap('Map:Cin');
%ua=struct('ColumnName',{{'name','value','tooltip'}}, ...
%    'table',{{'root','cur','name root substring'}}); 
%cinM.add=ua;  % Create a default CinCell structures to be used for vhandle.uo 
% projM.append(d_rail('nmap'));projM('Map:Cin')

cinM=projM('Map:Cin');% Create a default CinCell structures to be used for vhandle.uo 
cinM.add=['top(-700 300#%g#"contact refine start end (rail direction)")' ...
          'topCoarse(0#31#"do not refine coarse")' ...
          'gKC( #%g#"sleeper gradient properties")' ...
          'rail( #%s#"Rail Urn, eg U60")' ...
          'pad( #%s#"Pad Urn, eg PadFuSn{io4}")' ...
          'sleeper( #%s#"Sleeper Urn, eg mass{4}")' ...
          'sub( #%s#"Sub-structure Urn, eg spring{k,c}")' ...
          'sw(.6#%s#"slice width")' ...
          'quad(#31#"use quadratic elements")' ...
          'spos(#31#"position of sleeper in slice 1=half, 0=centered ")' ...
          'unit(SI#%s#"unit system SI/TM/MM, ... ")' ...
          'Tgrad(gc#%s#"type of gradient at the end xxx/ track ends")' ...
          'InterMPC(1#31#" use MCP for non-conform refinement")' ...
     ... % #Cin.contact   -3      
          'Kc(-23e9#%ug#"contact stiffness (>0 lin,<0 sqrt)")' ...
          'delta(20#%ug#"contact lowpass")' ...
          'Fl(20#%ug#"contact saturation")' ...
     ... %
          'Lc(50#%g#"characteristic length for rail/wheel mesh")' ...
          'CLc(2#%g#"characteristic length for rail/wheel contact")' ...
          'CDiam(12#%g#"width of potential contact")' ...
          'CThick(10#%g#"thickness of contact refined zone")' ...
     ... %
          'Wref(coarse1#%s#"wheel refinement strategy")' ...
     ... %
          'dt(1e-5#%g#"time step")' ...
          'Vel(20#%ug#"train velocity m/s")' ...
          'profile(9#31#"strategy for profiling")' ...
          'FirstDown(-.01#%g#" vertical wheel step to gain contact")' ...
          ];  


cinM.add={
'gr:Ctc','Contact parameters', ...
      {'Kc','delta'}
'gr:Wheel','Wheel refinement parameters', ...
        {'Wref'}
'gr:s1','Slice definition', ...
    {'rail','pad','sleeper','sub','unit'}
'gr:Track','Track definition parameters', ...
      {'top','topCoarse','sw','quad','Lc','CLc','CWide','CDeep','unit','TGrad'}
'gr:SimuTime','Time simulation parameters', ...
      {'dt','Vel','profile','FirstDown'}
'gr:All','Experiment parameters', ...
      {'gr:s1','gr:Track','gr:Wheel','gr:Ctc','gr:SimuTime'}
      };

% vhandle.uo('',C3.info,rail19('nmap.Map:Cin'))

%% #Map:Color -2
colM=vhandle.nmap({'Wheel','#FF6D00';'b','#FF9E00';'c','#00B4D8'; ...
    'Rail','#0077B6';'e','#023E8A'});
nmap('Map:Color')=colM;

%% #Map:OProp : wrong see rail u -2
nmap('Map:Oprop')=railu.getNmap('Map:OProp');
nmap('Map:Annot')=railu.getNmap('Map:Annot');

%% #GlobalFlags -2
%CNTC.un_cntc=1934;
nmap('Global_flags')=vhandle.uo([],{ ...
  'if_units',1934,'Unit system : Contact, SI, Simpack'
  'ic_config',0,'Configuration of the problem; 0 for left side, 1 for right side' 
  'ic_tang',3,'Tangential problem to be solved, T=3: steady state rolling --> G=0 GauSei'
  'ic_pvtime',2,['Relation of the current to the previous case or' ...
  ' previous time instance; P=2: no previous time'] 
  'ic_discns',2,'Type of contact area, D=2 planar contact surface'
  'if_wrtinp',1,' 0: no .inp-file needed'
  'ic_matfil',1,'subsurface stress input; A=0: no .mat-file needed'
  'ic_output',3,'subsurface stress output; O=1: min. output to .out-file'
  'ic_flow',2,['governs the extent of the flow trace to the screen and the ' ...
  'output-file; little progress output']
  'idebug',1,'1: just a bit of information from the library'
  'imodul', 1,'w/r contact'
  'ire', 1, 'experiment called iwhe in the example'
  'iwhe',1, '1/2 Left/right wheel'});

nmap('DataExchange')={'s_ws';'dywhl';'z_ws';'drollw';'dyaww';'pitch_ws';% pos'
       'vs';'vy';'vz';'vroll';'vyaw';'vpitch';        %vel     
       'dxwhl';'dywhl';'dzwhl';'drollw';'dyaww';'dpitchw'; %flex'
       'vxwhl';'vywhl';'vzwhl';'vrollw';'vyaww';'vpitchw';  
       'dyrail';'dzrail';'drollr';'vyrail';'vzrail';'vrollr'; %dev'
       }; 

%% #CNTCModel -2
% #Mod_CsteProf -2
nmap('Mod_CsteProf')={'Solver{GauSei default,maxgs 999,maxin 100, maxnr 30, maxout 1,eps 1e-5}'
    'Mat{Mater1 0, nu1 0.28, nu2 0.28, g1 82000,g2 82000}'
    'Friction{FrcLaw Coul, fstat 0.3, fkin 0.3}'
    'PotCntc{PosGrid WR, dx 0.2, ds 0.2, a_sep 90deg, d_sep 8.0, d_comb 4.0}'
    'Rolling{StepSize WRCn, dqrel 1}'
    ['Track{Design NewBoth, gaught 14, gaugsq 0, gaugwd 1435, cant 0.05' ...
   ', nomrad 0, dyrail 0, dzrail 0, drollr 0, vyrail 0, vzrail 0, vrollr 0']
    'setProfile{fname "MBench_UIC60_v3.prr",iswheel 0,mirrory 0, sclfac 1., smooth 0.}'
    'setProfile{fname "MBench_S1002_v3.prw",iswheel 1,mirrory 0, sclfac 1., smooth 0.}'
    'wheelsetDim{Ewheel NewAll, fbdist 1360, fbpos -70, nomrad 460}'
    };

% #Mod_WheelflatPBound penetration boundary condition -2
nmap('Mod_WheelflatPB')={'Solver{GauSei default,maxgs 999,maxin 100, maxnr 30, maxout 1,eps 1e-5}'
   'Mat{Mater1 0, nu1 0.28, nu2 0.28, g1 82000,g2 82000}'
   'Friction{FrcLaw Coul, fstat 0.3, fkin 0.3}'
   'PotCntc{PosGrid WR, dx 0.4, ds 0.4, a_sep 90deg, d_sep 8.0, d_comb 4.0}'
   'Rolling{StepSize WRCn, dqrel 1}'
    ['Track{Design NewBoth, gaught 14, gaugsq 0, gaugwd 1435, cant 0.05' ...
   ', nomrad 0, dyrail 0, dzrail 0, drollr 0, vyrail 0, vzrail 0, vrollr 0}']
    % the Gauge calculation require the track definition  
   'setProfile{fname "r300_wide.prr",iswheel 0,mirrory 0, mirrorz -1, sclfac 1, smooth 0}'
   'setProfile{fname "S1002_flat.slcw",iswheel 1,mirrory 0,mirrorz -1, sclfac 1, smooth 5}'
   'wheelsetDim{Ewheel NewAll, fbdist 1360, fbpos -70, nomrad 460}'
   };
% #Mod_WheelflatFBound force boundary condition -2
nmap('Mod_WheelflatFB')={'Solver{GauSei default,maxgs 999,maxin 100, maxnr 30, maxout 1,eps 1e-5}'
   'Mat{Mater1 0, nu1 0.28, nu2 0.28, g1 82000,g2 82000}'
   'Friction{FrcLaw Coul, fstat 0.3, fkin 0.3}'
   'PotCntc{PosGrid WR, dx 0.4, ds 0.4, a_sep 90deg, d_sep 8.0, d_comb 4.0}'
   'Bound{Cond Force, fz 125000}'
   'Rolling{StepSize WRCn, dqrel 1}'
    ['Track{Design NewBoth, gaught 14, gaugsq 0, gaugwd 1435, cant 0.05' ...
   ', nomrad 0, dyrail 0, dzrail 0, drollr 0, vyrail 0, vzrail 0, vrollr 0}']
   'wheelsetDim{Ewheel NewDimProfPosVel, fbdist 1360, fbpos -70, nomrad 460}'
   'setProfile{fname "r300_wide.prr",iswheel 0,mirrory 0, mirrorz -1, sclfac 1, smooth 0}'
   'setProfile{fname "S1002_flat.slcw",iswheel 1,mirrory 0,mirrorz -1, sclfac 1, smooth 5}'};

% #Mod_VarProf -2
nmap('Mod_VarProf')={'Solver{GauSei default,maxgs 999,maxin 100, maxnr 30, maxout 1,eps 1e-5}'
   'Mat{Mater1 0, nu1 0.28, nu2 0.28, g1 82000,g2 82000}'
   'Friction{FrcLaw Coul, fstat 0.3, fkin 0.3}'
   'Bound{Cond Force, fz 125000}'
   'PotCntc{PosGrid WR, dx 0.4, ds 0.4, a_sep 90deg, d_sep 8.0, d_comb 4.0}'
   'Rolling{StepSize WRCn, dqrel 1}'
    ['Track{Design NewBoth, gaught 14, gaugsq 0, gaugwd 1435, cant 0.05' ...
   ', nomrad 0, dyrail 0, dzrail 0, drollr 0, vyrail 0, vzrail 0, vrollr 0']
   'setProfile{fname "var_rail.slcs",iswheel 0,mirrory 0, sclfac 1, smooth 0}'
   'setProfile{fname "var_wheel.slcw",iswheel 1,mirrory 0, sclfac 1, smooth 0}'
   % 'setProfile{fname "S1002_flat.slcw",iswheel 1,mirrory 0, sclfac 1, smooth 0}'
   'wheelsetDim{Ewheel NewDimProfPosVel, fbdist 1360, fbpos -70, nomrad 460}'};

%% #Trajectory -2
LI=cntc.call;

% #Traj_Sin -2
C1=struct('X',{{[],{'dywhl';'s_ws'}}},'Xlab',{{'Step','Comp'}},'Y', []);
C1.Y=zeros(100,2);
C1.Y(:,2)=linspace(0,1000,100);
C1.Y(:,1)=5*cos(pi/50*C1.Y(:,2));
C1.X{1}=(1:size(C1.Y,1))';
nmap('Traj_Sin')=C1;

% #Traj_Classic -2
C1=struct('X',{{[],{'dzwhl';'vpitch';'vs'}}},'Xlab',{{'Step','Comp'}},'Y', [-0.5 -4.34811810 2000]);
C1.X{1}=(1:size(C1.Y,1))';
nmap('Traj_Classic')=C1;
% #Traj_ODG -2
% Vs=89000mm/s ; Vpitch = 
C1=struct('X',{{[],{'dywhl';'dyaww';'drollw';'vpitch';'vs'}}},'Xlab',{{'Step','Comp'}},'Y', [ ...
   [ 0 ]; [ 0 ]; [  0.00000000]; [ -4.34811810]]');
C1.X{1}=(1:size(C1.Y,1))';
C1.Y(:,5)=2000; % set vs
nmap('Traj_ODG')=C1;

% #Traj_Benchmark -2
C1=struct('X',{{[],{'dywhl';'dyaww';'drollw';'vpitch';'vs';'z_ws'}}},'Xlab',{{'Step','Comp'}},'Y', [ ... 
   [ 0 : 0.5 : 10 ];
   [ 0 : 0.0012 : 0.024 ];
   [  0.00000000, -0.00002304, -0.00005049, -0.00008103, -0.00011280, ...
            -0.00014570, -0.00018030, -0.00021680, -0.00025570, -0.00029770, ...
            -0.00035540, -0.00047770, -0.00062720, -0.00437600, -0.00639300, ...
            -0.00764200, -0.00860600, -0.00940800, -0.01010113, -0.01071386, ...
            -0.01126431 ];
   [ -4.34811810, -4.34741340, -4.34657520, -4.34624400, -4.34591970, ...
            -4.34556270, -4.34515030, -4.34466880, -4.34409300, -4.34337150, ...
            -4.33536370, -4.33188640, -4.32937180, -4.27488340, -4.26356290, ...
            -4.25757470, -4.25348570, -4.25032450, -4.24775610, -4.24556650, ...
            -4.24363750 ]
   ]');
C1.X{1}=(1:size(C1.Y,1))';
C1.Y(1:end,5)=2000; % set vs
C1.Y(1:end,6)=0.1; %set penetration
nmap('Traj_BenchManchester')=C1;

%  #Traj_TestPen -2
C1=struct('X',{{[],{'z_ws'}}},'Xlab',{{'Step','Comp'}},'Y', [[1,1,2,2,3,3,6,6]]');
C1.X{1}=(1:size(C1.Y,1))';
nmap('Traj_TestPen')=C1;

%  #Traj_WheelflatW pen defined on the wheel -2
C1=struct('X',{{[],{'s_ws','mm';'pitch_ws','rad';'vs','mm/s';'vpitch','rad/s';'dzwhl','[mm]'}}}, ...
 'Xlab',{{'Step','Comp'}},'Y',[]);
 C1.Y(:,2)=-pi/180*linspace(0,50,100)'; % rotation angle [rad] (pdf page 36)
 % C1.Y(:,2)=-pi/180*linspace(0,20,10)'; % rotation angle [rad] (pdf page 36)
if isfield(LI,'wheelsetDim')
 C1.Y(:,1)=-C1.Y(:,2)*LI.wheelsetDim.nomrad;
end
C1.X{1}=(1:size(C1.Y,1))';
C1.Y(:,3)=2000; % set vs [mm/s]
C1.Y(:,4)=-4.08190679; % set vpitch wheel rotation speed [rad/s]
% z_ws curve found with load boundary conditions fz=125000 N
C1.Y(:,5)=-[0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,...
 0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,...
 0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,...
 0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,...
 0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,...
 0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,...
 0.4436,0.4426,0.4432,0.4460,0.4518,0.4617,0.4726,0.4853,0.4979,...
 0.5151,0.5394,0.5705,0.6089,0.6533,0.7009,0.7485,0.7924,0.8305,...
 0.8593,0.8805,0.8938,0.9005,0.9002,0.8905,0.8743,0.8466,0.8018,...
 0.7528,0.7042,0.6598,0.6215,0.5882,0.5599,0.5353,0.5156,0.4993,...
 0.4845,0.4702,0.4590,0.4515,0.4461,0.4445,0.4436,0.4432,0.4432,...
 0.4433,0.4434];
nmap('Traj_WheelFlatW')=C1;

%  #Traj_Simp_WheelFlat pen defined on the wheel -2
C1=struct('X',{{[],{'s_ws','[mm]';'dzwhl','[mm]';'pitch_ws','[rad]'}}}, ...
 'Xlab',{{'Step','Comp'}},'Y',[]);
 C1.Y(:,3)=-pi/180*linspace(0,50,100)'; % rotation angle [rad] (pdf page 36)
 % C1.Y(:,2)=-pi/180*linspace(0,20,10)'; % rotation angle [rad] (pdf page 36)
if isfield(LI,'wheelsetDim')
 C1.Y(:,1)=-C1.Y(:,1)*LI.wheelsetDim.nomrad;
end
C1.X{1}=(1:size(C1.Y,1))';
% z_ws curve found with load boundary conditions fz=125000 N
C1.Y(:,2)=-[0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,...
 0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,...
 0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,...
 0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,...
 0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,...
 0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,0.4436,...
 0.4436,0.4426,0.4432,0.4460,0.4518,0.4617,0.4726,0.4853,0.4979,...
 0.5151,0.5394,0.5705,0.6089,0.6533,0.7009,0.7485,0.7924,0.8305,...
 0.8593,0.8805,0.8938,0.9005,0.9002,0.8905,0.8743,0.8466,0.8018,...
 0.7528,0.7042,0.6598,0.6215,0.5882,0.5599,0.5353,0.5156,0.4993,...
 0.4845,0.4702,0.4590,0.4515,0.4461,0.4445,0.4436,0.4432,0.4432,...
 0.4433,0.4434];
nmap('Traj_Simp_WheelFlat')=C1;

%  #Traj_Amplified
C1=struct('X',{{[],{'s_ws','mm';'pitch_ws','rad';'vs','mm/s';'vpitch','rad/s';'dzwhl','[mm]';'dywhl','[mm]'}}}, ...
 'Xlab',{{'Step','Comp'}},'Y',[]);
 C1.Y(:,2)=-pi/180*linspace(0,50,3)'; % rotation angle [rad] (pdf page 36)
 % C1.Y(:,2)=-pi/180*linspace(0,20,10)'; % rotation angle [rad] (pdf page 36)
if isfield(LI,'wheelsetDim')
 C1.Y(:,1)=-C1.Y(:,2)*LI.wheelsetDim.nomrad;
end
C1.X{1}=(1:size(C1.Y,1))';
C1.Y(:,3)=2000; % set vs [mm/s]
C1.Y(:,4)=-4.08190679; % set vpitch wheel rotation speed [rad/s]
C1.Y(:,5)=-[1 2 4];
C1.Y(:,6)=-[1 2 4];
nmap('Traj_Amplified')=C1;

%  #Traj_rail_disp
C1=struct('X',{{[],{'s_ws','mm';'pitch_ws','rad';'vs','mm/s';'vpitch','rad/s';'dyrail','[mm]';'dzrail','[mm]';'drollr','[rad]'}}}, ...%
 'Xlab',{{'Step','Comp'}},'Y',[]);
 C1.Y(:,2)=-pi/180*linspace(0,50,3)'; % rotation angle [rad] (pdf page 36)
 % C1.Y(:,2)=-pi/180*linspace(0,20,10)'; % rotation angle [rad] (pdf page 36)
if isfield(LI,'wheelsetDim')
 C1.Y(:,1)=-C1.Y(:,2)*LI.wheelsetDim.nomrad;
end
C1.X{1}=(1:size(C1.Y,1))';
C1.Y(:,3)=2000; % set vs [mm/s]
C1.Y(:,4)=-4.08190679; % set vpitch wheel rotation speed [rad/s]
C1.Y(:,5)=[1 2 4];
C1.Y(:,6)=[1 2 4];
C1.Y(:,7)=-pi/180*[1 5 10];
nmap('Traj_rail_disp')=C1;

%  #Traj_WheelflatR pen defined on the rail -2
% C1=struct('X',{{[],{'pitch_ws';'vs';'vpitch'}}},'Xlab',{{'Step','Comp'}},'Y', ...
%  [-pi/180*[25 : 1 : 50]']); % rotation angle [rad] (pdf page 36)
C1=struct('X',{{[],{'s_ws','mm';'pitch_ws','rad';'vs','mm/s';'vpitch','rad/s';'dzrail','[mm]'}}}, ...
 'Xlab',{{'Step','Comp'}},'Y',[]);
 C1.Y(:,2)=-pi/180*linspace(0,50,100)'; % rotation angle [rad] (pdf page 36)
 % C1.Y(:,2)=-pi/180*linspace(0,20,10)'; % rotation angle [rad] (pdf page 36)
if isfield(LI,'wheelsetDim')
 C1.Y(:,1)=-C1.Y(:,2)*LI.wheelsetDim.nomrad;
end
C1.X{1}=(1:size(C1.Y,1))';
C1.Y(:,3)=2000; % set vs [mm/s]
C1.Y(:,4)=-4.08190679; % set vpitch wheel rotation speed [rad/s]
% C1.Y(:,5)=0*ones(size(C1.Y,1),1);
% C1.Y(:,5)=linspace(0.1,1,10);
% C1.Y(:,5)= 0.5*ones(size(C1.Y,1),1);
% z_ws curve found with load boundary conditions fz=125000 N
C1.Y(:,5)=[-0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-...
 0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-...
 0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-...
 0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-...
 0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-...
 0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-0.4436,-...
 0.4436,-0.4426,-0.4432,-0.4460,-0.4518,-0.4617,-0.4726,-0.4853,-0.4979,-...
 0.5151,-0.5394,-0.5705,-0.6089,-0.6533,-0.7009,-0.7485,-0.7924,-0.8305,-...
 0.8593,-0.8805,-0.8938,-0.9005,-0.9002,-0.8905,-0.8743,-0.8466,-0.8018,-...
 0.7528,-0.7042,-0.6598,-0.6215,-0.5882,-0.5599,-0.5353,-0.5156,-0.4993,-...
 0.4845,-0.4702,-0.4590,-0.4515,-0.4461,-0.4445,-0.4436,-0.4432,-0.4432,-...
 0.4433,-0.4434];
nmap('Traj_WheelFlatR')=C1;

%  #Traj_WheelflatPen pen defined on the wheel -2
% C1=struct('X',{{[],{'pitch_ws';'vs';'vpitch'}}},'Xlab',{{'Step','Comp'}},'Y', ...
%  [-pi/180*[25 : 1 : 50]']); % rotation angle [rad] (pdf page 36)
C1=struct('X',{{[],{'s_ws','mm';'pitch_ws','rad';'vs','mm/s';'vpitch','rad/s';'dzwhl','[mm]'}}}, ...
 'Xlab',{{'Step','Comp'}},'Y',[]);
 C1.Y(:,2)=-pi/180*linspace(0,6,6)'; % rotation angle [rad] (pdf page 36)
if isfield(LI,'wheelsetDim')
 C1.Y(:,1)=-C1.Y(:,2)*LI.wheelsetDim.nomrad;
end
C1.X{1}=(1:size(C1.Y,1))';
C1.Y(:,3)=2000; % set vs [mm/s]
C1.Y(:,4)=-4.08190679; % set vpitch wheel rotation speed [rad/s]
C1.Y(:,5)=[0.1:0.1:0.6]';
nmap('Traj_WheelFlatPen')=C1;

%  #Traj_BasisDemo wheelflat with positiv wheel rotation-2
% C1=struct('X',{{[],{'pitch_ws';'vs';'vpitch'}}},'Xlab',{{'Step','Comp'}},'Y', ...
%  [-pi/180*[25 : 1 : 50]']); % rotation angle [rad] (pdf page 36)
C1=struct('X',{{[],{'s_ws','mm',[];'pitch_ws','rad',[];'vs','mm/s',[];'vpitch','rad/s',[]}}}, ...
 'Xlab',{{'Step','Comp'}},'Y',[]);
 C1.Y(:,2)=pi/180*linspace(20,-50,100)'; % rotation angle [rad] (pdf page 36)
if isfield(LI,'wheelsetDim')
 C1.Y(:,1)=C1.Y(:,2)*LI.wheelsetDim.nomrad;
end
C1.X{1}=(1:size(C1.Y,1))';
C1.Y(:,3)=2000; % set vs [mm/s]
C1.Y(:,4)=-4.08190679; % set vpitch wheel rotation speed [rad/s]
nmap('Traj_BasisDemo')=C1;

%  #Traj_TimeTest -2
C1=struct('X',{{[],{'pitch_ws','rad',[];'s_ws','mm',[];'vs','mm/s',[];'vpitch','rad/s',[]}}}, ...
 'Xlab',{{'Step','Comp'}},'Y', ...
 -pi/180*[0:720]'); % rotation angle [rad] (pdf page 36)
C1.X{1}=(1:size(C1.Y,1))';
if isfield(LI,'wheelsetDim')
 C1.Y(:,2)=-C1.Y(:,1)*LI.wheelsetDim.nomrad;
end
C1.Y(:,3)=2000; % set vs [mm/s]
C1.Y(:,4)=-4.08190679; % set vpitch wheel rotation speed [rad/s]
nmap('Traj_TimeTest')=C1;

%  #Traj_TestParam -2
%  #Traj_TestParam1 -3
l1={'s_ws' ;'dywhl' ;'z_ws' ;'drollw';'dyaww';'pitch_ws';...
    'vs'   ;'vy'   ;'vz'   ;'vroll'  ;'vyaw'  ;'vpitch';...
    'dxwhl';'dywhl';'dzwhl';'drollw' ;'dyaww' ;'dpitchw'; ...
    'vxwhl';'vywhl';'vzwhl';'vrollw' ;'vyaww' ;'vpitchw'};
C1=struct('X',{{[],l1}},'Xlab',{{'Step','Comp'}}, ...
 'Y',[]); % rotation angle [rad] (pdf page 36)
C1.Y(1,:)=zeros(1,24);C1.Y(1,3)=2;
C1.Y(2,:)=[0,3,0.4,0,0,0,...
           0,0,0,0,0,0,...
           1,2,1,0,0,0,...
           0,0,0,0,0,0];
C1.X{1}=(1:size(C1.Y,1))';
nmap('Traj_TestParam1')=C1;

%  #Traj_TestParamFlex1  -3
l1={'z_ws' ;'dywhl'};
C1=struct('X',{{[],l1}},'Xlab',{{'Step','Comp'}}, ...
 'Y',[[0.4,0.4,0.4,0.4];[-10,-5,0,5]]'); % rotation angle [rad] (pdf page 36)
C1.X{1}=(1:size(C1.Y,1))';
nmap('Traj_TestParamFlex1')=C1;

%  #Traj_TestParamFlex2  -3
l1={'z_ws';'dxwhl';'dywhl';'dzwhl';'drollw';'dyaww'};
C1=struct('X',{{[],l1}},'Xlab',{{'Step','Comp'}}, ...
 'Y',[[0.4 0.4 0.4]; [0 1 1]; [0 1 1]; [0 1 1]; [0 0 0.0175]; [0 0 0.0175]]'); % rotation angle [rad] (pdf page 36)
C1.X{1}=(1:size(C1.Y,1))';
nmap('Traj_TestParamFlex2')=C1;

%  #Traj_TestParamFlex3  -3
l1={'z_ws';'drollw' ;'dyaww' ;'dpitchw'};
C1=struct('X',{{[],l1}},'Xlab',{{'Step','Comp'}}, ...
 'Y',[[0.4 0.4]; [0 0.00175]; [0 0.00175]; [0 0.00175]]'); % rotation angle [rad] (pdf page 36)
C1.X{1}=(1:size(C1.Y,1))';
nmap('Traj_TestParamFlex3')=C1;

%  #Traj_TestParam3 -3
l1={ 'z_ws';'dywhl';'dyaww';'dywhl';'dyaww' ;};
C1=struct('X',{{[],l1}},'Xlab',{{'Step','Comp'}}, ...
 'Y',[]); % rotation angle [rad] (pdf page 36)
C1.Y=zeros(8,5);
C1.Y(:,2)=[-1,0,1,5,0,0,0,0] ;
C1.Y(:,3)=[-0.0175,0,0.0175,0.035,0,0,0,0] ;
C1.Y(:,4)=[0,0,0,0,-1,0,1,5] ;
C1.Y(:,5)=[0,0,0,0,-0.0175,0,0.0175,0.035] ;
C1.Y(:,1)=0.4; %set penetration
C1.X{1}=(1:size(C1.Y,1))';
nmap('Traj_TestParam3')=C1;

%  #Traj_TestParam2  -3
l1={'s_ws' ;'dywhl' ;'z_ws' ;'drollw';'dyaww';'pitch_ws';...
    'vs'   ;'vy'   ;'vz'   ;'vroll'  ;'vyaw'  ;'vpitch';...
    'dxwhl';'dywhl';'dzwhl';'drollw' ;'dyaww' ;'dpitchw'; ...
    'vxwhl';'vywhl';'vzwhl';'vrollw' ;'vyaww' ;'vpitchw'
    'dyrail';'dzrail';'drollr';'vyrail';'vzrail';'vrollr'};
C1=struct('X',{{[],l1}},'Xlab',{{'Step','Comp'}}, ...
 'Y',[]); % rotation angle [rad] (pdf page 36)
C1.Y=zeros(1,30); C1.Y(2,14)=-10; 
C1.Y(3,14)=10; C1.Y(:,3)=2; %set penetration
C1.X{1}=(1:size(C1.Y,1))';
nmap('Traj_TestParam2')=C1;

%% #nmap.Cst : constants to be reused -3
nmap('Ctc21')=struct('delta',20,'Fl',1e6,'Kc',-2.3e9);
nmap('CMsh21')=struct('CRailEdge',.1,'Vel',20,'ToolTip','xxx');

 if nargout==1;out=sdtm.stdNmapOut(nmap,key,nargout,CAM);
 elseif nargout>1;[out,out1,out2]=sdtm.stdNmapOut(nmap,key,nargout,CAM);
 else; sdtm.stdNmapOut(nmap,key,nargout,CAM);
 end


%% #end
elseif comstr(Cam,'wd')
 %% #wd : Resolve fullfile from ProjectWd ----------------------------2
 out=comgui('cd',fullfile(base_wd,varargin{carg:end}));
elseif comstr(Cam,'pwd')
 %% #pwd : Resolve fullfile from ProjectWd ----------------------------2
 out=comgui('cd',fullfile(sdtroot('PARAM.Project.ProjectWd'),varargin{carg:end}));
 if nargout>1; out1=exist(out,'file');end 

elseif comstr(Cam,'jup')
%% #jup : jupyter notebook building -2
wd=fullfile(fileparts(which('d_rail')),'../jup');

RO=struct('BuildDir',sdtm.safeFileName('@tempdir\_jup'), ...
      'reset',1,'book',{{'rail'}}); 
RO.helpDir=sdtu.f.safe('@sdt/helpj');
RO.endCopy='jup{rail}';
RO=sdtm.jup('build',RO);



%% #Latex #Tex #Hevea -------------------------------------------------------2
elseif comstr(Cam,'svg')
  cd(fullfile(fileparts(which('rail19')),'../tex/plots'));
  lat('tex2svg',{'nl_maxwell.tex';'nl_maxwell_m.tex';'nl_GapCyl.tex'; 
      'nl_sts.tex'})
    
elseif comstr(Cam,'latex')||comstr(Cam,'tex')||comstr(Cam,'hevea')

pw0=pwd;
%cd(fullfile(HOME,'tex'));
HOME=fileparts(fileparts(which('rail19')));
wd=sdth.fileutil('firstdir',{fullfile(HOME,'tex'),fullfile(HOME,'../tex')});
cd(wd);
wd=pwd;
st=fullfile(sdtcheck('SDTRootDir'),'tex');
if exist(st,'dir')
   % keywords, - between, but not at beginning
   setenv('sdtfun',['rail19-hbm_solve-hbmui-hbm_post-nl_solve-nl_mesh' ...
      '-nl_spring-nl_inout-chandle-mkl_utils-d_hbm'  ]); 
   wd1=sdtu.f.safe('@sdt/help/nlsim');
   if comstr(Cam,'hevea'); 
     wdh=lat('hevea -module SNCF-DyRail','dyrail');
     st=sprintf('lat(''tex2svg -eq'',''%s'')',pwd);sdtweb('_link',st);
     cd(wdh);!cp index.html dyrail.html 
   else;
      lat('colortex');lat('_dyrail.tex-shortlog');
      !bibtex _nlsim
      eval(sprintf('!mv _nlsim.pdf "%s"',fullfile(wd1,'nlsim.pdf')));
      sdtweb('_link',sprintf('sdtweb(''%s'')',...
              fullfile(fileparts(wd),'help','nlsim.pdf')),'PDF DOC')
   end
else
 st=sprintf('%s nlsim',sd_pref('get','OpenFEM','LatexFcn','!pdflatex'));
 eval(st);
end
cd(pw0);
sdtweb('_link','hbm_utils(''helpputall'')');
sdtweb('_link','sd helpcursync helphbm');


elseif comstr(Cam,'cvs'); out=sdtcheck('revision'); 
elseif comstr(Cam,'@');out=eval(CAM);
elseif ~isempty(CAM); error('%s',CAM);
end
end

%% #SubFunc





