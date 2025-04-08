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

else;error('Script%s',CAM)
end
sdtweb('_link','sdtweb(''t_exp19'',''ref'')');

elseif comstr(Cam,'solve');[CAM,Cam]=comstr(CAM,6);


elseif comstr(Cam,'load');[CAM,Cam]=comstr(CAM,5);
%% #Load 



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

%% Global flags 
%CNTC.un_cntc=1934;
nmap('BenchManchester_flags')=vhandle.uo([],{ ...
  'if_units',1934,'Unit system : Contact, SI, Simpack'
  'ic_config',0,'Configuration of the problem; 0 for left side, 1 for right side'
  'ic_tang',3,'Tangential problem to be solved, T=3: steady state rolling --> G=0 GauSei'
  'ic_pvtime',2,['Relation of the current to the previous case or' ...
  ' previous time instance; P=2: no previous time'] 
  'ic_discns',2,'Type of contact area, D=2 planar contact surface'
  'if_wrtinp',0,' 0: no .inp-file needed'
  'ic_matfil',0,'subsurface stress input; A=0: no .mat-file needed'
  'ic_output',3,'subsurface stress output; O=1: min. output to .out-file'
  'ic_flow',2,['governs the extent of the flow trace to the screen and the ' ...
  'output-file; little progress output']
  });

%% Trajectory 
%Classic
C1=struct('X',{{[],{'y';'yaw';'roll';'vpitch';'vx'}}},'Xlab',{{'Step','Comp'}},'Y', [ ...
   [ 0 ]; [ 0 ]; [  0.00000000]; [ -4.34811810]]');
C1.X{1}=(1:size(C1.Y,1))';
C1.Y(:,5)=2000; % set vx
nmap('Classic_traj')=C1;

%  Benchmark
C1=struct('X',{{[],{'y';'yaw';'roll';'vpitch';'vx'}}},'Xlab',{{'Step','Comp'}},'Y', [ ...
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
C1.Y(:,5)=2000; % set vx
nmap('BenchManchester_Traj')=C1;

%% #nmap.Cst : constants to be reused
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
RO.helpDir=sdtu.f.cffile(fileparts(which('feplot')),'helpj');
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
   eval(sprintf('!cp %s .',fullfile(st,'etienne.bib')));
   eval(sprintf('!cp %s .',fullfile(st,'macros_tex.sty')));
   eval(sprintf('!cp %s .',fullfile(st,'macros_hevea.hva')));
   eval(sprintf('!cp %s .',fullfile(st,'macros.tex')));
   eval(sprintf('!cp %s .',fullfile(st,'book.cls')));
   % keywords, - between, but not at beginning
   setenv('sdtfun',['rail19-hbm_solve-hbmui-hbm_post-nl_solve-nl_mesh' ...
      '-nl_spring-nl_inout-chandle-mkl_utils-d_hbm'  ]); 
   wd1=fullfile(fileparts(which('feplot')),'helpnlsim');
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





