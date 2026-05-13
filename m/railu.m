classdef railu
    
% railu : name space for methods related to RAIL
%
% <a href="matlab:methods railu">methods railu</a> for methods <a href="matlab:railu.pcin">railu.pcin</a> for parameters
%
% For revision information type <a href="matlab:railu.cvs">railu.cvs</a>

%       E. Balmes, G. Vermot des Roches, G. Guillet
%       Copyright (c) 1990-2025 by SDTools & SNCF, All Rights Reserved.

%#ok<*NOSEM,*PROP,*ASGLU,*STREMP>
properties (GetAccess = public , SetAccess = public)
% type
end

%{
```DocString  {module=base,src=railu.md#railu.f} 
 file handling methods 
%}

methods
function ob = railu(varargin)
     help railu
end
end % methods

%% #Static / generic methods  ----------------------------------------------
methods (Static)

function   out=asUo(r1);
%% #asUo : convert Dynavoie and rail inputs to UO format -2

 if isfield(r1,'ToolTip')&&isfield(r1,'data')
  r2=r1.data;
  DefBut=feval(dyn_ui('@genDefBut'));
  [table,level,loc]=sdtu.ivec.In2Tree(DefBut);
  cinM=containers.Map(table(:,1),table(:,2));
  out=sdtm.toCell(r2{2},'flatcell');
  for j1=1:size(out,1)
   st2=[r2{1},'.',out{j1,1}];
   if ~isKey(cinM,st2);
       st2=[r2{1},'.',regexprep(out{j1,1},'\.','P.','once')];
   end
   if isKey(cinM,st2);
     r3=cinM(st2);
     if isfield(r3,'ToolTip');out{j1,3}=r3.ToolTip;end
   end
   out{j1,1}=[r2{1} '.' out{j1,1}];
  end
 end
 if nargout==0;sdtm.toString(out);clear out;end
end
function out=pcin(varargin)
 %% #pcin : add cinM nodes
 CAM='pcin';
 % sdtu.f.open('@dynavoie/m/dynavoie_ui_en-us.csv')
 preRO={'key','ToolTip','DoOpt';
  'railu.Mesh.BeamMass','beam mass track slice model',{'DoOpt'
   'ncell(15#%g#"number of sleepers")' 
   'Slice.Gen.half';'Slice.ArmP.SMass';% 'half(#31#"1 one rail, 0 symmetric") ' 'ms(130#%g#"sleeper mass")' 
   'gaug(1.575#%g#"gaug [m|>20mm]") '
   'Slice.ArmP.Ltr' % 'Ltr(2.32#%g#"width") '
   'kp(300e6#%g#"pad stiffness") ' % Slice.ArmP.k_pad
   'kpr(1e6#%g#"pad torsion stiffness") '
   'kps(1e6#%g#"pad transverse torsion stiffness") '
   'kb(20e7#%g#"ballast stiffness") ' % value of Arlaud 2016
   'rand(0#%g#"ballast dispersion") '
   'randp(0#%g#"pad dispersion") '
   'simpack(0#%g#"use simpack compatible elements") '
   'Ftr(0#3#" renumber for simpack") '
   'eta(0#31#"use eta") '
   'cyc(0#31#"apply periodicity") '
   'unit(SI#%s#"unit system")' 
   'lc(.15#%g#"refine length [m|<1mm]")' 
   'pk(75e3#%g#"pad stiffness")';'projM';'preCase'}  
  'dyn_mesh.arm','spleeper/rail/pad',{'key';'Slice.ArmP.ArmType';'Slice.ArmP.RailType'
         'Slice.ArmP.PadType';'Slice.ArmP.SleeperType'
   }
  'railu.MeshPad','beam mass track slice model',{'DoOpt'
   'Slice.ArmP.PadType';'Slice.ArmP.lsem';'Slice.ArmP.hs';'Slice.ArmP.hs';% 'Slice.ArmP.hs>hs'
   'Slice.ArmP.k_pad'
   }
  'railu.RailSection','properties for rail section',{'DoOpt'
   'Slice.ArmP.RailType>RailType';'Slice.ArmP.hr3>hr3';'Slice.Gen.half>half'
   'RailX(0:100:600#%g#"rail cuts [mm|<1=m]")'
   }
   };
 r1=sdtu.f.ppath;carg=1;
 if nargin>0&&isa(varargin{1},'vhandle.nmap');cinM=varargin{1};carg=2;
 else; cinM=sdtm.pcin;
 end
 if ~isKey(cinM,'prero.railu.Range')
    DefBut=sdtm.rmfield(feval(dyn_ui('@genDefBut')),'PTree','Tab','Methods');
    sdtm.pcin(cinM,'append',DefBut);
 end
 if ~isKey(cinM,'prero.dyn_mesh.paramSlice')
  feval(dyn_mesh('@getParamEdit'),[],[],'Slice')
 end
 if ~isKey(cinM,'prero.dyn_mesh.paramTrack')
  feval(dyn_mesh('@getParamEdit'),[],[],'Track')
 end
 preRO(end+(1:2),1:3)={'railu.MeshSlice','Slice parameters',cinM('prero.dyn_mesh.paramSlice')
    'railu.MeshTrack','Track parameters',cinM('prero.dyn_mesh.paramTrack')
     };

  % augment cinM/osM using preRO/preOs
  sdtm.pInitPre([nargout exist('preRO','var') exist('preOS','var')]);
end

function mt=addVehicle(mt,mv,RO,RunOpt)
 %% #addVehicle : combine track and vehicle model with placement 
 if RunOpt.Reset 
  % Remove vehicle to allow reset
  % ve : defined by PID in the 1-100 range, nodes < 1000
  i1=setdiff(unique(mt.Node(:,2)),0);
  if ~isempty(i1)&&isfield(mt,'bas')&&~isempty(mt.bas)
    mt.bas(ismember(mt.bas(:,1),i1),:)=[];
  end
  i1=mt.Node(mt.Node(:,2)>0,1);
  if ~isempty(i1);
      mt.Elt=feutil('removeelt withnode',mt,i1);
      mt.Node(mt.Node(:,2)>0,:)=[];
  end
 end
 Geo=RO.Geo;
 safeDistM=dyn_mesh('@safeDistM');
  % put the vehicle *on* the rails : 
  if RO.half>0
    mv.Node(:,6)=(mv.Node(:,6)-min(mv.Node(:,6))); % centers y
  elseif RO.half<0
    mv.Node(:,6)=(mv.Node(:,6)+max(mv.Node(:,6))); % centers y
  else
    mv.Node(:,6)=mv.Node(:,6)-...
      (min(mv.Node(:,6))+max(mv.Node(:,6)))/2; % centers y
  end
  if max(mv.Node(:,6))>0
    if RO.double
        mv.Node(:,6)=mv.Node(:,6)./max(mv.Node(:,6))*diff(Geo.y_rail)+min(Geo.y_rail);
    elseif isfield(Geo,'y_rail')
        mv.Node(:,6)=mv.Node(:,6)./max(mv.Node(:,6))*max(Geo.y_rail);
    end
  end
  n1=(mt.Node(:,1)>1000); % mt without any vehicle
  if nnz(n1)==0; n1=true(size(mt.Node,1),1);end % Not shifted non-dynavoie
  if isfield(mv,'meta')&&isfield(mv.meta,'nomrad')
   mv.Node(:,7)=mv.Node(:,7)+Geo.z_rail_up+safeDistM(mv.meta.nomrad); %z/ vehicle put on the rails
   mv.bas=[];
  else
   mv.Node(:,7)=mv.Node(:,7)+max(mt.Node(n1,7))-...
     min(mv.Node(:,7)); %z/ vehicle put on the rails
  end

  % basis : add basis  1 99 range
  if ~isfield(mt,'bas');basid=398;mt.bas=[];
  else; basid=max([mt.bas(mt.bas(:,1)<100,1);0])+1;
  end
  y=mt.Node(mt.Node(n1,7)==max(mt.Node(n1,7)), ...
      6); %Ys of nodes on the top of the slice (pad nodes if any)
  if ~RO.half;
      ybas=(min(y)+max(y))/2;
  else
      ybas=0;
  end
  mt.bas=[mt.bas;
    basid 1 0    0.0 ybas 0.0  ...
    1.0 0.0 0.0  0.0 1.0 0.0  0.0 0.0 1.0]; % basis of the vehicle
  trajectory=RunOpt.trajectory;
  if ~isempty(trajectory) % trajectory given as an argument
    mt=dyn_solve(sprintf('trajectory %i',basid),mt,trajectory.table);
  end
    
  trajectory=stack_get(mt,'info','usertrajectory','getdata');
  if isempty(trajectory); sdtw('_nb','can''t find trajectory in stack'); 
  else
   RO.ID=[trajectory(:).ID];
   if ~any(RO.ID==basid)
    sdtw('_nb','trajectory %i should be fulfilled in the stack of the global model',basid);
   end
  end
  if ~isempty(trajectory) && any(RO.ID==basid)
    trajectory=trajectory(RO.ID==basid);
    mt.bas(end,4+size(trajectory.table,2)-2)=trajectory.table(1,2:end); % origin of the trajectory
  end
  
  % - associate nodes of mv and this basis :
  mv.Node(:,2:3)=basid*ones(size(mv.Node,1),2);
  
  % nodes and elts :
  % checks whether the free range of nodes at the beginning of track model
  % is sufficient :
  
  NodeId0=max([mt.Node(mt.Node(:,1)<1000,1);0])+1;

  mv=feutil('renumber',mv,...
    NodeId0:NodeId0+size(mv.Node,1)-1); % begins vehicles node numbering from the last vehicle

  if ~isempty(find(mt.Node(:,1)>=min(mv.Node(:,1)) & mt.Node(:,1)<=max(mv.Node(:,1)), 1)) 
    sdtw('_nb','there is not enough space in the range 1-1000 to add the vehicle')
    out=mt;
    return
  end
 
  % material and element properties : .il, .pl
  if ~isfield(mt,'il'); mt.il=[]; end;
  if ~isfield(mt,'pl')||isempty(mt.pl); mt.pl=zeros(1,5); end;
  if isfield(mv,'pl')
    if min(mv.pl(:,1))<401 || max(mv.pl(:,1))>500
      sdtw('_nb','vehicle matid should be in the 401-500 range'); end
    % matid shifting in mv.Elt :
    mpid=feutil('mpid',mv.Elt);
    mpid(mpid(:,1)~=0,1)=mpid(mpid(:,1)~=0,1)...
      -min(mpid(mpid(:,1)~=0,1))...
      +max([mt.pl(mt.pl(:,1)>400 & mt.pl(:,1)<=500,1);400])+1;
    mv.Elt=feutil('mpid',mv.Elt,mpid);
    % matid shifting in mv.pl :
    mv.pl(:,1)=mv.pl(:,1)-min(mv.pl(:,1))...
      +max([mt.pl(mt.pl(:,1)>400 & mt.pl(:,1)<=500,1);400])+1;
  end
  if isfield(mv,'il')
    if min(mv.il(:,1))<401 || max(mv.il(:,1))>500
      sdtw('_nb','vehicle proid should be in the 401-500 range'); end
    % proid shifting in mv.Elt :
    mpid=feutil('mpid',mv.Elt);
    mpid(mpid(:,2)~=0,2)=mpid(mpid(:,2)~=0,2)...
      -min(mpid(mpid(:,2)~=0,2))...
      +max([mt.il(mt.il(:,1)>400 & mt.il(:,1)<=500,1);400])+1;
    mv.Elt=feutil('mpid',mv.Elt,mpid);
    % proid shifting in mv.il :
    mv.il(:,1)=mv.il(:,1)-min(mv.il(:,1))...
      +max([mt.il(mt.il(:,1)>400 & mt.il(:,1)<=500,1);400])+1;
  end
  mv=stack_rm(mv,'info','OrigNumbering');
  mt=feutil('addtest -noOri;',mt,mv);
 
  % add model of vehicle to the others in the stack of the global model
  mo1=stack_get(mt,'info','mv','getdata'); % previous vehicles
  if isempty(mo1)||RunOpt.Reset
    mo1=mv;
  else
    mo1=feutil('addtest',stack_rm(mo1,'info','OrigNumbering'),mv);
  end
  mo1.pl=mt.pl;mo1.pl=feutil('getpl',mo1);
  mo1.il=mt.il;mo1.il=feutil('getil',mo1);
  mo1.bas=mt.bas(ismember(mt.bas(:,1),unique(mo1.Node(:,2:3))),:);
  % FixDof to have a 2D motion without x motion for the vehicle :
  [Case,vDOF]=fe_mknl('init',mo1); i1=fe_c(vDOF,.01*[1 2 4 6]','dof');
  mt = fe_case(mt,'FixDof','2D',i1);
  mo1=fe_case(mo1,'FixDof','2D',i1);

   r1=stack_get(mt,'SE','#s\d');
   if ~isempty(r1)&&~isfield(r1{1,3},'K') 
    %% for just mesh display full model (typically rail)
    mt=stack_set(mt,'SE','mv',mv);
    mo1=fesuper('seadd-unique -name"mv" -owrite',mt,mv);
    sdtm.store('base{mo1>mo1}')
    dbstack;
    feplot(mo1)
    mo2=fesuper('fsemodel-join2',mo1);
    feplot(mo2);fecom('showfipro')
   else
    mt=stack_set(mt,'info','mv',mo1); % vehicles
    % stack the contact dofs trajectory
    mt=dyn_solve('cdoftrajectory',mv,mt,trajectory);
   end
  out=mt; if nargout==0; feplot(out);fecom('colordatamat -alpha.3');end

end
function   mv=WheelSection(RO);
%% #WheelSection : wheel section mesh -1
 mv=[];
 if ischar(RO.VehType); 
  nmap=d_rail('nmap');
  if isKey(nmap,RO.VehType);RA=nmap(RO.VehType);
  elseif contains(RO.VehType,'.mat#');
   % railu.RailSection('U30_Sections.mat#ra')
   fname=RO.VehType(1:find(RO.VehType==35,1,'first')-1);tag=RO.VehType(length(fname)+2:end);
   FileName=d_rail('wd',fname);
   if ~exist(FileName,'file');error('''%s''< d_rail(''wd'',''%s'') not found',FileName,fname);end
   sections=[];sdtm.load(FileName,'sections');
   if isa(sections,'containers.Map');sections=vhandle.nmap(sections);end
   if ~isKey(sections,tag); error('%s not found in %s',tag,FileName);end
   mv=useOrDefault(sections,tag,'','','getValue');
   if ~isfield(mv,'name');mv.name=tag;end
   [~,fname,ext]=fileparts(fname); 
   mv.name=[fname '_' mv.name];
  else;
   nmap=nmap('Map:Sections');
   if isKey(nmap,RO);RA=nmap(RO);end
  end
 end
 if isempty(mv);error('Not implemented');end
 % Otr Middle of Track is typically at y==0
 y=mv.Node(:,6); r1=[min(y) mean(y) max(y)];
 if abs(RO.Geo.y_rail_m)<1&&max(abs(r1))>5
   %% convert to units of track (typically meters)
   mv.Node(:,5:7)=mv.Node(:,5:7)/1000;
   mv.pl=fe_mat(['convert' mv.unit 'SI'],mv.pl);
   mv.il=fe_mat(['convert' mv.unit 'SI'],mv.il);
   mv.unit='SI';
   r1=r1/1000;
 end
 n1=mv.Node(mv.Node(:,7)==max(mv.Node(:,7)),:);
 if sign(r1(1)*r1(2))==-1
  %% wheel is meshed at x==0
  mv.Node(:,5)=mv.Node(:,5)-n1(5);
  % cf.sel='withnode{x<=-1e-5} & selface & facing >.9 1000 0 0 & innode{x>-.01} & seledge'

 end
 if RO.half==-1 % Negative should be left wheel (positive y in dv)
  mv.Node(:,6)=(mv.Node(:,6)-r1(3));
 end 

end




function   m_rail=RailSection(RO,RC);
%{
```DocString {module=rail} -1
RailSection:  rail section mesh
```EXAMPLE
railu.RailSection('60-E1')
railu.RailSection('U30_Sections.mat#ra')
railu.RailSection(m_rail) % cleanup
```keywords
{var{mdl,pl},fcn{drmesh}}
%}
%%
name='';
if isfield(RO,'RailType');name=RO.RailType;
elseif ischar(RO);name=RO;
end
RA=struct;m_rail=[];
if ~isempty(name); 
  if strcmpi(name,'_genUIC60')    
%% #_genUIC60 manual generation of UIC60 by morphing -2
% G. Guillet section
 RO=railu.getNmap('60-E1'); 
 mo2=fe_gmsh('read',d_rail('wd','UIC60_GAE.msh'));mo2.Elt=feutil('selelt eltnamebeam1',mo2);
 mo2.Node=feutil('getnodegroupall',mo2);
 mo2.Node(:,7)=mo2.Node(:,7)-max(mo2.Node(:,7)); % prrLine is on top
 RO.prrLine=mo2.Node(mo2.Node(:,6)>0,5:7);

 f2=d_rail('wd','U30_Sections.mat');z=load(f2{1});mo3=z.sections('ra');
 m_rail=railu.RailSection(mo3); % remove cant manually

 mo3=m_rail;mo3.Node(:,6:7)=mo3.Node(:,6:7).*[RO.lar1/m_rail.meta.lar1 RO.hr3/m_rail.meta.hr3];
 mo3.Node(:,7)=interp1([0 -43 -140 -172],[0 -29 -140  -172],mo3.Node(:,7));

 n3=feutil('getnode inelt{seledge}',mo3);
 n2=feutilb('addNode-nearest -epsl 10',mo2.Node(:,5:7),n3(:,5:7));
 fecom('shownodemark',mo2.Node(n2,5:7))

 RM=struct('base',mo3.Node(mo3.Node(:,7)<-171.9,1),'Morph', ...
     {[{'Profile','Coarse'};num2cell(mo2.Node(n2,1)) num2cell(n3(:,5:7),2)]});
 mo4=fe_shapeoptim('Morph',feutil('addtest',mo2,mo3),RM);
 mo4.meta=sdth.sfield('addselected',mo3.meta,RO,{'lar1','hr3'});mo4.name='60-E1';
 mo4.meta.ToolTip='UIC60 ';
 mo4.meta.Holes='{h76.25, x 60 230 400,d23}';

 f2=d_rail('wd','UIC60_Sections.mat');r1=load(f2{1});
 r1.sections('ra')=mo4;sdtm.save(f2{1},'-struct','r1');
 
 %% now do the RA+BE section
 mo3=railu.RailSection(r1.sections('ra+be')); 
 mo3.Node(:,6:7)=mo3.Node(:,6:7).*[RO.lar1/m_rail.meta.lar1 RO.hr3/m_rail.meta.hr3]+[0 0];
 mo3.Node(:,7)=interp1([0 -43 -63.3 -140 -172.1],[0 -29 -172+76.25 -140  -172.1],mo3.Node(:,7),'linear','extrap');
 feplot(mo3);axis on; fecom view4

 n3=feutil('getnode inelt{seledge}',mo3);
 n2=feutilb('addNode-nearest -epsl 10',mo2.Node(:,5:7),n3(:,5:7));
 fecom('shownodemark',mo2.Node(n2,5:7))

 RM=struct('base',mo3.Node(mo3.Node(:,7)<-171.9,1),'Morph', ...
     {[{'Profile','Coarse'};num2cell(mo2.Node(n2,1)) num2cell(n3(:,5:7),2)]});
 mo4=fe_shapeoptim('Morph',feutil('addtest',mo2,mo3),RM);
 mo4.meta=sdth.sfield('addselected',mo3.meta,RO,{'lar1','hr3'});mo4.name='60-E1';
 mo4.meta.ToolTip='UIC60 ';
 mo4.meta.Holes='{h76.25, x 60 230 400,d23}';

 %% #UIC60_U85 with meshing using GMSH and Gaetan's cut -2
 mo4=abaqus('read',d_rail('wd','UIC60_U85.inp'));
 mo4.Elt=feutil('removeelt eltname bar',mo4);
 mo4.meta.Holes='{h76.25, x 60 230 400,d23}';
 mo4=feutil('addtest -noori;',mo4,feutil('symsel 35  0 1 0',mo4));

 feplot(mo4);

  return
  end
  nmap=d_rail('nmap');
  RB=struct;
  DoOpt=sdtm.pcin('prero.railu.RailSection');
  if any(name=='{'); 
    %[RB,st,name]=cingui('paramedit -DoClean2',DoOpt,{RO,name}); Cam=lower(CAM); st='';
    [name,RB]=sdtm.urnPar(name,DoOpt);
  end

  sectionM=[];if isKey(nmap,'Map:Sections');sectionM=nmap('Map:Sections');end
  if isKey(nmap,name);RA=followAlias(nmap,name);
  elseif isKey(sectionM,name);RA=followAlias(sectionM,name);
  elseif isKey(sectionM,RO.RailName);name=RO.RailName;RA=followAlias(sectionM,name);
  end
  RA=sdth.sfield('addmissing',RB,RA);
  if contains(name,'.mat#')||isfield(RA,'section')
  %% #RailSection.readMat railu.RailSection('U30_Sections.mat#ra')
   if isfield(RA,'section');st=RA.section;elseif ischar(name);st=name;end
   fname=st(1:find(st==35,1,'first')-1);tag=st(length(fname)+2:end);
   FileName=d_rail('wd',fname);if iscell(FileName);FileName=FileName{1};end
   if ~exist(FileName,'file');error('''%s''< d_rail(''wd'',''%s'') not found',FileName,fname);end
   sections=[];sdtm.load(FileName,'sections');
   if isa(sections,'containers.Map');sections=vhandle.nmap(sections);end
   if ~isKey(sections,tag); error('%s not found in %s',tag,FileName);end
   m_rail=useOrDefault(sections,tag,'','','getValue');
   i1=feutil('findelt eltname quad9',m_rail);
   if ~isempty(i1)
     m_rail.Elt(i1(1)-1,1:6)=[Inf abs('quadb')];
     m_rail.Elt(i1,:)=[m_rail.Elt(i1,[1:8 10:end]) i1*0];
   end
   m_rail=feutil('quad2lin',m_rail);m_rail.Node=feutil('getnodegroupall',m_rail);
   if ~isfield(m_rail,'name');m_rail.name=tag;end
   [~,fname,ext]=fileparts(fname); m_rail.name=[fname '_' m_rail.name];
  else;
   nmap=nmap('Map:Sections');
   if isKey(nmap,RO);RA=nmap(RO);end
  end
end
if isfield(RO,'ArmType');RA=RO;end
if isfield(RO,'Elt'); 
    m_rail=RO; RO=RC; % actually extrude
elseif isfield(m_rail,'Elt')
elseif isfield(RA,'lar1')
     %% dynavoie coarse rail model
      m_rail=struct('Node',...
     [1   0 0 0      0 -RA.lar1/2 0      % 1   y+1575/2 ?
      2   0 0 0      0 -RA.lar1/2 RA.hr1  % 2
      3   0 0 0      0 -RA.lar2/2 RA.hr1  % 3
      4   0 0 0      0 -RA.lar1/4 0       % 4
      5   0 0 0      0 +RA.lar1/4 0      % 5
      6   0 0 0      0 +RA.lar2/2 RA.hr1  % 6
      7   0 0 0      0 +RA.lar1/2 RA.hr1  % 7
      8   0 0 0      0 +RA.lar1/2 0     % 8
      9   0 0 0      0 -RA.lar3/2 RA.hr3  % 9
      10  0 0 0      0 -RA.lar3/2 RA.hr3-RA.hr2 % 10
      11  0 0 0      0 -RA.lar2/2 RA.hr3-RA.hr2  % 11
      12  0 0 0      0 -RA.lar2/2 RA.hr3  % 12
      13  0 0 0      0 +RA.lar2/2 RA.hr3 % 13
      14  0 0 0      0 +RA.lar2/2 RA.hr3-RA.hr2  % 14
      15  0 0 0      0 +RA.lar3/2 RA.hr3-RA.hr2  % 15
      16  0 0 0       0 RA.lar3/2 RA.hr3],...  % 16
                  'Elt',...
     [Inf abs('quad4');
      1  2  3  4   301 301;
      5  6  7  8   301 301;
      3  4  5  6   301 301;
      9  10 11 12  301 301;
      11 12 13 14  301 301;
      13 14 15 16  301 301;
      3  11 14 6   301 301]);
    if isfield(RA,'hs') % sleeper height (obsolete)
     m_rail.Node(:,7)=m_rail.Node(:,7)+RA.hs;
    end
    m_rail=feutil('Divide 2 2',m_rail);    
    if exist('dyn_mesh','file') % xxx should use PrePl 
     m_rail.pl=dyn_mesh('matdb 301',{{'Rail','Standard'}});
    else;
     m_rail.pl=d_rail('nmap.MatDb.rail');m_rail.pl(1)=301;
    end
    m_rail.il = [301  fe_mat('p_solid','SI',1) 0 10002 0 1]; % new integration rule
else; 
  error('Not implemented')
end

%% #RailSection.remove_cant
n1=sortrows(m_rail.Node,[6 7]); % bottom left
i1=feutil('geolinetopo',m_rail,struct('starts',n1(1),'dir',[0 1 0],'cos',.99));
n1=[n1(1,:);n1(n1(:,1)==i1{1}(end),:)];
cant=atan2(diff(n1(:,7)),diff(n1(:,6)));rot=[cos(cant) -sin(cant);sin(cant) cos(cant)];
if abs(cant)>1e-4 % Reorient the mesh
 m_rail.Node(:,6:7)=m_rail.Node(:,6:7)*rot;
end
r2=n1(:,6:7)*rot;
m_rail.Node(:,6:7)=m_rail.Node(:,6:7)+[-mean(r2(:,1)) -max(m_rail.Node(:,7))];
m_rail.meta=struct('GaugPoint',[NaN NaN],'lar1',diff(r2(:,1)), ...
    'OrigCant',abs(cant),'hr3',-min(m_rail.Node(:,7)));
m_rail.meta.OrigCant=abs(cant); 
% Node with mass element on the middle-top of the rail section  
n1=feutil('findnode z== epsl1e-4',m_rail,0);
m_rail.meta.TopNodes=n1; 


%%

  % checks properties ranges :
  if ~isfield(m_rail,'pl');m_rail.pl=[];end
  if ~isempty(m_rail.pl)
   if min(m_rail.pl(:,1))<301 || max(m_rail.pl(:,1))>400
    sdtw('_nb','rail matid should be in the 301-400 range'); 
   end
   if min(m_rail.il(:,1))<301 || max(m_rail.il(:,1))>400
    sdtw('_nb','rail proid should be in the 301-400 range'); 
   end
  end
  if ~isfield(m_rail,'il');m_rail=p_solid('default;',m_rail);end
  % Make sure it uses SI
  if max(max(abs(m_rail.Node(:,5:7))))>10
    m_rail.Node(:,5:7)=m_rail.Node(:,5:7)/1000;m_rail.unit='SI';
  end

  % position of m_rail
  m_rail.Node(:,5)=m_rail.Node(:,5)-min(m_rail.Node(:,5)); % x0=0
  if  max(abs(m_rail.Node(:,5)))==0;RC.Extrude=1;end
  Geo=struct;if isfield(RO,'Geo');Geo=RO.Geo;end
  if isfield(RA,'RailX'); if ischar(RA.RailX);RA.RailX=comstr(RA.RailX,-1);end
  elseif ~isfield(Geo,'XtrudeRail')||isempty(Geo.XtrudeRail); return % early return of section
  elseif RC.Extrude
    %% #Gamma rail track extrusion should support curved line, section straight -3
    % now position based on Or (rail contact point at 0 0) 
    nb_slices=1; % here only one slice    
    XtrudeRail=Geo.XtrudeRail;
    r3=XtrudeRail.'*ones(1,nb_slices)+RO.sw/2+...
        RO.sw*ones(length(XtrudeRail),1)*(0:nb_slices-1);
    r3=unique([0 r3(:).' nb_slices*RO.sw]);
    RA.RailX=r3; 
  end
  if isfield(RA,'RailX')&&~isempty(RA.RailX)
    if mean(diff(RA.RailX))>1; RA.RailX=RA.RailX/1000; end % Force SI
    m_rail=feutil('extrude 0 1 0 0',m_rail,RA.RailX); % extrudes section
    if ~isempty(feutil('selelt eltnamehexa',m_rail));
        RA.RailType='Volume';
    end
  end
  m_rail.Elt=feutil('orient;',m_rail);
  if isfield(Geo,'z_pad_up') % dynavoie 
   m_rail.Node(:,7)=m_rail.Node(:,7)+Geo.z_pad_up-...
     min(m_rail.Node(:,7)); % z0=top of the pad
  end
  if min(m_rail.Node(:,7))==max(m_rail.Node(:,7)) && min(m_rail.Node(:,6))==max(m_rail.Node(:,6)) % rail=beam
     m_rail.Node(:,7)=m_rail.Node(:,7)+RA.hr3; % top of the rail at the good z
  end


end % RailSection

function m_rail=RailPrr(RA,m_rail,ms);
  % #RailPrr Add rigid links on the sides of the volume rail -2
  if nargin==3
   %% dynavoie the pad is on the sleeper
   padxyz=feutil('getnode InElt{MatId 101 & selface & facing> 0 0 0 1e6}',ms);
   padxyz(:,1:4)=[];
  else 
   % the pad is on the rail
   padxyz=feutil('getnode InElt{ProId 101 & selface & facing> 0 0 0 -1e6}',m_rail);
   padxyz(:,1:4)=[];
   padxyz=[];
  end
done=1;
  if min(m_rail.Node(:,7))==max(m_rail.Node(:,7)) && ...
       min(m_rail.Node(:,6))==max(m_rail.Node(:,6)) % rail=beam : rigid link
  elseif strncmpi(RA.RailType,'Vol',3) 
    %% #prrLine_for_contact / rigid edges and rail guide for volume -3
    [node,i1]=sortrows(round(m_rail.Node*1e6),[5 7 6]);
    node(:,1)=m_rail.Node(i1,1);r2=unique(node(:,5));r3=unique(node(:,7));
    elt=[]; 
    for j1=1:length(r2)% x positions
      i1=find(node(:,5)==r2(j1));
      i2=i1(node(i1,6)==0&abs(node(i1,7)-r3(end))<1e-2);
      if isempty(i2) % Need to add centerline
        n2=r2; n2(1,3)=0;
        [m_rail.Node,i2]=feutil('addnode',m_rail.Node,n2/1e6); 
        el1=feutil('objectbeamline',m_rail.Node(i2,1));el1(2:end,3:5)=399;
        m_rail.Elt=feutil('addelt',m_rail.Elt,el1);
        done=0; break;
      end
      if j1==1||j1==length(r2);i1=node(i1); % Edges all face nodes
      else;  % only top nodes
        i1=feutil('findnode inelt{withnode} & x==',m_rail,node(i2,1),node(i2,5)/1e6);
      end
      i2=node(i2,1);i1=setdiff(i1,i2)*[1 1];i1(:,1)=i2;
      elt=[elt;i1]; %#ok<AGROW> % connections to centerline 
    end
    if ~isempty(elt);elt(:,3)=123;m_rail=feutil('addelt',m_rail,'rigid',elt);end

  else;
    m_rail.Node(:,7)=m_rail.Node(:,7)-min(m_rail.Node(:,7))+max(padxyz(:,3));
    if ~isfield(ms,'nmap');ms.nmap=vhandle.nmap;end
    nameM=ms.nmap('Map:SetName');
    nameM.append({'Mat:301','Rail';'Mat:101','Pad';'Mat:201','Sleeper'})

  end

  if ~done&& ~isempty(feutil('selelt proid399',m_rail))
    %% #prrLine_for_contact_sections (centered on y=0 for now) -3
    [node,i1]=sortrows(round(m_rail.Node*1e6),[5 7 6]);
    node(:,1)=m_rail.Node(i1,1);r2=unique(node(:,5));r3=unique(node(:,7));
    elt=[];
    for j1=1:length(r2)% x positions
      i1=find(node(:,5)==r2(j1));
      i2=i1(node(i1,6)==0&node(i1,7)==r3(end));
      if 1==2%j1==1||j1==length(r2);i1=node(i1);
      else; 
       i1=node(i1(abs(node(i1,6)-node(i2,6))<15e3&node(i1,7)-r3(end)>-3e3),:); 
       i1=i1(:,1);      
      end
      i2=node(i2,1);i1=setdiff(i1,i2)*[1 1];i1(:,1)=i2;
      elt=[elt;i1]; %#ok<AGROW>
    end
    elt(:,3)=123;
    m_rail=feutil('addelt',m_rail,'rigid',elt);
  end

end


function   [model,RO]=Mesh(varargin);
 %% #Mesh implement forwarding of sdtsys.stepmesh
 CAM=regexprep(varargin{1},'^Mesh\s*','','ignorecase');Cam=lower(CAM);carg=2;
 if comstr(Cam,'beammass')
 %% #Mesh.BeamMass (was Pinault) -2
 % initial revision was PJE/dynavoie/dv16.m

 DoOpt=sdtm.pcin('prero.railu.Mesh.BeamMass');
 %  railu.pcin('prero',{'d_rail.meshBeamMass',DefString});

if ~isempty(strfind(Cam,'getedit'));out=DefString; return;end 
if carg<=nargin; RO=varargin{carg};carg=carg+1;
  if isequal(RO,struct)&&carg<=nargin&&isfield(varargin{carg},'nmap');
   model=struct; 
   RO=varargin{carg};carg=carg+1;
  end
else; RO=struct;
end
[RO,st,CAM]=cingui('paramedit -DoClean2',DoOpt,{RO,CAM});Cam=lower(CAM);
RO=vhandle.uo.safeLenUnit(RO,'mm/1000');
 if RO.cyc; RO.ncell=1;
 elseif length(RO.ncell)<2;RO.ncell(2)=0;
 end
 if ~isfield(RO,'xsens'); RO.xsens=[];end


%      ms=feutil('setmat',ms,'dyn_mesh(matdb{399,prrLine})');
% feutil('setpro',ms,'d_rail(nmap.ProDb.prrLine)');
RO=railu.MatProDb(RO);

if RO.half==0
 %% #Two_rails -3
 % Nodes definition
 n1=[1 0 0 0  0 0 0;2 0 0 0 -.3 0 0;3 0 0 0 .3 0 0;
    4 0 0 0  0 0 -.1]; n1(:,6)=RO.gaug/2;
 model=struct('Node',n1);
 n1(:,6)=-RO.gaug/2; n1(:,1)=(5:8)';model.Node=[model.Node;n1;
     9 0 0 0  0 -RO.Ltr/2 n1(end);10 0 0 0  0 RO.Ltr/2 n1(end)
     11 0 0 0 0 0 n1(end)];
 % Element definition
 model.Elt=feutil('addelt','beam1', ...
    [2 1 301 301; 1 3 301 301; 6 5 301 301; 5 7 301 301; %Rail
     9 8 201 201; 8 11 201 201; 11 4 201 201; 4 10 201 201]); %Sleepers
 % Material and element properties 
 RO.PreIl={'name','ProId';'60E1beam',301;'SleeperA',201};
 RO.PrePl={'name','MatId';'Rail',301;'RSleeper',201};
 if RO.simpack
  error('Need review')
  model.Node(end+(1:2),:)=[11 0 0 0  0 RO.gaug/2 -.2;
      12 0 0 0  0 -RO.gaug/2 -.2];
  model=feutil('addelt',model,'beam1', ...
     [1 4 151 151 ; 4 11 152 152;
      5 8 151 151 ; 8 12 152 152]);
  model.pl(end+(1:2),:)=[
     151 fe_mat('m_elastic','SI',1) 1e9 .3    1
     152 fe_mat('m_elastic','SI',1) 1e9 .3    1];
  il=[151 fe_mat('p_beam','SI',1) 1e-5  .1*[RO.kpr RO.kps]/1e9  .1*RO.kp/1e9
      152 fe_mat('p_beam','SI',1) 1e-5  .1*RO.kpr/1e9*[.01 .01] .1*RO.kb/1e9];
  model=feutil('setpro',model,il);
  RO.preCase={'fixdof','2DRail','matid 301 -DOF 1 2 4 6', ...
     'FixDof','2dSleeper','matid 201 -DOF 1 2 5 6',...
     'FixDof','clamped bottom','z==-.2'};
 else
  %RO.PreIl
  n2=[10 4 11 8 9]';          % ballast nodes
  r3=[1 4 -03 0  0 0 RO.kp ;
      5 8 -03 0  0 0 RO.kp ;  % pads nodes
      n2(:) ones(length(n2),1)*[ 0 -03 0  0 0 RO.kb/length(n2)]; % ballast
      %1 4 -04 0  0 0 RO.kpr; 5 8 -04 0  0 0 RO.kpr; % axial rotation
      %1 4 -05 0  0 0 RO.kps; 5 8 -05 0  0 0 RO.kps; % transverse
      ];
  if isfield(RO,'eta')&&RO.eta % Use a loss factor
   loss=RO.eta; r3(r3(:,7)==RO.kb,10)=loss*RO.kb;
   loss=.1;  r3(r3(:,7)==RO.kp,10)=loss*RO.kp;
  else % set visc based on equivalence at target frequency
   loss=.05; freq=100; % ballast
   r3(r3(:,7)==RO.kb/length(n2),9)=loss*RO.kb./(2*pi*freq)/length(n2);
   loss=.1;  freq=500; % pad
   r3(r3(:,7)==RO.kp,9)=loss*RO.kp./(2*pi*freq)/length(n2);
  end
  model=feutil('addelt',model,'celas',r3);
  % Add random oscil around 50 Hz with 10 % sleeper mass
  if isfield(RO,'randm') % Add random mass
    RO.mm=40; model.Node=[model.Node;
        11 0 0 0  0 RO.gaug/2 -.11;
        12 0 0 0  0 -RO.gaug/2 -.11];
    model=feutil('addelt',model,'mass1',[11 RO.mm*[1 1 1];12 RO.mm*[1 1 1]]);
    RO.km=(40*2*pi)^2*RO.mm;
    r3=[4 11 -03;8 12 -03];r3(:,7)=RO.km;r3(:,10)=RO.km*.1;
    model=feutil('addelt',model,'celas',r3);
  end
  RO.preCase(end+(1:2),1:3)={'fixdof','2DRail','matid 301 -DOF 1 2 4 6';
     'FixDof','2dSleeper','matid 201 -DOF 1 2 5 6'};
 end
 
elseif RO.half==3
 %% #One_rail_and_volume_ballast -3
    mob=dv18('MeshCubeV2',struct('hetero','homog','lx',.6)); %ballast beam
    mob.name='BallastSoil';
    NodeB=feutil('GetNode x== & y== & z==',mob,mean(mob.Node(:,5)),...
        mean(mob.Node(:,6)),max(mob.Node(:,7))); %Upper middle node
    posnodeB=NodeB(5:7); indB=NodeB(1);
    indBmax=max(mob.Node(:,1)); 
    mo1=struct('Node',[indBmax+1 0 0 0 posnodeB+[0 0 .1];...
        indBmax+2 0 0 0 posnodeB+[-.3 0 .1];...
        indBmax+3 0 0 0 posnodeB+[.3 0 .1];...
        indB 0 0 0 posnodeB]);
    mo1.Elt=feutil('addelt','beam1',[indBmax+2 indBmax+1 301 301;...
        indBmax+1 indBmax+3 301 301;]); %rail beam
%     mo1=feutil('addelt',mo1,'mass1',[elt 0 0 80 0 0 0]); %sleeper mass
    mo1=feutil('addelt',mo1,'celas',[indBmax+1 indB 03 0 ...
        0 0 300e6*[1 0 0 .05]]); %railpad spring
    mo1.il=[301 fe_mat('p_beam','SI',1) 0 3.05e-5 3.05e-5 60/8000];
    mo1.pl=[301 fe_mat('m_elastic','SI',1) 210e9 .3 8000 0 .02];
    mo1.name='Track1';
    mo1=feutil('AddTest Merge',mo1,mob);
    
    %xxx_HP
    [NodeTc,IndTc]=feutil('GetNode x>= & x<= & z==',mob,...
        mean(mob.Node(:,5))-0.15,mean(mob.Node(:,5))+0.15,...
        max(mob.Node(:,7))); %Nodes connected to the sleeper
    [~,EltT]=feutil('FindElt ConnectedTo IndT');
    moT=struct('Node',NodeT,'Elt',EltT);
    moT=feutil('Extrude 2 0 0 0.17',moT); %Extrusion of 0.17cm
    moT=q16p('h8Toh125',moT);
 
elseif RO.half==4
 %% #One_rail_and_beam_ballast -3
 z=zeros(1,3);
 model=struct('Node',[1 z 0 0 0;2 z -.3 0 0;3 z .3 0 0;
    4 z  0 0 -.1;5 z 0 0 -.2;6 z -.3 0 -.2;7 z .3 0 -.2; ]);% Sleeper
 % Ballast
 model.Elt=feutil('addelt','beam1',[6 5 400 400 0  0 1 0;5 7 400 400 0   0 1 0]);
 model=feutil('refinebeam .05',model);
 model=feutil('addelt',model,'beam1',[2 1 301 301;1 3 301 301]);
 %model=feutil('addelt',model,'rigid',[5 4 12346]);
 model=feutil('addelt',model,'celas',[5 4 -3 0  0 0 10e9 0 10e9/(2*pi*30)*.5]);%ballast sleeper connect

 model=feutil('addelt',model,'mass1',[4 0 0 80 0 0 0]);
 model=feutil('addelt',model,'celas', ...
     [1 4 3 0  0 0 300e6*[1 0 0 .05]]); % This is the pad 
  
 % Add soil springs
 i1=feutil('findnode z==',model,min(model.Node(:,7)));
 i1(:,2)=0; i1(:,3)=-3; i1(:,7)=1e6;% xxx soil stiffness value
 model=feutil('addelt',model,'celas',i1);
 
 model=fe_case(model,'fixdof','2D',[.02;.04;.06]);
 model=fe_case(model,'fixdof','zSleeper','z==-.1 -DOF 2 4 5 6');
 model=fe_case(model,'fixdof','BallastCompression','z==-.2 -DOF 1');
 model.il=[301 fe_mat('p_beam','SI',1) 0 3.05e-5 3.05e-5 60/8000];error('revise')
 model.il=p_beam(model.il,'dbval 400 rectangle 1 .25');
 
 %model.Prepl={'name','MatId';'Rail',301;'Ballast',400}
 model.pl=[301 fe_mat('m_elastic','SI',1) 210e9 .3 8000 0 .25
     400 fe_mat('m_elastic','SI',1) 68e6 .40 1700 0 .02]; % Ballast
 model.name='BeamMass';

 %feplot(model);fecom('textnode');
 
 
elseif RO.half==5
 %% #3D FEM Track -4
 model=dyn_sem('ReadMeshV2');
else
 %% #OneRail -4
 model=struct('Node',[1 0 0 0  0 0 0;2 0 0 0 -.3 0 0;3 0 0 0 .3 0 0;
    4 0 0 0  0 0 -.1]);
 model.Elt=feutil('addelt','beam1',[2 1 301 301;1 3 301 301;]);
 model=feutil('addelt',model,'mass1',[4 0 0 80 0 0 0]);
 model=feutil('addelt',model,'celas', ...
     [1 4 03 0  0 0 300e6*[1 0 0 .05];
      4 0 03 0  0 0 RO.kb*[1 0 0 .05]]);
 model=fe_case(model,'fixdof','2D',[.01;.02;.04;.06;4.02;4.01;4.05]);
 model=fe_case(model,'fixdof','zSleeper','z==-.1 -DOF 1 2 4 5 6');
 RO.PreIl={'ProId','name';301,'RailA'};
 RO.PrePl={'MatId','name';[301 fe_mat('m_elastic','SI',1) 210e9 .3 8000 0 .02],'RailA'};
 model.name='OneRail';
end
model=stack_set(model,'info','MeshParam',RO);
if isfield(RO,'PreIl') % uses RO.projM('ProDb')
 model=feutil('setpro',model,RO);
end
if isfield(RO,'PrePl')% uses RO.projM('MatDb')
 model=feutil('setmat',model,RO);
end

n1=feutil('getnodematid301',model);
% Now add rail
x=[n1(n1(:,6)>-.1,5); RO.xsens(:)]; 
if isfield(RO,'addMass');x=[x;vertcat(RO.addMass{:,1})];end
x=sort(unique(round(x*1000)))/1000;
x=feutil(sprintf('refineline%.15g',RO.lc),x);%[x(1)-.3;x;x(end)+.3]);
model.Elt=feutil('removeelt matid 301',model);
for j1=1:2
 n2=x;n2(:,2)=abs(n1(1,6));n2(:,3)=n1(1,7); 
 if min(n1(:,6))==0;break;elseif j1==2; n2(:,2)=-n2(:,2);end
 [model.Node,i1]=feutil('addnode',model.Node,n2);i1=model.Node(i1,1);
 model=feutil('addelt',model,'beam1',[i1(1:end-1) i1(2:end) ones(length(i1)-1,2)*301]);
end

%ms=model; % slice model at this point
model.Elt(feutil('findelt eltname beam1',model),7)=1; % vz for vertical orient
model.unit='SI';

elseif comstr(Cam,'{')
 %% #MeshFromNomenclature d_rail('Mesh{U30{5,Gc,tc-350_400}:Ec41{Air}:W2{XaZa}}')
 RM=railu.nameToMeshRO(CAM(2:end-1));
 model=railu.MeshList(RM); % meshlist
elseif comstr(Cam,'trackli');
 %% #MeshTrackLi MeshRailDb : database of track configurations (includes before/after)
 
 if carg<=nargin;RM=varargin{carg};carg=carg+1;
 else; RM=struct; 
 end
 projM=RM.projM;% Expect projM of experiment

r1=struct('Nb',4,... % Number of bolt
       'Ns',5); % Number of sleepers (doubled by symmetry)
   % 504*2e3+.1 approximately 1 MN/m
RM=sdth.sfield('addmissing',RM,r1);
if RM.Ns==1 % Single slice for periodic computations
elseif isnumeric(RM.Nb); RM.Nb={num2str(RM.Nb),'Air'};
end

li=[];st1=RM.Nb{1};

li=d_rail(['MeshDb.' st1],projM);
if ischar(li{1})&&strcmpi(li{1},'seclab');li(1,:)=[];end
if isempty(li)||~iscell(li)
 error('''%s'' obsolete use of unknown Eclisse configuration',st1)
end
% switch RM.Nb{1} % Select configuration in rail19('nmap')
 % case {'61','xxx6bolt'};li=nmap('Ec61');
 % case {'62'};li=nmap('Ec62');  
 % case {'41','E120'};li=nmap('Ec41'); %%  Choice of 4 bolts
 % case {'42'};li=nmap('Ec42');
 % case {'0','no','Ec0'};li=nmap('Ec0');
 % otherwise
 %  error('''%s'', Not a known case see Ec entries in rail19(''nmap'')',RM.Nb{1})
 % end
li(cellfun(@(x)isequal(x,15),li(:,4)),4)={RM.Lc};

%% #MeshDb.add_sleeper before and after -4
Ns=2*RM.Ns; % Correspond of the real number of sleeper 
            % Permit to input also odd number at the beginning
if ~isfield(RM,'Lc');RM.Lc=15;end
if length(RM.Lc)==2
elseif isfield(RM,'Vel')
   RM.Lc(2)=RM.Vel*1e3/5e3*3; % advancement at 5 kHz 
   if isfield(RM,'topCoarse')&&any(isequal(RM.topCoarse,{1,'coar'}))
   else; RM.Lc(2)=RM.Lc(2)*3;
   end
end
if RM.Ns==1; % Single slice {seclab,xstart,xend,lc}
 li=[{'SecLab','xstart','xend','lc'};li];
elseif any(strcmpi(RM.Nb{1},{'Ec0','no'})); 
 % Sleeper width 255, cell length(600)
 r1=(0:Ns)'*600;r1=r1(:,[1 1 1])+[0 (600-255)/2 (600-255)/2+255];
 r1=r1-r1(end,1)/2;r1=[r1(1);reshape(r1(1:end-1,2:3)',[],1);r1(end,1)];
 li=[repmat({'ra'},length(r1),1) num2cell(r1)];li(:,4)={RM.Lc};
 li(diff(r1)==255,1)={'ra+s'};li{end,1}='';
 li=[{'SecLab','xstart','xend','lc'};li];
 
elseif ~any(strcmpi(RM.Nb{1},{'42'}));
 %% For all model except 4 bolts without sleeper
 li=[{'ra',li{1,2}-225,[],RM.Lc};li];

 for j1=1:(Ns-2)/2-1 % left side
  li=[{'ra+s',li{1,2}-255,[],RM.Lc};li];
  li=[{'ra',li{1,2}-345,[],RM.Lc};li];
 end
 li=[{'ra+s',li{1,2}-255,[],RM.Lc};li];% Add of the last term for the left side 
 li=[li;{'ra+s',li{end,2}+225,[],RM.Lc}];% Add the first term for the right side
 li=flipud(li); % We reverse the order to built the same way for the right

 for j1=1:(Ns-2)/2-1  % right side
  li=[{'ra',li{1,2}+255,[],RM.Lc};li];
  li=[{'ra+s',li{1,2}+345,[],RM.Lc};li];% Each time we add in a constant way 
 end
 li=[{'',li{1,2}+255,[],RM.Lc};li]; % Add the last term for the right side
 li=flipud(li);
 li=[{'SecLab','xstart','xend','lc'};li];

else
%% For 4 bolts without sleeper 
  li=[{'ra',li{1,2}-125,[],RM.Lc};li];
  
 for j1=1:(Ns-2)/2-1 % left side
  li=[{'ra+s',li{1,2}-255,[],RM.Lc};li];
  li=[{'ra',li{1,2}-345,[],RM.Lc};li];
 end
 li=[{'ra+s',li{1,2}-255,[],RM.Lc};li];% Add of the last term for the left side 
 %li=[{'ra',li{a,2}-225,[],RM.Lc};li];
 li=[li;{'ra+s',li{end,2}+125,[],RM.Lc}];% Add the first term for the right side
 li=flipud(li); % We reverse the order to built the same way for the right

 for j1=1:(Ns-2)/2-1  % right side
  li=[{'ra',li{1,2}+255,[],RM.Lc};li];
  li=[{'ra+s',li{1,2}+345,[],RM.Lc};li];% Each time we add in a constant way 
 end
 li=[{'',li{1,2}+255,[],RM.Lc};li]; % Add the last term for the right side
 li=flipud(li);
 li=[{'SecLab','xstart','xend','lc'};li];

end
if length(RM.Lc)==2
 r1=vertcat(li{2:end,2});
 li(2:end,4)={RM.Lc(2)};
 li(setdiff(find(r1<min(RM.top)),1),4)={RM.Lc(1)};
 li(find(r1(2:end)>max(RM.top))+2,4)={RM.Lc(1)};
end
%% #Mesh.TGrad : creation of gradient if not given in KC -4
if isfield(RM,'TGrad') % Selection of gradient types
 switch lower(RM.TGrad); 
 case 'ga';
  RM=sdth.sfield('AddMissing',RM,struct( ... 
       'Ngrad',3,... % From the 3rd sleeper, the gradient begins
       'C',32e3,... % Constant value of C later divided by surface
       'K',45e6,... % value for standard sleeper K (N/m) later divided by surface
       'Cedge',64e3,... % Maximum value of C for gradient
       'Kedge',90e6)); % Edge value K for gradient
     % Damping ratio xi=C/2sqrt(Km) with m=117.5kg the mass of half sleeper 
     % xi= 0.44 
  case 'gb';
  RM=sdth.sfield('AddMissing',RM,struct( ... 
       'Ngrad',3,... % From the 3rd sleeper, the gradient begins
       'C',32e3,... % Constant value of C later divided by surface
       'K',45e6,... % value for standard sleeper K (N/m) later divided by surface
       'Cedge',128e3,... % Maximum value of C for gradient
       'Kedge',180e6)); % Edge value K for gradient
     % Damping ratio xi=C/2sqrt(Km) with m=117.5kg the mass of half sleeper 
     % xi= 0.44 
 case 'gc';
 RM=sdth.sfield('AddMissing',RM,struct( ... 
       'Ngrad',3,... % From the 3rd sleeper, the gradient begins
       'C',32e3,... % Constant value of C later divided by surface
       'K',45e6,... % value for standard sleeper K (N/m) later divided by surface
       'Cedge',640e3,... % Maximum value of C for gradient
       'Kedge',900e6)); % Edge value K for gradient
     % Damping ratio xi=C/2sqrt(Km) with m=117.5kg the mass of half sleeper 
     % xi= 0.44
     case 'gd';
 RM=sdth.sfield('AddMissing',RM,struct( ... 
       'Ngrad',3,... % From the 3rd sleeper, the gradient begins
       'C',32e3,... % Constant value of C later divided by surface
       'K',45e6,... % value for standard sleeper K (N/m) later divided by surface
       'Cedge',32e3,... % Maximum value of C for gradient
       'Kedge',45e6)); % Edge value K for gradient
     % Damping ratio xi=C/2sqrt(Km) with m=117.5kg the mass of half sleeper 
     % xi= 0.44    
 otherwise
         error('''%s'' not a known gradient',RM.TGrad);
 end
else; error('Expecting TGrad');
end
if ~isfield(RM,'KC')
 K=[linspace(RM.Kedge,RM.K,RM.Ngrad)]'; % Built the gradient for K from input data
 C=[linspace(RM.Cedge,RM.C,RM.Ngrad)]'; % Built the gradient for C from input data
 RM.KC=struct('X',{{[],{'K';'C'}}},'Y',[K,C]);% define sleeper gradient
end
model=struct('KC',RM.KC,'rail',{li}); % sdtm.toString(RM.KC)
 elseif comstr(Cam,'joint')
%% #Mesh.Joint : meshing of rail joint from parts
RO=varargin{carg};carg=carg+1;
 error('Need implement')
 else; error('Need implement ''%s''',CAM)
 end
end


function   ms=MeshSlice(RO);
%% #MeshSlice : rail section mesh -2
% railu.MeshSlice('{UIC60{RailX0:50:600},M450{Slx-.6 -.105 .105 .6},SN,half0,cant.05}') %

 RO=sdtu.ui.cleanEntry(RO);
 if iscell(RO) 
  %% dynavoie list
  Range=dyn_solve('Range',struct('SliceCfg',{{'T',RO}}));
  dyn_solve('RangeLoopReset',Range,struct('ifFail','error'));  
  if sdtdef('isinteractive')
   PA=dyn_ui('paramvh'); cf=feplot(PA.ms);fecom(cf,'showfimat');
  end
 else % rail mesh 
  if isfield(RO,'nmap')&&isKey(RO.nmap,'CurSlice')
   RT=RO; RO=RT.nmap('CurSlice');
  else
   DoOpt=sdtm.pcin('prero.railu.MeshSlice');R1=RO;
   [~,RO]=sdtm.urnPar(RO,'{RailType%s,SleeperType%s,PadType%s}{}',struct('var',{DoOpt},'out','uo'));
  end
  RO=vhandle.uo.safeLenUnit(RO,'mm/1000');
  m_rail=railu.RailSection(RO);
  m_rail=railu.MeshPad(RO,m_rail);
  m_rail=railu.RailPrr(RO,m_rail); %% add top line
  % xxx robust vhandle.uo.safeLenUnit(RO,'mm/1000');
  if max(abs(m_rail.Node(:,7)))>1; m_rail.Node(:,5:7)=m_rail.Node(:,5:7)/1000;m_rail.unit='SI';end

  mo2=railu.MeshSleeper(RO,struct('rail',m_rail,'Slice',RO)); % combine rails
  sdtm.store(RT.nmap,sprintf('mo2>%s',RO.store));
 end

end
function   ms=MeshSliceCheck(ms,RC)
  %% #MeshSliceCheck

  if nargin<2||~isfield(RC,'Do');
    RC.Do={'EltId','plrange','railtop'}
  end
  for j0=1:length(RC.Do)
  switch lower(RC.Do{j0})
  case 'eltid'
   % #MeshSliceCheck.EltId fix eltids
   [u1,ms.Elt]=feutil('EltIdFix;',ms);
  case 'plrange'
   % #MeshSliceCheck.PlRange xxx sdtweb drmesh.matdb
   if min(ms.pl(:,1))<1 || max(ms.pl(:,1))>400
    sdtw('_nb','slice whith rail matid should be in the 1-400 range'); 
   end
   if min(ms.il(:,1))<1 || max(ms.il(:,1))>400
    sdtw('_nb','slice with rail proid should be in the 1-400 range'); 
   end
   case 'railtop'
   % Nodes of top rail
   r1=feutil('findnode pronameprrLine',ms);
   if isempty(r1)% old attempt to use line
    r1=fe_case(ms,'getdata','rail');
    if isempty(r1);error('Missing DofSet,Rail');end
    r1=unique(fix(r1.DOF));
   end
   RailSel=RC.projM('RailSel')
   dbstack; keyboard; 
   RailSel(j1*2-1,1:2)={RB.sname,feutil('addelt','mass1',r1)};
 case 'doperiod'
 %% #MeshSliceCheck.doPeriod
   %% xxx was DoPeriod
    if ~isfield(ms,'meta')||~isfield(ms.meta,'sw'); 
      r2=ms.Node(:,5); ms.meta.sw=max(r2)-min(r2);
    end
  data=fe_case(ms,'stack_get','cyclic','Symmetry','get'); 
  if isempty(data)||~isfield(data,'IntNodes')
  % periodicity condition :
  ms=fe_cyclic(sprintf('build -noslave -1 %g 0.0 0.0;',ms.meta.sw),ms); % builds periodicity
  % interfaces node renumbering :
  c1=fe_case(ms,'getdata','Symmetry');c1.IntNodes=sortrows(c1.IntNodes);  
  i1=[c1.IntNodes(:,1);setdiff(ms.Node(:,1),c1.IntNodes(:));
       c1.IntNodes(:,2)]; i1(:,2)=(1:length(i1))';
  if ~isfield(RC,'NoRenum')
   r2=fe_case(ms,'stack_get');
   m2=feutil('renumber-noOri;',stack_rm(ms,'info','OrigNumbering'),i1); 
   if size(r2,1)~=size(fe_case(ms,'stack_get'),1)
       error('Renumbering issue');
   else; ms=m2; 
   end
  end
 R1=stack_get(ms,'','MeshParam','get');
 [Case,name]=fe_case(ms,'getcase');
 r1=stack_get(Case,'','HalfTrack','g');
 r3=stack_get(Case,'','PmlDispMpc','g');
 if ischar(r1)&&~isempty(r3); [st,un1,r2]=comstr('-dof',[-25 1],r1,lower(r1));
    r1=feutil('getdof',feutil(['findnode' st],ms),r2/100);
    r1(fe_c(r1,r3.DOF(r3.slave),'ind'))=[];% Don't fix slaves
    ms=fe_case(ms,'FixDof','HalfTrack',r1);
    i1=fe_c(r3.DOF,r1,'ind'); r1=r3.DOF(r3.slave); 
    r3.c(:,i1)=[];r3.DOF(i1)=[];r3.slave=fe_c(r3.DOF,r1,'ind');
    ms=fe_case(ms,'mpc','PmlDispMpc',r3);
        
 end

  end % actually build

      otherwise 
          error(RC.Do{j0})
  end
  end % Loop on things to do
   



end
function   mt=MeshTrack(RT);
%% #MeshTrack : concatenate mesh slices as a track -2

RO=RT.nmap('CurTrack');projM=RT.nmap;
RO=vhandle.uo.safeLenUnit(RO,'mm/1000');

if ~isfield(RO,'Do'); RO.Do={'GenTrack','presens'}
end

for j0=1:length(RO.Do)
switch lower(RO.Do{j0})
case 'gentrack'  % generate track
% sdtweb dyn_mesh GenTrack
RB=struct('projM',projM); projM('RailSel')={};
RB.Do={'doPeriod','plrange','eltid'};
RB.bas=101; RB.NodeId0=100000;
RB.XinfoCol={'x','SE','bas'};RB.Xinfo=[];
RB.PreSleeper={'index','x'};
mt=struct('Node',[],'Elt',[],'bas',[]);
   if ~isfield(mt,'bas')||isempty(mt.bas)
   elseif max(mt.bas(:,1))<101;error('Expecting BasId<100 for train');
   end
   % - free the 1-100 basid range for the vehicle model 
   % - basid=100+index_slice for the slices : last row end of track with no s_
   % the s(j)i(k) use the j+100 basid
   % - NodeId of s(j) s(j)i(j) overlap by data.IntNodes 

  for j1=2:size(RO.preTrack,1) % Deal with multiple slices
   ev1=vhandle.tab.getValEvt(RO.preTrack,j1);RB.sname=ev1.SE;
   ms=projM(ev1.SE);
   ms=railu.MeshSliceCheck(ms,RB);   %RS=getParamInStack(ms);
   RB.xmin=min(ms.Node(:,5));
   data=fe_case(ms,'stack_get','cyclic','Symmetry','get');
   mt=stack_rm(mt,'SE',ev1.SE);
   i4=length(ev1.slNum); 
   %if length(RO.nb_slices)>1&&j1<length(RO.nb_slices);i4=i4-RO.nsi; end
   st1=sprintf('SEadd -trans %i %g 0 0 %i %i %i %s -bas %i',...
    i4,ms.meta.sw,size(data.IntNodes,1), ... % Shift interface nodes (should be full not reduced)
    RB.NodeId0,100000,RB.sname, ...
    RB.bas);
   RA=struct('Type','translated','opt',[i4 ms.meta.sw 0 0], ...
       'offset',[0 0 0],'bas',RB.bas,'MatId0',(1000+j1-1),'ProId0',(1000+j1-1), ...
       'NodeShift',size(data.IntNodes,1),'NodeId0',RB.NodeId0-1,'EltId0',1e5-1,'name',ev1.SE);
   if ~isempty(RB.Xinfo) % Position ms horizontally
    RA.bas(2:15)=[1 0  RB.Xinfo(end,1)-RB.xmin 0 0 reshape(eye(3),1,9)];
   else
    RA.bas(2:15)=[1 0  -RB.xmin 0 0 reshape(eye(3),1,9)];
   end
   ms.name=RA.name; 
   if ~isfield(mt,'nmap');mt.nmap=ms.nmap;end

   mt=fesuper(st1,mt,ms,RA); % add with translations and node range shifting
   mt=feutil('setpro',mt,[1000+j1-1 fe_mat('p_super','SI',1) 1. 1. 1.]);

   r1=mt.bas(sdtu.fe.NNode(mt.bas(:,1),RA.bas(1)-1+(1:i4)'),[4 1 1]);
   r1(1:end-1,2)=fesuper(['s_' RB.sname]); 
   r1(:,1)=r1(:,1)+min(ms.Node(:,5)); r1(end+1,1)=r1(end,1)+ms.meta.sw;
   if j1==2; RB.Xinfo=r1; % First slice
   else;
       RB.Xinfo=[RB.Xinfo(1:end-1,:);r1]; 
   end
   r3=ms.meta.SlX-RB.xmin+r1(1);
   if isscalar(r3);r3=r3+(0:RA.opt(1)-1)*ms.meta.sw;end 
   r3=r3(:);r3(:,2)=size(RB.PreSleeper,1)+(1:length(r3))-1;
   RB.PreSleeper(end+r3(:,2)-r3(1,2)+1,1:2)=num2cell(r3);
   if isequal(ev1.slNum,RO.slNum0) % Force origin 
     RB.sGamma0=r1(1)-RB.xmin; 
   end

   RB.NodeId0=mt.Elt(end,3)-size(data.IntNodes,1)+1;
   RB.bas=max(mt.bas(:,1))+1; % Leave space for an interface slice
  end % j1 for multi-slice 
  st1='mt>mt'; sdtm.store(RT.nmap,st1);
case 'presens'
  %% #MeshTrack.presens
  
  r3=cell2mat(RB.PreSleeper(2:end,:));r3=[r3 round((r3-[RB.sGamma0 0])./[.6 1],2)];
  %RB.preNode={'NodeId','lab'}
  st1=cellfun(@(x)sprintf('S\\%+.0f([+-]*)',x),num2cell(r3(:,3)),'uni',0);
  st1(:,2)=cellfun(@(x)sprintf(' %g $1',x),num2cell(r3(:,1)),'uni',0);
  colM=sdtu.ivec('ColList',RO.preSens(1,:));
  i1=colM('X');RB.iXYZ=colM({'X','Y','Z'});
  ta=RO.preSens(2:end,:);
  st2=[ta(:,i1) regexprep(ta(:,i1),st1(:,1),st1(:,2))];
  for j1=1:size(st2,1)
    try;r4=eval(st2{j1,2});catch;r4=0;end
    ta{j1,i1}=r4;RB.marker=ta{j1,colM('M')};
    pos=[ta{j1,RB.iXYZ}];
    if strcmpi(RB.marker,'mgl'); 
    elseif strcmpi(RB.marker,'mrl')
     bas=ms.bas(ms.bas(:,1)==399,:);
     pos=bas(4:6)'+reshape(bas(:,7:15),3,3)*pos(:);
     ta(j1,RB.iXYZ)=num2cell(pos);
    elseif strcmpi(RB.marker,'mrr')
     bas=ms.bas(ms.bas(:,1)==398,:);
     pos=bas(4:6)'+reshape(bas(:,7:15),3,3)*pos(:);
     ta(j1,RB.iXYZ)=num2cell(pos);
    end
  end
  ta(:,colM('NodeId'))=num2cell(-(1:size(ta,1))');
  % place sensors 
  RS=struct('tdof',{[colM.prop.name';ta]});
  RS.tdof(:,colM('M'))=[];
  mt=fe_case(mt,'SensDof','Out',RS);

otherwise; error('%s',RO.Do{j0})
end
end % Do Loop

end

function RM=MeshSleeper(name,evt);
 %% #MeshSleeper 
 % mo1=railu.MeshSleeper('m450pi')
 if isa(name,'vhandle.uo'); RS=name;
     name=RS.SleeperName;
 else; RS=struct; 
 end
 if nargin==1;evt=struct;end
 if strncmpi(name,'m450pi',4)
 %  #MeshSleeper.M450pi    sdtu.f.open('@onedrive/*/sncf*/e*/21*/Sleeper_IN00213.pdf#page=116') -3
 RO=struct;
 preMesh=struct('cuts', ... % mm
  {{[0     0  ;...
     290/2 0  ;...
     220/2 190;...
     0     190]...
    [0     0  ;...
     290/2 0  ;...
     220/2 220;...
     0     220]...
    [0     0  ;... 
     290/2 0  ;...
     220/2 220-350/20;...
     0     220-350/20]...  % pente 1/20      
    [0     0  ;...
     220/2 0  ;...3
     180/2 220-850/20;...
     0     220-850/20]...
    [0     0  ;...
     220/2 0  ;...
     180/2 220-850/20;...
     0     220-850/20]}}, ...
     'y',...
     [0 150 500 1000 1200]-1200);

 for j1=1:length(preMesh.cuts)
  node=[preMesh.cuts{j1}(:,1) repmat(preMesh.y(j1),length(preMesh.cuts{1}),1) preMesh.cuts{j1}(:,2)];
  preMesh.cuts{j1}=node;
 end
 mo2=feutil('objecthexa',[0 0 0;0 0 1;0 1 0;1 0 0],1,1,preMesh.y);mo2.Node(:,5:7)=vertcat(preMesh.cuts{:});
 mo2.Elt(2:end,[3 4 7 8])= mo2.Elt(2:end,[4 3 8 7]);
 mo2=feutil('addtest -noori;',mo2,feutil('symsel 17  1 0 0',mo2));
 mo2=feutil('addtest -noori;',mo2,feutil('symsel 17  0 1 0',mo2));
 mo2.Elt=feutil('orient;',mo2);
 mo2.Node(:,5:7)=mo2.Node(:,5:7)/1000;mo2.unit='SI';
 RO=railu.MatProDb(RO);
 mo2.Elt=feutil('setgroupall matid 201 proid 201',mo2);
 RO.PrePl={'name','MatId';'Sleeper',201};
 RO.PreIl={'name','il';'Sleeper',p_solid('dbval 201 d3 -3')};
 mo2=feutil('setpro',mo2,RO);
 mo2=feutil('setmat',mo2,RO);
 % adjust mass
 mo2.meta=struct('half',0,'mass',300,'cant',.05);
 mo2.Node(:,7)=mo2.Node(:,7)-.3655-2e-4; 'xxx zOff'
 r2=feutilb('geomrbmass',mo2); mo2.pl(1,5)=mo2.pl(1,5)/r2.mass*mo2.meta.mass;
 if isa(RS,'vhandle.uo'); % no need to return RS
     RS=sdth.sfield('addselected',RS,mo2.meta,{'half','cant'});
 end
 RM=mo2;
 if nargout==0; feplot(mo2);end

 end % type of sleeper
 if max(abs(mo2.Node(:,7)))>1; mo2.Node(:,5:7)=mo2.Node(:,5:7)/1000;mo2.unit='SI';end
 if isfield(evt,'rail')
  %% #MeshSleeper.Rail append rail to sleeper -3
  % reimplementation of sdtweb dyn_mesh Arm_Add_Rail 
  Slice=evt.Slice;
  if ~Slice.half;evt.LR=[1 -1]; % dynavoie has forward, zUp thus left is y>0
  elseif Slice.half>0; evt.LR=1; else;evt.LR=-1;
  end
  if isfield(Slice,'SlX')&&~isempty(Slice.SlX)
   %% replicate sleepers as needed
   for j1=1:length(Slice.SlX)
    mo3=mo2;mo3.Node(:,5)=mo3.Node(:,5)+Slice.SlX(j1);
    if j1==1; mo1=mo3;
    else; mo1=feutil('addtest-noori;',mo1,mo3);
    end
   end
   mo2=mo1;mo2.meta.SlX=Slice.SlX; 
  end
  nameM=mo2.nmap('Map:SetName');if ~isfield(mo2,'bas');mo2.bas=[];end
  for j1=1:length(evt.LR)   
   m_rail=evt.rail; Slr=evt.LR(j1); if Slr>0;st1='MrL';else;st1='MrR';end 
   cant=-Slr*Slice.cant; rot=[cos(cant) -sin(cant);sin(cant) cos(cant)];
   if ~isfield(Slice,'Or_tr'); Slice.Or_tr=[0 Slice.gaug/2 0];end
   r3=[0 Slice.Or_tr(2)*Slr -abs(Slice.Or_tr(3))];
   r2=eye(3);r2(2:3,2:3)=rot; rot=r2';
   m_rail.Node(:,5:7)= (m_rail.Node(:,5:7)*rot')+r3;
   mo2=feutil('addtest-noori;',mo2,m_rail);
   rot=rot*diag([1 -1 -1]);
   bas=[400-j1 1 0    r3  rot(:)'];
   nameM(sprintf('bas:%i',bas(1)))=st1;
   mo2.bas(end+1,1:length(bas))=bas;
  end
  RM=mo2; 
  %feplot(mo2); fecom view4; fecom ShowBas
 end

end
function m_rail=MeshPad(evt,m_rail);
 %% #MeshPad

 switch lower(evt.PadType)
 case 'sn'
  % 
  warning('need pad mesh cleaning SlX')
  %padxyz=feutil('getnode InElt{selface & facing> 0 0 0 -1e6}',m_rail)
  mo2=m_rail;mo2.Elt=feutil('selelt selface & innode{ z==}',m_rail,min(m_rail.Node(:,7)));
  mo2=feutil(sprintf('extrude 1 0 0 -%g',evt.hs),mo2);
  mo2.Elt=feutil('set groupall proid 101 101',mo2);
  m_rail.Node=mo2.Node;m_rail.Elt=feutil('addelt',m_rail.Elt,mo2.Elt);

  %dyn_mesh('pad',struct('PadType','volume','lsem',220,'hs',9,'ws'))
 end

end

function RO=MatProDb(RO);
 %% #MatProDb initialized database 
 r2=d_rail('nmap');
 if isfield(RO,'nmap')&&isa(RO.nmap,'vhandle.nmap');RO.projM=RO.nmap; 
 elseif ~isfield(RO,'projM');RO.projM=r2;
 end
 if ~isKey(RO.projM,'MatDb');RO.projM('MatDb')=r2('MatDb');end
 if ~isKey(RO.projM,'ProDb');RO.projM('ProDb')=r2('ProDb');end
end

function out=MeshList(RO)
%% #MeshList : joint meshing strategy from list -2

li=RO.Track; projM=RO.projM; 
%%

if ~isfield(RO,'gap');RO.gap=d_rail('MeshGap');end
 mo1=[]; 
 for j1=2:size(li,1) 
  %% #Mesh.Generate_segments from sections -3
  if isempty(li{j1}); continue;end
  if isempty(li{j1,3});li{j1,3}=li{j1+1,2};end % xend
  % Refine along x using lc 
  r1=feutil(sprintf('refineline %g',li{j1,4}),[li{j1,2:3}]);
  if max(abs(r1))>10; RO.unit='MM';else;RO.unit='SI';end
  st=sdth.findobj('_sub:',li{j1,1}); st={st.subs};% Split at :
  mo2=d_rail(sprintf('MeshDb.%s.%s',RO.rail,st{1}),projM);%U30.e ...

  if length(st)>1&&RO.gap.isKey(st{2})
   % change properties for the gap 'ra+e:Air'
   r2=RO.gap(st{2});
   mo2.Elt=feutil('setsel mat 104 pro 104',mo2,'proid7');
   mo1=feutil('setmat',mo1,fe_mat(['convert SI' RO.unit], ...
       [104 fe_mat('m_elastic','SI',1) r2.E,r2.nu,r2.rho]));
   mo1.info.xgap=[li{j1,2:3}];mo1.info.GapMat=st{2};
  end
  if all(mo2.Node(:,5)==0);mo2=feutil('extrude 0  1 0 0',mo2,r1);
  else;mo2.Node(:,5)=mo2.Node(:,5)+li{j1,2};
  end
  if isempty(mo1);% First segment initialized model
      mo1=mo2;
      mo1=sdth.sfield('addselected',mo1,RO,{'pl','il','unit','Stack'});
      mo1.nmap=vhandle.nmap; % should use material names here 
      mo1.nmap{'Map:MatName'}(5)='Bolt';
      mo1=sdth.urn('nmap.mat.set',mo1,struct('value',int32([1 4 5 7 8 9 10 11 ])', ...
          'ID',{{'Shaft','Eclisse','Bolt','Rail','Wheel','Pad','Sleeper','Screw'}}));
      mo1.nmap{'Map:ProName'}(5)='Bolt'; %xxx should not be needed 
      mo1=sdth.urn('nmap.pro.set',mo1,struct('value',int32([1 2 3 4 5 7 8 9 10 11 13])', ...
          'ID',{{'Shaft','CtcRail','CtcWheel','Eclisse','Bolt','Rail','Wheel','Pad','Sleeper','Screw','RailEdge'}}));

  else; mo1=feutil('addtest Merge -noori;',mo1,mo2);
  end
 end
 mo1.Node=feutil('getnode groupall',mo1);

 %%  #MeshList.Contacts -3
 mo1=feutil('joinhexa8',mo1); if ~isfield(mo1,'info');mo1.info=struct;end
 mo1=feutil('joinpenta6',mo1);
 if ~isfield(RO,'top')
 elseif ~isfield(RO,'Wheel')||isequal(RO.Wheel,{'W0'})
   i1=feutil('findnode z>82 & x> & x<',mo1,RO.top(1),RO.top(2));
   if isempty(i1); warning('No selected node t%s',sdtm.toString(RO.top));
   else
    mo1=fe_case(mo1,'DofLoad','Impact',struct('def',ones(size(i1)),'DOF',i1+.03));
   end
   RO.Wheel={'W0'};
 else
   d_rail('MeshRailDist{back}',projM); % Define rail profile no display
   if isnumeric(RO.top)
     RO.topsel=sprintf('selface & withnode {z>82 & x>%.15g & x<%.15g}',RO.top);
     mo1.info.top=RO.top;
   end
   if isfield(RO,'topsel')
    r2=struct('sel',RO.topsel);r2=sdth.sfield('addselected',r2,RO,{'topCoarse','projM'});
    mo1=d_rail('MeshRailTop',mo1,r2);
   end
 end

 %% #MeshList.Wheel -3
 if isfield(RO,'Wheel')
  RW=struct('wfrom','WheelCut','wref','coarse1','Lc',RO.Lc,'integ',-1, ...
      'top',RO.top);
  RW.projM=RO.projM; 

  if length(RO.Wheel)>1
   i1=find(strncmpi(RO.Wheel{2},'in',2));%in104
   if ~isempty(i1);RW.integ=comstr(RO.Wheel{2}{i1}(3:end),-1);end
   i1=find(strncmpi(RO.Wheel{2},'wref',2));%wref
   if ~isempty(i1);RW.wref=comstr(RO.Wheel{2}{i1}(5:end),1);end
  end
  dbM=projM('Map:Sections');
  if strcmpi(RO.Wheel{1},'wa');RO.Wheel{1}='Wheel';
      error('eb check')
  elseif strcmpi(RO.Wheel{1},'w1');
    %% W1 single point contact
    mo2=d_rail('MeshWheel',RW);
    dbM(RO.Wheel{1})=mo2;

  elseif strcmpi(RO.Wheel{1},'w2');
    %% W2 Change density to add brake disc
    if RW.integ==-1; RW.integ=-3; end
    mo2=d_rail('MeshWheel',RW);
    %mo1=feutil(sprintf('setpro 3 Integ=%i Match -2',RW.integ),mo1);
    %mo1.pl(1,5)=mo1.pl(1,5)+1.995e-5;
    dbM(RO.Wheel{1})=mo2;
  elseif ~isKey(dbM,RO.Wheel{1})&&~strcmpi(RO.Wheel{1},'w0')
    warning('''%s'' not known using ''Wheel''',RO.Wheel{1})
    RO.Wheel{1}='Wheel';
  end
 else; RO.Wheel={'Wheel','XuZu'};
 end
 if ~strcmpi(RO.Wheel{1},'w0')
  %% there is a wheel
  mo2=dbM(RO.Wheel{1});
  % place wheel at 0 
  data=mo2.info; %stack_get(mo1,'','wheel_geom','g');
  mo2.Node(:,5)=mo2.Node(:,5)-data.Origin(1);data.Origin(1)=-2;
  mo1=stack_set(mo1,'info','wheel_geom',data);
  %-(min(mo2.Node(:,5))+max(mo2.Node(:,5)))/2;
  mo1=feutil('addtest -noori -keeptest;',mo1,mo2);
  try; 
   mo1.nmap(sprintf('Wheel:%s',RO.Wheel{1}))=mo2.info;
  end 
  % renumber contact xxx check
  %n1=feutil('findnode matid 2;',mo1); 
  %if ~isempty(n1)
  % n1(:,2)=max(mo1.Node(:,1))+10e3+n1(:,1); 
  % mo1=feutil('renumber -noOri ;',mo1,n1);
  % mo1.Node=sortrows(mo1.Node);
  %end
  mo1=feutil('addsetNodeId',mo1,'Wheel','matid 1 8');
  mo1=feutil('addsetNodeId',mo1,'Rail','matid 7');
  mo1=d_rail(['MeshCase' RO.Wheel{2}{1}],mo1);%XaZu, ...
 end
 mo1.info=sdth.sfield('addselected',mo1.info, ...
     sdtm.rmfield(RO,'Stack','il','pl','sections','unit','topsel','gap'));
 %% #MeshList.KC_Under sleeper springs ProId 12  -3
  
 mo2=mo1;mo2.Elt=feutil('selelt selface & innode {z==}',mo2,min(mo1.Node(:,7)));
 mo2.Elt=feutil('divideingroups',mo2);
 [EGroup,nGroup]=getegroup(mo2.Elt);elt=[];
 NNode=sparse(mo2.Node(:,1),1,1:size(mo2.Node,1));
 for jGroup=1:nGroup
     cEGI=EGroup(jGroup)+1:EGroup(jGroup+1)-1;
     i1=unique(mo2.Elt(cEGI,1:4));
     n1=mo2.Node(NNode(i1),5);
     mo1.info.xgap(jGroup+1,1:2)=[min(n1) max(n1)];
     if ~isfield(RO,'KC'); continue;end
     if jGroup==1 %% Gradient
      r2=[RO.KC.Y; repmat(RO.KC.Y(end,:),nGroup-2*size(RO.KC.Y,1),1);flipud(RO.KC.Y)];
     end
     i1(:,3)=-3; 
     % Value given for K / C is total divide uniformly
     i1(:,7)=r2(jGroup,1)/size(i1,1); i1(:,9)=r2(jGroup,2)/size(i1,1);
     elt=[elt;i1];  %#ok<AGROW>
 end
 if ~isempty(elt)
   mo1=feutil('addelt',mo1,'celas',elt);
   if ~isempty(mo1.il);mo1.il(mo1.il(:,2)==12,:)=[];end
 else % uniform spring
  i1=feutil('findnode z==',mo1,min(mo1.Node(:,7)));
  i1(:,3)=-3; i1(:,5)=12; mo1.nmap{'Map:ProName'}(12)='UniSoil'; 
  mo1=feutil('addelt',mo1,'celas',i1);
 end
 
 %% #MeshList.Rail_end dampers (proid 13) -3
 if RO.Ns==1&&~isfield(RO,'CRailEdge') 
  %% use the periodic solution for a single slide
  mo1=fe_cyclic('build -1 600 0 0',mo1);
 else
  i1=feutil('findnode x== | x== & matid7',mo1,min(mo1.Node(:,5)),max(mo1.Node(:,5)));
  i1(:,3)=-3; i1(:,5)=13; mo1.nmap{'Map:ProName'}(13)='RailEnd'; 
  mo1=feutil('addelt',mo1,'celas',i1);
  if isfield(RO,'CRailEdge');
   mo1.il(mo1.il(:,1)==13,5)=RO.CRailEdge;
  end
  % boundary condition for multiple slices
  st=sprintf('x==%g|x==%g& proid 7 -DOF 1 2',min(mo1.Node(:,5)),max(mo1.Node(:,5)));
  mo1=fe_case(mo1,'fixdof','RailEdge',st);
 end

 %% Boundary conditions 
 mo1=fe_case(mo1,'fixdof','Ground',sprintf('z==%g -DOF 1 2',min(mo1.Node(:,7))));
 mo1=fe_case(mo1,'fixdof','SymY',sprintf('y==%g -DOF 2',max(mo1.Node(:,6))));
 
 %% how prepare output
 mo1=p_solid('default;',mo1);
 mo1.il(mo1.il(:,1)==3,:)=[];% Dont' duplicate contact
 if isfield(RO,'name');mo1.name=RO.name;end
  %% #MeshList.end -3
 if nargout==0;
    cf=feplot(2);cf=feplot(mo1);cf.sel='eltname~=celas';fecom('showfipro');
 else; out=mo1; 
 end
 
end % end MeshList


function RM=nameToMeshRO(name,evt);
%% #nameToMeshRO naming nomenclature to parameter 
%  'U30{5,Gc,tc-700_300,PadFuSN{io4}}:no:W2{XaZa}'
if name(1)=='{';name(1)=[];name(end)='';end
if nargin==1;evt=struct;end
r1=sdth.findobj('_sub:',name);st={r1.subs};
RM=struct; 
if isfield(evt,'projM')
 if isKey(evt.projM,'CMsh');RM=evt.projM('CMsh'); end
 RM.projM=evt.projM;
else; RM.projM=vhandle.nmap;
end
% Rail U30{3,Ga}
switch lower(st{1})
case 'wc2'; 
  [opts,model]=rail19('LoadForContact RailWCoarse2');
otherwise; 
  if isempty(st);error('Not found %s',name);end
  RM.rail=st{1};st(1)=[];st1=st{1};st(1)=[];
  RM.Ns=str2double(st1{1});
  i2=strncmpi(st1,'pad',3);if any(i2);RM.Pad=st1{i2};st1(i2)=[];end
  i2=strncmpi(st1,'G',1);if any(i2);RM.TGrad=st1{i2};end
  i2=strncmpi(st1,'t',1);if any(i2);RM.top=st1{i2};else;RM.top='';end
  i2=strncmpi(st1,'Lc',1);
  if any(i2);RM.Lc=comstr(st1{i2}(3:end),-1);else;RM.Lc=15;end
  if ~isempty(RM.top) % tc
   RM.topCoarse='fine';
   if strncmpi(RM.top,'tc',2); RM.topCoarse='coar';RM.top(2)='';end
   name=strrep(name,[',' RM.top],'');
   st2=regexp(RM.top,'t([^_]*)_(.*)','tokens');
   RM.top=[str2double(st2{1})];
  end
end
% Eclisse/Bracket Ec41{Air}
if ~isempty(st); i2=min(2,length(st)); if strncmpi(st{2},'w',1);i2=1;end
    RM.Nb=st(1:i2); st(1:i2)=[];
end
% Wheel Wa{XaZu}
if ~isempty(st);RM.Wheel=st(1:min(end,2)); end
RM.name=name;
r2=railu.Mesh('TrackLi',RM); 
RM.Track=r2.rail;RM.KC=r2.KC;
if isempty(evt);evt=struct;end
st=intersect(fieldnames(evt),fieldnames(RM));
for j1=1:length(st);RM.(st{j1})=evt.(st{j1});end
try
    cinM=d_rail('nmap.Map:Cin');%evt.projM('Map:Cin');
    r2=sdtm.pcin(cinM,'gr:All','asTabUo');
    st1=fieldnames(RM);st1(ismember(st1,r2.table(:,1)))=[];
    st1(sdtm.regContains(st1,'(gr.|ToolTip|projM|name'))=[];
    fprintf('%s not in sdtweb d_rail cinM, add ? \n',sdtm.toString(st1(:)'))
end

end % nameToMeshRO


function  [model,RO]=Case(varargin);
 %% #Case implement forwarding of sdtsys.stepmesh case handling
 [CAM,Cam]=comstr(varargin{1},1);model=varargin{2};RO=varargin{3};
 carg=4; if ~isfield(RO,'projM')&&isfield(RO,'nmap');RO.projM=RO.nmap;end
 for j1=1:length(RO.preCase)
  [st1,RC]=sdtm.urnPar(RO.preCase{j1},'{}{}');
  switch lower(st1)
  case 'rep' 
%% #Case.Rep Full model with slice repetitions
  [st1,RC]=sdtm.urnPar(RO.preCase{j1},'{}{ncell%g,nb_slices%g,fixEnd%g}');
  RC=sdth.sfield('rename;',RC,{'ncell','nb_slices'});
  %% Build track model in mt
 ms=model;
 RP=stack_get(ms,'info','MeshParam','g');RP=sdth.sfield('addmissing',RP,RC);
 mt=feutil(sprintf('repeatsel %i .6 0 0',RC.nb_slices),ms);
 mt=feutil('joinall',mt);
 mt.Elt=feutilb('SeparatebyProp',mt.Elt);

 if isKey(RO.projM,'NeedFix'); RO.fix=RO.projM('NeedFix');
    %mt=fe_case(model,RO.fix{:});
    mt=fe_case(mt,RO.fix{:});
 end
 if isfield(RP,'fixEnd')
  mt=fe_case(mt,'fixdof','ends',sprintf('x==-.3|x==%g',.6*nb_slices+.3));
 end
 mt=stack_set(mt,'info','MeshParam',RP);
 model=mt;


  case 'track' % combination of superelements
 dbstack; keyboard;
  case 'prrline' 
  %% #Case.prrline add rail surface reference line
  pro=feutil('getpro 301 -struct',mt);
  if isfield(pro,'I1') % beam
    n1=sortrows(feutil('getnode proid 301',mt),[6 5]);
    n2=n1(n1(:,6)<0,:);  % should allow curved rails
    if ~isempty(n2)
     el1=[n2(1:end-1,1) n2(2:end,1)];el1(:,3:4)=399;
     mt=feutil('addelt',mt,'bar1',el1);
    end
    n2=n1(n1(:,6)>0,:); 
    if ~isempty(n2)
     el1=[n2(1:end-1,1) n2(2:end,1)];el1(:,3:4)=399;
     mt=feutil('addelt',mt,'bar1',el1);
    end
    mt=feutil('setmat',mt,'dyn_mesh(matdb{399,prrLine})');
    mt=feutil('setpro',mt,'d_rail(nmap.ProDb.prrLine)');
  end
  model=mt;
  case 'vehicle' 
  %% #Case.vehicle add vehicle
  RV=RO.projM('Vehicle');
  mt=dyn_mesh('vehicleadd MovingLoad',mt,RV);
  model=mt;

 otherwise
  if exist(st1,'file')
   %% #Case.per #Case.cyc apply periodicity condition
   model=feval(st1,RC.Other{1},model,RC.Other{2:end});
  else % attempt to use fe_case
    model=fe_case(model,RO.preCase{j1,1:3});
  end
  end
 end

end


function out=getNmap(opt)
%% #getNMap : implementation of nmap (legacy compatibility)
% railu.getNmap('reset');
persistent gnmap
if nargin==0; opt=[];end
if isempty(gnmap)||isequal(opt,'reset')
  gnmap=vhandle.nmap; 
  % rails
  gnmap.append({ ...
   '60-E1', struct('hr1',11.5,'hr2',51,'hr3',172, ...
      'lar1',150,'lar2',16.5,'lar3',72,'ToolTip','Coarse DV Rail', ...
      'section','UIC60_Sections.mat#ra','cite','rail_IN10208.pdf#page=49') 
   'SN',struct('PadType','Volume','hs',9,'lsem',220,'k_pad',180e6)
   'UIC60',struct('alias','60-E1')
    '50-E6',struct('hr1',10.2,'hr2',49,'hr3',153, ...
      'lar1',140,'lar2',15.5,'lar3',65,'ToolTip','Coarse DV Rail')  
    '46-E2',struct('hr1',10.5,'hr2',47,'hr3',145, ...
       'lar1',134,'lar2',15,'lar3',62,'ToolTip','Coarse DV Rail');
      })

  %RO=struct;disp(comstr(RO,-30,struct('NoClip',2)))
  % Ward : @sdtdata/bug/ward/B70_RP_UIC60E1_mesh.bdf
  % M41_RP_UIC60E1_mesh.bdf
  % OEBB_RP_UIC60E1_mesh.bdf 
  % OEBB_RP_UIC60E1_mesh_fine.bdf
  gnmap.append({ ...
  'M450',struct('SleeperType','Monoblock','hb',230,'hsb',115,'Ltr',2411,'lab',290,'lsem',220); 
  'B450',struct('SleeperType','Biblock','hb',240,'hsb',120,'Ltr',2411,'lab',290,'lsem',220,'lob',840,'nent',3);
  'M240',struct('SleeperType','Monoblock','hb',187,'hsb',93.5,'Ltr',2243,'lab',300,'lsem',170);
  'Wooden',struct('SleeperType','Monoblock','hb',140,'hsb',70,'Ltr',2600,'lab',240, ...
      'lsem',240,'Mv201E',12000000000,'Mv201nu',0.45,'Mv201rho',800,'Mv201eta',0.01);
      })

  gnmap('Map:Arm')=vhandle.nmap( ...
      {'MainMono',struct('RailName','60-E1','SleeperName','M450') 
      'MainBi',struct('RailName','60-E1','SleeperName','B450') 
      'Secondary',struct('RailName','50-E6','SleeperName','M240') 
      'WoodenTrack',struct('RailName','50-E6','SleeperName','Wooden') 
      'SlabTrack',struct('RailName','60-E1','SleeperName','None') 
      'None',struct
      });
  % r1=d_cntc('nmap.Map:Oprop.prrEulSurf')

  propM=vhandle.nmap;gnmap('Map:OProp')=propM;
  propM.append({'prrLagSurf', struct('ToolTip','Rail surface in Larangian/global/iSys frame', ...
       'value',{{'EdgeColor',railu.color('EdgeGrey'),'EdgeAlpha',0.3,'FaceColor',railu.color('RailL'),'FaceAlpha',0.3}});
   'prrEulSurf',struct('ToolTip','Rail surface in Eulerian/tr frame', ...
       'value',{{'EdgeColor',railu.color('EdgeGrey'),'EdgeAlpha',0.3,'FaceColor',railu.color('RailE'),'FaceAlpha',0.7}});
   'prwEulSurf',struct('ToolTip','Wheel surface in Eulerian/w frame WheelW', ...
       'value',{{'EdgeColor','k','EdgeAlpha',0.5,'FaceColor',railu.color('WheelE'),'FaceAlpha',0.7}});
   'prwLine',struct('ToolTip','Ow(theta) line', ...
       'value',{{'Color',railu.color('WheelE'),'LineWidth',2}});
   'prrLine',struct('ToolTip','Or(x) line', ...
       'value',{{'Color',railu.color('RailE'),'LineWidth',2}});
   'prwLagSurf',struct('ToolTip','Wheel surface in Lagrangian/Wc frame', ...
       'value',{{'EdgeColor',railu.color('EdgeGrey'),'EdgeAlpha',0.3,'FaceColor',railu.color('WheelL'),'FaceAlpha',0.3}});
   'ViewYZ',struct('ToolTip','xxx -z vertical, y xxx', ...
      'value','view(90,0);iimouse(''viewh+180'')')
   'ViewXZ',struct('ToolTip','xxx -z vertical, y xxx', ...
      'value','view(0,0);iimouse(''views+180'')')
       })
  propM('line')=struct('ToolTip','dashed line','value',{{'Color',railu.color('EdgeGrey'),'linestyle','--'}});
  % propM('prrEulSurf')=struct('alias','prrEulSurf');
  % propM('prrLagSurf')=struct('alias','prrLagSurf');
  % propM('prrSlice')=struct('alias','prrLine');
  % propM('prwLagSurf')=struct('alias','WheelWc');
  % propM('prwEulSurf')=struct('alias','prwEulSurf');
  % propM('prwSlice')=struct('alias','prwLine');
  annotM=vhandle.nmap;gnmap('Map:Annot')=annotM;
  annotM.append({'prwCut',struct('ToolTip','Eulerian wheel profile slice', ...
      'from','prwEulSurf','tcL',1,'scL',1:10:300,'bas','w') 
    'prwLineFull',struct('ToolTip','Eulerian trajectory of OwL', 'from','prwLine','tcL',-180:2:180,'bas','w')
     });
  
  % xxx alias
  %railu.prop({'WheelW',''prwEulSurf'')
  %propM.append(d_cntc('nmap.Map:OProp'))
  if isequal(opt,'reset');opt=[];end


end
if ~isempty(opt); 
    out=gnmap(opt);
else
 out=gnmap;
end
end
function out=prop(key)
  %% prop : named properties
  propM=railu.getNmap('Map:OProp');
  if nargin==0; asTab(propM)
  else
   out=useOrDefault(propM,key,'','','getValue');
   if isfield(out,'alias');
       out=useOrDefault(propM,out.alias,'','','getValue');
   end
  end
end

function col=color(tag)
  %% #color tag based properties
  col={'WheelE','#C75204';'WheelL','#FF9E00';'RailE','#0077B6';'RailL','#00B4D8'
   'EdgeGrey','#808080'};
  if nargin==1;
      i1=find(strncmpi(col,tag,length(tag)),1,'first');
      if isempty(i1);i1=size(col,1);end
      col=col{i1,2};
  end

end


end % Static
end
