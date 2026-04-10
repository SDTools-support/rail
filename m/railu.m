classdef railu
    
% railu : name space for methods related to RAIL
%
% <a href="matlab:methods railu">methods railu</a> for methods
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
%% #asUo : convert Dynavoie and rail inputs to UO format -3

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
function cinM=pcin(varargin)
 %% #pcin : add cinM nodes
 CAM='pcin';
 li={'key','ToolTip','DoOpt';
  'railu.Mesh.BeamMass','beam mass track slice model',[ ...
   'ncell(15#%g#"number of sleepers")' ...
   'half(#31#"1 one rail, 0 symmetric") '... %xxx_HP check
   'gaug(1.575#%g#"gaug") '...
   'Ltr(2.32#%g#"width") '...
   'kp(300e6#%g#"pad stiffness") '...
   'kb(20e7#%g#"ballast stiffness") '... % value of Arlaud 2016
   'kpr(1e6#%g#"pad torsion stiffness") '...
   'kps(1e6#%g#"pad transverse torsion stiffness") '...
   'rand(0#%g#"ballast dispersion") '...
   'randp(0#%g#"pad dispersion") '...
   'simpack(0#%g#"use simpack compatible elements") '...
   'Ftr(0#3#" renumber for simpack") '...
   'eta(0#31#"use eta") '...
   'cyc(0#31#"apply periodicity") '...
   'ms(130#%g#"sleeper mass")' ...
   'unit(SI#%s#"unit system")' ...
   'lc(.15#%g#"refine length")' ...
   'pk(75e3#%g#"pad stiffness")']  
   'dyn_mesh.arm','spleeper/rail/pad',{'key';'Slice.ArmP.ArmType';'Slice.ArmP.RailType'
         'Slice.ArmP.PadType';'Slice.ArmP.SleeperType'
   }
   };
  sdtm.pcin(['prero',comstr(CAM,5)],li);

 r1=sdtu.f.ppath;carg=1;
 if nargin>0&&isa(varargin{1},'vhandle.nmap');cinM=varargin{1};carg=2;
 else; cinM=sdtm.pcin;
 end
 if ~isKey(cinM,'prero.railu.Range')
    DefBut=sdtm.rmfield(feval(dyn_ui('@genDefBut')),'PTree','Tab');
    sdtm.pcin(cinM,'append',DefBut);
 end
 if nargout>0; out=cinM;else; sdtm.pedit('{disp}',li);end
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
%% #WheelSection : wheel section mesh -3
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




function   m_rail=RailSection(RO);
%% #RailSection : rail section mesh -3

RA=struct;m_rail=[];
if ischar(RO); 
  nmap=d_rail('nmap');
  if isKey(nmap,RO);RA=nmap(RO);
  elseif contains(RO,'.mat#');
   % railu.RailSection('U30_Sections.mat#ra')
   fname=RO(1:find(RO==35,1,'first')-1);tag=RO(length(fname)+2:end);
   FileName=d_rail('wd',fname);
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
if isfield(RO,'Arm')
     RA=RO.Arm;
end
if isfield(m_rail,'Elt')
elseif isfield(RA,'lar1')
     % dynavoie coarse rail model
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
    if isfield(RA,'hs')
     m_rail.Node(:,7)=m_rail.Node(:,7)+RA.hs;
    end
     
    m_rail=feutil('Divide 2 2',m_rail);
    
    if exist('dyn_mesh','file')
     m_rail.pl=dyn_mesh('matdb 301',{{'Rail','Standard'}});
    else;
     m_rail.pl=d_rail('nmap.MatDb.rail');m_rail.pl(1)=301;
    end
    m_rail.il = [301  fe_mat('p_solid','SI',1) 0 10002 0 1]; % new integration rule
else; error('Not implemented')
end

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
n1=feutil('findnode z==',m_rail,0);
m_rail.meta.TopNodes=n1; 

end % RailSection

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
[RO,st,CAM]=cingui('paramedit -DoClean',DoOpt,{RO,CAM});Cam=lower(CAM);
RO.isFullModel=1;
 if RO.cyc; RO.ncell=1;
 elseif length(RO.ncell)<2;RO.ncell(2)=0;
 end
 if ~isfield(RO,'xsens'); RO.xsens=[];end


%      ms=feutil('setmat',ms,'dyn_mesh(matdb{399,prrLine})');
% feutil('setpro',ms,'d_rail(nmap.ProDb.prrLine)');
r2=d_rail('nmap');
if isfield(RO,'nmap')&&isa(RO.nmap,'vhandle.nmap');RO.projM=RO.nmap; 
 if ~isKey(RO.projM,'MatDb');RO.projM('MatDb')=r2('MatDb');end
 if ~isKey(RO.projM,'ProDb');RO.projM('ProDb')=r2('ProDb');end
elseif ~isfield(RO,'projM');RO.projM=r2;
end

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
 RO.PreIl={'name','ProId';'UIC60beam',301;'SleeperA',201};
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
  RO.fix={'fixdof','2DRail','matid 301 -DOF 1 2 4 6', ...
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
  RO.projM('NeedFix')={'fixdof','2DRail','matid 301 -DOF 1 2 4 6', ...
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


 else; error('Need implement')
 end
end
function  [model,RO]=Case(varargin);
 %% #Case implement forwarding of sdtsys.stepmesh case handling
 [CAM,Cam]=comstr(varargin{1},1);model=varargin{2};RO=varargin{3};
 carg=4; if ~isfield(RO,'projM')&&isfield(RO,'nmap');RO.projM=RO.nmap;end
 for j1=1:length(RO.Case)
  [st1,RC]=sdtm.urnPar(RO.Case{j1},'{}{}');
  switch lower(st1)
  case 'rep' 
%% #Case.Rep Full model with slice repetitions
  [st1,RC]=sdtm.urnPar(RO.Case{j1},'{}{ncell%g,nb_slices%g,fixEnd%g}');
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

 case {'cyc','per'}
%% #Case.per #Case.cyc apply periodicity condition
  model=fe_cyclic('build -1 .6 0 0;',model);

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
          keyboard; 
  end
 end

end

function   ms=MeshSlice(RO);
%% #MeshSlice : rail section mesh -3
 RO=sdtu.ui.cleanEntry(RO);
 if iscell(RO) 
  %% dynavoie list
  Range=dyn_solve('Range',struct('SliceCfg',{{'T',RO}}));
  dyn_solve('RangeLoopReset',Range,struct('ifFail','error'));  
  if sdtdef('isinteractive')
   PA=dyn_ui('paramvh'); cf=feplot(PA.ms);fecom(cf,'showfimat');
  end
 else % rail mesh 
  [~,R1]=sdtm.urnPar(RO,'{}{Sleeper%s,Rail%s,Pad%s,Sub%s}');
  error('Need extend implement slice meshing')
 end

end


function out=getNmap(opt)
%% #getNMap : implementation of nmap (legacy compatibility)
persistent gnmap
if nargin==0; opt=[];end
if isempty(gnmap)||isequal(opt,'reset')
  gnmap=vhandle.nmap; 
  % rails
  gnmap.append({ ...
   '60-E1', struct('hr1',11.5,'hr2',51,'hr3',172, ...
      'lar1',150,'lar2',16.5,'lar3',72,'ToolTip','Coarse DV Rail') 
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
