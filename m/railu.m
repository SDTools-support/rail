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
  if isfield(mv,'meta')&&isfield(mv.meta,'nomrad')
   mv.Node(:,7)=mv.Node(:,7)+Geo.z_rail_up+safeDistM(mv.meta.nomrad); %z/ vehicle put on the rails
   mv.bas=[];
  else
   mv.Node(:,7)=mv.Node(:,7)+max(mt.Node(n1,7))-...
     min(mv.Node(:,7)); %z/ vehicle put on the rails
  end

  % basis : add basis  1 99 range
  basid=max([mt.bas(mt.bas(:,1)<100,1);0])+1;
  n1=find(mt.Node(:,1)>1000); % mt without any vehicle
  n1=mt.Node(mt.Node(n1,7)==max(mt.Node(n1,7)), ...
      6); %Ys of nodes on the top of the slice (pad nodes if any)
  if ~RO.half;
      ybas=(min(n1)+max(n1))/2;
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
  mo1.bas=mt.bas(ismember(mt.bas(:,1),unique(mo1.Node(:,1))),:);
  % FixDof to have a 2D motion without x motion for the vehicle :
  [Case,vDOF]=fe_mknl('init',mo1); i1=fe_c(vDOF,.01*[1 2 4 6]','dof');
  mt = fe_case(mt,'FixDof','2D',i1);
  mo1=fe_case(mo1,'FixDof','2D',i1);

   r1=stack_get(mt,'SE','#s\d');
   if ~isfield(r1{1,3},'K') 
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
     error('Not implemented')
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
