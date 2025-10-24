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
   load(FileName,'sections');
   if isa(sections,'containers.Map');sections=vhandle.map(sections);end
   if ~isKey(sections,tag); error('%s not found in %s',tag,FileName);end
   m_rail=useOrDefault(sections,tag,'','','getValue');
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
 if all(ismember(fieldnames(RO),{'ToolTip','data'}))
   RO=RO.data; 
 end
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
  propM=vhandle.nmap;
  gnmap('Map:OProp')=propM;
  propM.append({'RailGl', struct('ToolTip','Rail surface in Larangian/global/iSys frame', ...
       'value',{{'EdgeColor','k','EdgeAlpha',0.3,'FaceColor',railu.color('iSys'),'FaceAlpha',0.3}});
   'RailTR',struct('ToolTip','Rail surface in Eulerian/tr frame', ...
       'value',{{'EdgeColor','k','EdgeAlpha',0.3,'FaceColor',railu.color('r'),'FaceAlpha',0.5}});
   'prwEulSurf',struct('ToolTip','Wheel surface in Eulerian/w frame WheelW', ...
       'value',{{'EdgeColor','k','EdgeAlpha',0.3,'FaceColor',railu.color('w'),'FaceAlpha',0.5}});
   'WheelWc',struct('ToolTip','Wheel surface in Lagrangian/Wc frame', ...
       'value',{{'EdgeColor','k','EdgeAlpha',0.3,'FaceColor',railu.color('wc'),'FaceAlpha',0.5}});
   'ViewYZ',struct('ToolTip','xxx -z vertical, y xxx', ...
      'value','view(90,0);iimouse(''viewh+180'')')
       })

  propM('prrEulSurf')=struct('alias','RailTR');
  propM('prrLagSurf')=struct('alias','RailGl');
  propM('prwLagSurf')=struct('alias','WheelWc');
  propM('WheelW')=struct('alias','prwEulSurf');
  % xxx alias
  %railu.prop({'WheelW',''prwEulSurf'')
  propM.append(d_cntc('nmap.Map:OProp'))
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
  end
end

function col=color(tag)
  %% #color tag based properties
  col={'WheelE','#FF6D00';'w','#FF6D00';
       'wc','#FF9E00';'wc','#FF9E00';
       'isys','#00B4D8';'gl','#00B4D8';
       'tr','#0077B6';'r','#0077B6';
       'a','#023E8A'};
  if nargin==1;
      i1=find(strncmpi(col,tag,length(tag)),1,'first');
      if isempty(i1);i1=size(col,1);end
      col=col{i1,2};
  end

end

end % Static
end
