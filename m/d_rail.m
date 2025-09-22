function [out,out1]=d_rail(varargin)

% d_rail: meshing utilities for rail applications
%
% <a href="matlab:d_rail('nmap')">List configurations (nmap)</a>
% <a href="matlab:d_rail('MeshDb')">Mesh database</a>

%  Etienne Balmes, Guillaume Vermot des Roches
% Copyright SDTools & SNCF I&R, 2019-2025

%#ok<*ASGLU,*NASGU,*NOSEM,*NBRAK>
obj=[];evt=[];
if nargin==0; help d_rail
    return; 
else; [CAM,Cam]=comstr(varargin{1},1);carg=2;
end

if comstr(Cam,'script');[CAM,Cam]=comstr(CAM,7);
%% 

if comstr(Cam,'dv1')
 %% #TutoDvbeam: d_rail('ScriptDv1'); 
%% Step1: initialize interface
wd=fullfile(sdtdef('tempdir'),'TutoBallast'); % change the working directory if necessary
dyn_ui('SetWD',wd);ofact('mklserv_utils -silent');

%% #TodoEb : debug

% Step2 : mesh 
RA=struct('ArmType','MainMono'); % standard monoblock parameters
RS=struct('SubType','Ballast'); % just ballast
RR=struct('RedType','Modal'); % default reduction parameters
SliceCfg={'Slice',struct('Arm',RA,'Sub',RS,'half',1);
          'Reduce',RR}; % format the SliceCfg for a half track
      
% define track parameters and place sensors
Sens={'s1(10:11):Bal:Sleeper:z'}; % define sensors
TrackCfg=struct('type','GenTrack','nb_slices',40,'Sens',{{Sens{:}}});

% Define the simulation to be done on the track model
SimuCfg={'Impact',struct('xload',7.2,'weight',10e-3);... % 10 kg at x=7.2m
         'Time',struct('tf',.8)}; % impact simulation
Range=struct('SliceCfg',{{'S1',SliceCfg}}, ...
             'TrackCfg',{{'T1',TrackCfg}}, ...
             'SimuCfg',{{'S1',SimuCfg}}); % format the range
dyn_solve('RangeLoop -reset',Range); % build the model and store in the interface
dyn_post('PostRecept'); % calculate transfer at sensors and plot

%%

l1={'Vehicle',struct('VehType','Axle','v0',80,'x0','$x_start-2','xf','$x_start-1');... % define vehicle speed
   'Time',struct('NeedV',1,'NeedA',1,'dt',1e-4)}; % time integration par.
Range=struct('SimuCfg',{{'S1',l1}});
dyn_solve('RangeLoop',Range);

dyn_post('ViewRestAll -light'); % display deformation animation


%% EndTuto

else;
        sdtweb('_link','sdtweb(''t_exp19'',''ref'')');
end

elseif comstr(Cam,'solve');[CAM,Cam]=comstr(CAM,6);


elseif comstr(Cam,'load');[CAM,Cam]=comstr(CAM,5);
%% #Load 

elseif comstr(Cam,'case');[CAM,Cam]=comstr(CAM,5);
%% #Case : DOE commands to define case (boundary cond, ...)
model=varargin{carg};carg=carg+1;
RO=varargin{carg};carg=carg+1; if isempty(CAM);[CAM,Cam]=comstr(RO.Case,1);end

if comstr(Cam,'21ref')
 %% #Case.21ref : 22 port of reference that was in  sdtweb abr21 v21testab 

 model=feutil('setmat',model,rail19('nmap.MatDb.PadAniso'));
 nmap=model.nmap; expM=RO.nmap; 
 RT=struct('V',[], ...
  'v',linspace(2,20,5829),'dt',1e-5, ...
 'Tdamp',{{2e-3 [0 0;3000 .05/5]; ...
    6000e-3 [0 0;3000 .05]}}, ...  % 6s is actually no change
 'TdampDown',[0 0;3000/.05*.1 1], ...% higher damping first
 'xstart',-300,'xend',2e4,...%'maxNout',200, ... %20/(80e-3)-300
 'CtcForm',13, ...
 'SaveFcn',@ChangeRayleigh, ... % where parameter shifts are implemented
 'profile',0,'stra',1);
 if ~isfield(model.info,'Vel');model.info.Vel=20;end
 RT.V=struct('X',{{[-1 .02 .047 1]'}},'Y',[2 2 model.info.Vel*[1 1]]'*[1 0]); %tdown and accel end
 if ~isKey(expM,'Ctc');expM('Ctc')=rail19('nmap.Ctc21'); end
 RT=sdth.sfield('addmissing',RT,expM('Ctc'));
 %RT.delta=20;RT.Fl=1e6; RT.Kc=-2.3e9;
 if isfield(model.info,'top')
    RT.xstart=model.info.top(1)+50;
    RT.xend=min(RT.xend,model.info.top(2)-30);
 end
'xxx tdamp2 3e3 .01'

 % Contact surface handled in wheel 
 if any(model.pl(:,1)==104)
  model=feutil('setmat 104 rho=0 nu=0',model); % Interrail space
 end
 %if ~isfield(model.info,'xgap');RT.xgap=[-7.55 7.55];end
 RT.FirstDown=-.01;
 %RT2.OutInd='cg=feplot('';'');opt.OutInd=fe_c(Case.DOF,cg.sel.Node,''ind'');';
 %RT2.OutputFcn='tout=t(1:10:end);';
 model.info=sdth.sfield('addmissing',RT,model.info); 
 out=model; 

elseif comstr(Cam,'pada')
 %% #Case.pada : enforced pad compression
 model=fe_case(model,'fixdof','base',sprintf('z==%.15g',min(model.Node(:,7))));
 n2=feutil('findnode z==',model,max(model.Node(:,7)));
 d2=struct('def',[zeros(size(n2,1)*2,1);ones(size(n2,1),1)], ...
     'DOF',[n2+.01;n2+.02;n2+.03],'curve',{{'Input'}},'KeepDof',1);
 model=fe_case(model,'DofSet','Top',d2);
 if strncmpi(model.name,'pada2',5)
   %% 
   model=fe_case(model,'fixdof','xsym','x==0 &z>0&z<.09 -DOF 1');
   model=fe_case(model,'fixdof','ysym','y==0 &z>0&z<.09 -DOF 2');
 end
 out=model;
elseif comstr(Cam,'dofsetsine')
 %% #CaseDofSetSine : port of Rafael's dofset see sdtweb dfr_ident case
 % Should be moved to sdtweb d_fetime Case.TopZ

 model.Node(:,5:7)=fix(model.Node(:,5:7)*1e5)/1e5;
   
   if ~isfield(RO,'Top'); RO.Top=sprintf('z==%.15g',max(model.Node(:,7)));end
   n1=feutil(['getnode' RO.Top],model);
   try;n1=sortrows(n1,[5 6 7],'descend');catch;n1=flipud(sortrows(n1,[5 6 7]));end
   tcell={'lab',    'SensId',  'X','Y','Z',      'DirSpec','unit';
       'Corner:z',n1(1)+.03,n1(1,5),n1(1,6),n1(1,7),'z','*1=mm'; % 
       'Top:RFz', n1(1)+.40,n1(1,5),n1(1,6),n1(1,7), ... % Driving resultant 
        'Fz{Resultant,"matid1 -fieldOut4","inelt{selface & facing>.8 0 0 1000}","[0 0 1]"}','*1e-0=N'; 
       };
   %     'resultant.{"matid1 -fieldOut4","inelt{selface & facing>.8 0 0 1000}","[0 0 1]"}','*1e-0=N'; 
   if ~isfield(RO,'Case');RO.Case='DofSet.Sine';end
   if contains(lower(RO.Case),'dofset')  
     %% #LoadCase.DofSet : enforced displacement using a given shape -3
     if strncmpi(RO.name,'unis',4) % traction
      data=struct('DOF',feutil('getdof',(1:3)'/100,n1(:,1)), ...
          'def',[zeros(2*size(n1,1),1);ones(size(n1,1),1)]);         
     else, error('define load'); 
     end
     data.type='trans';data.KeepDof=1;
     data.NLdata=struct('type','nl_inout','MexCb',{{d_fetime('@MBBryan')}});
     data.curve={'Curve'};
     model=fe_case(model,'DofSet','Input',data);
     r2=fe_case(model,'getdata','SimplePeriodicity');
     if ~isempty(r2) % Cant use MPC on DofSet
      i2=fe_c(r2.DOF(r2.slave),data.DOF,'ind');
      r2.slave(i2)=[]; r2.c(i2,:)=[];
      model=fe_case(model,'mpc','SimplePeriodicity',r2);
     end  
     if any(strcmpi(RO.CaseVal,'st3')) % StoreType3 
      i1=strcmpi(model.Stack(:,1),'pro');
      model.Stack{i1,3}.NLdata.StoreType=3; 
     end
   else, error('case %s',RO.Case);
   end
   %%
   % now define sensors for output
   model=fe_case(model,'SensDof -InFEM','Sensors',tcell);
   model=fe_case(model,'SensMatch radius 10','Sensors','groupall');
   NL=struct('type','nl_inout','slab',{{'Sensors:.*'}},'StoreType',1);
   model=stack_set(model,'pro','Observe',struct('il',999,'type','p_null', ...
       'NLdata',NL));  

   %% #LoadCase.curve Now add the time loading definition -3
   if contains(lower(RO.Case),'sine')||contains(lower(RO.Case),'tri')
     [~,r3]=sdtm.urnPar(['{' sprintf('%s,',RO.CaseVal{:}) '}'],'{}{A%ug,Nper%g,freq%ug}');
     st=sprintf('testeval %f*sin(2*pi*%f*t)',r3.A(1),r3.freq(1));%Strain % not used
   end
   model=fe_curve(model,'set','Curve',st); %define loadramp
   RO.nmap('CurModel')=model;out=[];
elseif comstr(Cam,'period')
  %% #CasePeriod : 2023 computation of periodic response
  model=fe_case(model,'remove','RailEdge');
  model.Elt=feutil('removeElt eltnamecelas & innode{x==300|x==-300}',model);
  model=fe_cyclic('build -1 600 0 0',model);
  RO.nmap('CurModel')=model;out=[];

elseif comstr(Cam,'anatraj')
  %% #CaseAnaTraj : analytical trajectories for validation
  [~,RO]=sdtm.urnPar(CAM,'{ty%s}:{t%g,x%g,z%g}');
  %cf=feplot;mo1=cf.mdl.GetData;

   % sdtweb rail19 nmap.mat.set
  %cf.sel='matid 7  11 104 4 5'
  i1=feutil('findnode matid 7 11 104 4 5',model);
  d1=feutilb('geomrb',model,[0 0 0],feutil('getdof',i1,(1:3)'/100));
  d1=feutilb('placeindof',model.DOF,d1);
  if isfield(RO,'z')
    d1.def=d1.def(:,3)*RO.z(:)'; d1.data=RO.t(:);
  end
  out=d1;
elseif comstr(Cam,'none');out=model;
else; error('Case%s',CAM);
end

elseif comstr(Cam,'mesh');[CAM,Cam]=comstr(CAM,5);
%% #Mesh : mesh loading/generation 
        
 
if comstr(Cam,'db');
 %% #MeshDb : database configurations

 PA=sdtroot('PARAMVh'); 
 if nargin==2&&isa(varargin{2},'vhandle.nmap');PA.nmap=varargin{2}; end
 if ~isfield(PA,'nmap');PA.nmap=vhandle.nmap;end;projM=PA.nmap;
 dbM=useOrDefault(projM,'Map:Sections',[],'d_rail');projM('Map:Sections')=dbM;
 sliceM=useOrDefault(projM,'Map:Slice',[],'d_rail');projM('Map:Slice')=sliceM;
 if ~isKey(dbM,'U30.ra')
   %% rail sections init
   f2=d_rail('wd','U30_Sections.mat'); 
   r1=load(f2);
   if isfield(r1,'sections')
    r2=[reshape(r1.sections.keys,[],1) reshape(r1.sections.values,[],1)];
    if strcmpi(r2{1},'Wheel');dbM('Wheel')=r2{1,2};r2(1,:)=[];end
    if ~strncmpi(r2{1},'u30',3)
      r2(:,1)=cellfun(@(x)['U30.' x],r2(:,1),'uni',0);
    end
    dbM.append(r2)
   end
 end
 if ~isKey(dbM,'WheelCut') % Initialize DB
   d_rail('MeshWheelCut',dbM);
 end

 if nargout>0; 
  st=CAM(4:end);
  if ~isempty(st)
     out=useOrDefault(dbM,st,[],sliceM,'getvalue');
  else
     out=projM;
  end
 else; % Display database contents 
   [~,RM]=sdtm.urnPar(CAM,'{}{}');
   if isfield(RM,'Other')&&any(sdtm.regContains(RM.Other,'imgen','i'))
     %% generate images % d_rail meshdb{ImGen}
     li=[dbM.keys' dbM.values'];
     sdtroot('SetProject',struct('PlotWd',fullfile(sdtdef('tempdir'),'sdtdemos/rail'), ...
         'Report','RailImages.docx'))
     
     for j1=1:size(li,1)
      r1=li{j1,2};
      if isfield(r1,'Elt')
       mo1=r1; mo1.name=li{j1,1};
       cf=feplot(mo1);iimouse('resetview')
       sdth.os(cf,'d.',{'FiPro','LgMdlName'},'p.',{'FnM'});
       feplot(cf);r2=comgui('imwrite',cf);
       [~,fname,ext]=fileparts(r2.Children{1}.FileName);
       li{j1,3}=sprintf('<img src="./_images/%s" alt="%s" width="200"/>',[fname ext],fname);
      end
     end
     fprintf('%s\n',li{:,3});
   else
     asTab(vhandle.tab(dbM));
   end
 end

elseif comstr(Cam,'pad')
%% #MeshPad
if comstr(Cam,'pada2')
 %% #MeshPadA2 from STL geom quater mesh for faster validation -3
 mo1=fe_gmsh(['read' d_rail('wd','vv16191_2D3.msh')]);
 mo1.Elt=feutil('selelt eltname==tria|eltname==quad',mo1);
 mo1.Node(:,5)=0;mo1.Node(:,5:7)=mo1.Node(:,[5 7 6]);
 mo1=feutil('extrude 3 .03 0 0',mo1);
 mo1=feutil('lin2quad',mo1); mo1.name='PadA2';
 mo1.Elt=feutil('set groupall matid 1 proid 1',mo1);
 mo1.pl=m_elastic('dbval steel');mo1=p_solid('default;',mo1);
 mo1.Node(:,2:3)=0;
 out=mo1; 
elseif comstr(Cam,'pada')
 %% #MeshPadA from STL geom -3
 mo1=fe_gmsh(['read' d_rail('wd','vv16191_2D.msh')]);
 mo1.Elt=feutil('selelt eltname==tria|eltname==quad',mo1);
 mo1.Node(:,5)=0;mo1.Node(:,5:7)=mo1.Node(:,[5 7 6]);
 mo1=feutil('extrude 18 .01 0 0',mo1);
 mo1=feutil('lin2quad',mo1);mo1.name='PadA';
 mo1.Elt=feutil('set groupall matid 1 proid 1',mo1);
 mo1.pl=m_elastic('dbval steel');mo1=p_solid('default;',mo1);
 mo1.Node(:,2:3)=0;
 out=mo1; 
elseif comstr(Cam,'pads')
 %% #MeshPadS very simple coarse volume -3
 [~,RM]=sdtm.urnPar(CAM,'{}{pad%s,lc%s}');
 if ~isfield(RM,'lc');RO.lc='20';end
 mo1=feutil('objecthexa 9 9',[0 0 0;eye(3)],feutil(['refineline' RM.lc],[0 140]), ...
     feutil(['refineline' RM.lc],[0 181]),[0 9]);
 mo1.pl=fe_mat('convertSIMM',rail19('nmap.MatDb.Pad'));
 mo1=p_solid('default',mo1); mo1.unit='MM';
 [a,b]=moe23(['Pad' RM.pad]);

 NL=struct('type','nl_inout','lab','uMaxw','MatTyp',{{1}},'keepLin',1, ...
   'StoreType',3, ... % Just along z
  'adofi',[],'MexCb',{{nlutil('@uMaxw')}}, ...
  'Fu',b);
mo1=feutil('setpro 9 in-3',mo1,'NLdata',NL);
out=mo1;

else
 %  sdtweb cbi20b SncfMeshRahssta
 model=varargin{carg};carg=carg+1;
 mo2=model;mo2.Elt=feutil('selelt proid 9 & seledgeall',model);
 i1=mo2.Elt(2:end,1:2);
 NNode=sparse(mo2.Node(:,1),1,1:size(mo2.Node,1));i2=full(NNode(i1));
 i1(abs(mo2.Node(i2(:,1),7)-mo2.Node(i2(:,2),7))<1,:)=[];i2=full(NNode(i1));
 r1=reshape(mo2.Node(i2,7),[],2);i3=r1(:,2)>r1(:,1); i2(i3,:)=i2(i3,[2 1]);
 r1=reshape(mo2.Node(i2,7),[],2);
 %cf=feplot(model);cf.sel='matid9';fecom('shownodemark',model.Node(i2(:,1),1));
 model.Node(i2(:,2),7)=model.Node(i2(:,1),7)-9; % pad thickness ?
 out=model;
end
 %% #MeshWheel : commands -2

elseif comstr(Cam,'wheelcut')
 %% #MeshWheelCut : append 2D section in x-z plane -3
   dbM=varargin{carg};carg=carg+1;
   % get a ray line on AXE
   model=dbM('Wheel');mo1=model; 
   mo1.Elt=feutil('selelt proid1',mo1); mo1.Elt=feutil('selelt selface',mo1);
   mo1.Node=feutil('getnodegroupall',mo1);
   mo1.Elt=feutil('selelt innode{y==max(y)}',mo1);
   mo1.Node=feutil('getnodegroupall',mo1); mo1.Elt=feutil('selelt seledge',mo1);
   mo1=feutil('joinall',mo1);
   mo1.Elt=feutil('divideingroups',mo1);   mo1.Elt=feutil('selelt group1',mo1);
   mo1.Node=feutil('getnodegroupall',mo1);
   r1=feutilb('geofindcircle',mo1,struct('nodes',mo1.Node(:,1)')); r1=r1{1};
   n1=mo1.Node(1,:);
      
   % ok now gen a node set with angular threshold and then exact line position
   n2=feutil('getnode proid1 8',model); % ROUE
   r2=n2(:,5:7)-reshape(r1.Origin,1,3);
   r0=n1(1,5:7)-r1.Origin;
   bas=basis(r1.normal,r0);
   r4=(r2-r0)*bas;   
   
   n3=n2(r4(:,3)>0&r4(:,2)>=0,1);
   n4=n2(abs(r4(:,3))<1e-3,1);
   
   % do selection to recover profile section elements
   mo1=model; %mo1.Node=feutil('getnodegroupall',mo1);
   mo1=feutil('addsetnodeid',mo1,'sectors',n3);
   mo1=feutil('addsetnodeid',mo1,'profile',n4);
   mo1.Elt=feutil('selelt proid1 8 & withnode{setname sectors} & selface & innode{setname profile}',mo1);
   mo1.Node=feutil('getnodegroupall',mo1);
   if ~isfield(r1,'profile')
     load(d_rail('wd','ProfWheel_U30_24.mat'),'n1');
     n1(:,3)=0;n1=[n1(:,2) n1(:,3) -n1(:,1)];figure(1);plot(n1(:,1),n1(:,3))
     r1.profile=n1;% [r theta z] 
   end
   mo1.info=r1;
   dbM('WheelCut')=mo1; 
   
elseif comstr(Cam,'wheel')
 %% #MeshWheel_DO -3
 RO=varargin{carg};carg=carg+1;
 mo1=useOrDefault(RO.projM,RO.wfrom,[],{'Map:Sections'});
 r1=mo1.info; % Origin normal 
 if ~isKey(RO.projM,'CbStickWheel')
  % CbStickWeel is used after refinement 
  RS=struct('sel',[ 'inElt{ProNameWheel & selface&innode{z<60}&facing >.9 0 0 -5000}'],'distFcn',[]);
  RS.distFcn=[];%lsutil('gen',[],{struct('shape','sphere','rc',l*5,'xc',l*1.5,'yc',l*.5,'zc',(5+1)*l)});
  RO.projM('CbStickWheel')={@lsutil,'SurfStick','$projM',RS};
 end
 % ProNameRail&selface
 RO.Ctc=struct('slave','inNode{proid2}', ...
   'master','proid3','ProId',3, ...
   'InitUnl0','rail19@InitCqWheel', ...% 'd_contact@surfStick'
   'StoreType',3,'MAP','tan1dx', ...
   'Fu','Kc 1e+15 stepx 0', 'Integ',[2], 'subtype',[2],'distFcn',[]);
 RO.Ctc.slave=['ProNameRail&selface&facing >.9 0 0 1e4&innode {distfcn"{box{' ...
     sprintf('%g %g %g,%g %g %g',[mean(RO.top) 7.7312   84.505],[diff(RO.top)/2 25 5]) ...
     ',1 0 0,0 1 0,0 0 1}}"}'];
 % xxx mo2=feutil(sprintf('setpro 3 Integ=%i Match -2',RW.integ),mo2);

 if length(RO.wref)>6; [u1,u2,RO.rLen]=comstr('coarse',[-25 2],RO.wref,lower(RO.wref));
 else; RO.rLen=[];
 end
%  cinM=RO.projM('Map:Cin');cinM('gr:Wheel')={'wfrom','wref','rLen','Lc','ctcSurf'}
   
 %% great now regen with base angular step   
  %    mo2=feutil(sprintf('rev %i o %.15g %.15g %.15g 360 %.15g %.15g %.15g',...
  %     length(r1.line),r1.Origin,r1.normal),mo1);
   % starts at horizontal
   if comstr(RO.wref,'coarse1')
    %% #MeshWeelCoarse1 Very coarse with few contact element -3
    if isfield(RO,'Lc')&&RO.Lc(1)>50
      r5=unique([setdiff(0:30:360,90) 89:92])/360;
    else
      r5=unique([setdiff(0:10:360,90) 89:.5:92])/360;
    end
    mo1.Node(:,7)=r1.Origin(3);
    mo2=feutil(sprintf('rev 0 o %.15g %.15g %.15g 360 %.15g %.15g %.15g',...
     r1.Origin,r1.normal),mo1,r5);

    elt=feutil('selelt selface',mo2);
    n1=feutil('getcg',mo2.Node,elt);
    i1=n1(:,1)>r1.Origin(1)-20& n1(:,1)<r1.Origin(1)+30&n1(:,2)>-10&n1(:,2)<10&n1(:,3)<100;
    %fecom('shownodemark',elt(i1,1:4))

    %[~,i1]=min(abs(n1(:,1)-r1.Origin(1))+abs(n1(:,2))+abs(n1(:,3)-85.02));
    mo2=feutil('addelt',mo2,'quad4',[elt(i1,1:4) ones(nnz(i1),1)*[3 3]]);
    %feplot(mo2);fecom('showfipro');
    mo2.Stack={};
    mo2.pl(~ismember(mo2.pl(:,1),[1 3 8]),:)=[];
    mo2.il(~ismember(mo2.il(:,1),[1 3 8]),:)=[];
    RO.projM('CbCtcGen')={@ctc_utils,'generatecontactpair','$projM',RO.Ctc};

   RD=struct('shape','Poly','Type','RevLine','isClosed',1, ...
     'Node',[(1:size(r1.profile,1))'*[1 0 0 0],r1.profile(:,[2 3 1])+[r1.Origin(1) 0 r1.Origin(3)]],'orig',r1.Origin,'axis',r1.normal);
   RD.orig(2)=0;
   li=feval(lsutil('@dToPoly'),'init',RD);
   if isfield(RO,'Debug')&&sdtm.Contains(RO.Debug,'ShowProfWheel')
    cf=feplot(mo2);cf.sel='innode{x>-180}';
    cf=feplot;lsutil('viewls',cf,li);fecom showCbTR; set(cf.ga,'clim',[-10 10])
    fecom('shownodemark',RD.Node(:,5:7),'marker','o')
   end
   % avoid wheel rotation 
   n1=sortrows(feutil('getnode z== & y==',mo2,mo2.info.Origin(3),mo2.info.Origin(2)),5);
   mo2=fe_case(mo2,'mpc','FixRot',struct('c',[1 -1],'DOF',n1([1 end],1)+.03,'slave',2));
   mo2.info.distFcn=li.distFcn;mo2.info=sdtm.rmfield(mo2.info,'line','profile','radius');
    out=mo2;
    return;
   elseif strncmpi(RO.wref,'fine',4) 
    %% #MeshWheelFine fairly large number of elements -3
    r5=[linspace(0,60,65*60/360) linspace(60,120,2.7*65*60/360) linspace(120,360,65*280/360)]/360;
    mo2=feutil(sprintf('rev 0 o %.15g %.15g %.15g 360 %.15g %.15g %.15g',...
     r1.Origin,r1.normal),mo1,unique(r5));
    y1=[-33.75 39.5];y1=[-10 10];
    x1=r1.Origin(1)+[-20 20];
    st1=sprintf('selface & withnode{y>%.15g & y<%.15g & x>%.15g&x<%.15g&z<2*min(z)} ',y1,x1);
    %st1=sprintf('withnode{y>%.15g & y<%.15g & x>%.15g&x<%.15g&z<124} ',y1,x1);
    mo2=fe_shapeoptim('refinehexasurf',mo2,st1);
    elt=feutil(['selelt' st1],mo2);
    mo2=feutil('addelt',mo2,feutil('set groupall matid 3 proid 3',elt));
    'xxx check radius liRad'
    out=mo2;
   else
   % #MeshWheelCoarse -3
   % @GV : detect sector and verify that generation by repetition will be correct
   % @GV : O:\sdtdata\dynavoie\rail19\profil_roue_castem.mat 
   %       subfunction that sets the radius of any point on the exact profile

 % n1=feutil('getNode inelt{setname C_RAIL}',model);y1=[min(n1(:,6)) max(n1(:,6))];
   % cannot coarsen profile easily
%    cf.sel=sprintf('withnode{y>%.15g & y<%.15g & x<%.15g}',y1,-1.05*max(n1(:,5)))
%    mo3=mo1;
%    [mo1.Elt,mo3.Elt]=feutil(sprintf('removeelt withnode{y>%.15g & y<%.15g & x<%.15g}',y1,-1.05*max(n1(:,5))),mo1);
%    mo3.Elt=feutil('selelt seledge',mo3)
%    mo3.Node=feutil('getnodegroupall',mo3);
    %% Obsolete version replace wheel
    %cf.sel=sprintf('selface & withnode{y>%.15g & y<%.15g & z<2*min(z)} ',y1)
    for j1=1:RO.rLen-1
     st1=sprintf('selface & withnode{y>%.15g & y<%.15g & z<2*min(z)} ',y1);
     mo2=fe_shapeoptim('refinehexasurf',mo2,st1);
    end
    mo1=abaqus('resolveset',model); % command that generates explicit sets
    %  mo3=feutilb('submodel',mo1,'proid1 8'); % removed
    %mo1=feutilb('submodel',mo1,'withoutnode{proid1 8}'); % kept
   
    % only recast sets ROUE AXE, C_ROUE
    mo2.Stack={};   mo2=feutil('removematpro-unused',mo2);
    mo2=feutil('addset eltid',mo2,'AXE','proid1');
    mo2=feutil('addset eltid',mo2,'ROUE','proid8');
   
    % regen contact set for wheel
    n1=feutil('getNode inelt{setname C_RAIL}',model);
    y1=[min(n1(:,6)) max(n1(:,6))];
  
    n2=feutil(sprintf('getnode inelt{selface & withnode{y>=%.15g & y<=%.15g}}',y1),mo2);
    n2=n2(n2(:,7)<min(n2(:,7))*1.1,:);
    [u1,mo2.Elt]=feutil('eltidfix;',mo2);
    mo2=feutil('addsetfaceid',mo2,'C_ROUE',sprintf('selface & withnode{nodeid %s}',num2str(n2(:,1)')));
   
    % add to base model
    mo1=feutilb('combinemodel-compatmatpro',mo1,mo2);
    mo1=stack_set(mo1,'info','wheel_geom',r1);
   
    %% now map profile for accurate gap
    f1=sdtu.f.safe('@sdtdata/rail19/INP/profil_roue_castem.mat');
    mo1=feval(t_exp19('@liRad'),'init',mo1,f1);  
    out=mo1;

   end
      
   
        
 
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

%% #MeshDb.add_sleeper before and after -3
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
%% #Mesh.TGrad : creation of gradient if not given in KC -3
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
out=struct('KC',RM.KC,'rail',{li}); % sdtm.toString(RM.KC)

  %% #MeshList.wc4 parameters see sdtweb rail19 MeshRailDb -3
  % sdtw('_ewt','should be done MeshRailDb')
  % RM=d_rail('MeshRailDb',struct('Nb',41, ... % Eclisse 41/42/61/62
  %      'Ns',3,...     % N sleeper before eclisse/brackets
  %      'Ngrad',1,...  % xxx
  %      'C',1e-1,'K',2e5,'Cmax',200,'Kmax',1.26e6));
  % li=RM.rail;
  % 
  % RO.name='WC4';


elseif comstr(Cam,'railfromcontour')
 %% #d_rail('MeshRailFromContour');
 %% ToDo25  sdtweb d_visco 'nl_mesh contour'
 if carg>nargin;f1='@d_rail.m/Rail_UIC60.stp';
 else; f1=varargin{carg};carg=carg+1;
 end
 if carg<=nargin;RO=varargin{carg};
 else
    %RO=struct('Run','-1 -order 2 -clmax30 -clim6 -v 1');
    RO=struct('quad',1,'lc',10,'spline',1,'ntol',1e-4) % sxxx algo
    %RO.algo='delquad'
    %RO.algo='quadqs'
    %RO.algo='front2d'
    %RO.algo='initial2d'
    RO.algo='meshadapt'
 end
% mo1=fe_gmsh('write',sdtu.f.safe(f1),RO);
 mo1=nl_mesh('contour clmax30',sdtu.f.safe(f1),RO); 

 % test: use lsutil cut to get a quad dominant mesh, pretty good, but robustness to handle
 mo2=nl_mesh('contour qcut qlc1',sdtu.f.safe(f1),RO);


 if nargout==0;
   feplot(mo1);
   fecom shownodemark groupall
 end


elseif comstr(Cam,'railtop')
%% #MeshRailTop_For_Contact
model=varargin{carg};carg=carg+1;
RO=varargin{carg};carg=carg+1;

if isfield(RO,'sel')

 if ~isfield(RO,'topCoarse')||isequal(RO.topCoarse,0)||isequal(RO.topCoarse,'fine')|| ...
         (ischar(RO.topCoarse)&&strncmpi(RO.topCoarse,'fin',3))% refine
     model=fe_shapeoptim('RefineHexaSurf',model,RO.sel);
 end
 elt=feutil(['selelt' RO.sel],model);
 %% if multi-layer of 104 adjust to 1 layer (expanding soft joint problem)
 i1=elt(:,5)==104;r1=elt(i1,:);i3=[];
 if ~isempty(r1) % If gap present fill it
  NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
  [r2,i2]=min(model.Node(NNode(r1(:,1)),[5 6])*[1;1]);
  mo2=model;mo2.Elt=[elt(1,:);r1]; 
  i2=feutil('geolinetopo',mo2,struct('starts',r1(i2),'dir',[0 1 0],'cos',.5));
  i2=feutil('geolinetopo',mo2,struct('starts',i2{1},'dir',[1 0 0],'cos',.5));
  i2=horzcat(i2{:})';i3=[];
  if 1==2
   % make flat over two edges (bad idea drop not seen)
   mo2.Elt=feutil('addelt','quad4',[i2(1:end-1,1) i2(2:end,1) i2(2:end,end) i2(1:end-1,end)]);
   mo2.Elt(1,end+1:size(elt,2))=0; 
   elt=[elt(~i1,:);mo2.Elt(2:end,:)];
  elseif 1==1 
   NNode=sparse(model.Node(:,1),1,1:size(model.Node,1));
   i3=NNode(i2(:,2:end-1));
  end
 end
 % now add and make sure that the orientation is correct
 elt(2:end,5:6)=2; 
 elt(1,8)=-1; % contact surface is display only'
 mo1=model;mo1.Elt=elt; mo1.Elt=feutil('orient 1 n 0 0 10000;',mo1);
 model=feutil('addelt',model,mo1.Elt);
 model=d_rail('MeshRailDist',RO.projM,model);
 if ~isempty(i3)
  % make a dip to that contact is lost
  model.Node(i3,7)=model.Node(i3,7)-3; %
  'xxx shapeoptim';
 end
 1;
 
else
%% early version
n1=feutil('getnode proid2',model);
sdtw('_ewt','obsolete ? ')
if min(n1(:,5))>min(RO.Cz.Y(:,1))
  % contact area does not cover trajectory
  st1=sprintf('matid 7 & selface & facing >.5 0 0 1e5 & innode{z>%g}',max(n1(:,7))-5);
  %cf.sel=st1;
  model.Elt=feutil('removeelt proid 2',model);
  elt=feutil(['selelt' st1],model);
  %% need to do a closure gap
  mo2=model;mo2.Elt=elt; mo2.Elt=feutil('selelt seledge',mo2);
  i1=mo2.Elt(isfinite(mo2.Elt(:,1)),1:2);
  n1=feutil('getnode',mo2,i1(:,1));
  n2=feutil('getnode',mo2,i1(:,2));
  [un1,i2]=max(abs(n1(:,5:7)-n2(:,5:7))');n1(i2==1,:)=[];n2(i2==1,:)=[];
  i2=abs(n1(:,5)-min(n1(:,5)))<50|abs(n1(:,5)-max(n1(:,5)))<50;
  n1(i2,:)=[];n2(i2,:)=[];
  r1=sortrows(round(unique([n1;n2],'rows')*1000),[5 6])/1000;
  i2=reshape(r1(:,1),[],2);
  elt(end+(1:size(i2,1)-1),[4 3 2 1])=[i2(1:end-1,1) i2(2:end,1:2) i2(1:end-1,2)];
  % Now add   
  i1=isfinite(elt(:,1));elt(i1,5:6)=2; elt(i1,7:end)=0;elt(1,8)=-1;
  model=feutil('addelt',model,elt);
  % cf=feplot(model);cf.sel='proid2';fecom showmap
end
end
 [eltid,model.Elt]=feutil('eltidfix;',model);
 r1=stack_get(model,'','Contact','g');
 if ~isempty(r1)
  model=feutil(sprintf('addset EltId -ID%i',r1.il(6)),model,'slave_Contact','proid 2');
 end
 out=model;
if nargout==0;
    cf=feplot(2);cf=feplot(out);fecom('showfipro');
end

elseif comstr(Cam,'gap');[CAM,Cam]=comstr(CAM,4);
%% #MeshGap

li={ ...
'CER',struct('h',5,'E',310e9,'nu',.25,'rho',3000,'Elim',9e12),'Alumine C799'
'EPC',struct('h',15,'E',30e9,'nu',.33,'rho',2100,'Elim',495e6),'Durostone'
'PA11',struct('h',30,'E',.18e9,'nu',.4,'rho',1049,'Elim',34e6),'Polyamide 11'
'Air',struct('h',30,'E',1e3,'nu',.4,'rho',1e-3,'Elim',9e12),'Air'
  };
if nargout==0
 disp(comstr(li,-30,struct('NoClip',1)))
end
if carg>nargin
 % Return material list
 for j1=1:size(li,1);li{j1,2}.name=li{j1,3};end
 out=containers.Map(li(:,1),li(:,2));
 return
end
model=varargin{carg};carg=carg+1;

% To mesh gap first find the two surface around the air
st='matid 7 & selface & facing >.9 10000 0 53 & innode {x>-8 & x<-7}';
%cf.sel=st;
mo1=model;mo1.Elt=feutil(['selelt ' st],mo1);
mo1=feutil('extrude 3 5 0 0',mo1);mo1.Elt=feutil('set group 1 matid 104 proid 104',mo1);

mo2=feutil('addelt',model,mo1.Elt);mo2.Node=mo1.Node;feplot(mo2); fecom showfipro
out=mo2;
elseif comstr(Cam,'build_sections');
%% #MeshList.Build_sections (init submodel map if not loaded from Sections.mat) -2 
 model=varargin{carg};carg=carg+1;
 RO=sdth.sfield('addselected',RO,model,{'pl','il','unit','Stack'});
 RO=stack_rm(RO,'','#G.*');RO=stack_rm(RO,'set','');
 dbM=containers.Map;
 % xxx u33, u60, u50 rail sections
 % eu102, 
 % sdtweb _textag ra+s % should contain more documentation
 st={'ra+s','selface & innode {x==-1369}','Rail U33 + sleeper'
     'ra','innode{x<=-1054}&selface&innode{x==-1054}','Rail U33 '
     'e','innode{x<=-2.5}&selface&innode{x==-2.5}&matid~=1 2 3 8','Eclisse U102'
     'ra+e','innode{x<=-7.5}&selface&innode{x==-7.5}','Rail + eclisse'
     'ra+se','innode{x<=-162.2}&selface&innode{x==-162.2}','Rail + sleeper+eclisse'
     'ra+sbe','innode{x>=-259 & x<=-211}&matid~=1 2 3 8','Rail + sleeper+eclisse bolt'
     'ra+be','innode{x>=-259 & x<=-211}&matid~=1 2 3 8 & proid~=9 10','Rail+eclisse bolt'
     };
 for j2=1:size(st,1)
  r2=struct('name',st{j2,1},'sel',st{j2,2},'Elt',[]);
  r2.Elt=feutil(['selelt',r2.sel],model);
  r2.Node=model.Node;r2.Node=feutil('getnode groupall',r2);
  r2.Node(:,5)=r2.Node(:,5)-min(r2.Node(:,5));
  if strcmpi(r2.name,'ra+se')||strcmpi(r2.name,'ra+e')
   % use the coarse mesh for top  
   r2.Elt=feutil('removeelt innode {z>58}',r2);
   r3=dbM('ra');r3.Elt=feutil('selelt innode {z>58}',r3);
   r3.Node=feutil('getnode groupall',r3);
   r2=feutil('addtest Merge -noori;',r2,r3);
  elseif strcmpi(r2.name,'ra+sbe')
   % use the coarse mesh for top  
   r2.Elt=feutil('removeelt innode {z>58}',r2);
   r3=feutil('extrude 0 1 0 0',r3,[0 12 24 36 48]);
   r2=feutil('addtest Merge -noori;',r2,r3);
   
  elseif strcmpi(r2.name,'ra+be')
   % use the coarse mesh for top  
   r2.Elt=feutil('removeelt innode {z>58}',r2);
   
   r2=feutil('addtest Merge -noori;',r2,r3);
  end
  
  dbM(r2.name)=r2;
 %% #MeshList.WheelSection : add wheel section to list -3
 mo2=model;mo2.Elt=feutil('selelt proid 1 8 3',mo2);mo2.Stack=[];
 mo2.Node=feutil('getnode groupall',mo2);
 dbM('Wheel')=mo2; 
 save(f2,'-struct','RO')
end
elseif comstr(Cam,'case');[CAM,Cam]=comstr(CAM,5);
%% #MeshCase : deal with loading strategy xa/u za/u
model=varargin{carg};carg=carg+1;
model=fe_case(model,'stack_rm','DofLoad');
model=fe_case(model,'stack_rm','DofSet');
if comstr(Cam,'xazu')
 %XAZU : horizontal acceleration load and vertical disp
  mo1=fe_case(model,'reset');mo1.Elt=feutil('selelt proid 1 8',model);
  mo1.name=''; mo1=fe_case(mo1,'assemble -matdes 2 NoT -SE');
  d1=struct('def',sum(mo1.K{1}(:,fe_c(mo1.DOF,.01,'ind')),2),'DOF',mo1.DOF);
  d1.name='Mx';
  model=fe_case(model,'DofLoad','Xa',d1);
  i1=feutil('findnode y<-110 & z>638',model);
  i2=i1; %i2=feutil('findnode inelt{ withnode{proid1}}',model);
  d1=struct('def',ones(size(i1)),'DOF',i1+.03);
  model=fe_case(model,'DofSet','Zu',d1);
elseif comstr(Cam,'xaza')
 %% horizontal and vertical load
  mo1=fe_case(model,'reset');mo1.Elt=feutil('selelt proid 1 8',model);
  mo1.name=''; mo1=fe_case(mo1,'assemble -matdes 2 NoT -SE');
  d1=struct('def',sum(mo1.K{1}(:,fe_c(mo1.DOF,.01,'ind')),2),'DOF',mo1.DOF,'name','Mx');
  model=fe_case(model,'DofLoad','Xa',d1);
  i1=fe_c(mo1.DOF,.03,'ind'); 
  d1=struct('def',sum(mo1.K{1}(:,i1),2),'DOF',mo1.DOF,'name','Mz(117 kN)');
  d1.def=d1.def/sum(d1.def(i1))*-117e6; %fe_mat('convertSIMM')
  model=fe_case(model,'DofLoad','Za',d1);
  
elseif comstr(Cam,'xczc')
 %% 2024 reimplement
elseif comstr(Cam,'xuzu')
 %% Enforce x and z motion of bearing 
 i1=feutil('findnode y<-110 & z>638',model);
 i2=i1; %i2=feutil('findnode inelt{ withnode{proid1}}',model);
 d1=struct('def',[ones(size(i2))*[1 0];ones(size(i1))*[0 1]],'DOF',[i2+.01;i1+.03]);
 model=fe_case(model,'remove','TrainWeight');
 model=fe_case(model,'DofSet','XuZu',d1);
end
out=model;

elseif isempty(Cam)&&isfield(varargin{2},'rail')
%% #MeshList : generic meshing strategy
 li=varargin{carg};carg=carg+1;
 if isfield(li,'rail');RO=li;li=li.Track;else; RO=struct;end
 projM=RO.projM; 

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
 
elseif comstr(Cam,'raildist');
 %% #MeshRailDist : distance to rail (force clean) 
 if isa(varargin{carg},'vhandle.nmap');projM=varargin{carg};carg=carg+1;end
 if ~isKey(projM,'RailDist')||isempty(projM('RailDist'))
  %d_rail('MeshRailDist',projM)   
  f2=d_rail('wd','Profil_RAIL.blk');mo2=nasread(f2);
  mo2.Node=feutil('getnodeGroupall',mo2);
  n1=mo2.Node(:,5:7); mo2.Elt(end+1,1:2)=[5094 5141];
  mo2.orig=min(n1);mo2.orig(2)=max(n1(:,2));n2=repmat(mo2.orig,size(n1,1),1);
  a=-atand(1/20);% 1/20 rail tilt
  n1=(mo2.Node(:,5:7)-n2)*[1 0 0;0 cosd(a) sind(a);0 -sind(a) cosd(a)]+n2;
  fecom('shownodemark',n1,'markersize',5)
  RD=struct('shape','Poly','Type','ExtLine','Node',[mo2.Node(:,1:4) n1],'axis',[1 0 0],'Elt',mo2.Elt);
  li=feval(lsutil('@dToPoly'),'init',RD); 
  projM('RailDist')=li;
 else;
    li=projM('RailDist');
 end
 if carg<=nargin
  mo1=varargin{carg};carg=carg+1;
  out=lsutil('surfstick',mo1,struct('sel','proid2','distFcn',li.distFcn));
 elseif ~sdtm.Contains(Cam,'back')% if isfield(RO,'Debug')&&sdtm.Contains(RO.Debug,'ShowProfWheel')
  cf=feplot; cf.sel='matid 7 104';
  lsutil('viewls',cf,li);set(cf.ga,'clim',[-30 30])
  fecom('shownodemark',li.Node(:,5:7),'markersize',5)
  mo1=cf.mdl;
  n1=feutil('getnode proid2',mo1);
  %d_rail('MeshRailDist')
 end
elseif comstr(Cam,'unis');
 %% #MeshUniS : load uniaxial mesh test (from sdtweb t_ident secmodtest)
 %     [mo1,RO]=dfr_ident('MatSimo',mo1,RO); %get mat properties set dt

 fname=d_rail('wd','UniS_sector.mat');%fname=d_rail('wd','mat/23_HR/UniS_sector.mat');
 load(fname,'mo1')

 mo1.Stack{end}.NLdata.f=[20;1000];
 mo1.Stack{end}.NLdata.g=[.1;.1];
 dbstack
 out=mo1;
 
elseif nargin==3&&isfield(varargin{3},'nmap') 
 %% #MeshDOE 2022 DOE call from d_mesh
 %mo1=rail19('load',CAM); out=mo1;
 if ~isempty(evt);
 elseif isfield(varargin{3},'RangeEvt');evt=varargin{3}.RangeEvt;
 else;evt=[];
 end
 if nargin>2&&~isfield(evt,'projM')&&isfield(varargin{3},'nmap');
     evt.projM=varargin{3}.nmap;
 end
 RM=nameToMeshRO(CAM,evt);
 out=d_rail('Mesh',RM); % sdtweb d_rail meshlist
 % if ishandle(2)
 %   try % Set feplot(2).nmap to experiment
 %    c2=get(2,'userdata');c2.data.nmap=struct('nmap',evt.projM);
 %    r3=evalin('caller','RO.RangeEvt');r4=nmap('CurDoLi'); 
 %    if isequal(r3,r4.RangeEvt);nmap('CurName')=r4.name{1};
 %    else;nmap('CurName')=r3.CurName;
 %    end
 %   end
 % end
else
   error('Mesh%s',CAM);
end
elseif comstr(Cam,'nl');[CAM,Cam]=comstr(CAM,3);
%% #NL DOE commands to define non-linearities 
model=varargin{carg};carg=carg+1;
RO=varargin{carg};carg=carg+1;
if strcmpi(RO.NL,'pada')
  %% #NL.PadA first attempt at hyperelastic material
  % Zhuravlev table 2.3 page 48 : {Yeoh,C10 2.5264, C20 -0.9177, C30 0.4711, g1 0.5688, tau 2.635}
  %  W=\sum_{p=1}^3 C_p (I-3)^p 
  sdtw('_ewt','implement NL:padA')
  %[st,r2]=sdtm.urnPar(RO.CaseVal,'{}{pro%s,mat%s,kappa%ug,v%s}');
  %model.pl=[ 1 fe_mat('m_elastic','SI',1) 6e6 .4999 1500]; % dd iso
  model.pl=[ 1 fe_mat('m_elastic','SI',1) 6e6 .49925 1500]; % dd iso *500 larger
  dd=feutil('getdd',[1 1 3 1],model);dd.ci_ts_egt=[1 5 9 8 7 4]';
  dd.iso=[ones(3) zeros(3);zeros(3,6)]*(ones(9,1)\reshape(dd.dd(1:3,1:3),[],1));
  dd.dev=dd.dd-dd.iso; dd.obs=[dd.iso(3,:);dd.dev(3,:)];
  model.info.dd=dd;% model.info.dd.iso./model.info.dd.dev

 model=stack_set(model,'info','oProp',mklserv_utils('oprop','RealNonSym'));
 % model=stack_set(model,'info','oProp',{'method','ldl'});
 %   [st,r2]=sdtm.urnPar(RO.mat,'{}{pro%s,mat%s,kappa%ug,NL%ug}');
 model.unit='SI';
  if any(strcmpi(RO.CaseVal,'vb')); 
   %% PadA{vb}  use UP formulation 
   model.il=p_solid('dbval 1 d3 20002');
   NLdata=struct('type','nl_inout', ...
     'opt',[0,0,0], ...
     'MexCb',{{m_hyper('@hyper_g'),struct}},'adofi',-.99*ones(29,1));
   model=feutil('setpro 1 isop100',model,'NLdata',NLdata);
  end

  % sdtweb d_mesh MatHyUp
  out=model;
 
else
   error('NL%s',CAM);
end
elseif comstr(Cam,'view');[CAM,Cam]=comstr(CAM,5);
%% #View commands
if comstr(Cam,'wheelctca')
%% #ViewWheelCtcA : wheel contact A 
 cf=get(2,'userdata');
 cf.sel={'withnode {proid 3}& proid8','colordataEnerK'};
 iimouse('view',cf.ga,[ -170.4 -6.865 -460.1 -170.4 -6.865 88.78 0.00 -1.00 0.00 9.79]);
elseif comstr(Cam,'railgz')
 cf.sel={'proid~=1 8 & eltname~=celas','showfiCevalz'};
 i1=intersect(feutil('findnode proid3',cf.mdl),cf.sel.Node);
 cf.SelF{1}.fvcs=sprintf('r1=r1(:,3)-r1(%i,3);',find(i1(1)==cf.sel.Node));
 fecom colorscaleinstanttight

elseif comstr(Cam,'railarrow')
 %% #viewRailArrowA % from fe_time iterNewton rail19('viewrailArrowA');
 if evalin('caller','~exist(''dq'',''var'')')
  d1=evalin('caller','struct(''def'',Case.T*[u r],''DOF'',Case.mDOF,''data'',[norm(u(:),''inf'');norm(r(:),''inf'')])');
 elseif evalin('caller','exist(''v'',''var'')')
  d1=evalin('caller','struct(''def'',Case.T*[u r dq v],''DOF'',Case.mDOF,''data'',[norm(u(:),''inf'');norm(r(:),''inf'');norm(dq(:),''inf'');norm(v(:),''inf'')])');
 else
  d1=evalin('caller','struct(''def'',Case.T*[u r dq],''DOF'',Case.mDOF,''data'',[norm(u(:),''inf'');norm(r(:),''inf'');norm(dq(:),''inf'')])');
 end
 cf=feplot;
 if Cam(end)=='a';
   st1='proid 2 3';if ~strcmpi(st1,cf.sel.selelt);cf.sel=st1;fecom('undefline');end
 end
 cf.def=d1;
 assignin('caller','cf',cf);fecom(cf,'animfreq');
 
elseif comstr(Cam,'defx1')
 %% #viewDefX1 % from fe_time iterNewton rail19('viewDefX1');
 % allow amplification of directions other than x 

 cf=feplot; 
 cna=cf.SelF{1}.cna; if isempty(cna); return;else; cna=cna{1};end
 if isnumeric(cna)
  cf.SelF{1}.cna{1}=@(def,coef)reshape(reshape(cna*double(def),[],3)*diag([1./coef 1 1]),[],1);
  feplot;
 end   
 
elseif comstr(Cam,'railctca')
%% #ViewRailCtcA : Rail contact A : stress in initial top section 
 % cf.def=fe_def('subdef',cf.def,'end');t_exp19('viewRailCtcA')
 if carg<=nargin; model=varargin{carg};carg=carg+1;else;model=[];end
 cf=sdth.urn('feplot(2)');if isempty(cf.mdl.Elt);cf.model=model;end
 %cf.sel={' proid7 &innode {x>-259 & z >83 & x<200}','colordataEnerK'};
 %cf.sel={' proid7 &withnode {x>-259 & z >83 & x<200 & y>-9}','colordataEnerK'};
 
 cf.sel={' proid7 &withnode {x>-170 & z >73 & x<-140 & y>-9.3 & y<23.3}','colordataStress -comp3'};
 fecom('colorbar',d_imw('get','CbTR','String','\sigma_{zz} kPa'));
 if 1==2
  cf.sel={' proid7 &withnode {x>-170 & z >73 & x<-140 & y>-9.3 & y<23.3}','colordataStress mises'};
  cf.os_('CbTR_Mises');%d_imw('get','CbTR','String','z (mu)'))     
 %     'FiMat', {'@PlotInfo',{'colordataMat-edgealpha.1-alpha1','','Set o1 def0',''}}
 end
 
 fecom('colormap parula(5)')
 r1=findall(2,'type','ColorBar');if length(r1)>1; delete(r1(2:end));end
 g=r1(1);
 set(g,'units','normalized','position',[.7 .5 .04 .4])
 
 fecom colorscaleinstanttight
 
 iimouse('view',cf.ga,[ -328.2 -220.7 241.7 -155 5 77.38 0.30 0.40 0.87 7.80]);
 
elseif comstr(Cam,'railctc')
%% #ViewRailCtc [BC...] : z disp of top rail surface
%  ctcC contact surfaces
 if carg<=nargin; model=varargin{carg};carg=carg+1;else;model=[];end
 if carg<=nargin; d2=varargin{carg};carg=carg+1;else;d2=[];end

cf=sdth.urn('feplot(2)');if isa(model,'sdth');cf=model;model=[];end
if ~isempty(cf.mdl.Elt);elseif ~isempty(model);cf.model=model;
else; error('do first :  cf=feplot(2);cf.model=mo1');end
if isempty(comgui('iminfo',cf))
  t_exp19('reset');cingui('plotwd',cf,'@OsDic(SDT Root)',{'ImToFigN','ImSw50','WrW49c'});  
end
 if ~isempty(strfind(Cam,'-woff'))
  %% attempt at correcting vertical wheel offset 
   i1=feutil('findnode proid 1',cf.mdl);
   i2=feutil('findnode proid 8',cf.mdl);
   i3=fe_c(d2.DOF,[i1;i2]+.03,'ind');
   d2.def(i3,:)=d2.def(i3,:)-mean(d2.def(fe_c(d2.DOF,i1,'ind'),1));
 end
st2='colordataEvalZ-edgealpha.1';st3=';colorscaleinstanttight';
st4='';

sdtroot('SetOsDic',{'FiTimeCEnd', ...
     {'@PlotInfo',{'colorscaleinstanttight','';'chend',''}}});
   
sdtroot('SetOsDic',{'FiCtcE', ...
     {'@ColorBar',{'units','normalized','position',[0.7 0.5 0.04 0.4], ...
     'YAxisLocation','left','FontSize',12,'TickLabelsMode','auto', ...
     '@xlabel',{'String',' w (mm)','FontSize',[14]}}, ...
       '@PlotInfo',{'colorscaleinstanttight','';'coloralpha',''}}});
sdtroot('SetOsDic',{'FiCtcENoAlpha', ...
     {'@ColorBar',{'units','normalized','position',[0.7 0.5 0.04 0.4], ...
     'YAxisLocation','left','FontSize',12,'TickLabelsMode','auto', ...
     '@xlabel',{'String',' w (mm)','FontSize',[14]}}, ...
       '@PlotInfo',{'colorscaleinstanttight','';'coloralphaoff','';
         'coloredgealpha.1-alpha.1',''}}});

if comstr(Cam,'railctcb')
 % #CtcB : just top of rail -3
 st='proid7 & selface & innode {z>81 & x > -270 & x < 600}';
elseif comstr(Cam,'railctcc')
 %% CtcC top rail and bottom wheel -3 
    % rail19('viewRailCtcC -woff',cf,d1);
 st='proid 2 3';st4='FiCtcENoAlpha';
 
elseif comstr(Cam,'railctcd')
 % rail19('ViewRailCtcD')% cut along the wheel
 % rail19('viewRailCtcD -woff',cf,d1);
 st='proid 7 8 & innode {x>-167 & x<-150 & z<100}';st4='FiCtcE';
elseif comstr(Cam,'railctce')% CtcE the wheel
 st='proid 7 8 & innode{y<30}';st4='FiCtcE';
elseif comstr(Cam,'railctcf')
 %% #ViewCtcF t_exp19('viewRailCtcF',model,d3) : edge
 st='seledge';st4='FiCtcE';
 
elseif comstr(Cam,'railctcg')
 %% #ViewCtcG rail19('viewRailCtcg',model,d3) : edge+surface
 % rail19('viewRailCtcg')
 if isfield(model,'def');cf.def=model;model=[];end
 if isempty(model);model=cf.mdl.GetData;end
 mo1=model;
 if contains(Cam,'nowheel')
  mo1.Elt=feutil('selelt matid ~=1 3 8 & seledge',model);
  mo1=feutil('addelt',mo1,feutil('selelt proid 2',model));
 else
  mo1.Elt=feutil('selelt seledge',model);
  mo1=feutil('addelt',mo1,feutil('selelt proid 2 3',model));
 end
 sel=feutil('getpatch new',mo1);
 cf.SelF{1}=sel;st='';st4='';
 cla(cf.ga);fecom(cf,'colorbarOff');fecom('colordataEvalZ-edgealpha.1');
 if sdtkey('cvsnum','feplot')>1.527
  cf.SelF{1}.ColorDataInit={{@rail19,'ViewRailCtcG-x1'},sdtdef};
  d1=cf.def; 
  if isempty(d1)||~isfield(d1,'data');
  elseif d1.data(2,1)<1e-4;cf.DefF{1}.data(:,1)=d1.data(:,1)*1000;
  end
 end
  
else; st='';
end
% cf.SelF{1}.fvcs='r1=r1(:,3)/norm(r1(:,3),''inf'');'
try; if ~isempty(st)&&isequal(cf.sel.selelt,st);st='';end;end
if ~isempty(st);cf.sel={st,st2};end
if ~isempty(d2);cf.def=d2; cf.os_('FiTimeCEnd');end
if contains(Cam,'x1');  % rail19('viewRailCtcf -x1 -time10');
    rail19('viewDefX1');
    fecom colorscaleinstanttight
end
if any(Cam=='{')
 [~,RA]=sdtm.urnPar(CAM,'{}{x%g,time%g}');
 if isfield(RA,'time');fecom(sprintf('animtime%i',RA.time));end%rail19('viewRailCtcf -x1 -time10');
end
cf.os_(st4);

elseif comstr(Cam,'viewminm')
 %% #ViewMinM
 cf=feplot;
 mo2=cf.mdl.GetData;
 [mo2,C2]=fe_case(mo2,'assemble -matdes 20 1 NoT -SE');
 r=diag(mo2.K{1});%r=r/r(1);%figure(1);semilogy(sort(r))
 r=-log10(abs(r));figure(1);plot(sort(r))
 i2=fe_c(mo2.DOF,.03,'ind');
 cf.def=struct('def',r(i2),'DOF',mo2.DOF(i2));
 cf.sel='reset';fecom('showfiCEvalz');
 set(cf.ga,'clim',[3.5 4],'alim',[3.5 4]);fecom coloralpha
 fecom('colorbaron');
 
elseif comstr(Cam,'esttime');
 %% #ViewEstTime
 
 RO=varargin{carg};carg=carg+1;
 RO.test=(36-31)/100;
 % fprintf('Estimated total time %s (h:m)\n',datestr(50e-3/8e-7*/3600/24,'HH:MM'))
 % fprintf('Estimated total time %s (h:m)\n',datestr(50e-3/1e-7*.12/3600/24,'HH:MM'))

 fprintf('Estimated total time %s (h:m) for 50 ms\n',datestr(50e-3/RO.dt*RO.tstep/3600/24,'HH:MM'))

elseif comstr(Cam,'cs');
 %% #ViewCs : extract information about plastified elements

 d2=varargin{carg};carg=carg+1;
 r2=varargin{carg};carg=carg+1;
 model=r2.model; NL=model.NL{2,3};
 ind=feutil('findelt proid 17',model);
 C1=fe_def('subdofind',d2.FNL,sprintf('%i+(1:%i)',NL.iopt(1),NL.iopt(5)*NL.iopt(4)));
 C1=struct('X',{{(1:6)',(1:8)',ind,[d2.FNL.data]}}, ...
    'Xlab',{{'comp','gauss','eltind','time'}}, ...
    'Y',reshape(C1.def,6,8,length(ind),size(d2.FNL.def,2)));
 a=squeeze(any(squeeze(any(C1.Y~=0,1)),1));
 ie=find(any(a,2));it=find(any(a,1));{it,ie}

 % Display the elements with a certain plasticity

 if nargout==0
  cf=feplot;if ~isequal(cf.mdl.Elt,model.Elt);cf.model=r2.model;end
  cf.sel(1)='proid 2';
  cf.sel(2)=['eltind' sprintf('%i ',C1.X{3}(ie))];
  cf.o(1)='sel 1 def 0 '; % mesh
  cf.o(2)='sel 2 def 0 '; % mesh
 end
 out=C1; 
 if 1==2
  squeeze(reshape(C1.Y(:,:,ie,it(1)),6*8,[]))
  ie=[];
  nnz(squeeze(C1.Y(:,1,:,:,:)))
 end
elseif comstr(Cam,'pad')
%% #ViewPad{} : display min vertical motion
 
[~,RO]=sdtm.urnPar(CAM,'{x%g}{jpar%g,tmin%g}');
if ~isfield(RO,'tmin')||isempty(RO.tmin);RO.tmin=.05;end
nmap=evalin('base','RT.nmap');
C4=nmap('RangeMacro');
if ~isfield(RO,'x')||isempty(RO.x);RO.x=[-300 300 900];end
st=cellfun(@(x)sprintf('Pad%i',x),num2cell(RO.x),'uni',0);
it=C4.X{1}>RO.tmin; 
gf=1;figure(gf);clf;
for j2=1:size(C4.Y,3)
for j1=1:length(st)
 i1=find(sdtm.regContains(C4.X{2}(:,1),st{j1}));
 r2=vhandle.cdm.urnVec(C4,sprintf('{y,l%i}{y,l%i}{cp}',i1(2),i1(1)));
 h=plot(C4.Y(it,i1(2),j2),C4.Y(it,i1(1),j2),'DisplayName',[C4.X{3}{j2}]); 
 xlabel(regexprep(r2.Xlab{1},'[-\d]',''))
 try;ylabel(regexprep(r2.X{2}{1},'[-\d]',''));end
 hold on
end
end
hold off
legend('location','best','interpreter','none')
cingui('plotwd',gf,'@OsDic(SDT Root)',{'ImToFigN','ImSw80','WrW49c'});
grid on; axis tight; 

%{y,l1}

elseif comstr(Cam,'dz')
%% #ViewDz : display min vertical motion
%rail19('viewdz',model,d1);
if carg>nargin; cf=feplot; mo1=cf.mdl; d1=cf.def;
else
 mo1=varargin{carg};carg=carg+1;
 d1=varargin{carg};carg=carg+1;
end
if comstr(Cam,'gap')
 i1=feutil('findnode proid2',mo1);
 i2=feutil('findnode proid3',mo1);
 figure(1);
 plot(-d1.data(:,4),min(fe_c(d1.DOF,i1+.03)*d1.def)-min(fe_c(d1.DOF,i2+.03)*d1.def));
elseif comstr(Cam,'dzrail')
 %% #ViewDzRail : display map 
 n1=feutil('getnode proid2',mo1);
 [n2,i2]=sortrows(round(n1*100),[5 6]);n1=n1(i2,:);
 i1=reshape(n1(:,1),length(unique(n2(:,6))),[])';
 r1=reshape(n1(:,5),size(i1,2),[])';
 r2=reshape(n1(:,6),size(i1,2),[]);
 C1=struct('X',{{d1.data*1000,r1(:,1),r2(:,1)}},'Xlab',{{{'Time','ms',[]},{'x','mm',[]},'y'}}, ...
     'Y',reshape((fe_c(d1.DOF,i1(:)+.03)*d1.def)',size(d1.def,2),size(r1,1),[]));
 C1.name='DzRail';C1.Ylab='Dz'; 
 C1.PlotInfo=ii_plp('plotinfo2D -type surface',C1); 
 C1=sdsetprop(C1,'PlotInfo','ua.YFcn','r3=-r3;');
 ci=iiplot(13,';');
 cingui('plotwd',ci,'@OsDic(SDT Root)',{'ImToFigN','ImSw80','WrW49c'});
 iicom('curveinit',C1);ci.osd_('CbTR{string,-z}','CmParulaB');iiplot(ci);
 
end
out=ci;

elseif comstr(Cam,'fz');[CAM,Cam]=comstr(CAM,3);
%% #ViewFz : display the vertical force

% 12e3*9.81*1e3*1e-6
if carg<=nargin; d2=varargin{carg};carg=carg+1;
else; cf=feplot;d2=cf.def;
end
if carg<=nargin&&isfield(varargin{carg},'Elt');mo1=varargin{carg};carg=carg+1;
end
% extract vertical resultant from FNL (not from K)
if ~isfield(d2.FNL,'data');d2.FNL.data=0;end
i3=find(strcmpi(d2.FNL.NL(:,2),'Contact'));NL=d2.FNL.NL{i3,3}; %#ok<FNDSB> 
NL.iNL=reshape(NL.iNL,[],length(NL.wjdet));
C2=struct('X',{{d2.FNL.data*diag([1000 1e-6 1 1]), ...
    {'Fx','kN';'Fy','kN';'Fz','kN';'Sc','mm^2'}}}, ...
    'Xlab',{{{'Time','ms',[];'Force','kN',[];'X','mm',[]},'Comp'}}, ...
    'Y',reshape(d2.FNL.def(NL.iNL,:),3,[],size(d2.FNL.def,2)),'Ylab',{{'Load','kN'}});
if ~isfield(d2,'name');C2.name='Forces (Fz target 117.72)';
else; C2.name=d2.name;
end
j1=strcmpi(d2.FNL.NL(:,2),'Contact');

 %out.data(:,2)=(NL.wjdet(:)'*out.FNL.def(NL.iNL(3,:),:))'; 

r3=NL.wjdet'*double(squeeze(C2.Y(3,:,:))~=0);% surface
r2=ones(size(NL.wjdet))'/1e6; % kN
C2.Y=reshape(r2*reshape(permute(C2.Y,[2 1 3]),size(C2.Y,2),[]),3,[])';
C2.Y(:,4)=r3;C2.X{2}(4,1:2)={'Surf','mm^2'};
[CAM,Cam,v]=comstr('v',[-25 2],CAM,Cam);
if comstr(Cam,'x');dbstack; keyboard;end
i1=~isfinite(C2.X{1}(:,2));C2.Y(i1,:)=[];C2.X{1}(i1,:)=[];
if isfield(d2,'info')
  C2.info=d2.info;
  if isfield(C2.info,'name');C2.name=C2.info.name;end
  if isfield(d2.info,'xgap')&&~isempty(d2.info.xgap)
   % Display band for gap
   try
   r2=interp1(C2.X{1}(:,3),C2.X{1}(:,1),d2.info.xgap);
   if all(isfinite(r2))
    C2.ID={struct('po',[r2],'marker','band', ...
        'tooltip',sprintf('gap %.1f %.1f',d2.info.xgap))};
    %c3.Stack{c3.ua.sList{1}}.ID=struct('po',c3.Stack{c3.ua.sList{1}}.info.xgap(2:end,:),'marker','band')
    % ii_plp(struct('po',C3.info.xgap(2:end,:),'marker','band'))
   end
   end
  end
end
if size(d2.data,2)>3
 coef=-208e8*1e-6; if isfield(d2.info,'Kpad');coef=-d2.info.Kpad(3);end
 C2.Y(:,end+(1:size(d2.data,2)-3))=coef*d2.data(:,4:end);
 C2.X{2}(end+(1:size(d2.data,2)-3),1)=cellfun(@(x)sprintf('Pad%i',x),num2cell(1:size(d2.data,2)-3),'uni',0);
end

if nargout>0;out=C2;
elseif size(C2.Y,1)<4; fprintf('Fx Fy Fz\n');disp(C2.Y);
elseif comstr(Cam,'fft')
  %% rail19('ViewFz fft -fmax 2000 tmin 52',d1);
  C2=ii_mmif([Cam  '-struct'],C2);
  c4=iiplot(4,';');iicom('curveinit','Fz',C2);
  cingui('plotwd',c4,'@OsDic(SDT Root)',{'ImToFigN','ImSw50','WrW99c'});
  iicom(c4,'ch',{'Load','Fz'});
  iicom('sub 1 1 1');
   
else
 ci=iiplot(3,';'); cingui('plotwd',ci,'@OsDic(SDT Root)',{'ImToFigN','ImSw50','WrW99c'});
 C2.PlotInfo={'sub 1 1 1 1 3','';'ua.axProp',{'@EndFcn',rail19('@trajText')}};
 iicom('curveinit',C2.name,C2);
% if size(ci.ax,1)~=1; iicom(ci,'sub 1 1 1  1 3');end
% ci.ua.axProp={'@EndFcn',rail19('@trajText')};
 iicom('chall');
end
 
% r1=[-170 -7.5 7.5]

elseif comstr(Cam,'def')
%% #ViewDef : do step and add trajectory
def=varargin{carg};carg=carg+1;
RO=varargin{carg};carg=carg+1;
def=fe_def('subdef',def,sprintf('nstep%i',RO.nstep));
i1=fe_c(def.DOF,feutil('findnode proid 1 8',RO.mdl)+.01,'ind');
def.def(i1,:)=def.def(i1,:)+RO.v*repmat(def.data(:,1)',length(i1),1);
out=def;

elseif comstr(Cam,'raild1')
%% #ViewRailD1

cf=feplot(3);d1=varargin{carg};carg=carg+1;

cf.sel={'proid 7 & withnode {x<50 & x>-370}'};fecom showficevalz
d2=d1;
d2.data(:,3)=of_time('lininterp',d1.FNL.NL{1,3}.table(:,1:2),d2.data(:,1),zeros(1,3))-175;
d2.LabFcn='sprintf(''t=%.0f ms x=%.1f mm F=%.2f kN'',def.data(ch,[1 3 2]).*[1000 1 1e-6])';

model=cf.mdl.GetData;
%r2=min(fe_c(d1.DOF,feutil('findnode proid8',model)+.03)*d1.def);
r2=(fe_c(d1.DOF,49343+.03)*d1.def);
d2.def=d2.def-ones(size(d2.def,1),1)*r2;

cf.def=d2; fecom('chend');
fecom colorscaleinstanttight;cf.os_('CmParulaB');



elseif comstr(Cam,'fcxmax')
%% #ViewFcxmax : extract max along x
% rail19('viewFcxmax')
cg=feplot(10); 
mo1=cg.mdl.GetData;d1=cg.def; 
if isempty(mo1.Node); rail19('ViewFcField');mo1=cg.mdl.GetData;d1=cg.def; end

[r1,i1]=min(mo1.Node(:,5:6)*[1;1]);i1=mo1.Node(i1,1);
i1=feutil('geolinetopo',mo1,struct('starts',i1,'dir',[0 1 0]));
i1=feutil('geolinetopo',mo1,struct('starts',i1{1},'dir',[1 0 0]));
i1=horzcat(i1{:}); % Contact surface topology
i1(:,~any(ismember(i1,fix(d1.DOF))))=[];
i1(~any(ismember(i1,fix(d1.DOF)),2),:)=[];
d1=feutilb('placeindof',i1(:)+.03,d1);
r2=squeeze(max(reshape(d1.def,size(i1,1),size(i1,2),[]),[],2));

C1=struct('X',{{mo1.Node(i1(:,1),5),d1.data(:,1)*1000}}, ...
    'Xlab',{{'x',{'Time','ms',[]}}}, ...
    'Y',r2); C1.DimPos=[2 1];
C1.Ylab={'Pressure','MPa'};
C1.PlotInfo=ii_plp('PlotInfo2D -type "surface"',C1);
C1=sdsetprop(C1,'PlotInfo','ua.YFcn','r3=abs(r3);');
ci=iiplot;iicom(ci,'curveinit','FcXMax',C1);
ci.osd_('CbTR{string,MPa}');

elseif comstr(Cam,'fcfield')
%% #ViewFcField : display contact field for dianostic
 
if carg>nargin; cf=get(2,'userdata'); d1=cf.def;
else;d1=varargin{carg};carg=carg+1;
end

r2=d1.FNL;
NL=r2.NL{end};n1=NL.GaussCoor; 
w=repmat(NL.wjdet',3,1);w=diag(sparse([1./w(:);0]));
r2.def=w*r2.def*1e-3; % go to MPa from milli-N/mm^2 
r2=fe_def('subdofind',r2,any(r2.def,2)); 

n1(:,1)=n1(:,1);
mo1=fe_fmesh('delaunay',n1);mo1.Elt=feutilb('t32q4',mo1);

cg=feplot(10);cg.model=mo1; cg.sel='withnode {x>-20 & x<20 &y>-10 & y<15}';
rail19('reset');cingui('plotwd',cg,'@OsDic(SDT Root)',{'ImToFigN','ImSw50','CrAll','WrW49c'});  
fecom showficevalz;

i2=sparse([51 52 50],1,[1 2 3]); 
r2.DOF=fix(r2.DOF)+i2(round(rem(r2.DOF,1)*100))/100; r2.fun=[0 4];
r2.data=d1.data; 
r2.data(:,3)=of_time('lininterp',d1.FNL.NL{1,3}.table(:,1:2),r2.data(:,1),zeros(1,3))-175;

r2.LabFcn='sprintf(''t=%.1f ms x=%.1f mm F=%.2f kN'',def.data(ch,[1 3 2]).*[1000 1 1e-6])';
cg.def=r2; fecom('chend');

r3=r2.def(:,min(80,size(r2.def,2))); [fe_c(r2.DOF(r3~=0)) num2cell(r3(r3~=0))];ans(1:20,:)
fecom colorscaleinstanttight
cg.os_('CmParulaB');cg.os_('CbTRNoLab');
r1=get(cg.ga,'colormap');r1(end,:)=1; set(cg.ga,'colormap',r1);
cg.SelF{1}.fvcs='r1=-r1(:,3);';feplot(cg);
 
%iimouse('view',gca,[ -0.3938 4.503 397.2 -0.3938 4.503 83.89 0.00 1.00 0.00 9.26]);
%iimouse('view',gca,[ -0.38 -0.1919 345.8 -0.38 -0.1919 84.15 0.00 1.00 0.00 8.17]);
iimouse('view',gca,[ 0.3539 2.522 332.1 0.3539 2.522 85.72 0.00 1.00 0.00 9.77]);
axis on

% i2=unique(fix(r2.DOF(any(r2.def(:,100:end-50),2))));
% figure(1);clf; plot(n1(i2,1),n1(i2,2),'+'); axis equal
% fecom('imwrite',struct('PlotWd',{{'@OsDic',{'FnI','ImMovie'}}}));

elseif comstr(Cam,'freqdef');[CAM,Cam]=comstr(CAM,8);
%% #ViewFreqDef : compute spectrum and display with shapes too

  cf=feplot(2,';');d1=varargin{carg};carg=carg+1; 
  d2=ii_mmif(['fft -struct' CAM],fe_def('subdof',d1,cf.sel.Node));
  d2=fe_def('subch curve',d2,{'DOF',find(any(d2.Y))});
  
  [r1,i1]=sort(mean(abs(d2.Y))); d2=fe_def('subch curve',d2,{'DOF',i1});
  c4=iiplot(4);iicom(c4,'curveinit','Def',d2);
  cf.def=struct('def',d2.Y.','DOF',d2.X{2},'data',d2.X{1}, ...
      'LabFcn','sprintf(''%i at %.4g Hz'',ch,def.data(ch))');
elseif comstr(Cam,'freq{');
 %% #ViewFreq{ renew frequency analysis 

 [~,RF]=sdtm.urnPar(CAM,'{}{tmin%g,spectro%ug}');
 if ~isfield(RF,'tmin')||isempty(RF.tmin);RF.tmin=[];
 else; RF.tmin=['tmin' sprintf(' %g',RF.tmin)];
 end



 if isfield(RF,'spectro')% rail19('viewFreq{spectro.02,tmin.05}');
  c21=sdth.urn('iiplot(3).Clone{21}');
  c3=get(3,'userdata');
  if carg<=nargin; C3=varargin{carg};carg=carg+1;
  else; C3=c3.Stack{c3.ua.sList{1}};
  end
  st=sprintf('spectro BufTime %g Overlap .95 fmin 200 2500 -window hanning %s',RF.spectro,RF.tmin);
  spec=ii_signal(st,C3);
  C2=spec.GetData; 
  C2.X{1}=[interp1(C3.X{1}(:,1),C3.X{1}(:,2),C2.X{1}) C2.X{1}];
  C2.Xlab{1}={'WheelPos','mm',[];'Time','s',[]};
  iicom(c21,'curveinit',C2)
  c21.os_('CbTR');iicom('ylin');

 elseif sdtm.Contains(CAM,'def') % rail19('viewFreq{def}') 
   cf=feplot(2,';'); nmap=cf.data.nmap.nmap;
   cf.sel='matid ~=1 3 8 & seledge';fecom showficevalz
   d4=nmap('CurTime');
   d3=ii_mmif(['fft  -fmin 1 5e3 -struct', RF.tmin],d4,struct('Out','def'));
   %d3=ii_mmif(['fft  -fmin 1 5e3 -struct', RF.tmin],d4,struct('Out','def'));
   d3.TR=fe_def('subdof',d4.TR,cf.sel.Node);
   cf.def=d3; 

 else
  c20=sdth.urn('iiplot(3).Clone{20}');
  c3=get(3,'userdata');
  if carg<=nargin; C3=varargin{carg};carg=carg+1;
  else; C3=c3.Stack{c3.ua.sList{1}};
  end
  ii_mmif(['fft -Display20 -fmin 1 5e3',RF.tmin],C3)

  iicom(c20,';sub 1 1 1 1 1;xlog')
  cingui('plotwd',3,'@OsDic(SDT Root)',{'ImGrid','ImToFigN','ImSw80','WrW49c'});
  cingui('plotwd',20,'@OsDic(SDT Root)',{'ImGrid','ImToFigN','ImSw80','WrW49c'});
 end
  % d3=ii_mmif('fft  -fmin 1 5e3 -struct',nmap('CurTime'))


elseif comstr(Cam,'stress');[CAM,Cam]=comstr(CAM,8);
%% #ViewStress : compute spectrum and display with shapes too

 RO=varargin{carg};carg=carg+1; 
 cf=feplot(2,';');d1=cf.def;
 
 mo1=fe_case(cf.mdl.GetData,'reset');
 mo1.Elt=feutil('selelt withnode',mo1,RO.at);
 cut=fe_caseg('StressCut-SelOut',struct('type','Gauss'),mo1);
 C2=fe_stress('stressMises',mo1,d1);

elseif comstr(Cam,'cz')
%% #ViewCz view wheel trajectory
 RT=varargin{carg};carg=carg+1; 
figure(10);clf % xxxAB pas figure 2 car feplot
subplot(1,2,1);
plot(RT.Cz.X{1}*1000,RT.Cz.Y(:,2)/1000,'r');%plot  the acceleration-time;
grid on
xlabel('Time (ms)'); ylabel('Acceleration (m/s^2)'); axis tight;
subplot(1,2,2); 
plot([RT.Cz.Y(:,1)'],RT.Cz.Y(:,4)/1000,'r');
legend('off')
hold on;xline(-7.5,'b');legend('speed','gap');grid on
xlabel('Displacement (mm)'); ylabel('Speed (m/s)'); axis tight;

elseif comstr(Cam,'hr');[CAM,Cam]=comstr(CAM,3);
%% #ViewHr : standard views for HR
if comstr(Cam,'full/red')
 r2=varargin{2};r3=varargin{3};carg=4;
 RB=struct;RB.Range=feutil('rmfield',r3,'Res','VectFcn','nmap','Project','level'); 
 RB.InterpFcn=sdtm.urnCb('interp1(linear,extrap)');
 C3=fe_def('curvejoin',RB,r3.Res);C3.name='Reduced';
 C3=fe_def('subchcurve',C3,{'Time',@(x)x(:,1)>.2});
 if isempty(r2); 
  C4=C3;
  RO.ch={'Comp',2;'jPar','all'};
  RO.li=sprintf('Li{L{%i-,%inone},M{%inone,%i.},C{%i}}',size(C3.Range.val,1)*ones(1,5));
 else
  C2=fe_def('curvejoin',struct('InterpFcn',{sdtm.urnCb('interp1(linear,extrap)')}),r2.Res);C2.name='Full';
  C4=fe_def('curvejoin',struct('InterpFcn',{sdtm.urnCb('interp1(linear,extrap)')}),C2,C3);
  C4.Xlab(3:4)={'A','Ty'};C4.name='Compare';
  RO.ch={'Comp',2;'Ty',1:2;'A',1:2}; RO.li='Li{L{2-,2none},M{2none,2.},C{2}}';
 end
 c12=iiplot(12,';');
 iicom('curveinit',C4);
 iicom('polar.Comp(2)');iicom('ch',RO.ch)
 d_imw(RO.li,12)
elseif comstr(Cam,'chfreq')
 r3=varargin{2},r4=varargin{3};carg=4;
 C3=fe_def('curvejoin',r3.Res);C3.name='Reduced';
 C3=fe_def('subchcurve',C3,{'Time',@(x)x(:,1)>r4(2)});
 %C4=fe_def('curvejoin',struct('InterpFcn',{sdtm.urnCb('interp1(linear,extrap)')}),C2,C3);
 %C4.Xlab(3:4)={'A','Ty'};C4.name='Compare';
 c12=iiplot(r4(1),';');
 iicom('curveinit',C3);
 iicom('polar.Comp(1)');
 iicom('ch',{'Comp',2;'Z join','all'})
 %d_imw('Li{L{2-,2none},M{2none,2.},C{2}}',np)
else;error('ViewHR%s',CAM);
end

elseif comstr(Cam,'traj');[CAM,Cam]=comstr(CAM,3);
%% #ViewTraj : view trajectory 

figure(1);subplot(121);
plot(C3.X{1}(:,1),gradient(C3.X{1}(:,3))/diff(C3.X{1}(1:2,1))/1000)
subplot(122); vhandle.cdm.urnVec(C3,'{x,1}{x,3}{gf1,tight,os{ImGrid}}');
cingui('plotwd',1,'@OsDic',{'ImLw40','WrW90c'});comgui('imwrite',1)

elseif comstr(Cam,'entry');[CAM,Cam]=comstr(CAM,6);
%% #ViewEntry
projM=d_rail('nmap');
while carg<=nargin
 tag=varargin{carg};val=[];
 if isKey(projM,tag); val=projM(tag);end
 if isempty(val);cM=projM('Map:Slice');if isKey(cM,tag);val=cM(tag);end;end
 
end

else; error('View%s',CAM);
end  
elseif comstr(Cam,'traj')
 %% #traj : build trajectory curve ----------------------------2
[~,RT]=sdtm.urnPar(CAM,'{}{v%ug,Fn%ug,zu%ug}');
dbstack; keyboard; 

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

%% #UniS : hyper reduction uniaxial test cases -2

RT=struct('nmap',vhandle.nmap);
li={'CamMesh';';'
    'SimuCfg{RO{NperPer12k,Nper8 5,snap100,freq20},SigA,Exp{4u,.4,uva110,rt-1e-3}}'
      'RunCfg{rail19($mo1,$RO,SolveHrRun{run,hr{1m,100u,chall,kinStra0},runhr{to0}})}'};
RT.nmap('FullRed')=li;

%% #UniS.vt : applied load virtual test -3
RT.nmap('dt')=1e-5; 
RT.nmap('Reduce')='nl_solve(ReducFree -Float2 -SetDiag -SE -Tgt -HR)';
RT.nmap('SigA')='Sig{cnInput,Sin(1.5,20,TR1)}';RT.nmap('JacCoef')=1;
li={'CamMesh';';' % Model Gart (2 DofLoad) sdtweb d_fetime MeshGart 
     'SimuCfg{SigA,EigOpt{2 2 0},ModalNewmark{$dt$,10,fc,chandle1,BetaR1.5e-4}}';';'
     'RunCfg{Reduce}'};
RT.nmap('InitVT')=li;

%% #UniS.mesh variants -3
RT.nmap('MeshA')='MeshCfg{rail19(UniS{SimoA}),DofSetSine{A15,freq20,st3}}';
RT.nmap('MeshB')='MeshCfg{d_fetime(VolTrack{matMatCur,st3}),TopZ{}}';
RT.nmap('MeshC')='MeshCfg{d_fetime(VolTrack{matMatCur,st3}),TopFZ{}}';
RT.nmap('MatCur')='m_hyper(urnSimoA{1,1,0,30,3,f5 20,g .33 .33,rho2.33n,unSI,isop100})';
%RT.nmap('MatZhu')='m_hyper(urnZhu{tyYeoh,C10 2.5264, C20 -0.9177, C30 0.4711, g1 0.5688, tau 2.635,rho2.33n,unSI,isop100})';
% m_hyper('urnzhu{Yeoh,C10 2.5264, C20 -0.9177, C30 0.4711,kappa 1000,rho2.33n,unSI,isop100}')
RT.nmap('MatZhu')='m_hyper(urnZhu{Yeoh,C10 2.5264, C20 -0.9177, C30 0.4711,kappa 1000,rho2.33n,unSI,isop100})';
RT.nmap('MatEZhu')='m_elastic(urnZhu{E15.2,nu.48,G5.05,rho2.33n,unSI,tyuMaxw,cells 8 .1 290 0 8 .1 390 0})'; % Elastic
RT.nmap('MatRafaelA')='m_hyper(urnSimoA{2,2,0,500,fv5,f50 200,g .033 .033,rho2.33n,unTM,isop100,tydfr_ident(matsimoa)})';


RT.nmap('SigA')='Sig{cnInput,Sin(15,20,TR32)}';

RT.nmap('CamMesh')=RT.nmap('MeshA');

RT.nmap('RedA')='rail19($mo1,$RO,SolveHrRun{runhr{dt10u,to0,zea},vcompa})';
RT.nmap('RedB')={'RunCfg{rail19($mo1,$RO,SolveHrRun{runhr{dt$dt$u,to0,ze$ze$,outResult}})}'};
RT.nmap('RangeLoopOpt')=struct('RangeLoopResKeys',{{'ObsT'}});

nmap('UniS')=RT; 

%% #nmap.Cst : constants to be reused
nmap('Ctc21')=struct('delta',20,'Fl',1e6,'Kc',-2.3e9);
nmap('CMsh21')=struct('CRailEdge',.1,'Vel',20,'ToolTip','xxx');

 %% Without eclisse
 if exist('rail19','file') % Append
   r2=rail19('nmap',struct('nmap',nmap));
 end
 ecM=useOrDefault(nmap,'Map:Slice');
 li= {'SecLab','xstart','xend','lc';
     'ra',-300,[],15
     'ra+s',-127.5,[],15
     'ra',127.5,[],15
     '',300,[],15};
 ecM('Ec0')=li;

%% #nmap.Trk configurations 
  nmap('Trk21ref')=vhandle.nmap.RT({ ...
    'Ref21','U30{5,Gc,tc-350_400}:Ec41{Air}:W2{XaZa}';
    'Coarse21','U30{5,Gc,tc-500_400,lc100}:Ec41{Air}:W2{XaZa,rigid}';
    'Coarse','U30{5,Gc,tc-500_400,lc100}:E0{Air}:W2{XaZa,rigid}';
    'SnapExp',{'MeshCfg{rail19(Ref21):21ref}', ...
     'RunCfg{rail19(SolvePadIO),rail19(SolveTimeB),rail19(SolveSnap{testA,1000,200})}'}
    'CoarseExp',{'MeshCfg{rail19(Coarse21):21ref}', ...
     'RunCfg{rail19(SolvePadIO),rail19(SolveTimeB)}'}
    'CurExp','SnapExp'
    'Map:Cin',nmap('Map:Cin');
    'Ctc',nmap('Ctc21');'CMsh',nmap('CMsh21');
      },'ToolTip','Reference transient');
  nmap('Trk24res')=vhandle.nmap.RT({ ...
    'Res24','U30{1}:W2{XaZa}';
    'DownExp',{'MeshCfg{rail19(Res24):21ref}', ...
         'RunCfg{}'}
    'CurExp','DownExp'
    'Ctc',nmap('Ctc21');'CMsh',nmap('CMsh21');
      },'ToolTip','Wheel down test');


  nmap('Trk21HakPos')=vhandle.nmap.RT({'HakPos','U30{7,Gc,tc-3500_1500}:Ec41{Air}:W2{XaZa}';
    'CurExp',{'MeshCfg{rail19(HakPos):21ref}', ...
     'RunCfg{rail19(SolvePadIO),rail19(SolveTimeB),rail19(SolveSnap{testA,1000,200})}'}
   },'ToolTip','xxx');

 nmap('Trk23PadSN')=vhandle.nmap.RT({ ...
  'FineExp',{'MeshCfg{rail19(U30{5,Gc,tc-700_300,$Pad$}:Ec41{Air}:W2{XaZa}):21ref}';
   'RunCfg{rail19(SolvePadIO),rail19(SolveTimeB)}'};
  'CoarseExp',{'MeshCfg{rail19(U30{5,Gc,tc-700_300,lc100,$Pad$}:Ec41{Air}:W2{XaZa}):21ref}';
   'RunCfg{rail19(SolvePadIO),rail19(SolveTimeB)}'};
  'CoarseExpS',{'MeshCfg{rail19(U30{5,Gc,tc-700_-615,lc100,$Pad$}:Ec41{Air}:W2{XaZa}):21ref}';
   'RunCfg{rail19(SolvePadIO),rail19(SolveTimeB)}'};
  'CoarseExpN',{'MeshCfg{rail19(U30{5,Gc,tc-1000_1000,lc100,$Pad$}:Ec0:W2{XaZa,in5}):21ref}';
   'RunCfg{rail19(SolvePadIO),rail19(SolveTimeB)}'};
  'CoarseExpNb',{'MeshCfg{rail19(U30{5,Gc,tc-1000_1000,lc100,$Pad$}:Ec0:W2{XaZa,in2}):21ref}';
   'RunCfg{rail19(SolvePeriod{rred2:9}),rail19(SolvePadIO),rail19(SolveTimeB)}'};
  'Map:Cin',nmap('Map:Cin');'Ctc',nmap('Ctc21');'CMsh',nmap('CMsh21');
  'RangeLoopOpt',struct('ifFail','Error','RangeLoopResKeys',{{}});
  'PadIo23',1;'DoSaveMacro',1},'ToolTip','DOE on pad properties');
 %RT=nmap('Trk23PadSN');

 %% #Trk24Noise model and range -3
 nmap('Trk24Noise')=vhandle.nmap.RT({ ...
  'Ref21','U30{2,Gc,tc-700_300,PadFuSN{io4}}:Ec0{}:W2{XcZc}';
  'Map:Cin',nmap('Map:Cin')
  'BcList',struct('ToolTip','Define contacts, edit BC', ...
    'li',{{ ...
   ...% 'Case{reset}'
    'Case{DofLoad,Top,"rb{setname Wheel &y<-110,dir 1 3,curveDownForward,KeepDof}"}'
    'Case{FixDof,FixYRailMinMax,"proid13 -DOF2"}'
    'Case{FixDof,YSim,"inelt{proid1 10&selface & facing<-.9 0 -2e5 0} -DOF2"}'
    'Case{Pcond,Scld,p_contact(''PcondScld'')}'
    'InitQ0'; % use the InitQ0 in RangeEvt
    'CbRefWheel';'CbRefRail';'CbStickWheel';'CbCtcGen'
    }});
  'CurExp',{'MeshCfg{rail19(Ref21):None{BcList}}';
  'SimuCfg{"Imp{.1m,.1,chandle1,BetaR7e-6}"}'
  'RunCfg{rail19(SolvePadIO),Time}'};
  'PadIo23',1;'DoSaveMacro',1;
  });
 %'DownForward',struct('X',{{[0 .01 1.01]*2/20, ... % T=dist2/vel 20
 %   {'x';'z'}}},'Y',linspace(0,1,length(t))'*[1 -.1], ...
 %   'LabFcn','sprintf(''x%.4f z%.4f'',def.data(ch,2:3))','fun',[0 4]);
 C2=struct('X',{{[0;.1],{'x';'z'}}},'Xlab',{{'Time','dir'}},'Y',[0 0;-1e3 -1e3]');
 TimeOpt=struct('InitAcceleration','feval(rail19(''@status24''));');
 Range=fe_range('buildGrid',struct('TimeOpt',{{'V1',TimeOpt}}, ...
    'InitQ0',{{'D1','elem0(vect{x0,y0,z-.16,"selwithnode{setNameWheel}"})'}}, ...
    'DownForward',{{'V2',C2}}));
 RT=nmap('Trk24Noise');RT.nmap('RangeA')=Range;

%% #PadA : pad compression test cases -2

RT=struct('nmap',vhandle.nmap);
RT.nmap('MeshA')='MeshCfg{rail19(PadA2),PadA{va},PadA}';
RT.nmap('MeshB')='MeshCfg{rail19(PadA2),PadA{vb},PadA}';
RT.nmap('CamMesh')=RT.nmap('MeshA');

li={'CamMesh';';'
    'SimuCfg{RO{NperPer12k,Nper8 5,snap100,freq20},SigA,Exp{4u,.4,uva110,rt-1e-3}}'
     };

RT.nmap('FullRed')=li;
RT.nmap('InT')='Sig{cnInput,Table(@lin(0,4,5),@lin(-.1e-4,-2e-3,5),ilin)}';
li={'MeshA';';'
    'SimuCfg{Static{-1m,chandle1,jcallstep1,fcrail19@vhe,MaxIter2},InT}';';'
   'RunCfg{Time}'};
RT.nmap('Static')=li;

nmap('PadA')=RT; 
RO.CAM=CAM;

li={'Shaft',[1 fe_mat('m_elastic','SI',1)  200e9 .3 7829];
    'Eclisse',[4 fe_mat('m_elastic','SI',1)  200e9 .3 7829] % Eclisse
    'Bolt',[5 fe_mat('m_elastic','SI',1)  200e9 .3 7829] % Ecrou
    'Rail',[7 fe_mat('m_elastic','SI',1)  200e9 .3 7829] % Rail Yield 800 MPa
    'Wheel',[8 fe_mat('m_elastic','SI',1)  200e9 .3 7829] % Roue Yield 800 MPa
    'Pad',[9 fe_mat('m_elastic','SI',1)  14e6 .45 1000 ] % Semelle, damping 0.1 kN/(m/s)
    'Sleeper',[10 fe_mat('m_elastic','SI',1)  20e9 .2 2500] % Traverse
    'Screw',[11 fe_mat('m_elastic','SI',1)  200e9 .3 7829] % Vis
    'PadAniso',[9  fe_mat('m_elastic','MM',6) 14e3 14e3 14e3 0.45 .45 .45 50e3 50e3 14e3]
    };
r2=vhandle.nmap(li,'Map:MatDb');
nmap('MatDb')=r2;

%% IO 
 %% #IO23SCNF: sensor configuration for pad testing with SDT ----2
 R2=struct('InCh',...
  struct('ChassisName','','ModID',...
  {'Mod1','Mod1','Mod1','Mod1','Mod2','Mod2','Mod2','Mod2'},...
  'ChannelID',{'ai0','ai1','ai2','ai3','ai0','ai1','ai2','ai3'},...
  'SensType',[{'Hammer'},repmat({'Accelerometer'},1,6),{'LoadCell'}], ...
  'ElecUnit',repmat({'mV'},1,8),...
  'PhyUnit',[{'N'},repmat({'g'},1,6),{'N'}],...
  'Sensitivity',{0.2191 10 10 10 10 10 10 1},... %
  'Point',{'Hammer{in}','Tri1','Tri1','Tri1','Mono1','Tri2','Tri2','LoadCell'},...
  'Dir',{'N/A','+X','+Y','+Z','+X','+Y','+Z','N/A'},...
  'IEPE',[repmat({1},1,7) {0}]));
 nmap('IO23pad')=R2;

 % sdtm.stdNmapOut('call')
 if nargout==1;out=sdtm.stdNmapOut(nmap,key,nargout,CAM);
 elseif nargout>1;[out,out1,out2]=sdtm.stdNmapOut(nmap,key,nargout,CAM);
 else; sdtm.stdNmapOut(nmap,key,nargout,CAM);
 end

%% #end
elseif comstr(Cam,'wd')
 %% #wd : Resolve fullfile from ProjectWd ----------------------------2
 wd={'@d_rail.m','@sdtdata/rail19/mat/21_sections'};
 if nargin==2&&ischar(varargin{2})
     out=sdtu.f.find(wd,varargin{2});
     if iscell(out)&&isscalar(out);out=out{1};end
 else
  base_wd=sdtu.f.firstdir(wd);
  out=comgui('cd',fullfile(base_wd,varargin{carg:end}));
 end
elseif comstr(Cam,'pwd')
 %% #pwd : Resolve fullfile from ProjectWd ----------------------------2
 out=comgui('cd',fullfile(sdtroot('PARAM.Project.ProjectWd'),varargin{carg:end}));
 if nargout>1; out1=exist(out,'file');end 

elseif comstr(Cam,'jup')
%% #jup : jupyter notebook building 
wd=fullfile(fileparts(which('d_rail')),'../jup');

RO=struct('BuildDir',sdtm.safeFileName('@tempdir\_jup'), ...
      'reset',1,'book','rail'); 
RO.helpDir=sdtu.f.cffile('@sdt/helpj');
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
   setenv('sdtfun',['rail19-d_rail'  ]); 
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

function RM=nameToMeshRO(name,evt);
%% #nameToMeshRO analyze naming nomenclature 
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
case 'wc2'; [opts,model]=rail19('LoadForContact RailWCoarse2');
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
r2=d_rail('MeshTrackLi',RM); 
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


function RT=WheelTraj(varargin)
%% #WeelTraj computation of wheel trajectory
RT=varargin{1}; if isfield(RT,'Node');model=RT;RT=varargin{2};else;model=[];end

if 1==2 
 %% Obsolete 1
 RT.Tdown=.02;
 H=[RT.xstart linspace(RT.xstart+(RT.v(1)*1000*RT.Tdown),-50,length(RT.v))]; % The space steps
 T=[0 [diff(H(:,1:end))]./[RT.v*1000]]; % The time steps
 [a,b]=size(H);T1=cumsum(T);
 [a2,b2]=size(RT.v);
 RT.v=RT.v*1000;% To calculate acceleration need the first term
 %B=[0 RT.v(end)]'; %Boundary conditions to solve the equation of seconde degree
 B=[0 0 (RT.v(end)-RT.v(1))]';
 %A=[T1(end)^2 T1(end); T1(end)^3/3 T1(end)^2/2]; % The matrix 2x2 
 A=[RT.Tdown^2 RT.Tdown 1 ;T1(end)^2 T1(end) 1; (T1(end)^3/3 -RT.Tdown^3/3) (T1(end)^2/2 -RT.Tdown^2/2) T1(end)-RT.Tdown]; % The matrix 3x3 
 X=A\B; % Find the unknows of the equation
 f=@(t) X(1)*t.^2+X(2)*t+X(3); % The equation of the acceleration
 V1=@(t) (X(1)*t.^3)/3+(X(2)*t.^2)/2+X(3)*t+RT.v(end)-(X(1)*t(end).^3)/3-(X(2)*t(end).^2)/2-X(3)*t(end);
 v2=V1(T1(1,2:end));

 F=f(T1(1,2:end)); % discretization
 V1=@(t) (X(1)*t.^3)/3+(X(2)*t.^2)/2+X(3)*t+RT.v(end)-(X(1)*T1(end).^3)/3-(X(2)*T1(end).^2)/2-X(3)*T1(end);
 v2=V1(T1(1,2:end));
 T1=[T1 (T1(end)+(RT.xend-H(end))/RT.v(end))];
 F(1)=0;
 H=[H +[linspace(H(end)+0.001,RT.xend,length(RT.v))]];
 v2=[v2 [v2(end)*ones(1,length(RT.v))]];
 RT.v=[RT.v [RT.v(end)*ones(1,length(RT.v))]];
 T2=[0 [diff(H(:,1:end))]./[RT.v]];
 T3=cumsum(T2);
 RT.Cz2=struct('X',{{T3',{'Xu';'Xa';'Zu'}}}, ...
 'Y',[H; ... % position
 [[0 F [zeros(1,b2)]]]; % Acceleration on mm/ms^-2
 [-0  -ones(1,(length(H)-1))]*0.075]'); % Vertical traj
 RT.vz0=diff(RT.Cz2.Y(1:2,3))/diff(RT.Cz2.X{1}(1:2))
elseif 1==3
 %% Obsolete 2 
 t=(0:RT.dt:RT.V.X{1}(end))';
 v=interp1(RT.V.X{1},RT.V.Y,t,'linea','extrap');
 x=RT.xstart+cumsum(v*RT.dt);
 RT.v=RT.v*1000;% To calculate acceleration need the first term
 V=RT.V; V.Y=V.Y*1000;%V.X{1}(2)=RT.Tdown; % V in m/s mesh in mm
 Fu=feval(nlutil('@cubeInterp'),'fu',V);
 dFu=feval(nlutil('@cubeInterp'),'dfu',V);
 H=cumsum(Fu(t))*RT.dt;H=H-H(1)+RT.xstart;
 H(H>RT.xend)=[];t(length(H)+1:end)=[];
 RT.Cz=struct('X',{{t,{'Xu';'Xa';'Zu'}}}, ...
 'Y',[H ... % position
  dFu(t) ... % Acceleration on mm/ms^-2
  [-0;-ones((length(H)-1),1)]*0.075]); % Vertical traj
 RT.vz0=-3.75;
%figure(1);plot(RT.Cz2.X{1},RT.Cz2.Y(:,2),'--',RT.Cz.X{1},RT.Cz.Y(:,2),'-')
%figure(1);plot(RT.Cz2.X{1},RT.Cz2.Y(:,1),'--',RT.Cz.X{1},RT.Cz.Y(:,1),'-')

 %figure(1);plot(t,dFu(t));
  
end

if ~isempty(model)
 Cz=RT.Cz; % possibly set opt.OutputFcn to select output time steps
 %% set curves based on Cz
 r1=fe_case(model,'stack_get','dofset');
 if isempty(r1) % xaza
  i1=strcmpi(Cz.X{2},'Xa');
  if nnz(i1)==0;error('Expecting Xa');
  else;
   r2=Cz; r2.Y(:,~i1)=[]; r2.X{2}(~i1)=[];
   model=fe_case(model,'SetCurve','Xa','Xacurve',r2);
  end
  r2=struct('X',{{[0;Cz.X{1}(end)],{'Fz coef'}}},'Y',[1;1]);
  model=fe_case(model,'SetCurve','Za','Zacurve',r2);
 elseif strcmpi(r1{1,2},'XuZu')
  if size(Cz.Y,2)~=2; error('Mismatch');end
  model=fe_case(model,'SetCurve','XuZu','XZ',Cz);
 elseif strcmpi(r1{1,2},'XaZu')
  %% xa, zu define curves for horizontal acceleration and vertical disp
  %  model=fe_case(model,'SetCurve','Xa','Xacurve',Cz);model.info.Cz=Cz;
  i1=strcmpi(Cz.X{2},'Xa');
  if nnz(i1)==0;error('Expecting Xa');
  else;
   r2=Cz; r2.Y(:,~i1)=[]; r2.X{2}(~i1)=[];
   model=fe_case(model,'SetCurve','Xa','Xacurve',r2);
  end
  i1=strcmpi(Cz.X{2},'Zu');
  if nnz(i1)==0;error('Expecting Zu');
  else;
   r2=Cz; r2.Y(:,~i1)=[]; r2.X{2}(~i1)=[];
   model=fe_case(model,'SetCurve','Zu','Zucurve',r2);
  end  
 elseif strcmpi(r1{1,2},'XaZa')
  %% xa, za define curves for horizontal and vertical acceleration

  i1=strcmpi(Cz.X{2},'Xa');
  if nnz(i1)==0;error('Expecting Xa');
  else;
   r2=Cz; r2.Y(:,~i1)=[]; r2.X{2}(~i1)=[];
   model=fe_case(model,'SetCurve','Xa','Xacurve',r2);
  end
  i1=strcmpi(Cz.X{2},'Zu');
  if nnz(i1)==0;error('Expecting Zu');
  else;
   r2=Cz; r2.Y(:,~i1)=[]; r2.X{2}(~i1)=[];
   model=fe_case(model,'SetCurve','Zu','Zucurve',r2);
  end  
  dbstack; keyboard;
 else;error('Not a valid case')
 end
 q0=stack_get(model,'','q0','g');
 if isempty(q0) 
  %% Initialize at trajectory
  mo2=model;mo2.Elt=feutil('removeelt proid2 3',mo2);[C1,DOF]=fe_mknl('initNoCon',mo2);
  q0=struct('def',C1.DOF*0,'DOF',C1.DOF);
  i1=fe_c(q0.DOF,feutil('findnode proid 1 8',model)+.01,'ind');
  i2=strcmpi(Cz.X{2},'Xu');
  q0.def(i1,1)=Cz.Y(1,i2); % Initial wheel position
  q0.def(i1,2)=diff(Cz.Y(1:2,i2))/diff(Cz.X{1}(1:2));% Initial X velocity
  i1=fe_c(q0.DOF,feutil('findnode proid 1 8',model)+.03,'ind');  
  q0.def(i1,2)=RT.vz0;% Initial vertical velocity
  model=stack_set(model,'info','q0',q0);
 end
 if isfield(RT,'q0')
  RT.q0=feutilb('placeindof',q0.DOF,RT.q0);
  if size(RT.q0.def,2)<2; RT.q0.def(:,2)=q0.def(:,2);end
  model=stack_set(model,'info','q0',RT.q0);
 end
 RT=feutil('rmfield',RT,'v'); 
 model.info=RT; RT=model;
end % Fill model 

end


