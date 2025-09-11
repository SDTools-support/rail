function [out,out1,out2,out3]=railui(varargin)

% RailUI / dynavoie callbacks for GUI
%  RB=railui('CbGeoInfo',mt)
 
% Etienne Balmes, Jean-Philippe Bianchi
% Copyright (c) 2004-2025 by SNCF/DRT & SDTools
% For revision information use railui('cvs')
 
persistent Verbose_Mode UI defaultSet

if isempty(Verbose_Mode); 
  defaultSet=sdtroot('@defaultSet');
    % GuiGF DefBut treeF
  if isempty(UI); 
      UI=struct('gf',[],'treeF',cinguj('@treeF'),'DefBut',[], ...
          'MainFcn',@railui,'toCinCell',cinguj('@toCinCell'), ...
          'tag','railui','DicFcn',{{'d_dynavoie','OsDic'}});
  end
  Verbose_Mode=0; 
  %dyn_utils('pathhelp') % Check that help is in path
end
% $CheckLic

%#ok<*NOSEM>
%#ok<*NASGU,*ASGLU,*TRYINC,*TRYNC,*CTCH>

%% OpenUI figure
if isempty(UI.DefBut); UI.DefBut=genDefBut; end

obj=[];evt=[];
if nargin>=1;
    if ischar(varargin{1});
     [CAM,Cam]=comstr(varargin{1},1);carg=2;
    else; % Callback
      [CAM,Cam]=comstr(varargin{3},1);carg=4;obj=varargin{1};evt=varargin{2};
    end
else; UI=initFig(UI); return
end

if sp_util('diag')>1; fprintf('railui(''%s'')\n',Cam); end

%% #Reset : UI reset -----------------------------------
if comstr(Cam,'reset'); [CAM,Cam]=comstr(CAM,6);
clear sdtroot; clear railui; railui

%% #Set : UI set parameters with dependencies -----------------------------------
elseif comstr(Cam,'set'); [CAM,Cam]=comstr(CAM,4);
 
[CAM,Cam,RunOpt.nodepend]=comstr('-nodepend',[-25 3],CAM,Cam); % no dependency mode (for load process for instance)
 % the nodepend mode must be handle when real dependencies exist (setcompute for exmaple)
UI=initFig(UI);PARAM=v_handle('uo',UI.gf);uo=[];

if carg<=nargin;RO=varargin{carg};carg=carg+1;
else; 
 if isempty(CAM) % generic Setfunction getting 
        [RO,uo,CAM,Cam]=clean_get_uf('getuo',['SetStruct' CAM],obj,evt);
 else RO=struct; % uo remain empty, CAM is specify outside
 end
end
if isempty(RO); RO=struct; end % RO need to be structure for fieldnames
RB.refresh={};RB.nodepend=0;
% if ~isfield(PARAM,'Project');r1j=[];else;r1j=PARAM.Project;end
% if isempty(r1j);railui('initproject');end
toCinCell=cinguj('@toCinCell');
  
%% #SetSlice : define slice parameters - - - - - - - - - - - - - - - - - - - -
if comstr(Cam,'slice');
    
 [r1j,r1,st]=defaultSet(UI,PARAM,RO,'Slice',{});
 j1=0;
 while j1<length(st); j1=j1+1;
  try 
   val=RO.(st{j1});
  catch ME
      return;
  end% xxx to clean

% railui('SetSlice',struct('half',12,'e1',0.5));
  if strcmpi(st{j1},'Elt') % InitFrom MeshParam in model.Stack
     RO=stack_get(RO,'','MeshParam','get');
     RO=feutil('rmfield',RO,'Save','Track_dim','Coarse','Slice');
     
  elseif strcmpi(st{j1},'InitFrom') % Possible reinit from MeshCfg
   Range=railui('PARAM.Range -safe');
   r1=fe_range('ValCfg',Range,struct('MeshCfg',RO.InitFrom));
   RO=sdth.sfield('AddMissing',RO,r1.MeshCfg);
   RO=feutil('rmfield',RO,'InitFrom');
   
  elseif strcmpi(st{j1},'Cfg') % Display current configuration
    R1=railui('PARAM.Range -safe'); 
    R1=fe_range('valCfg',R1,  ...
     struct('MeshCfg',R1.val(1,strcmpi(R1.lab,'MeshCfg'))));R1=R1.MeshCfg;
    if iscell(R1);R1=R1{2};end
    if ~isempty(fieldnames(R1));RO=R1;end

%% #Substructure types -3
  elseif strcmpi(st{j1},'SubType') 
    li={'e1','e2','e3','e4','e5','em','Lbs','Lsc','Lso','a1','a2'};
    PARAM.Slice=delPar(r1,li); % xxx not clean

    DB=UI.DefBut;
    rSub=sdth.sfield('MergeI',DB.Slice.SubP.Gen,DB.Slice.SubP.(['Sub' val]));
    
    [r1j,r1]=sdt_locale('ButReplace', ...
         struct('Remove','none','Root','Slice','SortCol',sprintf('{''%s''}','Sub'),...
                'List',{{'Sub',rSub}}),PARAM.Slice);   
    PARAM.Slice=r1j;
    sdcedit(r1j,st{j1},'value',val); 
    railui('InitSlice');
    
%% #Armement types -3
  elseif strcmpi(st{j1},'ArmType') 
    RO=feval(dyn_mesh('@ArmTypes'),RO);
    RO=rmfield(RO,'ArmType');
    sdcedit(r1j,st{j1},'value',val);
    railui('SetSlice',RO);

%% #Sleeper types -4
  elseif strcmpi(st{j1},'SleeperType') 
    li={'lob','nent',};
    PARAM.Slice=delPar(r1,li);
    switch val
        case 'Monoblock'
          RF=struct;
        case 'Biblock'
          R1=struct('type','double','callback',{{'railui','Set'}},...
               'name','Slice.ArmP.lob','value',840,'level',[3 0],...
               'ToolTip','Block length [mm]');
          R2=struct('type','double','callback',{{'railui','Set'}},...
               'name','Slice.ArmP.nent','value',3,'level',[3 0],...
               'ToolTip','Number of beams for the tie-bar');
          RF=struct('lob',R1,'nent',R2);
    end
    [r1j,r1]=sdt_locale('ButReplace', ...
         struct('Remove','none','Root','Slice','SortCol',sprintf('{''%s''}','Arm'),...
                'List',{{'Arm',RF}}),PARAM.Slice);   
    sdcedit(r1j,st{j1},'value',val);
    PARAM.Slice=r1j;
    railui('InitSlice');
    
%% #Rail Names -4
  elseif strcmpi(st{j1},'RailName') 
    RO=feval(dyn_mesh('@ArmTypes'),RO);
    RO=rmfield(RO,'RailName');
    sdcedit(r1j,st{j1},'value',val);
    railui('SetSlice',RO);

%% #Sleeper Names -4
  elseif strcmpi(st{j1},'SleeperName') 
    RO=feval(dyn_mesh('@ArmTypes'),RO);
    RO=rmfield(RO,'SleeperName');
    sdcedit(r1j,st{j1},'value',val);
    railui('SetSlice',RO);
    
  
  elseif strcmpi(st{j1},'MatTable')
  %% #SetSliceMatTable : display material table -3
  if strcmpi(val,'close') % close table and maybe save
      st=questdlg('Save modifications to the material database?','Saving material properties?', ...
              'Close without saving','Save and close','Cancel close', ...
              'Close without saving');
      switch st
      case 'Close without saving'; 
        return;
      case 'Save and close';  
        pl=evt.Source.Data; pl(:,[2 5])=pl(:,[2 5])*1e6;
        pl=[pl(:,1) ones(size(pl,1),1)*fe_mat('m_elastic','SI',1) pl(:,2:end)];
        PARAM=stack_set(PARAM,'info','mattab',pl); % save tab
        return;
      case 'Cancel close'; 
      end
  end
  [pl,PARAM]=dyn_mesh('MatTab',PARAM);
  pl2=pl(:,[1 3:end]); % suppress type
  pl2(:,[2 5])=pl2(:,[2 5])/1e6; % MPa
  
  uif=uifigure('Name','Material table');
  uit=uitable(uif,'Data',pl2,'DeleteFcn',{@railui,'SetSlice',struct('MatTable','close')},...
      'ColumnEditable',[false true true true true true],'Position',[0 0 5*70+45+173 400],...
      'ColumnName',{'MatID','E [MPa]','nu [-]','rho [kg/m3]','G [MPa]','eta [-]'},...
      'ColumnWidth',{40, 70, 70, 70, 70, 70},...
      'RowName',{'Ballast';'Ballast shoulder';'UnderLayer';'FormLayer';'Embankment';'Soil';'Sleeper block';'Sleeper bar';'Rail'});
    
    
  elseif strcmpi(st{j1},'MatVary')
  %% #SetSliceMatVary: matvary -3
  R1=fe_def('cleanEntry',r1j);
  st=['Mv' num2str(R1.MatID) R1.MatType];
  r2=R1.MatVal;
  if r2<0 % clear
    r1j=delPar(r1,st);
    PARAM.Slice=r1j;
    railui('InitSlice');
  else
    railui('SetSlice',struct(st,r2));
  end
    
  elseif comstr(lower(st{j1}),'mv')
  %% #MV fields -3
  st2=setdiff(fieldnames(RO),{'TrackType','Save','Track_dim','Cfg'});
  i1=find(~cellfun(@isempty,regexpi(st2,'^Mv\d+','once')));
  if ~isempty(i1) % Missing generic
   R2=sdth.sfield('GetSpecified',RO,st2(i1));
   R3=struct;
   for j1=1:length(i1);
       R3.(R2{j1,1})= ...
           struct('type','string','value',R2{j1,2}, ...
           'name',sprintf('Slice.SubP.%s',R2{j1,1}),'level',[1 0]);
   end 
   [r1j,r1]=sdt_locale('ButReplace', ...
      struct('Remove','none','Root','Slice','SortCol','{''Mv''}', ...
      'List',{{'Mv',R3}}),r1j);   PARAM.Slice=r1j;
   st2(i1)=[];
   railui('InitSlice');
   
  end
  st=st2;
    
  elseif strcmpi(st{j1},'Coarse')
  %% #SetSliceCoarse-3
    fprintf('Coarse mesh in progress...\n');
    R1=cleanSliceCfg(r1j,UI);
    R1.FileName='Visu'; R1.getcoarse=1;
    m_coarse=dyn_mesh('Slice-coarse',R1);
    dyn_post('ViewCoarse',m_coarse,R1);
    fprintf('...done.\n');
    
  elseif strcmpi(st{j1},'View')
  %% #SetSliceView-3
    fprintf('Slice view in progress...\n');
    R1=cleanSliceCfg(r1j,UI);
    if R1.PML;  R1=addPML(R1);  end
    ms=dyn_mesh('Slice',R1);
    dyn_post('ViewSlice',ms);
    fprintf('...done.\n');

    
  elseif strcmpi(st{j1},'SliceCfg')
  %% #SetSliceSliceCfg: add slicecfg to stack -3
    if ~isfield(PARAM,'Range'); railui('InitRange'); end;
    RO=fe_def('cleanEntry',r1j);
    R1=cleanSliceCfg(r1j,UI);
    if R1.PML;  R1=addPML(R1);  end
    li={'Slice', R1;
        'Reduce', struct('RedType',RO.RedType)};
    Range=stack_get(PARAM,'info','Range','get');
    if isempty(Range); Range=struct('SliceCfg',{{}}); end;
    if ~isfield(Range,'SliceCfg'); Range.SliceCfg={}; end;
    if isempty(Range.SliceCfg); ind=0;
    else; ind=any(cellfun(@(x) comstr(x,RO.Name),Range.SliceCfg(:,1)));
    end
    if ~ind
      Range.SliceCfg=[Range.SliceCfg; {R1.Name, li}];
      PARAM=stack_set(PARAM,'info','Range',Range);
    
      [r2j,r2,st]=defaultSet(UI,PARAM,struct(),'Range',{});
      R3=dyn_solve('Range',Range); R3=R3.param.SliceCfg; 
      R3=struct('type','push','callback',{{'railui','Set'}},...
              'name',sprintf('Range.Slices.%s',R1.Name),...
              'value','Delete','level',[2 1],'ToolTip','Delete this configuration');
      R3=struct(R1.Name,R3);
      [r2j,r2]=sdt_locale('ButReplace', ...
              struct('Remove','none','Root','Range','SortCol','{''SliceCfg''}',...
                'List',{{'SliceCfg',R3}}),PARAM.Range);   
       PARAM.Range=r2j;
       railui('InitRange');
     else
        warning('%s already in stack.',RO.Name);
     end

  else % default
      sdcedit(r1j,st{j1},'value',val);
  end
 end
%% #SetTrack (combination of slices) - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'track');
  if ~isfield(PARAM,'Track');PARAM.Track=[];end
  
 [r1j,r1,st]=defaultSet(UI,PARAM,RO,'Track',{});
 j1=0;
 while j1<length(st); j1=j1+1;
   val=RO.(st{j1});
%% #SetTrackTrackCfg: save trackcfg -3
  if strcmpi(st{j1},'TrackCfg')
    
    if ~isfield(PARAM,'Range'); railui('InitRange'); end;  
    RO=fe_def('cleanEntry',r1j);
    RO=rmfield(RO,{'TrackCfg','View'});
    Range=stack_get(PARAM,'info','Range','get');
    if isempty(Range); Range=struct('TrackCfg',{{}}); end;
    if ~isfield(Range,'TrackCfg'); Range.TrackCfg={}; end;
    
    if isempty(Range.TrackCfg); ind=0;
    else; ind=any(cellfun(@(x) comstr(x,RO.Name),Range.TrackCfg(:,1)));
    end
    if ~ind
      Range.TrackCfg=[Range.TrackCfg; {RO.Name, RO}];
      PARAM=stack_set(PARAM,'info','Range',Range);
    
      [r2j,r2,st]=defaultSet(UI,PARAM,struct(),'Range',{});
      R3=dyn_solve('Range',Range); R3=R3.param.TrackCfg; 
      R3=struct('type','push','callback',{{'railui','Set'}},...
              'name',sprintf('Range.Tracks.%s',RO.Name),...
              'value','Delete','level',[2 1],'ToolTip','Delete this configuration');
      R3=struct(RO.Name,R3);
      [r2j,r2]=sdt_locale('ButReplace', ...
              struct('Remove','none','Root','Range','SortCol','{''TrackCfg''}',...
                'List',{{'TrackCfg',R3}}),PARAM.Range);   
       PARAM.Range=r2j;
       railui('InitRange');
    else
        warning('%s already in stack.',RO.Name);
    end
    
%% #SetTrackView: view track -3
  elseif strcmpi(st{j1},'View')
    R1=fe_def('cleanEntry',r1j);
    fprintf('Track view in progress...\n');
    PA=railui('paramvh');
    Range=stack_get(PA,'info','Range','get');
    for j1=1:size(Range.SliceCfg,1)
       R2=Range.SliceCfg{1,2}{1,2};
       ms=dyn_mesh('Slice',R2);
       dyn_post('ViewTrack',ms,struct('nb_slices',R1.nb_slices,'Sub',R2.Sub));
    end
    
    fprintf('...done.\n');
    
%% #SetTrackIsHet: heterogeneity handle in track -3
  elseif strcmpi(st{j1},'isHet')
      warning('Not implemented yet');
  end
  
 end
  
PARAM.Track=r1j;
  

%% #SetRange : define range parameters - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'range');
 [r1j,r1,st]=defaultSet(UI,PARAM,RO,'Range',{});
 j1=0;
 while j1<length(st); j1=j1+1;
   val=RO.(st{j1});
%% #SetRangeRun: save trackcfg -3
  if strcmpi(st{j1},'Run')
     if sdtdef('isinteractive')
      st=questdlg('Run on all stored configurations?','DYNAVOIE Run call', ...
              'Yes','No','No');
     else; st='Yes';
     end
     switch st
     case 'Yes';
       Range=stack_get(PARAM,'info','Range','get');
       if ~isempty(Range)
          dyn_solve('RangeLoop -reset',Range); 
          railui('InitPost');
       else
          warning('No configuration saved in Range');
       end
       return;
     case 'No';
       return;
     end
     
%% #SetRangeSlices: View stored SliceCfg -3
  elseif strcmpi(st{j1},'Slices')
   Range=stack_get(PARAM,'info','Range','get');
   disp('      STORED SLICE CONFIGURATIONS     ');

   Cfg=Range.SliceCfg;
   for j1=1:size(Cfg,1)
       fprintf('--------------%s--------------------------- \n',Cfg{j1,1});
       for j2=1:size(Cfg{j1,2})
          fprintf('------%s------- \n',Cfg{j1,2}{j2,1});
          RO=Cfg{j1,2}{j2,2};
          if isfield(RO,'Arm');
              disp(rmfield(RO,{'Arm','Sub'}));
              disp('      --Arm--'); disp(RO.Arm);
              disp('      --Sub--'); disp(RO.Sub);
          else
              disp(RO);
          end
       end
   end
      
%% #SetRangeTracks: View stored TrackCfg -3
  elseif strcmpi(st{j1},'Tracks')
   Range=stack_get(PARAM,'info','Range','get');
   disp('      STORED TRACK CONFIGURATIONS     ');
   Cfg=Range.TrackCfg;
   for j1=1:size(Cfg,1)
       fprintf('--------------%s--------------------------- \n',Cfg{j1,1});
       disp(Cfg{j1,2});
   end

%% #SetRangeSlices: View stored SimuCfg -3
  elseif strcmpi(st{j1},'Simulations')
  Range=stack_get(PARAM,'info','Range','get');
  disp('      STORED SIMU CONFIGURATIONS     ');

   Cfg=Range.SimuCfg;
   for j1=1:size(Cfg,1)
       fprintf('--------------%s--------------------------- \n',Cfg{j1,1});
       for j2=1:size(Cfg{j1,2})
          fprintf('------%s------- \n',Cfg{j1,2}{j2,1});
          RO=Cfg{j1,2}{j2,2};
          disp(RO);
       end
   end

%% #SetRangeDelete: Delete a configuration -3
  elseif strcmpi(val,'Delete')
    st3=questdlg(sprintf('Are you sure to delete %s?',st{j1}),'Delete configuration', ...
              'Yes','No','No');
     switch st3
     case 'Yes'; 
       Range=stack_get(PARAM,'info','Range','get');
       li2=fieldnames(Range);
        for j2=1:length(li2)
           st2=li2{j2};
           for j3=1:size(Range.(st2),1)
              try
               if strcmpi(st{j1},Range.(st2){j3,1})
                 Range.(st2)(j3,:)=[];
               end
              end
           end
           if size(Range.(st2),1)==0; Range=rmfield(Range,st2); end;
        end
        
        PARAM=stack_set(PARAM,'info','Range',Range);
        
        PARAM.Range=delPar(r1,st{j1});
        railui('InitRange');
     case 'No';
       return;
     end
  end
 end

%% #SetSimu : define simulation parameters - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'simu');
 [r1j,r1,st]=defaultSet(UI,PARAM,RO,'Simu',{});
 j1=0;
 while j1<length(st); j1=j1+1;
   val=RO.(st{j1});
%% #SetSimuCfg: save Simucfg -3
  if strcmpi(st{j1},'SimuCfg')
    if ~isfield(PARAM,'Range'); railui('InitRange'); end;  
    RO=fe_def('cleanEntry',r1j);
    RO=rmfield(RO,{'SimuCfg'});
    if comstr(RO.SimuType,'Vehicle')
       R1=struct('kc',RO.kc,'v0',RO.v0,'VehType',RO.VehicleType);
       R2=rmfield(RO,{'SimuType','VehicleType','kc','v0'});
       li={'Vehicle',R1;
           'Time',R2}; 
    else
       li={'Impact',rmfield(RO,{'SimuType','tf','dt','NeedV','NeedA'});
           'Time',rmfield(RO,{'SimuType','xload','weight','tlen'})};
    end
    Range=stack_get(PARAM,'info','Range','get');
    if isempty(Range); Range=struct('SimuCfg',{{}}); end;
    if ~isfield(Range,'SimuCfg'); Range.SimuCfg={}; end;
    
    if isempty(Range.SimuCfg); ind=0;
    else; ind=any(cellfun(@(x) comstr(x,RO.Name),Range.SimuCfg(:,1)));
    end
    if ~ind
      Range.SimuCfg=[Range.SimuCfg; {RO.Name, li}];
      PARAM=stack_set(PARAM,'info','Range',Range);
    
      [r2j,r2,st]=defaultSet(UI,PARAM,struct(),'Range',{});
      R3=dyn_solve('Range',Range); R3=R3.param.SimuCfg; 
      R3=struct('type','push','callback',{{'railui','Set'}},...
              'name',sprintf('Range.Tracks.%s',RO.Name),...
              'value','Delete','level',[2 1],'ToolTip','Delete this configuration');
      R3=struct(RO.Name,R3);
      [r2j,r2]=sdt_locale('ButReplace', ...
              struct('Remove','none','Root','Range','SortCol','{''SimuCfg''}',...
                'List',{{'SimuCfg',R3}}),PARAM.Range);   
       PARAM.Range=r2j;
       railui('InitRange');
    else
        warning('%s already in stack.',RO.Name);
    end
    
    
  %% #SimuType -3
  elseif strcmpi(st{j1},'SimuType')
    sdcedit(r1j,st{j1},'value',val);
    li={'xload','weight','tlen'}; PARAM.Simu=delPar(r1,li); % reset optional par
    switch val % add optional par
    case 'None'
      RF=struct;
      
    case 'Static'
      R1=struct('type','double','callback',{{'railui','Set'}},...
              'name','Simu.Static.xload','value',10,'level',[2 1],...
              'ToolTip','Position of the load [m]');
      R2=struct('type','double','callback',{{'railui','Set'}},...
              'name','Simu.Static.weight','value',10,'level',[2 1],...
              'ToolTip','Weight of the load [t]');
      RF=struct('xload',R1,'weight',R2);
      
     case 'Impact'
      R1=struct('type','double','callback',{{'railui','Set'}},...
              'name','Simu.Impact.xload','value',10,'level',[2 1],...
              'ToolTip','Position of the load [m]');
      R2=struct('type','double','callback',{{'railui','Set'}},...
              'name','Simu.Impact.weight','value',10e-3,'level',[2 1],...
              'ToolTip','Weight of the load [t]');
      R3=struct('type','double','callback',{{'railui','Set'}},...
              'name','Simu.Impact.tlen','value',1e-2,'level',[2 1],...
              'ToolTip','Duration of the impact [s]');
      RF=struct('xload',R1,'weight',R2,'tlen',R3);
     case 'Vehicle'        
      R1=struct('type','pop','callback',{{'railui','Set'}},...
              'name','Simu.Vehicle.VehicleType','value',3,'level',[2 1],...
              'choices',{{'MovingLoad','MCK3','Axle','Bogie'}},...
              'ToolTip','Available vehicles');
      R2=struct('type','double','callback',{{'railui','Set'}},...
              'name','Simu.Vehicle.kc','value',1e9,'level',[2 1],...
              'ToolTip','Contact stiffness');
      R3=struct('type','double','callback',{{'railui','Set'}},...
              'name','Simu.Vehicle.v0','value',80,'level',[2 1],...
              'ToolTip','Vehicle speed [km/h]');
      RF=struct('VehicleType',R1,'kc',R2,'v0',R3);
    end
    
    [r1j,r1]=sdt_locale('ButReplace', ...
         struct('Remove','none','Root','Simu','SortCol',sprintf('{''%s''}',val),...
                'List',{{val,RF}}),PARAM.Simu);   
    PARAM.Simu=r1j;
    railui('InitSimu');

  end
  
 end
 
 
%% #SetPost : define post parameters - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'post');
 [r1j,r1,st]=defaultSet(UI,PARAM,RO,'Post',{});
 j1=0;
 while j1<length(st); j1=j1+1;
   val=RO.(st{j1});
%% #SetPostViewSlice: View Slice -3
   if strcmpi(st{j1},'ViewSlice')
       if ~isfield(PARAM,'ms'); error('No slice model computed'); end;
       dyn_post('ViewSlice',PARAM.ms);  
       
   %% #SetPostViewTrack: View Track -3
   elseif strcmpi(st{j1},'ViewTrack')
       if ~isfield(PARAM,'mt'); error('No track model computed'); end;
       dyn_post('ViewTrack',PARAM.mt); 
   end
 end

%% #SetGeo - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'geo');
  if ~isfield(PARAM,'Geo');PARAM.Geo=[];end
  UI=initFig(UI);PARAM=v_handle('uo',UI.gf); 
  mt=PARAM.mt; 
  if isempty(mt)&&isfield(PARAM,'m_slice')
   mt=PARAM.ms;
  end
  if isempty(mt)
   PARAM.ms=dyn_mesh('track',fe_def('cleanentry',PARAM.Slice));
   mt=PARAM.ms;
  end
  PARAM.Geo=railui('CbGeo',mt);
  
  if 1==1
  elseif nargout==0; 
   table=sdt_dialogs('uatable -end0','info','References',ref); 
   ua=struct('table',{table},'name','Geo', ...
            ...'ExploName','PTree.Simu',...
             'ToolTip', 'mesh references','NeedClose',1);  
   ua.ColWidth=[90 100 -1];
   [u0,gf]=cinguj('tabbedpaneAdd','railui',ua);
  else;out=RO;
  end
  
%% #SetSimu : SimuCfg (time, vehicle) - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'simu');

[r1j,r1,st]=defaultSet(UI,PARAM,RO,'Simu',{'none'}); 
% xxx RB.refresh='railui(''InitSimu'')';

if 1==2
end
if any(strcmpi(st,'cfg'))
 R1=PARAM.Range; 
 R1=fe_range('valCfg',R1,  ...
     struct('SimuCfg',R1.val(1,strcmpi(R1.lab,'SimuCfg'))));
 R1=R1.SimuCfg;
 % TreeTable view Veh.type='MovingLoad','traj'
 % Set with veh_m=1e5 Time_NeedV=1
 %  TimeOpt.NeedA, NeedV, dt, tend
 dbstack; keyboard    
end
%'x0','xf','dt','dx','v0','PKi','PKf'
  
PARAM.Simu=r1j;

%% #SetCompute - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  elseif comstr(Cam,'compute');
      
   error('Obsolete should be based on range');
   if ~isfield(PARAM,'Compute');PARAM.Compute=[];end
   if isempty(PARAM.Compute)
    [r1j,r1]=cinguj('objEditJ',UI.DefBut.Compute,UI.gf);
    RB.refresh='railui(''InitCompute'')';
   else; r1j=PARAM.Compute;r1=cinguj('objToStruct',r1j); 
   end
  
  st1=fieldnames(RO);
  for j1=1:length(st1)
   st3=st1{j1};
   sdcedit(r1j,st3,RO.(st3));
   if ~RunOpt.nodepend % dependencies
    if comstr(st3,'All_')
     st2=fieldnames(r1); st5=st1{j1}(5:end); 
     for j2=1:length(st2)
      if comstr(fliplr(st2{j2}),fliplr(st5))&&~comstr(st2,'All_')
       sdcedit(r1j,st2{j2},RO.(st1{j1}));
      end
     end
     if comstr(st1{j1},'All_save') && RO.All_save==1; %compute also
      railui('SetCompute',struct('All_comp',1))
     end
     if comstr(st1{j1},'All_repo') && RO.All_repo==1; %compute also
      railui('SetCompute',struct('All_comp',1))
     end
     if comstr(st1{j1},'All_comp') && RO.All_comp==0; %no save and no repo
      railui('SetCompute',struct('All_save',0,'All_repo',0))
     end
    elseif (comstr(fliplr(st3),'evas') || comstr(fliplr(st3),'oper')) && RO.(st3)==1 
     % impose computation if asked for saving or report
     st5=[st1{j1}(1:end-4) 'comp'];
     sdcedit(r1j,st5,1);
    elseif comstr(fliplr(st3),'pmoc_') && RO.(st3)==0 
     % no save if not computed
     st5=[st1{j1}(1:end-4) 'save']; 
     sdcedit(r1j,st5,0);
     % no report if not computed
     st5=[st1{j1}(1:end-4) 'repo'];
     sdcedit(r1j,st5,0);
    end
   end
  end
  
  PARAM.Compute=r1j;

%% #SetResult : contextual tab associated to results - - - - - - - - - - - - -
  elseif comstr(Cam,'result'); [CAM,Cam]=comstr(CAM,7);
   if ~isfield(PARAM,'Result');PARAM.Result=[];end
   if isempty(PARAM.Result)
    [r1j,r1]=cinguj('objEditJ',UI.DefBut.Result,UI.gf);
   else; r1j=PARAM.Result;r1=cinguj('objToStruct',r1j); 
   end
   RB.refresh='railui(''InitResult'')';   
   
   if ~isempty(CAM) % set name that is used in callbacks of the tab
    sdcedit(r1j,'name','value',CAM); % 4 args : low level put for enabled off
   end
   
   st1=fieldnames(RO);
   for j1=1:length(st1)
    st3=st1{j1};
    sdcedit(r1j,st3,RO.(st3));
   end
  
   PARAM.Result=r1j;
%    if isfield(RB,'refresh') && ~isempty(RB.refresh) % factorized below !
%     eval(RB.refresh)
%    end

elseif comstr(Cam,'ui')
 %% #SetUI: Update UI of railui with data generated externally -2
 UIc=RO;
 if isequal(UIc.gf,UI.gf); UI=UIc;
 elseif isempty(UIc.gf) % Case for project loading (no merging of defbut yet)
 else; error('Attempt at SetUI with bad UI.gf')
 end
 return;
 
elseif comstr(Cam,'wd')
 %% #SetWD: eventually mkdir and point interface to wd -2
 wd=comgui('cd',sdtu.f.safe(RO));
 if ~exist(wd,'dir');mkdir(wd); end
 if ~exist(fullfile(wd,'plots'),'dir');mkdir(fullfile(wd,'plots'));end
 railui('SetProject',struct('ProjectWd',wd,'PlotWd',fullfile(wd,'plots'),...
                            'Report',fullfile(wd,'Report_V1.docx'),...
                            'root','V1'));
 if nargout>0;out=wd;end

elseif comstr(Cam,'rail')
 sdtroot('setMesh',RO,UI)
%% #SetEnd : propagate to sdtroot 
else;  %sdtroot(['set' CAM],RO,UI)
 try;
  if ~isempty(obj)&&isfield(uo,'EditFcn'); % composite call
   sdtroot(varargin{:}); % transparent forward to keep full context
  else; sdtroot(sprintf('Set%s',CAM),RO,UI); % standard Set to a field in PARAM through UI
  end
 catch ERR; rethrow(ERR);
 end
end


% Fill RB.Refresh in the Set command if reinit of tab is needed
if ~isempty(RB.refresh)&&~RB.nodepend;
 if ischar(RB.refresh);eval(RB.refresh);
 else;
  for j1=1:size(RB.refresh,1); % obj=UI, evt='refresh'
   feval(RB.refresh{j1,1},UI,'refresh',RB.refresh{j1,2:end});
  end
 end
elseif ~isempty(uo);cingui('resize');
end

%% #init : intialisation of GUI ----------------------------------------------
elseif comstr(Cam,'init'); [CAM,Cam]=comstr(CAM,5);

if isempty(CAM);UI=initFig(UI,varargin{carg});return;end% First init
UI=initFig(UI);PARAM=v_handle('uo',UI.gf);

%% #InitProject : initialize project tab - - - - - - - - - - - - - - - - - - - 
if comstr(Cam,'project'); [CAM,Cam]=comstr(CAM,8);
%   sdtroot('InitProject');

 if ~isfield(PARAM,'Project')||isempty(PARAM.Project);
   feval(UI.MainFcn,'setproject',[],UI);
 end
 % generate tab
 % simple 3 columns table
 ua=struct('table', ...
  {sdt_dialogs('uatable-end0','info','Project',PARAM.Project)}, ...
  'name','Project','ExploName','PTree.Project','ToolTip', ...
  'Project properties','NeedClose',1, ...
  'ColumnName',{{'Param','Value','info'}}, ...
  'ColWidth',[70 -1 180],'RowHeight',20);
 ua.jProp={'jt_setTableHeader',[]}; 

 % Keep ordering of base fields
 [i1,i2]=ismember({'ProjectWd','PlotWd','Report','LastWd',...
  'root','name','Description'},ua.table(:,1));
 ua.table(sort(i2(i1)),:)=ua.table(i2(i1),:);
 % Description box is multiline:
 i1=find(strcmpi('Description',ua.table(:,1)));
 if ~isempty(i1);ua.table{i1,2}.put('type','stringa'); end%Multi-line
 [u0,gf]=cinguj('tabbedpaneAdd',UI.gf,ua);
 if ~isempty(i1); u0.JTable.setRowHeight(i1-1,100);end
 
 if comstr(Cam,'-reset')
    wd=fileparts(which('t_dynavoie')); 
    railui('SetProject;',struct('ProjectWd',wd, ... % Root file location
      'PlotWd',fullfile(wd,'plots'), ...          % Plot directory
      'Report',fullfile(wd,'Dynavoie_report.docx')));
  railui('InitProject');
 end

%% #InitRange : Range configuration handling - - - - - - - - - - - - - - - - - - - 
elseif comstr(Cam,'range'); [CAM,Cam]=comstr(CAM,6);

 if ~isfield(PARAM,'Range')||isempty(PARAM.Range)||comstr(Cam,'reset') % Init if necessary
   DB=UI.DefBut;

   r1=sdt_locale('butMerge',struct('Root','Range','List',...
       {{'SliceCfg',DB.Range.SliceCfg;
         'TrackCfg',DB.Range.TrackCfg;
         'SimuCfg',DB.Range.SimuCfg;
         'Gen',DB.Range.Gen;}}));
   PARAM.Range=cinguj('ObjEditJ',r1,UI.gf);
 end
 r1j=PARAM.Range;
 ua=struct('name','Range','IntegerHandle','off','ToolTip','Range configuration','NeedClose',1, ...
      'ParentPannel',UI.gf,'ExploName','PTree.Range');
 ua=sdt_dialogs('uatable-end0 -tipAsCol',ua,r1j);
 ua.setSort=4; ua=rmfield(ua,'groupable');
 % Sort names
 ua=sortUInames(ua,{'Slices','Tracks','Simulations','Gen'},...
     {'SliceCfg','TrackCfg','SimuCfg','Gen'});
 
 % suppress SortCols
 ua.table(:,3:end)=[];
 % tab data implementation
 % store data in java pointer, then display the tab
 [ub,ga]=cinguj('TabbedPaneAdd',UI.gf,ua); 
 ub.JTable.setRowSelectionAllowed(false);
 ub.JTable.setColumnSelectionAllowed(false);
 ub.JTable.setCellSelectionEnabled(false);
 ub.JTable.setRowHeight(20)


        
%% #InitSlice : geometry & reduction of a slice - - - - - - - - - - - - - - - - - - - 
elseif comstr(Cam,'slice'); [CAM,Cam]=comstr(CAM,6);

  if ~isfield(PARAM,'Slice')||isempty(PARAM.Slice)||comstr(Cam,'reset') % Init if necessary
   DB=UI.DefBut;
   
   curSub=DB.Slice.SubP.Gen.SubType.choices{DB.Slice.SubP.Gen.SubType.value};
   rSub=sdth.sfield('MergeI',DB.Slice.SubP.Gen,DB.Slice.SubP.(['Sub' curSub]));

   r1=sdt_locale('butMerge',struct('Root','Slice','List',{ ...
     {'Cfg',DB.Slice.Cfg;...
      'Arm',DB.Slice.ArmP;...
      'Sub',rSub;...
      'Reduction',DB.Slice.Red;...
      'Generic',DB.Slice.Gen}}));
   PARAM.Slice=cinguj('ObjEditJ',r1,UI.gf);
  end
  r1j=PARAM.Slice;%r1=cinguj('objToStruct',r1j)
  ua=struct('name','Slice','IntegerHandle','off','ToolTip','Slice configuration','NeedClose',1, ...
      'ParentPannel',UI.gf,'ExploName','PTree.Slice');
  ua=sdt_dialogs('uatable-end0 -tipAsCol',ua,r1j);
  %ua.groupable=3;
  ua.setSort=4; ua=rmfield(ua,'groupable');
  
  ua=sortUInames(ua,{'Name','ArmType','SubType','RedType','Generic'},...
     {'Cfg','Arm','Sub','Reduction','Generic'});

  ua.table(:,3:end)=[];
  
  % tab data implementation
  % store data in java pointer, then display the tab
  [ub,ga]=cinguj('TabbedPaneAdd',UI.gf,ua); 
   ub.JTable.setRowSelectionAllowed(false);
   ub.JTable.setColumnSelectionAllowed(false);
   ub.JTable.setCellSelectionEnabled(false);
   ub.JTable.setRowHeight(20)
%    % Set default parameters
%    railui('SetSlice',struct('SubType','Simple'));
%    railui('SetSlice',struct('SleeperType','Monoblock'));

%% #InitTrack : track generation - - - - - - - - - - - - - - - - - - - 
elseif comstr(Cam,'track')

  if ~isfield(PARAM,'Track')||isempty(PARAM.Track)||contains(Cam,'reset') % Init if necessary
   DB=UI.DefBut;

   r1=sdt_locale('butMerge',struct('Root','Track','List',{ ...
     {'Cfg',DB.Track.Cfg;...
      'Generic',DB.Track.Gen}}));
   PARAM.Track=cinguj('ObjEditJ',r1,UI.gf);
  end
  r1j=PARAM.Track;%r1=cinguj('objToStruct',r1j)
  ua=struct('name','Track','IntegerHandle','off','ToolTip','Track configuration','NeedClose',1, ...
      'ParentPannel',UI.gf,'ExploName','PTree.Track');
  ua=sdt_dialogs('uatable-end0 -tipAsCol',ua,r1j);
  %ua.groupable=3;
  ua.setSort=4; ua=rmfield(ua,'groupable');
  ua.table(:,3:end)=[];
  % tab data implementation
  % store data in java pointer, then display the tab
  [ub,ga]=cinguj('TabbedPaneAdd',UI.gf,ua); 
   ub.JTable.setRowSelectionAllowed(false);
   ub.JTable.setColumnSelectionAllowed(false);
   ub.JTable.setCellSelectionEnabled(false);
   ub.JTable.setRowHeight(20)
   
%% #InitSimu : simulation generation - - - - - - - - - - - - - - - - - - - 
elseif comstr(Cam,'simu')

  if ~isfield(PARAM,'Simu')||isempty(PARAM.Simu)||contains(Cam,'reset') % Init if necessary
   DB=UI.DefBut;

   r1=sdt_locale('butMerge',struct('Root','Simu','List',{ ...
     {'Cfg',DB.Simu.Cfg}}));
   PARAM.Simu=cinguj('ObjEditJ',r1,UI.gf);
  end
  r1j=PARAM.Simu;%r1=cinguj('objToStruct',r1j)
  ua=struct('name','Simu','IntegerHandle','off','ToolTip','Simulation configuration','NeedClose',1, ...
      'ParentPannel',UI.gf,'ExploName','PTree.Simu');
  ua=sdt_dialogs('uatable-end0 -tipAsCol',ua,r1j);
  %ua.groupable=3;
  ua.setSort=4; ua=rmfield(ua,'groupable');
  ua.table(:,3:end)=[];
  % tab data implementation
  % store data in java pointer, then display the tab
  [ub,ga]=cinguj('TabbedPaneAdd',UI.gf,ua); 
   ub.JTable.setRowSelectionAllowed(false);
   ub.JTable.setColumnSelectionAllowed(false);
   ub.JTable.setCellSelectionEnabled(false);
   ub.JTable.setRowHeight(20)
   
  
%% #InitPost : initialize track tab - - - - - - - - - - - - - - - - - - - 
elseif comstr(Cam,'post')

  ua=struct('table',{{}}, ...
      'name','Post','ExploName','PTree.Post','ToolTip', ...
      'Available post-processing','NeedClose',1, ...
      'ColumnName',{{'','','',''}}, ...
      'ColWidth',[100 -1 40 40],'RowHeight',20,'setSort',4);
  % model views
  if isfield(PARAM,'mt')
   mt=PARAM.mt;
   ua.table(end+1,1:2)={'mt','Track model'}; ua.level(size(ua.table,1),1)=1;
   ua.table{end,3}=struct('name','ViewTrack','type','push',...
         'value','Track','callback','dyn_post(''ViewTrack'')');
     
   SE=stack_get(mt,'SE');SE=SE(:,[2 2 2 3]);
   r1=struct('name','ViewCoarse','type','push','value','Coarse', ...
       'callback','');
   r2=struct('name','ViewSlice','type','push','value','Slice', ...
       'callback','');
      for j1=1:size(SE,1)
    r1.callback=sprintf('dyn_post(''ViewCoarse %s'')',SE{j1,1});SE{j1,2}=r1;
    r2.callback=sprintf('dyn_post(''ViewSlice %s'')',SE{j1,1});SE{j1,3}=r2;
   end
   r2=fesuper('s_',mt);
   if ~isempty(r2)
     for j1=1:size(SE,1)
       i1=vertcat(r2{strcmpi(r2(:,1),SE{j1,1}),3})-100;
       if ~isempty(i1);SE{j1,5}=sprintf('slices %i:%i',i1([1 end]));end
     end
   end
   r1=SE(:,[1 5:end 2 3]);
   ua.table(end+(1:size(SE,1)),1:size(r1,2))=r1;
   ua.level(end+1:size(ua.table,1),1)=2;
  end
  sel=stack_get(PARAM,'sel');
  if ~isempty(sel)
   ua.table{end+1,1}='Sel';   ua.level(size(ua.table,1),1)=1;
   r1=sel(:,[2 1]);
   if isfield(PARAM,'dstat')
    r2=struct('type','push','value','feplot', ...
        'callback','dyn_post(''ViewRest'')');
    r1(:,3)={r2};
    r2=struct('type','push','value','iiplot', ...
        'callback','dyn_post(''PostSens'')');
    r1(:,4)={r2};
   end

   ua.table(end+(1:size(sel,1)),1:size(r1,2))=r1;
   ua.level(end+1:size(ua.table,1),1)=2;
  end
  % Static
  if isfield(PARAM,'dstat')&&~isempty(PARAM.dstat);
   ua.table{end+1,1}='Static';   ua.level(size(ua.table,1),1)=1;
   ua.table{end+1,1}='Animation';  ua.level(size(ua.table,1),1)=2;
     ua.table{end,2}='Animated 3D view of static deformation'; 
   ua.table{end,3}=struct('name','ViewRestStat','type','push',...
         'value','Anim','callback','dyn_post(''ViewRestAll -dstat'')');
  end
  
  if isfield(PARAM,'def')&&~isempty(PARAM.def);
   if size(PARAM.def.Xlab{2},2)==2; % vehicle
     ua.table{end+1,1}='Vehicle';   ua.level(size(ua.table,1),1)=1;
     ua.table{end+1,1}='Animation';  ua.level(size(ua.table,1),1)=2;
     ua.table{end,2}='Animated 3D view of the pass-by'; 
     ua.table{end,3}=struct('name','ViewRest','type','push',...
         'value','Anim','callback','dyn_post(''ViewRestAll -light'')');
   else % impact
     ua.table{end+1,1}='Impact';   ua.level(size(ua.table,1),1)=1;
     ua.table{end+1,1}='Receptance';  ua.level(size(ua.table,1),1)=2;
     ua.table{end,2}='Plot receptance'; 
     ua.table{end,3}=struct('name','ViewRecept','type','push',...
         'value','Plot','callback','dyn_post(''PostRecept'')');
     ua.table{end+1,1}='Animation';  ua.level(size(ua.table,1),1)=2;
     ua.table{end,2}='Animated 3D view of the impact'; 
     ua.table{end,3}=struct('name','ViewRest','type','push',...
         'value','Anim','callback','dyn_post(''ViewRestAll -light'')');
   end
  end
  % periodic
  if isfield(PARAM,'perResp')&&~isempty(PARAM.perResp);
   ua.table{end+1,1}='Periodic response';   ua.level(size(ua.table,1),1)=1;
   ua.table{end+1,1}='2D-DFT';  ua.level(size(ua.table,1),1)=2;
     ua.table{end,2}='Plot space-time Fourier Transform'; 
     ua.table{end,3}=struct('name','View2DDFT','type','push',...
         'value','Plot','callback','dyn_post(''Post2DDFT'')');
  end
  
  [u0,gf]=cinguj('tabbedpaneAdd','railui',ua);
  sdtTreeTable('expandrow',{UI.gf,'Post'},0,true);

%% #InitMeshCfg : generic display of MesCfg - - - - - - - - - - - - - - - - - - - 
elseif comstr(Cam,'meshcfg')

if ~isfield(PARAM,'Range')|| ...
        ~isfield(PARAM.Range,'param')|| ~isfield(PARAM.Range.param,'MeshCfg')
  warning('PA.Range.param.MeshCfg expected and not found');
  return
end

RO=PARAM.Range.param.MeshCfg;
ua=RO.choices(:)';
ua(2,:)={struct('RootRow', cinguj('ObjToCinCell', ...
      struct('type','push','value','ViewCoarse',...
    'tab',{{'railui','MeshCfg'}},'callback',{{'dyn_post','view'}})))};
for j1=1:size(ua,2);ua{2,j1}.data=RO.data{j1};end
ua=struct(ua{:}); ua=fe_range('uiRO2Table',ua);
ua.name='MeshCfg';ua.ColWidth=[100 80 -1*ones(1,size(ua.table,2)-2)];
ua=sdsetprop(ua,'jProp','setSort',4,'NeedClose',1);
[u0,gf]=cinguj('tabbedpaneAdd','railui',ua);

%fe_range('UICfgTree')

%% #InitGeo : geometry & reduction of a slice - - - - - - - - - - - - - - - - - - - 
elseif comstr(Cam,'geo')

 if ~isfield(PARAM,'Geo')||isempty(PARAM.Geo); railui('SetGeo'); end
 
 ua=struct('name','Geo','ToolTip', 'mesh references','NeedClose',1);  
 ua.ColWidth=[90 100 -1];
 ua=sdt_dialogs('uatable -end0 -tipascol 1',ua,PARAM.Geo); 
 [u0,gf]=cinguj('tabbedpaneAdd','railui',ua);
   
%% #InitViewTable : initialize ViewTable tab - - - - - - - - - - - - - - - - - 
elseif comstr(Cam,'viewtable') % simulation tab
% railui('InitViewTable',table)
% visualize tab, adding Buttons to save/copy the table

 table=varargin{carg}; carg=carg+1;
 
 % Add button for save/copy operation
 b1=struct('type','push','value','Copy as Excel', ...
      'ToolTip','Copy as excel','name','xxx', ...
      'callback','railui(''CBCopyExcel'')');  
 b2=struct('type','push','value','Copy as Tex', ...
      'ToolTip','Copy as Tex','name','xxx', ...
      'callback','railui(''CBCopyLatex'')');  
 r1=struct('b1',b1,'b2',b2);
 r1j=cinguj('ObjEditJ',r1,Peer);
 table{end+1,1}=b1;  table{end,2}=b2; table=table([end 1:end-1],:);

 ua=struct('table',{table},'name','xxx','ToolTip', ...
           'Statfc visu','NeedClose',1);
 [u0,gf]=cinguj('tabbedpaneAdd','railui',ua);
  
%% #InitPTree : initialize project tree - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'ptree'); [Cam,CAM]=comstr(Cam,6);
    
  if carg<=nargin;RunOpt=varargin{carg};carg=carg+1;
  else;RunOpt=struct;
  end
  [CAM,Cam,RunOpt.NoInit]=comstr('-noinit',[-25 3],CAM,Cam);  

  % init empty tree entries :
  if ~isfield(PARAM,'Project');PARAM.Project=[];end
  if ~isfield(PARAM,'Track');PARAM.Track=[];end
    
  % project info
  rp=fe_def('cleanentry',PARAM.Project); 
  r4={}; 
  %if ~isempty(rp) && comstr(rp.Mode,'expert'); r4={UI.DefBut.PTree.DOE}; end
  
  % TREE TABLE:
  ua=struct('table',{ ...
      sdt_dialogs('uatable-end0','info','Info',UI.DefBut.PTree)}, ...
      'level',[],'name','ProjectTree','ToolTip', ...
           'ProjectTree','NeedClose',2);
  %ua.table={
  %    UI.DefBut.PTree.Project;      UI.DefBut.PTree.Slice;
  %    UI.DefBut.PTree.Simu;        UI.DefBut.PTree.Compute};
  ua.table=ua.table(:,2);ua.level=ones(size(ua.table,1),2);
  for j1=1:size(ua.table,1);ua.table{j1}.Peer=UI.gf;end
 
  % Deal with active level in tree:
  if isempty(UI.treeF); RunOpt.lastname='';
  elseif ~isfield(RunOpt,'lastname')||isempty(RunOpt.lastname)
     [RunOpt.lastname,tree]=UI.treeF('explolastname',UI.gf);
  end
  [tree,gf]=cinguj('tabbedpaneAddTree','railui',ua);
  tree.getSelectionModel.setSelectionMode( ...
   javax.swing.tree.TreeSelectionModel.SINGLE_TREE_SELECTION)
  if ~isempty(RunOpt.lastname)&&~RunOpt.NoInit
     node=UI.treeF('scrollToNameSelect',tree,RunOpt.lastname); 
  end
  % 'TList.Peer needs to be set'  u0.table{2,2}.get('parent')
  if nargout==0;
  else; out=ua;
  end
  
%% #InitMt : table with mt content - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'mt')
    
mt=PARAM.mt;
dbstack; keyboard
   
%% #InitMain : open main window - - - - - - - - - - - - - - - - - - - - - -
elseif comstr(Cam,'main')
  railui('InitPTree'); %  Project tree must be done first
  railui('tbrinit') % Initialize the toolboar
  railui('InitProject -reset'); % project field must be filled before any action
%% InitSimu is handled by sdtroot - - - - - - - - - - - - - - - - - - - 
%% #InitEnd  propagate to sdtroot (Project,Range)
elseif ~isempty(CAM);  sdtroot(['init' CAM],UI);
end
   
%% #tbr : Toolbar building commands ------------------------------------------
elseif comstr(Cam,'tbr');[CAM,Cam]=comstr(CAM,4);

if comstr(Cam,'init')
   gf=UI.gf;
   
   %% #MenuBar : deal with the main menu bar - - - -
   delete(findall(gf,'type','uimenu')); % delete existing menu then rebuild
   % File - - -
   r2=uimenu('parent',gf,'label','&File','callback','');
   uimenu('parent',r2,'label','Open project',...
             'callback',{'railui','CbProjectLoad'})
   r2=uimenu('parent',gf,'label','&Help','callback','');
   uimenu('parent',r2,'label','Table of contents',...
             'callback','sdtweb(''dynavoie'')');
   %% #ToolBar : deal with the toolbar icons - - - -
   delete(findall(gf,'tag','FigureToolBar'));
   h=uitoolbar(gf,'handlevisibility','off','tag','taglist_toolbar');
  
   % Refresh - - -
    uipushtool('parent',h, ...
     'HandleVisibility','off','BusyAction','Cancel', ...
     'TooltipString','refresh', ...
     'cdata',sdticons('refresh'), ...
     'clickedcallback','cingui(''Resize'')');
   % New - - -
    uipushtool('parent',h, ...
     'HandleVisibility','off','BusyAction','Cancel', ...
     'TooltipString','New project', ...
     'cdata',sdticons('new'), ...
     'clickedcallback',{'railui','CbProjectClose'}); % cf File/New
    % Load - - -
    uipushtool('parent',h, ...
     'HandleVisibility','off','BusyAction','Cancel', ...
     'TooltipString','Load an existing project', ...
     'cdata',sdticons('open'), ...
     'clickedcallback',{'railui','CbProjectLoad'}); % cf File/Load
   % Save - - -
    uipushtool('parent',h, ...
     'HandleVisibility','off','BusyAction','Cancel', ...
     'TooltipString','Save current project', ...
     'cdata',sdticons('save'), ...
     'clickedcallback',{'railui','CbProjectSave'}); % cf File/Save

   % delete tab
   if isunix
    uipushtool('parent',h,...
    'HandleVisibility','off','BusyAction','Cancel', ...
    'TooltipString','Close current Tab', ...
    'cdata',sdticons('exit'), ...
    'clickedcallback',['uf=clean_get_uf(gcbo);if size(uf.tStack,1)>0;'...
     'i1=uf.tab;obj=uf.JPeer;uf.tStack(i1,:)=[];uf.tab=max(1,uf.tab-1);'...
     'set(uf.JPeer.getClientProperty(''peer''),''userdata'',[],''userdata'',uf);'...
     'obj.remove(i1-1);end'],...
     'separator','on');
   end
  
else;error('Tbr%s unknown',CAM);
end

%% #_vh : obtain v_handle to standard objects --------------------------------
elseif comstr(Cam,'vh');[CAM,Cam]=comstr(CAM,3);

UI=initFig(UI);
PARAM=v_handle('uo',UI.gf);
% if carg<=nargin; r1=varargin{carg}; else; r1=[]; end 

%% #VHRemove - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if comstr(Cam,'remove');[CAM,Cam]=comstr(CAM,7);% remove properly handle
 if isfield(PARAM,CAM)
  %PARAM.(CAM)=[];
  PARAM.(CAM)=railui(sprintf('vh%s',CAM),[]);
  PARAM=feutil('rmfield',PARAM,CAM);
 end

%% #VH - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
elseif ~isempty(Cam)
%   comstr(Cam,'mo')|| comstr(Cam,'def')|| comstr(Cam,'fc')|| ...
%   comstr(Cam,'statfc')|| comstr(Cam,'uplift')|| comstr(Cam,'r2') || ... 
%   comstr(Cam,'elasticity')|| comstr(Cam,'sag')|| comstr(Cam,'tension') 
  % mo(i)
  h=findall(UI.gf,'tag',CAM);
  if isempty(h);
      h=uicontrol('parent',UI.gf,'visible','off','tag',CAM, ...
          'DeleteFcn',{@railui,'close'});
  end 
  if carg<=nargin; r1=varargin{carg}; 
    if isa(r1,'v_handle')&&~ishandle(r1.GetVHandleHandle);
     out=v_handle('uo',h);
    else
     out=v_handle('uo',h,r1);
    end
  else
     out=v_handle('uo',h);
  end

%%
else; out=PARAM;
end

%% #Rep : Utilities for report generation -----------------------------
elseif comstr(Cam,'rep'); [CAM,Cam]=comstr(CAM,4);
 %% #RepWrite : write report line in fprint or word file -----------------
 if comstr(Cam,'write');[CAM,Cam]=comstr(CAM,6);
    li=varargin{carg}; carg=carg+1;
    switch li{1}
        case {'subsection','text'}
          if 1==1
             fprintf(li{2});
          else
            sdtacx('WordRecInsert',struct('type',li{1},...
              'value',li{2},'Project','railui'));
          end
    end
 end

%% #CB : Callback for UI -----------------------------------------------------
elseif comstr(Cam,'cb');[CAM,Cam]=comstr(CAM,3);
 UI=initFig(UI);PARAM=v_handle('uo',UI.gf);
 
 %% #CBGeo : get geometry references ----------------------------------------------------
 if comstr(Cam,'geo');[CAM,Cam]=comstr(CAM,4);
 if nargin==1
     lab={'x_slice_l','Beginning of the slice';
          'x_slice_m','Middle of the slice';
          'x_slice_r','End of the slice';
          'y_track','Middle of the track';
          'y_slice','Middle of the slice';
          'x_scorners','X-corners of the sleeper';
          'x_sleeper_m','Middle of sleeper';
          'y_slice','Middle of the slice';
          'y_scorners','Y-corners of the sleeper';
          'y_rail_m','Middle of the rails';
          'y_rail_ext','Exterior of the rail';
          'y_sleeper_ext','Exterior of the sleeper';
          'z_rail_web','Rail web';
          'z_rail_foot','Rail foot';
          'z_pad_up','Upper face of pad';
          'z_pad_mid','Middle pad';
          'z_sleeper_up','Upper face of sleeper';
          'z_sleeper_mid','Middle of sleeper';
          'z_sleeper_down','Below sleeper'
          'z_soil_up','Upper face of soil'};
      disp(lab);
  else
  RO=varargin{carg};carg=carg+1;
  if ~isfield(RO,'sw'); RO.sw=.6; end;
  if ~isfield(RO,'targL'); RO.targL=0; end;
  if ~isfield(RO,'double'); RO.double=0; end;
  ref={}; % references
  if isfield(RO,'Elt');mt=RO;RO=feval(dyn_mesh('@getParamInStack'),mt) ;end
  %% Substructure geometry parameters
  if isfield(RO,'Arm'); % is a slice
  RA=RO.Arm; RS=RO.Sub; 
 
 if isfield(RS,'m_coarse')
 %% custom mesh
  x_slice_l=min(RS.m_coarse.Node(:,5));
  ref(end+1,1:3)={'x_slice_l',x_slice_l,'Beginning of the slice'};
  ref(end+1,1:3)={'x_slice_m',x_slice_l+RO.sw/2,'Middle of the slice'};
  ref(end+1,1:3)={'x_slice_r',x_slice_l+RO.sw,'End of the slice'};
  if isfield(RO,'double') && RO.double % double track
    y_slice=RO.Sub.Lso;
    y_track=0;
  else % simple track
    y_track=0;
    y_slice=0;
  end
  ref(end+1,1:3)={'y_slice',y_slice,'Middle of the slice'};
  z_slice=max(RS.m_coarse.Node(:,7));
  ref(end+1,1:3)={'z_slice',z_slice,'Top of the slice'};
      
 elseif comstr(lower(RA.ArmType),'slabtrack')
 %% slab track
   warning('Geometry parameters of slab track to be implemented');
   
 else
  %% ballasted track
  % X references
  ref(end+1,1:3)={'x_slice_l',0,'Beginning of the slice'};
  ref(end+1,1:3)={'x_slice_m',RO.sw/2,'Middle of the slice'};
  ref(end+1,1:3)={'x_slice_r',RO.sw,'End of the slice'};
  
  % Y references

  if isfield(RO,'double') && RO.double % double track
    y_slice=RO.Sub.Lso;
    y_track=0;
  else % simple track
    y_track=0;
    y_slice=0;
  end
  ref(end+1,1:3)={'y_track',y_track,'Middle of the track'};
  ref(end+1,1:3)={'y_slice',y_slice,'Middle of the slice'};
    
  % Z references:
   z_slice=0;
   ref(end+1,1:3)={'z_slice',z_slice,'Top of the slice'};
   ref(end+1,1:3)={'z_ballast_up',z_slice,'Upper face of ballast'};
   ref(end+1,1:3)={'z_ballast_mid',z_slice-(RS.e1+RA.hsb)/2,'Middle of ballast'};
   ref(end+1,1:3)={'z_ballast_down',z_slice-(RS.e1+RA.hsb),'Below ballast'};
   if isfield(RS,'SubType')&&comstr(lower(RS.SubType),'simple')
      z_soil_up=-(RS.e1+RA.hsb);
   elseif isfield(RS,'e2')&&isfield(RS,'e3');
      z_soil_up=-(RS.e1+RA.hsb)-RS.e2-RS.e3;
   else; z_soil_up=NaN; 
   end
   ref(end+1,1:3)={'z_soil_up',z_soil_up,'Upper face of soil'};
   
   %   SE=stack_get(mt,'SE','#s.*');
%   if ~isempty(SE);SE(~cellfun('isempty',regexp(SE(:,2),'r$')),:)=[];end % No r
%   if ~isempty(SE);SE(~cellfun('isempty',regexp(SE(:,2),'i\d')),:)=[];end % No i(j)
%   if isempty(SE);
%       if strcmpi(mt.name,'s1');SE={'SE','s1',mt};% Single slice
%       else; SE={'SE','s1',mt};
%       end; 
%   end %No
% 
%   for jElt=1:size(SE,1)
%    % layers Z:
%    MatId=unique(feutil('MatId',SE{jElt,3})); MatId(MatId==0)=[];
%    i1=min(MatId)-1; % MatId Shift
%    MatId(MatId>i1+99)=[]; 
%    st=SE{jElt,2};if strcmpi(st,'slice');st='s1';end
%    for j1=1:length(MatId)
%     n1=feutil(sprintf('GetNode MatId%i',MatId(j1)),SE{jElt,3});
%     ref(end+1,1:3)={sprintf('z_%s_%i_up',st,MatId(j1)),max(n1(:,7)),...
%                        sprintf('Top of layer %i in %s',MatId(j1),st)};
%     ref(end+1,1:3)={sprintf('z_%s_%i_mid',st,MatId(j1)),0.5*(min(n1(:,7))+max(n1(:,7))),...
%                        sprintf('Middle of layer %i in %s',MatId(j1),st)};
%     ref(end+1,1:3)={sprintf('z_%s_%i_down',st,MatId(j1)),min(n1(:,7)),...
%                        sprintf('Bottom of layer %i in %s',MatId(j1),st)};
%    end
%   end
 end
 
 %% Sleepers geometry parameters
 if ~comstr(lower(RA.SleeperType),'none')
  if ~RO.spos
    x_scorners=RA.lab*[-.5 .5]+RO.sw/2;
    x_smid=RO.sw/2;
  else
    x_scorners=[RA.lab/2, RO.sw-RA.lab/2];
    x_smid=0; % arbitrary left reference
  end
  ref(end+1,1:3)={'x_scorners',x_scorners,'Sleeper corners'};
  ref(end+1,1:3)={'x_sleeper_m',x_smid,'Middle of sleeper'}; 
  if RO.half&&~RO.double % half single track
      if comstr(lower(RA.SleeperType),'monobloc')
         y_scorners=RA.Ltr/2;
      elseif comstr(lower(RA.SleeperType),'bibloc')
         y_scorners=[RA.Ltr/2-RA.lob RA.Ltr/2];
      end
  elseif ~RO.half&&~RO.double||RO.double % full single track
    if comstr(lower(RA.SleeperType),'monobloc')
         y_scorners=RA.Ltr/2*[-1 1];
    elseif comstr(lower(RA.SleeperType),'bibloc')
         y_scorners=[-RA.Ltr/2 -RA.Ltr/2+RA.lob RA.Ltr/2-RA.lob RA.Ltr/2];
    end
  end
   y_scorners=y_scorners+y_slice(end);
   ref(end+1,1:3)={'y_scorners',y_scorners,'Y-corners of the sleepers'};
   ref(end+1,1:3)={'y_sleeper_ext',y_scorners(end),'Exterior of the sleeper'};
   ref(end+1,1:3)={'z_sleeper_up',RA.hb-RA.hsb,'Upper face of sleeper'};
   ref(end+1,1:3)={'z_sleeper_mid',RA.hb/2-RA.hsb,'Middle of sleeper'};
   ref(end+1,1:3)={'z_sleeper_down',-RA.hsb,'Bottom face of sleeper'};
 end
  
 %% Rail and pads geometry parameters
  if RO.half&&~RO.double % halftrack
     ref(end+1,1:3)={'y_rail_m',y_slice(end)+RA.gaug/2,'Middle of the rails'};
  elseif ~RO.half&&~RO.double
     ref(end+1,1:3)={'y_rail_m',y_slice(end)+RA.gaug/2*[-1 1],'Middle of the rails'};
  elseif RO.half&&RO.double
     ref(end+1,1:3)={'y_rail_m',y_slice(end)+RA.gaug/2*[-1 1],'Middle of the rails'};
  end
  ref(end+1,1:3)={'y_rail_ext',ref{end,2}(end)+RA.lar2/2,'Exterior of the rail'};

  if comstr(lower(RA.ArmType),'slabtrack')
    z_pad_down=z_slice-RA.hre;
  else
    z_pad_down=RA.hb-RA.hsb;
  end
  ref(end+1,1:3)={'z_pad_down',z_pad_down,'Bottom face of pad'};  
  ref(end+1,1:3)={'z_pad_mid',z_pad_down+RA.hs/2,'Middle pad'};
  ref(end+1,1:3)={'z_pad_up',z_pad_down+RA.hs,'Upper face of pad'};
  ref(end+1,1:3)={'z_rail_up',z_pad_down+RA.hs+RA.hr3,'Rail up'};
  ref(end+1,1:3)={'z_rail_web',z_pad_down+RA.hs+RA.hr1+(RA.hr3-RA.hr2-RA.hr1)/2,'Rail web'};
  ref(end+1,1:3)={'z_rail_foot',z_pad_down+RA.hs+RA.hr1,'Rail foot'};
  

  
  %% Track geometry parameters
  else; % This is a track  
   RT=RO; mt=PARAM.mt; ms=stack_get(mt,'SE'); 
   R1=stack_get(ms{1,3},'info','MeshParam','get');
   R2=stack_get(ms{end,3},'info','MeshParam','get');
  % X references:
  if isfield(RT,'xoffset');
  elseif isfield(RO,'bas');xoffset=min(mt.bas(:,4)); mt=rmfield(mt,'bas');
  else;xoffset=0;
  end
  ref(end+1,1:3)={'x_offset',xoffset,'First acceptable point'};
  ref(end+1,1:3)={'x_start',xoffset+(RT.nb_edge+4)*R1.sw,'First acceptable point'};
  ref(end+1,1:3)={'x_end',xoffset+(sum(RT.nb_slices)-RT.nb_edge-4)*R2.sw,'Last point'};
  if ref{end,2}-ref{end-1,2}<0; error('Track too short');end
  % Y references
  ms=PARAM.ms; RO=stack_get(ms,'info','MeshParam','get');
  if ~isempty(RO)&&isfield(RO,'Geo')
     ref(end+1,1:3)={'y_track',RO.Geo.y_slice,'Middle of the track'};
  else
     ref(end+1,1:3)={'y_track',0,'Middle of the track'};
  end
  
  end
  
  %% Affect parameters
  if exist('RT','var'); % is a track
      if comstr(Cam,'info') %
       st=ref(:,1:2)';RT.Geo=struct(st{:});
       out=RT;
      elseif comstr(Cam,'mapstack') %
       map=containers.Map(ref(:,1),ref(:,2),'UniformValues',false);
       mt=stack_set(mt,'info','FormalGeo',map);out=mt;
      elseif comstr(Cam,'cell');out=ref; %
      else
       R1=struct;
       for j1=1:size(ref,1); 
          R1.(ref{j1,1})=struct('type','double', ...
              'value',ref{j1,2},'ToolTip',ref{j1,3});
       end
          RT.Geo=R1;
          out=RT;
      end
  else % is a slice
      if comstr(Cam,'info') %CbGeoInfo : return struct
       st=ref(:,1:2)';RO.Geo=struct(st{:});
       out=RO;
      elseif comstr(Cam,'mapstack') %
       map=containers.Map(ref(:,1),ref(:,2),'UniformValues',false);
       mt=stack_set(mt,'info','FormalGeo',map);out=mt;
      elseif comstr(Cam,'cell');out=ref; %
      else
       R1=struct;
       for j1=1:size(ref,1); 
          R1.(ref{j1,1})=struct('type','double', ...
              'value',ref{j1,2},'ToolTip',ref{j1,3});
       end
          RO.Geo=R1;
          out=RO;
      end
      
   end
 end

 %% #CBPost : Callbacks for Post tab - - - - - - - - - - - - - - - - - - - - -  
 elseif comstr(Cam,'post');[CAM,Cam]=comstr(CAM,5);
  error('THIS IS NOW OBSOLETE : uses compute post instead ! xxxjp: remove when moved')    
  
 %% #CBShow : display results - - - - - - - - - - - - - - - - - - - - - - - -2
 elseif comstr(Cam,'show');[CAM,Cam]=comstr(CAM,5);
  [CAM,Cam,RunOpt.report]=comstr('-report',[-25 3],CAM,Cam); % do something to capture..
  if isfield(PARAM,'Project'); rp=fe_def('cleanentry',PARAM.Project); 
  else; rp=[]; 
  end
  
 %% #CBResult : Callbacks for contextual menu on result - - - - - - - - - - -
 elseif comstr(Cam,'result');[CAM,Cam]=comstr(CAM,7);
  
   % get calling result (contextual menu should be associated to one result, 
      % so is the tab) :
   r1j=PARAM.Result;
   r1=cinguj('ObjToStruct',r1j);
   name=fe_def('cleanentry',r1); name=name.name;
  
  %% #CBResultSave : Save in .mat file  - - - - - - - - - - - - - - - - - - -3 
  if comstr(Cam,'save');[CAM,Cam]=comstr(CAM,5);
 
   % ask for mat file to save
   [fname, wd, filterindex] = uiputfile( ...
      {'*.mat', 'MAT-files (*.mat)'},sprintf('Save %s',name));
   if ~ischar(fname); return;end % Cancel
   % save variable
   fname=fullfile(wd,fname);
   eval(sprintf('%s=PARAM.%s.GetData;',name,name)) % get variable in current workspace
   save(fname,name,'-v7.3') % save it
   
  %% #CBResultCapture : Capture (call to show-report) image - - - - - - - - -3 
  elseif comstr(Cam,'capture');[CAM,Cam]=comstr(CAM,8);
 
   % ask for fileroot to save
   st1='*.png'; st2='png-files (*.png)';
   if comstr(name,'def'); st1='*.avi'; st2='avi-files (*.avi)'; end
   [fname, wd, filterindex] = uiputfile( ...
      {st1, st2},sprintf('Capture %s',name));
   if ~ischar(fname); return;end % Cancel
   % call to show-report:
   [wd,root,ext]=fileparts(fullfile(wd,fname));
   railui(sprintf('CBShow%s-report-wdplot"%s"-root"%s"',name,wd,root))
   
  %% #CBResultExport : Export current result in Matlab base workspace - - - -3 
  elseif comstr(Cam,'export');[CAM,Cam]=comstr(CAM,7); 
    
   eval(sprintf('%s=PARAM.%s.GetData;',name,name))% get variable in current workspace
   eval(iigui({name},'SetInBaseC')) % set in base workspace
   
  else; error('CBResult%s unknown',CAM);
  end
  
  
 %% #CBMenu : Callbacks for menu entries - - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'menu');[CAM,Cam]=comstr(CAM,5);
  %% #CBMenuAbout : About - - - - - - - - - - - - - - - - - - - - - - - - - -3 
  if comstr(Cam,'about');[CAM,Cam]=comstr(CAM,6); 
   % About box
   %h=dialog('WindowStyle', 'normal', 'Name', 'About OSCAR');
   
   st=sdtkey('version');
   
   xt=160; yt=220;
   h=findall(0,'type','figure','tag','AboutOSCAR');
   if isempty(h); h=figure('visible','off','integerhandle','off'); 
   else figure(h);
   end
 
   if exist('AboutOSCAR.jpg','file')
    a=imread('AboutOSCAR.jpg');
    h=image(a,'buttondownfcn','web(''http://www.sncf.fr'')');
    axis image;
   end
   axis off;set(gca,'position',[0 0 1 1]);
   i1=get(gcf,'position');

   h=text('units','pixels','position',[xt yt 0],'string',sprintf('SDT ver : %s',st));
   i2=1;
   while get(h,'extent')*[0;0;1;0]<300 && i2<10
    set(h,'fontsize',get(h,'fontsize')+1); i2=i2+1;
   end

   try; st=sdtdef('key'); catch; st='Undefined key'; end
   h=text('units','pixels','position',[xt yt-30 0],'string',sprintf('SDT Key : %s',st));

   st5=railui('cvs');
   st5=strrep(st5,'$',''); st5=strrep(st5,'::',':'); st5=strrep(st5,'#','');
   h=text('units','pixels','position',[xt yt-60 0],'string',sprintf('GUI : %s',st5));
  
   i2=1;
   while get(h,'extent')*[0;0;1;0]<250 && i2<10
    set(h,'fontsize',get(h,'fontsize')+1); i2=i2+1;
   end

   set(gcf,'position',[i1(1:2) 908/2 474/2],'menubar','none','tag','AboutOSCAR', ...
           'visible','on','NumberTitle','off','Name','About OSCAR');
   axis off

%   else
%    msgbox('Copyright SNCF 2012','About OSCAR')
%   end
%    

  %% #CBMenuResultsReportAll :  - - - - - - - - - - - - - - - - - - - - - - -3 
  elseif comstr(Cam,'resultsreportall');[CAM,Cam]=comstr(CAM,17); 

   if comstr(Cam,'as')
    rp=fe_def('cleanentry',PARAM.Project);
    [fname,wd,filterindex] = uiputfile( ...
     {'*.html', 'HTML-files (*.html)'}, 'Save report as', ...
      sprintf('%s.html',rp.root));
    if ~ischar(fname); return;end % Cancel
    [un1,root,ext]=fileparts(fname);
    wdplot=fullfile(wd,'plots');
    railui(sprintf('ComputeReport-FromInt-AllComputed-wdplot"%s"-root"%s"',wdplot,root))
   else
    railui('ComputeReport-FromInt-AllComputed')
   end
   
  else; error('CBMenu%s unknown',CAM);
  end
 %% #CBProjectSave  - - - - - - - - - - - - - - - - - - - - - - -3 
 elseif comstr(Cam,'projectsave');[CAM,Cam]=comstr(CAM,12); 

  R1=struct('Project',railui('PARAM.Project -safe'));
  if isfield(PARAM,'Range');R1.Range=PARAM.Range;end
  if nargin>3&&isfield(varargin{4},'Save');RO=varargin{4};carg=5;
  elseif carg<=nargin; RO=varargin{carg};carg=carg+1;
  else;RO=struct;
  end
  if ischar(RO);RO=struct('Save',RO);
  elseif ~isfield(RO,'Save');RO.Save='do';
  end
  if ~isfield(RO,'KeepField');
      RO.KeepField={'mt','ms','dstat','def','Range','Stack','perMode','perHist','perResp'};
  end
  sdtroot('cbProjectSave',UI,R1,RO);
  
   %% #CBProjectClose  - - - - - - - - - - - - - - - - - - - - - - -3 
 elseif comstr(Cam,'projectclose');[CAM,Cam]=comstr(CAM,13); 
  if sdtdef('isinteractive')
  st=questdlg('Close DYNAVOIE?','Closing request', ...
              'Close without saving','Save then close','Cancel close', ...
              'Close without saving');
  switch st
  case 'Close without saving'; 
   i1=1;
  case 'Save then close';  
   railui('CBProjectSave-justsave') % SAVE
   i1=1;
  otherwise; 
   return; % out=0 in this case
  end
  else; i1=1;
  end
  if i1==1
     feval(cinguj('@tabChange'),'closeall',UI.gf) 
     set(UI.gf,'userdata',[]); % Clean UI
     delete(UI.gf);
     return 
  end
  
     
 %% #CBGetModel struct('MeshCfg',val,'step',1)  - - - - - - - - - - - -3 
 elseif comstr(Cam,'getmodel');
  RO=varargin{carg};carg=carg+1;
  data=PARAM.Range; 
  data=data.param.MeshCfg;data=data.data{strcmpi(data.choices,RO.MeshCfg)};
  if ~isempty(strfind(lower(RO.step),'coarse'))
   out=dyn_mesh('build',data);
   r1=RO.MeshCfg;if isjava(r1);r1=char(r1);end
   PARAM=stack_set(PARAM,'SE',[r1 '_mesh'],out);
  else;error('Not implemented');
  end
  
 %% #CB : reroute to sdtroot  - - - - - - - - - - - - - - - - - - - - - - -3 
 else
  if carg<=nargin&&isfield(varargin{carg},'MainFcn');carg=carg+1;end   
  eval(iigui({'sdtroot([''cb'' CAM],UI,varargin{carg:end})',nargout}, ...
      'OutReDir'))
  %sdtroot(['cb' CAM],UI,varargin{carg:end});
 end

%% #PARAM : RunOpt to/from PARAM ---------------------------------------------
elseif comstr(Cam,'param');[CAM,Cam]=comstr(CAM,6);
 
 
 %% #PARAM2RO : Build RunOpt from PARAM - - - - - - - - - - - - - - - - - - -
 %% - - - - - - - - - - - - -
 if comstr(Cam,'2ro');[CAM,Cam]=comstr(CAM,4);
 
  UI=initFig(UI);PARAM=v_handle('uo',UI.gf);
  
  % First : pass through all tabs (since PARAM is filled only when tab is
  % displayed):
  % Init all tabs: xxx
  st={'Project','Slice','Track','Veh','Simu'};
  for j1=1:length(st);
   railui(['Set' st{j1}],[]);
  end
  
  RO=fe_def('cleanentry',PARAM.GetData,struct('sub',{st}));

%   % Existing results
RO=feutil('rmfield',RO,'mt');
%   list=resultslist; % get normalized results
%   for j1=1:length(list)
%    RO.(list{j1})=railui(sprintf('vh%s',list{j1}));
%    r5=RO.(list{j1});
%    RO.(list{j1})=r5.GetData; % Generate holder for mo1
%   end

  if nargout==0; 
   i1=cellfun('isclass',struct2cell(RO),'v_handle');
   st=fieldnames(RO);st=st(i1);list={'xxx'};
   disp(comstr(feutil('rmfield',RO,st{:},list{:}),-30));
  else;
   out=RO;
  end
 
 %% #PARAMFromRO : Build PARAM from RunOpt - - - - - - - - - - - - - - - - - -
 elseif comstr(Cam,'fromro');[CAM,Cam]=comstr(CAM,7);
  % just add results of computed steps to PARAM 
  RO=varargin{carg}; carg=carg+1;
  if isfield(RO,'Elt')&&isfield(RO,'Node') % model is given
   RO=stack_get(RO,'info','MeshParam','getdata');
  end
    
  UI=initFig(UI);set(UI.gf,'userdata',[]); % Empty former PARAMs
  railui('SetProject',[])
  PARAM=v_handle('uo',UI.gf);
  toCinCell=cinguj('@toCinCell');
  rset={}; % parameters that has been used to set ...
  
  % Init all tabs: xxx
   st={'Project','Slice','Track','Vehicle','Simulation'};
   for j1=1:length(st);
    railui(['Set' st{j1}],[]);
   end
   
  % Get all parameters defined in UI, and assign corresponding par from RO
  [obj,objName]=sdcget(UI.gf,'',''); % objname contain tab name
  for j1=1:length(obj)
   name=sdcedit(obj{j1},'_get','name');
   i1=strfind(name,'.'); if ~isempty(i1); name=name(i1(end)+1:end); end
   if ~isfield(RO,name)
    disp(sprintf('Missing parameter %s in RO, keep UI default',name))
   else
    rset{end+1,1}=name;
    railui(sprintf('Set%s',objName{j1}),struct(name,RO.(name)))
   end
  end
  
  r3=setdiff(fieldnames(RO),rset); % parameterss that have not been set
  disp('Following RO param havent been used to set UI : (some are unused param due to the fact that RO contains ALL possible param.')
  disp(r3)
      
elseif comstr(Cam,'idlayer'); 
%% #PARAMIdLayer : information about layers
 if carg<=nargin; mo1=varargin{carg};carg=carg+1; % should be robust for slice
 elseif ~isempty(UI.gf);   PARAM=v_handle('uo',UI.gf);
 else; error('Model needed');
 end
 pl=feutil('getpl',mo1);
 out={'MatId','ymean','zmin','zmax','zmean','rho','mass','E','Nu'};
 out(end+size(pl,1),1)={''};
 r1=mo1;r1.il(r1.il(:,4)==10002,4)=-3;r1=feutilb('geomrb Mass -bymat',r1);
 
 for j1=1:size(pl,1)
  n1=feutil(sprintf('getnode matid%i',pl(j1,1)),mo1);
  n1=[min(n1(:,6)) max(n1(:,6)) min(n1(:,7)) max(n1(:,7))];
  mat=feutil(sprintf('getpl%i-struct',pl(j1,1)),mo1);
  % missing robust to anisotropic material
  if isfield(mat,'E')
   out(j1+1,1:9)={pl(j1,1) mean(abs(n1(1:2))) n1(3) n1(4) mean(n1(3:4))  ...
      mat.Rho r1.mass.Y(r1.mass.X{1}==pl(j1,1),1), ...
      mat.E mat.Nu};
  else
   out(j1+1,1:9)={pl(j1,1) mean(abs(n1(1:2))) n1(3) n1(4) mean(n1(3:4))  ...
      mat.Rho r1.mass.Y(r1.mass.X{1}==pl(j1,1),1), ...
      'ortho' 'Ortho'};
  end
 end
 if nargout==0
   if isempty(UI.gf); railui('init',mo1);PARAM=railui('vh');end
   comstr(out(2:end,:),-17,'tab', ...
       struct('ColumnName',{out(1,:)}, ...
         'name','Info','FigTag','railui'))  
 end
%% #railui('ParamVH')
elseif comstr(Cam,'vh');
    UI=initFig(UI); out=v_handle('uo',UI.gf);
    
elseif comstr(Cam,'.'); 
%% #PARAM.get parameters
 UI=initFig(UI);out=sdtroot(['PARAM' CAM],UI);
%  UI.DefBut=genDefBut;
%  [CAM,Cam,st]=comstr('.',[-25 4 1],CAM,Cam);
%  out=eval(sprintf('DefBut.%s;',st));
%  if nargout==0
%    % could be 
%    % sdt_dialogs('uatable -end0 -tipascol 1','info','Ok',a);fe_def('cleanentry',ans)
%    out=sdt_locale(['defbutTreeName' CAM],out);
%    r1=out.table(:,1:3);i1=cellfun('isclass',r1,'struct');
%    r1(i1)=cellfun(@(x)x.value,r1(i1),'uni',0);
%    if ~isempty(strfind(Cam,'tex'));
%        comstr(r1,-17,'tex');clear out
%    else;comstr(r1,-17,'text');clear out;
%    end
%  end
 
%% #PARAM_end - - - -
elseif isempty(Cam); % Renew & return v_handle to param
 DefBut=genDefBut;
 if nargout==0; 
   sdt_locale('DefButTreeName',DefBut)
   gf=findall(0,'tag','DefButTree'); uo=[];
   uf=clean_get_uf(gf);
   ua=get(uf.Explo,'userdata');ua.OutFmt='tex';set(uf.Explo,'userdata',ua);
 else;out=DefBut;
 end
elseif comstr(Cam,'.') % GetValue oscarui('param.Project.WdCat')
 UI=initFig(UI);out=sdtroot(['PARAM' CAM],UI);
 %RO=oscarui('PARAM2RO'); % Need optimze xxxeb
 %out=eval(sprintf('RO%s',CAM));
else; 
  UI=initFig(UI);
  if nargout==0; sdtroot(sprintf('PARAM%s',CAM),UI,varargin{carg:end});
  else; out=sdtroot(sprintf('PARAM%s',CAM),UI,varargin{carg:end});
  end
end

  
%% #End
%% #icon : icon loading, calls sdticons -------------------------------------2
elseif comstr(Cam,'icon');[CAM,Cam]=comstr(CAM,5);
   % 'at sdtools'   load from m/icons : ask eb to clean sdticons
   % when deployed use 'railui.jar'
   
%% #wd ----------------------------------------------------------------------2
elseif comstr(Cam,'wd');[CAM,Cam]=comstr(CAM,3);

if comstr(Cam,'_')
 %wd_  : subcommands    
 % In particular distinguish : Models,Results,plots,Meshsource
else
 UI=initFig(UI);PARAM=v_handle('uo',UI.gf);
end

%% #Verbose_mode ------------------------------------------------------------2
elseif comstr(Cam,'verbose_mode');
 if carg<=nargin;  Verbose_Mode=varargin{carg}; carg=carg+1;
 else Verbose_Mode=0;
 end

elseif comstr(Cam,'getpar')
 out=getPar(varargin{carg:end});
elseif comstr(Cam,'cvs')
 out='$Revision: 1201 $  $Date:: 2020-07-17 08:09:59#$';
elseif comstr(Cam,'@');out=eval(CAM);
else
 try
  %if ~isempty(obj)
  % eval(iigui({'sdtroot(varargin{1:3},UI,varargin{4:end})',nargout},'OutReDir'))
  %else
   eval(iigui({'sdtroot(varargin{1},UI,varargin{2:end})',nargout},'OutReDir'))
  %end
 catch err
  error('%s unknown',CAM);
 end
end
 
 
end

%% #SubFunc


%% #initFig
function UI=initFig(UI,varargin);

if isempty(UI.gf)||~ishandle(UI.gf);
 if exist('com.sdtools.cinguj.EditT','class')
     cinguj('InitSwing');
 end
 tag='railui';
 gf=findall(0,'type','figure','tag',tag);
 if nargin>1&&strncmpi(varargin{1},'init',4)&&~isempty(gf);
   GuiGF=gf;%PARAM=v_handle('uo',gf);
 elseif isempty(gf);
   gf=figure('IntegerHandle','off','NumberTitle','off','name','DYNAVOIE', ...
        'tag',tag,'CloseRequestFcn','railui(''CbProjectClose'')');
   r6=get(0,'screensize'); % screen dimensions
   r5=get(gf,'outerposition');
   dx=640; % OSCAR UI width
   dy=480; % OSCAR UI heigh
   r5(3:4)=[dx dy]; % OSCAR windows size : 640x480
   r5(1)=(r6(3)-dx)/2; % x position : center
   r5(2)=(r6(4)-dy-50); % y position : top
   set(gf,'outerposition',r5);
   GuiGF=gf;
   PARAM=v_handle('uo',gf);
   if ~isempty(varargin)&&isfield(varargin{1},'Elt');
       PARAM.mt=varargin{1};
   else; PARAM.mt=[];
   end
   assignin('caller','GuiGF',GuiGF);
   railui('initMain');
   [uf,gf]=clean_get_uf(GuiGF);uf.MainFcn=@railui;
   set(gf,'userdata',[],'userdata',uf);

 else;GuiGF=gf;%PARAM=v_handle('uo',gf);
 end
 UI.gf=GuiGF;
 if usejava('awt')==0; UI.treeF=[]; end
 
 jFrame=get(UI.gf,'javaframe');
 wd=fileparts(which('railui'));
 try
  jicon=javax.swing.ImageIcon(fullfile(wd,'LogoDynavoie_32.png'));
  jFrame.setFigureIcon(jicon);
 end
 
%else;PARAM=v_handle('uo',GuiGF);
end
end


%% #check_empty_project : check if current project is empty
function r1=check_empty_project
 
 PARAM=railui('vh'); % get param
 r1=1; % empty
 dbstack;
 'xxxjp emptyparam'
 if isempty(PARAM); return; end
 
 % check if panto/section filled 
 if isfield(PARAM,'Panto') && ~isempty(PARAM.Panto); r1=0; end
 if isfield(PARAM,'Section') && ~isempty(PARAM.Section); r1=0; end
 
 % xxxjp : something to done to compare simu to default ?
 
 % check if results exist
 if 1==2

  list=resultslist;
  r2=intersect(list,fieldnames(PARAM));
  if ~isempty(r2); r1=0; end % there is a result
 else
     'xxxjp check'
 end
end
 


%% #getPar -----------------------------------------------------    
% value=feval(railui('@getPar'),model,name);
function [out]=getPar(varargin); 

if nargin==0 
    [DefBut,ParString]=feval(railui('@genDefBut'));
else; model=varargin{1};carg=2;
end
if carg<=nargin;name=varargin{carg};carg=carg+1;
else; name='';
end
if isfield(model,'Elt'); RO=feval(dyn_mesh('@getParamInStack'),model);
else;RO=model;
end
if isempty(name); out=RO;
elseif isfield(RO,name); out=RO.(name);
elseif iscell(name)
  RO=railui('CbGeoInfo',model);
elseif strncmpi(name,'TrackG',6)
  RO=railui('CbGeoInfo',model); out=RO.TrackG;
  [CAM,Cam]=comstr(name,8);
  if isempty(CAM)
  elseif isfield(out,CAM);out=out.(CAM);
  else; out=str2double(CAM);
  end
elseif strncmpi(name,'gravity',6); % Default gravity
  out=-9.80665;
elseif isfield(RO,'Slice')
  out=RO.Slice{1,3};
  [CAM,Cam]=comstr(name,8);if isempty(CAM);[CAM,Cam]=comstr(name,1);end
  if isempty(CAM)
  elseif isfield(out,CAM);out=out.(CAM);
  elseif isfield(out,name);out=out.(name);
  elseif isfield(out,'Geo')&&isfield(out.Geo,name);out=out.Geo.(name);
  else; out=str2double(CAM);
  end
  
else
  DefBut=feval(railui('@genDefBut'));
  error('Sub indexing not done yet');
end
end
%% #setPar -----------------------------------------------------    
% model=feval(railui('@setPar'),model,name,value);
function [out]=setPar(varargin); %#ok<DEFNU>

model=varargin{1};carg=2;
name=varargin{carg};carg=carg+1;
error('not implemented');
if isfield(model,'Elt')
   RO=stack_get(m_global,'info','MeshParam','getdata');
end
if isfield(RO,name); out=RO.(name);
else; error('Sub   indexing not done yet');
end
end

%% #delPar -----------------------------------------------------    
function out=delPar(varargin); %#ok<DEFNU>
  r1=varargin{1}; st=varargin{2};
  r2j=com.sdtools.cinguj.EditT;
  st2=fieldnames(r1);
  for j2=1:length(st2)
    if ~strcmpi(st2{j2},st)
        r2j.setField(st2{j2},r1.(st2{j2}))
    end
  end
  out=r2j;
end

%% #genDefBut -----------------------------------------------------    
% [DefBut,ParString]=feval(railui('@genDefBut'));
function r1=genDefBut
 %setpref('SDT','locale',{'en-us','fr-fr'})
 % [DefBut,ParString]=feval(railui('@genDefBut'));
 try
   r1=sdt_locale('mergelang',dyn_utils('wd','dynavoie_ui_'),'',struct('FixCb',1));
 catch err; 
   disp(err)
   sdtw('_err','dynavoie_ui csv is missing or corrupted')   
   r1=[]; return; % Bypass for old SDT with new oscar / I&R case..
 end
 
 %% cingui('paramprint',struct(st2{j2},r2.(st2{j2})))
 
 r3=sdtroot('paramUI');r2=r3.DefBut.PTree;
 r1.PTree=sdth.sfield('AddMissing',r1.PTree,r2);
 r1.Tab=r3.DefBut.Tab; 
 
end

%% #setUtil
function out=setUtil(varargin)
% These parameters are defined for each slice. Other parameters are
% assumed to be read in defbut for corresponding slice, or from coarse model
% meshing script

[CAM,Cam]=comstr(varargin{1},1);
if comstr(Cam,'emptyslice')
 error('xxx to be updated')
  % Possibly remove some fields of panto setUtil('emptypanto','name','oldname')
  st=varargin{2};oldst=varargin{3};
  if strcmp(oldst,st); return;end
  eval(iigui({'GuiGF','r1','r1j','RB'},'GetInCaller'))
  r1=feutil('rmfield',r1,setdiff(fieldnames(r1),PtoGenField)); %#ok<NODEF>
  r1j=cinguj('ObjEditJ',r1,GuiGF);
  RB.refresh='oscarui(''InitPanto'')';
  eval(iigui({'r1','r1j','RB'},'SetInCaller'))
elseif comstr(Cam,'init'); 
%% Extract fields from existing RO
  if isempty(strfind(Cam,'field'))
    r1=evalin('caller','UI.DefBut');st=fieldnames(r1.Track);
  else; st=setUtil(CAM(5:end));
  end
  RO=varargin{2}; mt=varargin{3};
  if isempty(mt); out=RO; return;end
  r2=getPar(mt);
  for j1=1:length(st)
   if isfield(r2,st{j1});RO.(st{j1})=r2.(st{j1});end 
  end
  out=RO;
else;error('%s',CAM);
end
end


%% #addPML
function RO=addPML(RO)
  if RO.half
     RO.Sub.PML=struct('Lp',[0 0 1 0 1 0]);
  else
     RO.Sub.PML=struct('Lp',[0 1 1 0 1 0]);
  end
  RO=rmfield(RO,'PML');
  RO.PMLblock=1; % block pmls
end


