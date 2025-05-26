classdef cntc

 % cntc : demonstration files for CONTACT
 %
 % <a href="matlab:methods cntc">methods cntc</a> for methods
 % <a href="matlab:sdtweb cntc">sdtweb cntc</a> for HTML documentation
 % <a href="matlab:sdtweb _taglist cntc">sdtweb _taglist cntc</a> for taglist
 %
 %
 % For revision information type <a href="matlab:cntc.cvs">cntc.cvs</a>

 %  SDT wrapper E. Balmes, G. Guillet
 %  Copyright (c) 1990-2025 by SDTools, All Rights Reserved.
 % CONTACT wrapper functions
 %  Copyright 2008-2023 by Vtech CMCC.
 %  Licensed under Apache License v2.0.  See the file "LICENSE_CONTACT.txt" for more information.

 %#ok<*NOSEM,*PROP,*ASGLU,*STREMP,*OR2,*AND2,*INUSD,*NBRAK1,*NBRAK2,*AGROW,*NASGU>
 properties (GetAccess = public , SetAccess = public)
  % type
  libname
 end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %% #SDT/CONTACT initialization

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %% #idx.refresh{cntc.help(pdf),cntc.help(md)} -2
 %{
```DocString  {module=rail,src=cntc.md} -2
cntc : interface between SDT and CONTACT
 %}

 methods
  function ob = cntc(varargin)
   help cntc
  end

 end % methods

 % #Static / generic methods  -------------------------------------------- -2
 methods (Static)

  %------------------------------------------------------------------------------------------------------------
  function []=testlibrary()
   [CNTC, ifcver, ierror] = cntc_initlibrary;
  end
  %------------------------------------------------------------------------------------------------------------

  %------------------------------------------------------------------------------------------------------------
  function []=initializeflags()
   % #initializeflags-2
   if nargin==0
    LI=cntc.call;
    if ~isfield(LI,'CNTC');cntc.init;end
    CNTC=LI.CNTC;
    li=LI.global;
    cntc.setglobalflags(CNTC.if_idebug, li.ire);
    [ifcver, ierror] = cntc.initialize(li.ire,li.imodul);
   else
    disp('ERROR cntc.initializeflags : no parameter is needed');
   end
  end
  %------------------------------------------------------------------------------------------------------------

  %------------------------------------------------------------------------------------------------------------
  function []=clearLI()
   % #clearLI -2
   LI=cntc.call;
   if isfield(LI,'libname')
    if ~isempty(LI.libname)||libisloaded(LI.libname)
     cntc.closelibrary;
    end
    evalin('base','clear("LI")')
   else
    disp('Error : LI is already empty')
   end
  end
  %------------------------------------------------------------------------------------------------------------

  %------------------------------------------------------------------------------------------------------------
  function []=setglobalflags(params, values)
   % [ ] = cntc.setglobalflags(params, values)
   %
   % used for configuring flags that are the same for all contact problems
   %
   %  params - codes of the parameters to be communicated to CONTACT
   %  values - values of the parameters to be communicated to CONTACT
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 0: m=*, glob   - no icp needed
   % #setglobalflags -2



   if (nargin<2 | isempty(params) | isempty(values))
    disp('ERROR in cntc.setglobalflags: params and values are mandatory.');
    return
   end
   if (length(params) ~= length(values))
    disp('ERROR in cntc.setglobalflags: params and values must provide same number of items.');
    return
   end

   cntc.call('cntc.setglobalflags', length(params), params, values);

  end % cntc.setglobalflags

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [ ifcver, ierror]=initialize(ire, imodul, c_outdir, idebug)
   % [ ifcver, ierror ] = cntc.initialize(ire, imodul, [outdir], [idebug])
   %
   % upon first call: initialize the addon internal data and initialize output channels,
   %                  print version information;
   % for each ire:   initialize and return the addon version number.
   %
   %  in:  integer    ire          - result element ID
   %       integer    imodul       - module number 1=w/r contact, 3=basic contact
   %       character  outdir(*)    - output folder
   %       integer    idebug       - show (1) or hide (0) error messages
   %  out: integer    ifcver       - version of the CONTACT add-on
   %       integer    ierror       - error flag
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 0: m=*, glob   - no icp needed
   % #initialize -2


   cntc_err_allow = -12;

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(imodul))
    disp(sprintf('cntc.initialize: please select module 1 or 3'));
    return;
   end
   if (nargin<3 | isempty(c_outdir))
    c_outdir = ' ';
   end
   if (nargin<4 | isempty(idebug))
    idebug = 1;
   end

   p_ierr = libpointer('int32Ptr',-1);
   p_ver  = libpointer('int32Ptr',-1);

   p_ierr.value = 0;
   cntc.call('cntc.initialize', ire, imodul, p_ver, p_ierr, c_outdir, length(c_outdir));
   % disp(sprintf('test_caddon: obtained ver=%d, ierr=%d', p_ver.value, p_ierr.value));

   ifcver = double(p_ver.value);
   ierror = double(p_ierr.value);

   if (idebug>=1 & ierror==cntc_err_allow)
    disp(sprintf('cntc.initialize: no license found or license invalid, check output-file (%d).',ierror));
   elseif (idebug>=1 & ierror<0)
    disp(sprintf('cntc.initialize: an error occurred in the CONTACT library (%d).',ierror));
   end

  end % cntc.initialize

  %------------------------------------------------------------------------------------------------------------

  %------------------------------------------------------------------------------------------------------------
  function [ CNTC, ifcver, ierror ] = initlibrary(c_wrkdir, c_outdir, c_expnam, idebug);
   % [ CNTC, ifcver, ierror ] = cntc.initlibrary(wrkdir, outdir, expnam, idebug);
   %
   % load the library into Matlab, initialize its internal data and output channels
   %
   %  in:  character  wrkdir(*)    - [optional] effective working folder
   %       character  outdir(*)    - [optional] output folder
   %       character  expnam(*)    - [optional] experiment name, default 'contact_addon'
   %       integer    idebug       - [optional] show (1) or hide (0) error messages
   %  out: integer    CNTC         - struct with 'magic numbers' for configuring CONTACT
   %       integer    ifcver       - version of the CONTACT add-on
   %       integer    ierror       - error flag


   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 0: m=*, glob   - no icp needed
   % #initlibrary -2

   % retrieve current directory, form library name
   [pathstr, name, ext] = fileparts(which('cntc.initlibrary'));

   if (isunix())
    libname='contact_addon_linux64';
    lib_ext='.so';
   else
    libname=['contact_addon_', computer('arch')];
    lib_ext='.dll';
   end

   % unload the library if it was loaded before
   %xxxgae cant unload the library when it crash and the library is still
   %on
   if libisloaded(libname)
    LI=cntc.call;LI.libname=libname;
    cntc.closelibrary;
   end
   if nargin>0
    if isequal(c_wrkdir,'close');
     for iwhe= 1:2
      cntc.finalize(iwhe);
     end
     cntc.closelibrary;
     return
    end %
   end
   % load the library into Matlab
   if (~libisloaded(libname))
    pathstr = [deblank(pathstr) filesep '..' filesep 'bin'];
    if ispref('SDT','CONTACT_Path')
     pathstr=getpref('SDT','CONTACT_Path');
    end
    fullname = fullfile(pathstr, [libname lib_ext]);

    if ~exist(fullname,'file')
     disp(['ERROR: cant find library: ', fullname]);
     return
    else
     loadlibrary(fullname, 'contact_addon.h');
    end
   end

   % initialize the internal data of the library, open its output streams
   if (nargin<1 | isempty(c_wrkdir))
    c_wrkdir = ' ';
   end
   if (nargin<2 | isempty(c_outdir))
    c_outdir = ' ';
   end
   if (nargin<3 | isempty(c_expnam))
    c_expnam = ' ';
   end
   if (nargin<4 | isempty(idebug))
    idebug = 1;
   end
   len_wrkdir = length(c_wrkdir);
   len_outdir = length(c_outdir);
   len_expnam = length(c_expnam);
   ioutput = 0;

   p_ierr = libpointer('int32Ptr',-1);
   p_ver  = libpointer('int32Ptr',-1);

   p_ierr.value = 0;
   cntc.call(libname,'cntc.initializefirst_new', p_ver, p_ierr, ioutput, c_wrkdir, c_outdir, c_expnam, ...
    len_wrkdir, len_outdir, len_expnam);
   % disp(sprintf('test_caddon: obtained ver=%d, ierr=%d', p_ver.value, p_ierr.value));

   ifcver = double(p_ver.value);
   ierror = double(p_ierr.value);

   % return a struct with 'magic numbers' for setting flags later on

   CNTC = cntc.getmagicnumbers();

   if (idebug>=1 & ierror==CNTC.err_allow)
    % error('cntc.initlibrary: no license found or license invalid, check output-file (%d).',ierror);
   elseif (idebug>=1 & ierror<0)
    % error('cntc.initlibrary: an error occurred in the CONTACT library (%d).',ierror);
   end

  end % cntc.initlibrary

  %------------------------------------------------------------------------------------------------------------

  function [ ierror]=managelicense(lic_command, ints, string1, string2, string3)
   % [ ierror ] = cntc.managelicense(lic_command, ints, string1, string2, string3)
   %
   % perform license management actions:
   %  - 'activate' <license_id> <password>
   %  - 'refresh'
   %  - 'print'
   %  - 'deactivate'
   %
   %  in:  character  lic_command(*) - desired action
   %       integer    ints(*)        - 'activate': ints(1) == <license_id>
   %       character  string1(*)     - 'activate': string1 == <password>
   %       character  string2(*)     - reserved for future use
   %       character  string3(*)     - reserved for future use
   %  out: integer    ierror         - error flag

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 0: m=*, glob   - no icp needed
   % #managelicense -2


   if (isempty(libname) | ~libisloaded(libname))
    % load library when needed, using default outpath and experiment name
    [ CNTC, ifcver, ierror ] = cntc.initlibrary();
   end

   if (nargin<1 | isempty(lic_command))
    lic_command = 'print';
   end
   disp('cntc.managelicense: output is printed to contact_addon.out.');

   ierror = 0;
   p_ierr = libpointer('int32Ptr',-1);

   % command == 'activate':

   if (strcmp(lic_command, 'activate'))
    if (nargin<2 | isempty(ints))
     disp('cntc.managelicense: activation needs the license id.');
     ierror = -1;
    else
     license_id = ints(1);
    end
    if (nargin<3 | isempty(string1))
     disp('cntc.managelicense: activation needs the license password.');
     ierror = -2;
    else
     password = string1;
    end

    if (ierror==0)
     p_ierr.value = 0;
     cntc.call('cntc.managelicense', lic_command, 1, license_id, password, ' ', ' ', ...
      length(lic_command), length(password), 1, 1, p_ierr);
     ierror = double(p_ierr.value);
    end

    % command == 'refresh', 'print' or 'deactivate':

   elseif (strcmp(lic_command, 'refresh') || strcmp(lic_command, 'print') || ...
     strcmp(lic_command, 'deactivate'))

    cntc.call('cntc.managelicense', lic_command, 0, 0, ' ', ' ', ' ', ...
     length(lic_command), 1, 1, 1, p_ierr);
    ierror = double(p_ierr.value);

   end

  end % cntc.managelicense

  %------------------------------------------------------------------------------------------------------------

  function [ ierror]=readinpfile(ire, inp_type, c_fname)
   % [ ierror ] = cntc.readinpfile(ire, inp_type, fname)
   %
   % read settings from inp-file
   %
   %  in:  integer    ire          - result element ID
   %       integer    inp_type     - type of inp-file: cntc.inp_spck, ...
   %       character  fname(*)     - filename
   %  out: integer    ierror       - error flag

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 2: m=1, wtd    - no icp needed
   % #readinpfile -2


   CNTC = cntc.getmagicnumbers();

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2)
    inp_type = [];
   end
   if (~any(inp_type==[CNTC.inp_spck]))
    disp('ERROR(readinpfile): inp_type must be CNTC.inp_spck');
    return;
   end
   if (nargin<3 | isempty(c_fname))
    disp('ERROR(readinpfile): fname is mandatory')
    return;
   end

   p_ierr = libpointer('int32Ptr',-1);

   p_ierr.value = 0;
   cntc.call('cntc.readinpfile', ire, inp_type, c_fname, length(c_fname), p_ierr);

   ierror = double(p_ierr.value);

   if (ierror~=0)
    disp(sprintf('cntc.readinpfile: an error occurred in the CONTACT library (%d).',ierror));
   end

  end % cntc.readinpfile

  %------------------------------------------------------------------------------------------------------------

  function []=resetcalculationtime(ire, icp)
   % [ ] = cntc.resetcalculationtime(ire, icp)
   %
   % reset the accumulated cpu-time and wall-clock-time used for a contact problem


   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 6: m=*, cp     - require icp>0, default 1
   % #resetcalculationtime -2



   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.resetcalculationtime: not available for icp=%d',icp));
    return
   end

   cntc.call('cntc.resetcalculationtime', ire, icp);

  end % cntc.resetcalculationtime

  %------------------------------------------------------------------------------------------------------------

  function []=closelibrary();

   % [ ] = cntc.closelibrary();
   %
   % clean-up, close files and unload the library from Matlab


   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 0: m=*, glob   - no icp needed
   % #closelibrary -2

   LI=cntc.call;
   cntc.call('cntc.finalizelast');
   %dbstack;keyboard;
   unloadlibrary(LI.libname);
   LI.libname='';

   %evalin('base','clear("LI")');
  end % cntc.closelibrary

  %------------------------------------------------------------------------------------------------------------

  function []=finalize(ire)
   %%------------------------------------------------------------------------------------------------------------
   % [ ] = cntc.finalize(ire)
   %
   % Finalize calculations and clean-up for a result element


   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 0: m=*, glob   - no icp needed
   % #finalize -2

   if (nargin<1 | isempty(ire))
    ire = 1;
   end

   cntc.call('cntc.finalize', ire);

  end % cntc.finalize

  %------------------------------------------------------------------------------------------------------------

  function out=call(varargin)
   %% #call catch library calls to set variables  ------------------------- -2

   persistent LI

   if isempty(LI)
    if evalin('base','exist(''LI'',''var'')')
     %% use base workspace version if exists
     LI=evalin('base','LI');
    else
     % Initialize if does not exist
     LI=vhandle.uo([]); assignin('base','LI',LI);
     LI.callLog=vhandle.nmap([],[],'calls and parameters');
    end
   end
   if nargin==1&&isa(varargin{1},'vhandle.uo')&&isfield(varargin{1},'callLog')
    LI=varargin{1};return
   end
   r1=struct;
   for j1=1:nargin
    %% store inputs by name
    st1=inputname(j1);
    if ~isempty(st1);r1.(st1)=varargin{j1};end
   end
   if nargin==0; out=LI;return;end
   CAM=varargin{1};
   if ~isprop(LI,'callLog');addprop(LI,'callLog');LI.callLog=vhandle.nmap;end
   if isa(LI.callLog,'vhandle.nmap')
    r1.now=now;r1.dbstack=dbstack; %#ok<TNOW1>
    LI.callLog(CAM)=r1; % Store calls to allow logging
    diary off;f1='C:\Users\0021221S.COMMUN\MATLAB\xxx.log';
    if exist(f1,'file');diary(f1);end
    % debug message
    % disp(CAM);sdtm.toString(r1)
    % disp(f1);
   end

   if strncmpi(CAM,'contact_addon',13)
    LI.libname=CAM;varargin(1)=[];CAM=varargin{1};
   end
   CAM=strrep(CAM,'cntc.','cntc_');
   if ~isprop(LI,'libname')||isempty(LI.libname);
    if strcmpi(CAM,'cntc_finalizelast');return
    else
     error('Library not initialized');
    end
   end
   if strcmp(CAM,'finalizelast') % already loaded
    try
     calllib(LI.libname,'cntc_finalizelast');
    end
   else
    calllib(LI.libname,CAM,varargin{2:end});
   end
   if nargout>0; out=LI;end

  end

  %------------------------------------------------------------------------------------------------------------

  function out=help(CAM)
   %% #help : configure help -2
   wd0=pwd;
   if ispref('SDT','CONTACT_Path');wd=getpref('SDT','CONTACT_Path');
   else
    wd=sdtu.f.firstdir({'D:\APP\win64\contact_v24.1\bin';
     'C:\Program Files\Vtech CMCC\contact_v24.1\bin';
     '/o/APP/contact_v24.1/bin'});
    setpref('SDT','CONTACT_Path',wd);
   end
   cd(wd0);
   f1=sdtu.f.cffile(fullfile(wd,'../doc/user-guide.pdf'));
   if exist(f1,'file')
    RO=sdtu.f.ppath(struct,{'sdtu.idx.tagI'});
    ev1=sdtu.f.WhichName(f1);
    if ~isKey(RO.tagI.fileM,ev1.FileName)
     evt=struct('FileName',f1,'list',{{...
      'setflags(section.2.3) setflags configuration of flags, control digits'
      'setsolverflags(section.4.6) - configuration of solver parameters'
      }});
     sdtu.idx.addFileMatch(evt);
    end
    % xxx should test if needed
   end
   if nargin==0
   elseif strcmpi(CAM,'@examples')
     out=sdtu.f.cffile(fullfile(wd),'../examples');
   elseif strcmpi(CAM,'pdf')
    out=f1;
   elseif strcmpi(CAM,'md')
    out=sdtu.f.cffile(sdtu.f.safe('@cntc.m/../../rail/jup/cntc.md'));
   end

  end

  % xxx should test if needed
  %  end
  %  if nargin==0
  %  elseif strcmpi(CAM,'pdf')
  %    out=f1;
  %  elseif strcmpi(CAM,'md')
  %    out=sdtu.f.cffile(sdtu.f.safe('@cntc.m/../../rail/jup/cntc.md'));
  %    if nargout==0; sdtu.idx.jup(out);clear out;end
  %  end
  %
  % end


  function LI=init
   %% #init -2

   % if (~exist('cntc_initlibrary.m','file'))
   %    % set location of CONTACT installation folder
   %    % contactdir = 'C:\Program Files\Vtech CMCC\contact_v24.1';
   %    contactdir = '..';
   %    addpath([contactdir, '\matlab_intfc']);
   %    addpath([contactdir, '\matlab']);
   % end

   %if (exist('cntc_initlibrary')~=2)
   %   disp('ERROR: cant find CONTACT library on the Matlab search path');
   %   return
   %end
   cntc.help;
   [ CNTC, ifcver, ierror ]=cntc.initlibrary;
   LI=cntc.call; LI.CNTC=CNTC;

  end
  function out=cvs(obj); % #cvs -2
   out=sdtcheck('revision','$Revision: b515365 $  $Date: 2025-05-23 15:38:11 +0200 $');
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% #Spline

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ v_out, x_out, y_out, z_out ] = eval_2dspline( spl2d, u_in, v_in, y_in, idebug )

   % [ v_out, x_out, y_out, z_out ] = eval_2dspline( spl2d, u_in, v_in, y_in, idebug )
   %
   % evaluate parametric 2d spline {tui, tvj, cij_x, cij_y, cij_z}
   %  - if v_in is given, determine [xyz]_out(u_in, v_in) for given (u_in, v_in)
   %  - if y_in is given, determine v_out and corresponding x_out,z_out(u_in, v_out)
   %
   % the inputs (u,v) or (u,y) are considered a list when of the same size, else a tensor grid

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % #eval_2dspline -2
   if (isstruct(spl2d) & isfield(spl2d, 'slc_file'))
    slcs = spl2d;
    spl2d = slcs.spl2d;
   end
   if (nargin<3)
    v_in = [];
   end
   if (nargin<4)
    y_in = [];
   end
   if (nargin<5 | isempty(idebug))
    idebug = 0;
   end
   if ( (isempty(v_in) & isempty(y_in)) | (~isempty(v_in) & ~isempty(y_in)) )
    disp('ERROR: either v_in or y_in must be given');
    return;
   end
   if (min(size(u_in))>1 || min(size(v_in))>1 || min(size(y_in))>1)
    disp('ERROR: u_in, v_in, y_in must be 1D arrays');
    return;
   end
   if (size(u_in,2)>size(u_in,1))
    u_in = u_in';     % make column vector
   end
   if (size(v_in,1)>size(v_in,2))
    v_in = v_in';     % make row vector
   end

   if (spl2d.use_cylindr)
    u_in = cntc.apply_wrap_around( spl2d, u_in );
   end

   if (~isempty(v_in))

    % forward calculation u_in, v_in --> x_out, y_out, z_out

    % disp('u_in specified, forward (u,v) --> (x,y,z)')
    v_out = v_in;
    [ x_out, y_out, z_out ] = cntc.eval_2dspline_forward(spl2d, u_in, v_in, idebug);

   else

    % inverse calculation u_in, y_in --> v_out

    % disp('y_in specified, inverse (u,y) --> (v)')
    v_out = cntc.eval_2dspline_inverse( spl2d, u_in, y_in, idebug );

    % for tensor grid nu x nv, reshape 2d v_out to 1d list
    nu = length(u_in); ny = length(y_in);
    if (nu~=ny)
     u_in  = repmat(u_in, ny, 1);
     v_out = reshape(v_out, nu*ny, 1);
    end

    % forward calculation v_out --> x_out, z_out

    % disp('y_in specified, forward (u,v) --> (x,z)')
    y_out = y_in;
    [ x_out, ~, z_out ] = cntc.eval_2dspline_forward(spl2d, u_in, v_out, idebug);
    if (nu~=ny)
     x_out = reshape(x_out, nu, ny);
     y_out = repmat(y_in, nu, 1);
     z_out = reshape(z_out, nu, ny);
    end

   end

  end % eval_2dspline

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [x_out, y_out, z_out] = eval_2dspline_forward(spl2d, u_in, v_in, idebug)

   %[x_out, y_out, z_out] = eval_2dspline_forward(spl2d, u_in, v_in, idebug)
   %
   % if u_in and v_in have the same size, evaluate 2d spline at pairs (u_in(i), v_in(i))
   %                                else, evaluate 2d spline at tensor grid u_in x v_in.
   % #eval_2dspline_forward -2
   has_xij = isfield(spl2d, 'cij_x');

   if (all(size(u_in)==size(v_in)))
    use_list = 1; typ_eval = 'list';
   else
    use_list = 0; typ_eval = 'tens.grid';
   end
   k = 4;

   if (idebug>=1)
    disp(sprintf('Forward spline evaluation (u,v) --> (x,y,z), u: %d x %d, v: %d x %d (%s)', ...
     size(u_in), size(v_in), typ_eval));
   end

   % determine collocation matrix for v-direction: evaluate each B-spline at each output location
   %  - check basic interval [tj(4),tj(nknot-k+1)]

   n0 = nnz(v_in < spl2d.tvj(4));
   n1 = nnz(v_in > spl2d.tvj(end-3));
   if (n0+n1>0)
    disp(sprintf('Warning: there are %d v-positions before and %d after spline v-range [%3.1f,%3.1f]', ...
     n0, n1, spl2d.tvj(4), spl2d.tvj(end-3)));
   end

   numnan = nnz(isnan(v_in));
   if (idebug>=3 | (idebug>=1 && numnan>0))
    disp(sprintf('array v_in has %d NaN-values', nnz(isnan(v_in))));
   end

   [~, ~, ~, Bmat] = cntc.eval_bspline_basisfnc( spl2d.tvj, v_in, k );
   % nknotu = length(spl2d.tui); nsplu  = nknotu - k;
   % nknotv = length(spl2d.tvj); nsplv  = nknotv - k;

   % matrix-multply to get 1d spline coefficients for u-direction

   % ci_y: [ nsplu, nout ], Bmat: [ nout, nsplv ], cij_y: [ nsplu, nsplv ], v_in: [ nout ]
   if (has_xij)
    ci_x = (Bmat * spl2d.cij_x')';
   end
   ci_y = (Bmat * spl2d.cij_y')';
   ci_z = (Bmat * spl2d.cij_z')';
   mask = (Bmat * spl2d.mask_j')';

   % determine collocation matrix for u-direction: evaluate each B-spline at each output location

   tiny  = 1e-10;
   n0 = nnz(u_in < spl2d.tui(4)-tiny);
   n1 = nnz(u_in > spl2d.tui(end-3)+tiny);
   if (n0+n1>0)
    disp(sprintf('Warning: there are %d u-positions before and %d after spline u-range [%3.1f,%3.1f]', ...
     n0, n1, spl2d.tui(4), spl2d.tui(end-3)));
   end

   [~, ~, ~, Bmat] = cntc.eval_bspline_basisfnc( spl2d.tui, u_in, k );

   % matrix-multply to get output values

   % y_out: [ noutu, noutv ], Bmat: [ noutu, nsplu ], ci_y: [ nsplu, noutv ], u_in: [ noutu ]
   if (has_xij)
    x_out = Bmat * ci_x;
   else
    x_out = u_in * ones(1,length(v_in));       % half-parametric spline: x==u
   end
   y_out = Bmat * ci_y;
   z_out = Bmat * ci_z;
   mask  = Bmat * mask;

   % remove points that don't have full support

   ix    = find(mask<1-tiny);
   x_out(ix) = NaN;
   y_out(ix) = NaN;
   z_out(ix) = NaN;

   % output a list when both inputs are of same size

   if (use_list)
    x_out = diag(x_out);
    y_out = diag(y_out);
    z_out = diag(z_out);
    if (any(size(y_out)~=size(u_in)))
     x_out = x_out';
     y_out = y_out';
     z_out = z_out';
    end
   end

  end %eval_2dspline_forward

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function [ v_out ] = eval_2dspline_inverse(spl2d, u_in, y_in, idebug)

   %[ v_out ] = eval_2dspline_inverse(spl2d, u_in, y_in, idebug)
   %
   % if u_in and y_in have the same size, invert 2d spline at pairs (u_in(i), y_in(i))
   %                                else, invert 2d spline at tensor grid u_in x y_in.
   % #eval_2dspline_inverse -2
   tiny   = 1e-10 * max(abs(spl2d.ui));
   k      = 4;
   nknotu = length(spl2d.tui);
   nknotv = length(spl2d.tvj);
   nsplu  = nknotu - k;
   nsplv  = nknotv - k;

   if (all(size(u_in)==size(y_in)))
    use_list = 1; typ_eval = 'list';
    nu    = length(u_in);
    ny    = 1;
    v_out = ones(size(u_in));
   else
    use_list = 0; typ_eval = 'tens.grid';
    nu    = length(u_in);
    ny    = length(y_in);
    v_out = ones(nu, ny);
   end

   if (idebug>=1)
    disp(sprintf('Inverse spline evaluation (u,y) --> (v),   u: %d x %d, y: %d x %d (%s)', ...
     size(u_in), size(y_in), typ_eval));
   end

   % check basic interval in u-direction, knots [t(4), t(nknot-k+1)]

   if (any(u_in<spl2d.tui(k)) | any(u_in>spl2d.tui(nknotu-k+1)))
    disp(sprintf('ERROR: all input u_in must lie in range of spline u=[%3.1e,%3.1e]', ...
     spl2d.tui(k), spl2d.tui(nknotu-k+1)));
    return;
   end

   % determine collocation matrix for computing coefficients cj_y

   [~, ~, ~, Bmatx] = cntc.eval_bspline_basisfnc( spl2d.tui, u_in );

   % loop over output positions u_in

   for iout = 1 : nu

    % determine spline coefficients cj_y at position u_in(iout)

    % cj_y: [ nsplv, 1 ], Bmat: [ nu, nsplu ], cij_y: [ nsplu, nsplv ], u_in: [ nu ]

    cj_y = (Bmatx(iout,:) * spl2d.cij_y)';

    % n+k knots gives n basisfunctions and n spline coefficients;
    % basisfunction j becomes nonzero at knot j until knot j+k.
    % values of y in segment j = [t_j,t_j+1) are defined by spline coefficients j-k+1:j
    % n+k knots gives n+k-1 segments. The first k-1 and last k-1 segments are outside the basic interval
    % this gives n-k+1 segments in the basic interval with numbers k to n

    % determine interval [ylow, yhig] for segments j, j = k to n

    tmp            = NaN * ones(4, nknotv);
    tmp(:,k:nsplv) = [ cj_y(1:nsplv-3)';  cj_y(2:nsplv-2)';  cj_y(3:nsplv-1)';  cj_y(4:nsplv)' ];

    ylow = min(tmp); yhig = max(tmp);

    % no inner loop when using list-input, size(u_in) == size(y_in)

    if (use_list)
     j0 = iout; j1 = iout;
    else
     j0 = 1; j1 = ny;
    end

    for jout = j0 : j1

     % determine segments jseg that may contain y(jout) at xout

     y    = y_in(jout);
     jseg = find( ylow(1:nsplv)<=y & y<=yhig(1:nsplv) );
     nseg = length(jseg);

     if (nseg<=0)
      % disp(sprintf('ERROR: no segments for (u,y) = (%6.1f,%6.1f)', u_in(iout), y));
      found = 0;
     else
      if (idebug>=2)
       disp(sprintf(['(i,j)=(%2d,%2d), (u,y)=(%6.1f,%6.1f): found %d possible segments jseg = ', ...
        '%d, %d, %d, %d, %d, %d'], iout, jout, u_in(iout), y_in(jout), nseg, jseg(1:min(6,nseg))));
       if (idebug>=-4)

        [~, ~, ~, Bmatv] = cntc.eval_bspline_basisfnc( spl2d.tvj, spl2d.vj );
        % ytmp: [ nv, 1 ], Bmatv: [ nv, nv ], cj_y: [ nv, 1 ]
        ytmp = Bmatv * cj_y;
        tvj  = spl2d.tvj;
        tgrev(1:nsplv) = (tvj(2:nsplv+1) + tvj(3:nsplv+2) + tvj(4:nsplv+3)) / 3;
        tm   = (tgrev(1:nsplv-1) + tgrev(2:end)) / 2;

        figure(7); clf; hold on;
        tmpt  = reshape( [1;1]*tvj(1:nsplv+1), 1, 2*nsplv+2 );
        tmpyl = reshape( [1;1]*ylow(1:nsplv), 1, 2*nsplv );
        tmpyh = reshape( [1;1]*yhig(1:nsplv), 1, 2*nsplv );
        plot(tmpt(2:end-1), tmpyl);
        plot(tmpt(2:end-1), tmpyh);
        plot(spl2d.vj([1,end]), y*[1 1], '--');
        plot(spl2d.vj, ytmp);
        plot([1;1]*(tvj(jseg)+tvj(jseg+1))/2, [ylow(jseg);yhig(jseg)], 'color',cntc.matlab_color(3))
        plot(spl2d.tvj(nsplv+1), cj_y(end), '*');
        for jj = jseg
         text( (tvj(jj)+tvj(jj+1))/2, ylow(jj), sprintf('j=%d',jj), ...
          'horizontalalignment','center', 'verticalalignment','top');
        end
        grid on
        xlabel('v_j'); ylabel('y_r');
        legend('y_{low}(v;u)', 'y_{hig}(v;u)', 'target y_{out}', 'actual y(v), u=u_{out}', ...
         'location','northwest')
       end
      end

      % try potential segments one-by-one until a solution is found

      jj = 1;
      found = 0;
      while(jj<=nseg & ~found)

       % spline segment: interval [v_j, v_{j+1}]

       vseg = spl2d.tvj( jseg(jj)+[0,1] );

       % PP-spline coefficients for segment:

       coef  = cntc.get_ppcoef_1seg_sparse(spl2d, cj_y, jseg(jj));
       if (0==1)
        coef0 = cntc.get_ppcoef_1seg_full(spl2d, cj_y, jseg(jj));
        coef1 = cntc.get_ppcoef_1seg_sparse(spl2d, cj_y, jseg(jj));
        if (max(abs(coef0-coef1)>1e-10))
         disp('coef:')
         disp([coef0; coef1])
        end
       end

       % Solve cubic equation for segment:

       [v, found] = cntc.solve_cubic_eq(coef, vseg, y, idebug, iout, jout, jseg(jj));
       if (found & idebug>=2)
        disp(sprintf('jj = %d: found v =%7.2f', jj, v));
       elseif (~found & idebug>=3)
        disp(sprintf('jj = %d: no solution',jj));
       end

       found = 1;
       if (~found), jj = jj + 1; end
      end % jj

     end % nseg>0

     % store solution

     if (~found)
      if (idebug>=2)
       disp(sprintf('No solution for (i,j)=(%d,%d), setting NaN', iout, jout));
      end
      v = NaN;
     end
     if (use_list)
      v_out(iout) = v;
     else
      v_out(iout,jout) = v;
     end

    end % for jout
   end % for iout

  end %eval_2dspline_inverse

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ th_ev ] = apply_wrap_around( spl2d, th_ev )

   % - apply wrap-around th_ev --> [-pi,pi)
   % - restrict th_ev \in spline range [th0,th1] to obtain constant extrapolation
   % #apply_wrap_around -3

   th_ev = mod(th_ev+pi, 2*pi) - pi;

   th0   = spl2d.tui(4);
   th1   = spl2d.tui(end-3);
   th_ev = max(th0, min(th1, th_ev));

  end %apply_wrap_around

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % #get_ppcoef -2
  function [ coef ] = get_ppcoef_1seg_full(spl2d, cj_y, jseg)

   %[ coef ] = get_ppcoef_1seg_full(spl2d, cj_y, jseg)
   %
   % determine PP-spline coefficients for segment jseg

   % average step sizes over 1/2/3 adjacent intervals
   % #get_ppcoef_1seg_full -3
   nknot  = length(spl2d.tvj);
   dtj_k1 = (spl2d.tvj(2:end) - spl2d.tvj(1:nknot-1))';
   dtj_k2 = (spl2d.tvj(3:end) - spl2d.tvj(1:nknot-2))' / 2;
   dtj_k3 = (spl2d.tvj(4:end) - spl2d.tvj(1:nknot-3))' / 3;

   % inverse average step sizes, zero at zero step size

   dtj_inv1 = zeros(size(dtj_k1)); ix  = find(dtj_k1>0); dtj_inv1(ix) = 1./ dtj_k1(ix);
   dtj_inv2 = zeros(size(dtj_k2)); ix  = find(dtj_k2>0); dtj_inv2(ix) = 1./ dtj_k2(ix);
   dtj_inv3 = zeros(size(dtj_k3)); ix  = find(dtj_k3>0); dtj_inv3(ix) = 1./ dtj_k3(ix);

   % backward difference matrices using average step lengths 1/2/3

   D1   = spdiags(dtj_inv1, 0, nknot-1, nknot-2) - spdiags(dtj_inv1(2:end), -1, nknot-1, nknot-2);
   D2   = spdiags(dtj_inv2, 0, nknot-2, nknot-3) - spdiags(dtj_inv2(2:end), -1, nknot-2, nknot-3);
   D3   = spdiags(dtj_inv3, 0, nknot-3, nknot-4) - spdiags(dtj_inv3(2:end), -1, nknot-3, nknot-4);

   % evaluate function and derivative values at position vj(jseg)

   [B1, B2, B3, B4] = cntc.eval_bspline_basisfnc( spl2d.tvj, spl2d.tvj(jseg)+1e-9 );

   ay0 = B4 * cj_y;
   ay1 = B3 * D3 * cj_y;
   ay2 = B2 * D2 * D3 * cj_y / 2;
   ay3 = B1 * D1 * D2 * D3 * cj_y / 6;
   coef = [ay0, ay1, ay2, ay3];

   if (0==1)
    tmp = [B4,0,0,0;  B3,0,0;  B2,0;  B1];
    disp(full(tmp(:,jseg:jseg+3)))

    figure(5); clf;
    spy(tmp);
    c = get(gca,'children');
    set(c,'marker','*');
    axis([jseg-1 jseg+4 0 5]);
    axis normal
    hold on
    plot((jseg+3)*[1 1], [0 5], '--');
    grid on;
    title(sprintf('non-zeros for j = %d',jseg))
   end

  end %get_ppcoef_1seg_full

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ coef ] = get_ppcoef_1seg_sparse(spl2d, cj_y, jseg)

   % [ coef ] = get_ppcoef_1seg_sparse(spl2d, cj_y, jseg)
   %
   % determine PP-spline coefficients for segment jseg

   % average step sizes over 1/2/3 adjacent intervals
   % #get_ppcoef_1seg_sparse -3

   nknotv = length(spl2d.tvj);
   tvj    = spl2d.tvj(jseg-3:jseg+4);            % dtj_k3 uses (jseg+1)+3 - (jseg+1)
   dtj_k1 = (tvj(2:end) - tvj(1:end-1))';        % used for jseg  :jseg+1, last one * 0.0
   dtj_k2 = (tvj(3:end) - tvj(1:end-2))' / 2;    % used for jseg-1:jseg+1, last one * 0.0
   dtj_k3 = (tvj(4:end) - tvj(1:end-3))' / 3;    % used for jseg-2:jseg+1, last one * 0.0

   % inverse average step sizes, zero at zero step size

   dtj_inv1 = zeros(9,1); ix  = find(dtj_k1>0); dtj_inv1(ix) = 1./ dtj_k1(ix);
   dtj_inv2 = zeros(8,1); ix  = find(dtj_k2>0); dtj_inv2(ix) = 1./ dtj_k2(ix);
   dtj_inv3 = zeros(7,1); ix  = find(dtj_k3>0); dtj_inv3(ix) = 1./ dtj_k3(ix);
   % disp(dtj_inv3')

   % evaluate function and derivative values at position tvj(jseg)

   B1   = zeros(1,5);
   B2   = zeros(1,5);
   B3   = zeros(1,5);
   B4   = zeros(1,5);

   u = tvj(4)+1e-9;
   for jj = 4 : 4       % B(j;1) is nonzero for jseg only
    B1(1,jj)    = (u >= tvj(jj)  & u < tvj(jj+1));
   end
   for jj = 3 : 4       % B(j;2) uses B(j+1;1), nonzero for jseg-1 and jseg,
    % using dt1(jseg:jseg), using 0 * dt1(jseg-1) and 0 * dt1(jseg+1)
    B2(1,jj) = ( (u-tvj(jj  )) .* B1(:,jj  ) * dtj_inv1(jj  )  +  ...
     (tvj(jj+2)-u) .* B1(:,jj+1) * dtj_inv1(jj+1) ) / 1;
   end
   for jj = 2 : 4       % B(j;3) uses B(j+1;2), nonzero for jseg-2 to jseg, using dt2(jseg-1:jseg)
    B3(1,jj) = ( (u-tvj(jj  )) .* B2(:,jj  ) * dtj_inv2(jj  )  +  ...
     (tvj(jj+3)-u) .* B2(:,jj+1) * dtj_inv2(jj+1) ) / 2;
   end
   for jj = 1 : 4       % B(j;4) uses B(j+1;3), nonzero for jseg-3 to jseg, using dt3(jseg-2:jseg)
    B4(1,jj) = ( (u-tvj(jj  )) .* B3(:,jj  ) * dtj_inv3(jj  )  +  ...
     (tvj(jj+4)-u) .* B3(:,jj+1) * dtj_inv3(jj+1) ) / 3;
   end
   % disp(B1); disp(B2); disp(B3); disp(B4)

   d1_cy = zeros(4,1);
   d2_cy = zeros(4,1);
   d3_cy = zeros(4,1);
   d1_cy(2:4) = diff( cj_y(jseg-3:jseg)) .* dtj_inv3(2:4);  % using dt3(jseg-2:jseg)
   d2_cy(3:4) = diff( d1_cy(    2:   4)) .* dtj_inv2(3:4);  % using dt2(jseg-1:jseg)
   d3_cy(4:4) = diff( d2_cy(    3:   4)) .* dtj_inv1(4:4);  % using dt1(jseg  :jseg)

   ay0 = B4(1,1:4) *  cj_y(jseg-3:jseg);
   ay1 = B3(1,2:4) * d1_cy(     2:   4);
   ay2 = B2(1,3:4) * d2_cy(     3:   4) / 2;
   ay3 = B1(1,4:4) * d3_cy(     4:   4) / 6;
   coef = [ay0, ay1, ay2, ay3];

  end % get_ppcoef_1seg_sparse

  function [ dx_out, dy_out, dz_out ] = eval_2dspline_deriv( spl2d, u_in, v_in, idir, idebug )

   % [ dx_out, dy_out, dz_out ] = eval_2dspline_deriv( spl2d, u_in, v_in, idir, idebug )
   %
   % evaluate 1st derivatives of parametric 2d spline {tui, tvj, cij_x, cij_y, cij_z}
   %    for parameter idir = 1 (u) or 2 (v)
   %
   % the inputs (u,v) are considered a list when of the same size, else a tensor grid

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % #eval_2dspline_deriv -2
   if (nargin<4 | isempty(idir))
    idir = 1;
   end
   if (nargin<5 | isempty(idebug))
    idebug = 0;
   end
   if (min(size(u_in))>1 | min(size(v_in))>1)
    disp('ERROR: u_in, v_in must be 1D arrays');
    return;
   end

   has_xij = isfield(spl2d, 'cij_x');

   if (all(size(u_in)==size(v_in)))
    use_list = 1; typ_eval = 'list';
   else
    use_list = 0; typ_eval = 'tens.grid';
   end
   if (size(u_in,2)>size(u_in,1))
    u_in = u_in';     % make column vector
   end
   if (size(v_in,1)>size(v_in,2))
    v_in = v_in';     % make row vector
   end
   k = 4;

   nknotu = length(spl2d.tui);
   nknotv = length(spl2d.tvj);
   nsplu  = nknotu - k;
   nsplv  = nknotv - k;

   if (idebug>=1)
    disp(sprintf('2d spline evaluation (u,v) --> (dx,dy,dz), u: %d x %d, v: %d x %d (%s), idir=%d', ...
     size(u_in), size(v_in), typ_eval, idir));
   end

   if (idir==1)

    % determine derivatives dx/du, dy/du, dz/du

    %  1. form noutv u-splines at output-locations v_in

    [~, ~, ~, B4] = cntc.eval_bspline_basisfnc( spl2d.tvj, v_in, k, idebug );

    % ci_y: [ nsplu, noutv ], B4: [ noutv, nsplv ], cij_y: [ nsplu, nsplv ], v_in: [ noutv ]
    if (has_xij)
     ci_x = (B4 * spl2d.cij_x')';
    end
    ci_y = (B4 * spl2d.cij_y')';
    ci_z = (B4 * spl2d.cij_z')';
    mask = (B4 * spl2d.mask_j')';

    %  2.a collocation matrix for 1st derivative of u-splines at u_in

    [~, ~, B3, B4] = cntc.eval_bspline_basisfnc( spl2d.tui, u_in, k, idebug );

    %  2.b backward difference matrix using average step lengths

    dtj_k3 = (spl2d.tui(4:end) - spl2d.tui(1:nknotu-3)) / 3;
    dtj_inv3 = zeros(size(dtj_k3)); ix  = find(dtj_k3>0); dtj_inv3(ix) = 1./ dtj_k3(ix);
    D3     = spdiags(dtj_inv3, 0, nknotu-3, nknotu-4) - spdiags(dtj_inv3(2:end), -1, nknotu-3, nknotu-4);

    %  2. evaluate 1st derivative of u-splines at u_in

    % y_out: [ noutu, noutv ], B3: [ noutu, nsplu+1 ], D3: [nsplu+1, nsplu],
    %                                                  ci_y: [ nsplu, noutv ], u_in: [ noutu ]

    dy_out = B3 * D3 * ci_y;
    if (has_xij)
     dx_out = B3 * D3 * ci_x;
    else
     dx_out = ones(size(dy_out));    % half-parametric spline: x==u
    end
    dz_out = B3 * D3 * ci_z;
    mask   = B4 * mask;

   else

    % determine derivatives dx/dv, dy/dv, dz/dv

    %  1. form noutu v-splines at output-locations u_in

    [~, ~, ~, B4] = cntc.eval_bspline_basisfnc( spl2d.tui, u_in, k, idebug );

    % cj_y: [ noutu, nsplv ], B4: [ noutu, nsplu ], cij_y: [ nsplu, nsplv ], u_in: [ noutu ]

    if (has_xij)
     cj_x = B4 * spl2d.cij_x;
    end
    cj_y = B4 * spl2d.cij_y;
    cj_z = B4 * spl2d.cij_z;
    mask = B4 * spl2d.mask_j;

    %  2.a collocation matrix for 1st derivative of v-splines at v_in

    [~, ~, B3, B4] = cntc.eval_bspline_basisfnc( spl2d.tvj, v_in, k, idebug );

    %  2.b backward difference matrix using average step lengths

    dtj_k3 = (spl2d.tvj(4:end) - spl2d.tvj(1:nknotv-3))' / 3;
    dtj_inv3 = zeros(size(dtj_k3)); ix  = find(dtj_k3>0); dtj_inv3(ix) = 1./ dtj_k3(ix);
    D3     = spdiags(dtj_inv3, 0, nknotv-3, nknotv-4) - spdiags(dtj_inv3(2:end), -1, nknotv-3, nknotv-4);

    %  2. evaluate 1st derivative of v-splines at v_in

    % y_out: [ noutu, noutv ], B3: [ noutv, nsplv+1 ], D3: [nsplv+1, nsplv],
    %                                                  cj_y: [ noutu, nsplv ], v_in: [ noutv ]

    dy_out = (B3 * D3 * cj_y')';
    if (has_xij)
     dx_out = (B3 * D3 * cj_x')';
    else
     dx_out = zeros(size(dy_out));    % half-parametric spline: x==u
    end
    dz_out = (B3 * D3 * cj_z')';
    mask   = (B4 * mask')';

   end

   % remove points that don't have full support

   tiny  = 1e-10;
   ix    = find(mask<1-tiny);
   dx_out(ix) = NaN;
   dy_out(ix) = NaN;
   dz_out(ix) = NaN;

   % output a list when both inputs are of same size

   if (use_list)
    dx_out = diag(dx_out);
    dy_out = diag(dy_out);
    dz_out = diag(dz_out);
    if (any(size(dy_out)~=size(u_in)))
     dx_out = dx_out';
     dy_out = dy_out';
     dz_out = dz_out';
    end
   end

  end % eval_2dspline_deriv

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ B1, B2, B3, B4, B5, B6, B7, B8 ] = eval_bspline_basisfnc( tj, si, k, idebug )

   % [ B1, B2, B3, B4, B5, B6, B7, B8 ] = eval_bspline_basisfnc( tj, si, [k], [idebug] )
   %
   % evaluate B-spline basis functions B_j,k for knot-vector tj at sample locations si
   %  k = spline order, default k=4 for cubic splines
   %
   % B1 = basis functions for k=1 using knot-vector tj  --  size (nout, nknot-1)
   % B2 = basis functions for k=2 using knot-vector tj  --  size (nout, nknot-2)
   % B3 = basis functions for k=3 using knot-vector tj  --  size (nout, nknot-3)
   % B4 = basis functions for k=4 using knot-vector tj  --  size (nout, nknot-4)

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % #eval_bspline_basisfnc -2
   if (min(size(si))>1)
    disp('ERROR(eval_bspline_basisfnc): si must be 1d array');
    return;
   end
   if (nargin<3 | isempty(k))
    k      = 4;
   end
   if (nargin<4 | isempty(idebug))
    idebug = 1;
   end
   tiny  = 1e-10;
   nknot = length(tj);
   npnt  = length(si);

   % determine value of each basis-function at each measurement location si

   if (size(si,2)>size(si,1)), si = si'; end

   % check for values si outside basic interval [tj(4),tj(nknot-k+1)]

   if (idebug>=1)
    n0 = nnz(si < tj(k)-tiny);
    n1 = nnz(si > tj(nknot-k+1)+tiny);
    if (n0+n1>0)
     disp(sprintf('Warning: there are %d u-positions before and %d after basic interval t=[%3.1f,%3.1f]', ...
      n0, n1, tj(k), tj(nknot-k+1)));
    end
   end

   numnan = nnz(isnan(si));
   if (idebug>=3 | (idebug>=1 & numnan>0))
    disp(sprintf('array si has %d NaN-values', numnan));
   end

   for ik = 1 : k

    % initialize sparse array

    Bk = sparse(npnt, nknot-ik);

    % average step sizes over ik-1 adjacent intervals

    dtj_k   = (tj(ik:end) - tj(1:nknot-ik+1)) / (ik-1);

    % inverse average step sizes, zero at zero step size

    dtj_inv = zeros(size(dtj_k)); ix  = find(dtj_k>0); dtj_inv(ix) = 1./ dtj_k(ix);

    % determine index j for start of last 'true' interval.
    % j_last == nknot-k in case the end-knot is added k-1 times

    j_last = find(tj<tj(end), 1, 'last');

    if (ik==1) % ik=1: piecewise constant B_j,1

     for j = 1 : nknot-ik
      if (ik==1 & j==j_last)
       Bk(:,j)  = (si >= tj(j)  & si <= tj(j+1));
      else
       Bk(:,j)  = (si >= tj(j)  & si <  tj(j+1));
      end
     end

    else

     % ik=2: piecewise linear B_j,2

     % disp([ik nknot-ik length(dtj_k) length(dtj_inv)])

     for j = 1 : nknot-ik
      % disp([ik, j, length(tj) length(dtj_inv) size(Bp)])
      Bk(:,j) = ( (si - tj(j)) .* Bp(:,j) * dtj_inv(j) + ...
       (tj(j+ik) - si) .* Bp(:,j+1) * dtj_inv(j+1) ) / (ik-1);
     end

    end

    % store result in output B1, B2, etc

    eval(sprintf('B%d = Bk;', ik));

    % copy to Bp for previous B_{ik-1} in next iteration

    Bp = Bk;

   end

  end % eval_bspline_basisfnc


  function [ s_out, y_out, z_out ] = eval_spline( spl, s_in, y_in, idebug )

   % [ s_out, y_out, z_out ] = eval_spline( spl, s_in, y_in, idebug )
   %
   % evaluate parametric spline {s, a0y-a3y, a0z-a3z}
   %  - if s_in is given, determine y(s_in), z(s_in)
   %  - if y_in is given, determine corresponding s(y) and z(s(y))

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % #eval_spline -2
   if (nargin<2)
    s_in = [];
   end
   if (nargin<3)
    y_in = [];
   end
   if (nargin<4 | isempty(idebug))
    idebug = 0;
   end
   if ( (isempty(s_in) & isempty(y_in)) | (~isempty(s_in) & ~isempty(y_in)) )
    disp('ERROR: either s_in or y_in must be given');
    return;
   end

   if (~isempty(s_in))
    s_out = s_in;
   else
    s_out = cntc.spline_get_s_at_y( spl, y_in, idebug );
   end

   if (~isempty(y_in))
    y_out = y_in;
   else
    y_out = cntc.eval_1d_spline(spl.s, spl.ay0, spl.ay1, spl.ay2, spl.ay3, s_out, idebug);
   end

   if (nargout>=3)
    z_out = cntc.eval_1d_spline(spl.s, spl.az0, spl.az1, spl.az2, spl.az3, s_out, idebug);
   end

  end % eval_spline

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ s_out ] = spline_get_s_at_y(spl, yev, idebug)
   % #spline_get_s_at_y -2
   if (nargin<3 | isempty(idebug))
    idebug = 0;
   end

   if (size(yev,2) > size(yev,1)), yev = yev'; end;

   % make input y-values ascending
   if (yev(end) > yev(1))
    is_ascending = 1;
   else
    yev = flipud(yev);
    is_ascending = 0;
   end
   if (any(sort(yev)~=yev))
    disp('ERROR: input y_ev should be sorted (ascending or descending)');
    return
   end

   s  = spl.s;
   a0 = spl.ay0;
   a1 = spl.ay1;
   a2 = spl.ay2;
   a3 = spl.ay3;
   np = length(a0);

   tiny_a2 = 1e-10;
   tiny_a3 = 1e-10;
   s_out = zeros(size(yev));

   for iout = 1 : length(yev)

    % find segment iseg [ a0(i), a0(i+1) ] containing yev

    if (0==1 | ~isfield(spl, 'top_ybrk'))
     iseg = cntc.locate_segment( np, a0, yev(iout) );
    else
     npart = length(spl.top_ybrk) - 1;
     itop  = cntc.locate_segment( npart+1, spl.top_ybrk, yev(iout) );
     if (itop<=0)
      iseg = 0;
      %disp(sprintf('yev(%2d)=%8.3f lies in top part %d, y=   (-inf,%7.3f]', iout, yev(iout), itop, ...
      %       spl.top_ybrk(itop+1)));
     elseif (itop>=npart+1)
      iseg = np+1;
      %disp(sprintf('yev(%2d)=%8.3f lies in top part %d, y=[%7.3f, +inf)', iout, yev(iout), itop, ...
      %       spl.top_ybrk(end) ));
     else
      %disp(sprintf('yev(%2d)=%8.3f lies in top part %d, y=[%7.3f,%7.3f]', iout, yev(iout), itop, ...
      %       spl.top_ybrk(itop), spl.top_ybrk(itop+1)));
      isec = spl.top_sec(itop);
      ip0 = spl.isec_uni_y(isec);
      ip1 = spl.isec_uni_y(isec+1);
      % disp(sprintf('yev(%2d)=%7.3f lies in section %d, y=[%7.3f,%7.3f]', iout, yev(iout), isec, ...
      %           a0(ip0), a0(ip1)));
      iseg   = ip0-1 + cntc.locate_segment( ip1-ip0+1, a0(ip0:ip1), yev(iout) );
     end
    end

    if (idebug>=1)
     if (iseg<=0)
      disp(sprintf('Searching y(%2d)=%8.3f: iseg=%4d, sseg= ( +/- inf,%8.3f], yseg=( +/- inf,%8.3f]',...
       iout, yev(iout), iseg, s(1), a0(1)));
     elseif (iseg>=np)
      disp(sprintf('Searching y(%2d)=%8.3f: iseg=%4d, sseg= [%8.3f, +/- inf), yseg=[%8.3f, +/- inf)',...
       iout, yev(iout), iseg, s(np), a0(np)));
     else
      disp(sprintf('Searching y(%2d)=%8.3f: iseg=%4d, sseg= [%8.3f,%8.3f], yseg=[%8.3f,%8.3f]',...
       iout, yev(iout), iseg, s(iseg), s(iseg+1), a0(iseg), a0(iseg+1)));
     end
    end

    if (isempty(iseg) | iseg<=0)   % before start of spline range: yev < a0(1)

     s_out(iout) = s(1) - 1e-9;  % set s slightly before s(1)

    elseif (iseg>=np)           % after end of spline range: yev > a0(end)

     s_out(iout) = s(end) + 1e-9; % set s slightly after s(end)

    else

     % Solve cubic equation for segment:

     ppcoef = [a0(i), a1(iseg), a2(iseg), a3(iseg)];
     sseg   = s(iseg+[0,1]);
     [sout1, found] = cntc.solve_cubic_eq(ppcoef, sseg, yev(iout), idebug, iout, 1, iseg);

     if (found & idebug>=2)
      disp(sprintf('iout = %d: found sl =%7.2f', iout, sout1));
     elseif (~found & idebug>=3)
      disp(sprintf('iout = %d: no solution',iout));
     end

     s_out(iout) = sout1;

    end % yev in interior
   end % iout

   if (~is_ascending)
    s_out = flipud(s_out);
   end

  end % spline_get_s_at_y

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ iseg ] = locate_segment( np, vp, vi )
   % #locate_segment -2
   ascending = (vp(end) >= vp(1));

   if (ascending)
    if (vi < vp(1))
     iseg = 0;
    elseif (vi > vp(end))
     iseg = np + 1;
    else
     % find segment iseg: vp(iseg) <= vi <= vp(iseg+1)
     iseg = find( vp<=vi, 1, 'last');
    end
   else
    if (vi > vp(1))
     iseg = 0;
    elseif (vi < vp(end))
     iseg = np + 1;
    else
     % find segment iseg: vp(iseg) >= vi >= vp(iseg+1)
     iseg = find( vp>=vi, 1, 'last');
    end
   end

  end % locate_segment

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ y ] = eval_1d_spline(xsp, a0, a1, a2, a3, xev, idebug)
   % #eval_1d_spline -2
   % evaluate plain (non-parametric) spline y=y(x)

   iev_debug = 1;

   y = zeros(size(xev));

   % copy NaN inputs to output

   ix = find(isnan(xev));
   y(ix) = NaN;

   % assume that xev is sorted

   n = length(xsp) - 1;

   % all x < x(1): constant extrapolation

   i = find(xev<xsp(1));
   y(i) = a0(1);

   if (idebug>=5 & any(i==iev_debug))
    disp(sprintf('xev(%d) lies before start xsp(1)', iev_debug));
   end

   % all x(1) <= x <= x(n+1) : cubic polynomials

   for isp = 1 : n
    i = find(xev>=xsp(isp) & xev<xsp(isp+1));
    xloc = xev(i) - xsp(isp);
    y(i) = a3(isp) * xloc.^3 + a2(isp) * xloc.^2 + a1(isp) * xloc + a0(isp);

    if (idebug>=5 & any(i==iev_debug))
     j = find(i==iev_debug);
     disp(sprintf('xev(%d) lies in segment isp %d, xloc=%8.3f', iev_debug, isp, xloc(j)));
    end
   end

   % all x > x(n+1): constant extrapolation

   i = find(xev>=xsp(n+1));
   xloc = xsp(n+1) - xsp(n);
   y(i) = a3(n) * xloc.^3 + a2(n) * xloc.^2 + a1(n) * xloc + a0(n);

  end % eval_1d_spline

  function [ dy_out, ddy_out, dz_out, ddz_out, kappa, rcurv ] = eval_spline_deriv( spl, s_out, idebug )

   % [ dy_out, ddy_out, dz_out, ddz_out, kappa, rcurv ] = ...
   %                               eval_spline_deriv( spl, s_out, idebug )
   %
   % evaluate derivatives for parametric spline {s, a0y-a3y, a0z-a3z}
   %  - determine dy/ds, d2y/ds2, dz/ds and d2z/ds2

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % #eval_spline_deriv -2
   if (nargin<2 | isempty(s_out))
    s_out = spl.s;
   end
   if (nargin<3 | isempty(idebug))
    idebug = 0;
   end

   [dy_out, ddy_out] = cntc.eval_1d_spline_deriv(spl.s, spl.ay0, spl.ay1, ...
    spl.ay2, spl.ay3, s_out);

   [dz_out, ddz_out] = cntc.eval_1d_spline_deriv(spl.s, spl.az0, spl.az1, ...
    spl.az2, spl.az3, s_out);

   kappa = (dy_out.*ddz_out - dz_out.*ddy_out) ./ (dy_out.^2 + dz_out.^2).^(3/2);

   k_min = 1e-6;
   rcurv = zeros(size(kappa));
   ix = find(abs(kappa)>=k_min); rcurv(ix) = 1 ./ kappa(ix);
   ix = find(abs(kappa)< k_min); rcurv(ix) = sign(kappa(ix)) / k_min;

  end % eval_spline_deriv

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ dy, ddy ] = eval_1d_spline_deriv(xsp, a0, a1, a2, a3, xev)
   % #eval_1d_spline_deriv -2
   % evaluate plain (non-parametric) spline y=y(x)

   dy  = zeros(size(xev));
   ddy = zeros(size(xev));

   % assume that xev is sorted

   n = length(xsp) - 1;

   % all x < x(1): constant extrapolation

   i = find(xev<xsp(1));
   dy(i)  = a1(1);
   ddy(i) = a2(1);

   % all x(1) <= x <= x(n+1) : cubic polynomials

   for isp = 1 : n
    i = find(xev>=xsp(isp) & xev<xsp(isp+1));
    xloc = xev(i) - xsp(isp);
    dy(i)  = 3 * a3(isp) * xloc.^2 + 2 * a2(isp) * xloc    + a1(isp);
    ddy(i) = 6 * a3(isp) * xloc    + 2 * a2(isp);
   end

   % all x > x(n+1): constant extrapolation

   i = find(xev>=xsp(n+1));
   xloc = xsp(n+1) - xsp(n);
   dy(i)  = 3 * a3(n) * xloc.^2 + 2 * a2(n) * xloc    + a1(n);
   ddy(i) = 6 * a3(n) * xloc    + 2 * a2(n);

  end % eval_1d_spline_deriv

  %% #Load
  function [ s1, s2, s3, s4, s5, s6, s7, s8, s9, s10 ] = loadcase( expnam, icase, ipatch )
   %
   % [ s1, s2, ... ] = loadcase( expnam, [icase], [ipatch] )
   %
   % Read CONTACT mat-file(s) for one case and return struct(s) with its(their) contents.
   %
   % expnam    = the experiment name used, i.e. the filename without .mat extension
   % icase     = the case-number to be read
   % ipatch    = optional contact patch number, <=0 means all available
   %               0 = all available, returned as an array of structures [s1]
   %              -1 = all available, returned using separate structures  s1, s2, ...
   %             default: 0 when nargin<=1, -1 when nargin>1
   % s1,...    = structure(s) with results of the computation

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % #loadcase
   max_patch = 10;
   s1 = []; s2 = []; s3 = []; s4 = []; s5 = [];
   s6 = []; s7 = []; s8 = []; s9 = []; s10 = [];

   if (nargin<2 | isempty(icase))
    icase = 1;
   end
   if (nargin<3 | isempty(ipatch))
    if (nargout<=1)
     ipatch = 0;
    else
     ipatch = -1;
    end
   end

   if (ipatch<0)
    use_array = 0;
   else
    use_array = 1;
   end

   % check for the presence of output-files "experiment.<icase>[a-z].mat", ...

   has_output = zeros(max_patch,1);
   for jpatch = 1:max_patch
    if (icase<=9999)
     fname = [expnam, sprintf('.%04d%c.mat', icase, 96+jpatch)];
    else
     fname = [expnam, sprintf('.%06d%c.mat', icase, 96+jpatch)];
    end
    if (exist(fname,'file'))
     has_output(jpatch) = 1;
    end
   end

   % set default number of patches: all available

   if (max(has_output)<=0)
    npatch = 1;
   else
    npatch = max(find(has_output));
   end

   % set requested range of patches

   if (ipatch<=0)
    ipatch = [1 : npatch];
   end

   % loop over the requested patches

   for jpatch = ipatch

    % npatch may be reduced to value obtained from actual mat-files

    if (jpatch > npatch)
     continue
    end

    tmp = [];
    if (jpatch>1 & ~has_output(jpatch))

     % short-cut for [b-z]: no data found

     disp(sprintf('No output found for contact patch %d.',jpatch));
     continue

    elseif (has_output(jpatch))

     % load file "<icase>[a-z]" if present

     if (icase<=9999)
      fname = [expnam, sprintf('.%04d%c.mat', icase, 96+jpatch)];
     else
      fname = [expnam, sprintf('.%06d%c.mat', icase, 96+jpatch)];
     end
     tmp = load(fname, '-ascii');

    else

     % else, check for presence of file "experiment.0001.mat" or "experiment.000001.mat"

     if (icase<=9999)
      fname = [expnam, sprintf('.%04d.mat', icase)];
     else
      fname = [expnam, sprintf('.%06d.mat', icase)];
     end
     fname_ic = fname;

     if (exist(fname,'file'))

      tmp = load(fname, '-ascii');

     else

      % else, check for presence of file "experiment.mat"

      fname = [expnam,'.mat'];
      if (exist(fname,'file'))
       tmp = load(fname, '-ascii');
      else
       disp(['ERROR: cannot find file ',fname,' nor ',fname_ic]);
      end
     end
    end

    % create solution struct

    if (isempty(tmp))
     strc = [];
    else
     strc = cntc.fill_struct(tmp, fname);
    end

    % Each patch tells the actual npatch for the run in which it was created
    % warn & update npatch as needed

    if (~isempty(strc) & strc.meta.npatch>0 & strc.meta.npatch < npatch)
     disp(sprintf('Warning: found %d mat-files for run with %d contact patches. Left-over from earlier run?', ...
      npatch, strc.meta.npatch));
     npatch = min(npatch, strc.meta.npatch);
    end

    % put solution struct into output argument

    eval(sprintf('s%d = strc;', jpatch));

   end

   if (use_array)
    s1 = [s1, s2, s3, s4, s5, s6, s7, s8, s9, s10];
   end

  end % loadcase

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ sol ] = fill_struct(tmp, fname)
   % #fill_struct -2
   % check the file format

   fmtmat = tmp(1,end);
   if (fmtmat>2410)
    disp('ERROR: the mat-file is created with a newer version of CONTACT.');
    disp('       Please use the corresponding version of this Matlab script.');
    sol=[];
    return
   end

   numcol = size(tmp,2);
   irow = 1;

   sol = struct('config',[], 'h_digit',[], 'meta',struct(), 'mater',struct(), 'fric',struct(), ...
    'kincns',struct('t_digit',[]) );

   % extract the metadata from the first three rows

   if (fmtmat>=2320)
    sol.meta.tim      = tmp(irow,1);     % time: may be provided by calling program
    sol.meta.s_ws     = tmp(irow,2);     % s-position along track curve, wrt. inertial reference
    sol.meta.th_ws    = tmp(irow,3);
    sol.meta.y_r      = tmp(irow,4);     % position of rail marker wrt track reference
    sol.meta.z_r      = tmp(irow,5);
    sol.meta.roll_r   = tmp(irow,6);     % orientation of rail marker wrt track reference
    sol.meta.rnom_rol = tmp(irow,7);
    irow = irow + 1;

    sol.meta.x_w      = tmp(irow,1);     % position of wheel marker wrt track reference
    sol.meta.y_w      = tmp(irow,2);
    sol.meta.z_w      = tmp(irow,3);
    sol.meta.roll_w   = tmp(irow,4);     % orientation of wheel marker wrt track reference
    sol.meta.yaw_w    = tmp(irow,5);
    sol.meta.rnom_whl = tmp(irow,6);
    irow = irow + 1;

    sol.meta.xcp_r    = tmp(irow,1);     % position of contact ref.point wrt rail marker
    sol.meta.ycp_r    = tmp(irow,2);
    sol.meta.zcp_r    = tmp(irow,3);
    sol.meta.deltcp_r = tmp(irow,4);     % orientation of contact ref.point wrt rail marker
    sol.meta.xcp_w    = tmp(irow,5);     % position of contact ref.point wrt wheel marker
    sol.meta.ycp_w    = tmp(irow,6);
    sol.meta.zcp_w    = tmp(irow,7);
    sol.meta.npatch   = tmp(irow,8);
    sol.meta.ipatch   = tmp(irow,9);
    sol.meta.spinxo   = tmp(irow,10);
    sol.meta.spinyo   = tmp(irow,11);
    irow = irow + 1;
   elseif (fmtmat>=1801)
    sol.meta.tim      = tmp(irow,1);     % time: may be provided by calling program
    sol.meta.s_ws     = tmp(irow,2);     % s-position along track curve, wrt. inertial reference
    sol.meta.x_w      = tmp(irow,3);     % position of wheel marker wrt track reference
    sol.meta.y_w      = tmp(irow,4);
    sol.meta.z_w      = tmp(irow,5);
    sol.meta.roll_w   = tmp(irow,6);     % orientation of wheel marker wrt track reference
    sol.meta.yaw_w    = tmp(irow,7);
    sol.meta.pitch_w  = 0;
    sol.meta.y_r      = tmp(irow,8);     % position of rail marker wrt track reference
    sol.meta.z_r      = tmp(irow,9);
    sol.meta.roll_r   = tmp(irow,10);    % orientation of rail marker wrt track reference
    irow = irow + 1;

    sol.meta.xcp_r    = tmp(irow,1);     % position of contact ref.point wrt rail marker
    sol.meta.ycp_r    = tmp(irow,2);
    sol.meta.zcp_r    = tmp(irow,3);
    sol.meta.deltcp_r = tmp(irow,4);     % orientation of contact ref.point wrt rail marker
    sol.meta.xcp_w    = tmp(irow,5);     % position of contact ref.point wrt wheel marker
    sol.meta.ycp_w    = tmp(irow,6);
    sol.meta.zcp_w    = tmp(irow,7);
    sol.meta.npatch   = tmp(irow,8);
    sol.meta.ipatch   = tmp(irow,9);
    sol.meta.rnom_whl = tmp(irow,10);
    sol.meta.rnom_rol = tmp(irow,11);
    irow = irow + 1;
   else
    sol.meta.tim      = 0;
    sol.meta.s_ws     = 0;
    sol.meta.x_w      = 0;
    sol.meta.y_w      = 0;
    sol.meta.z_w      = 0;
    sol.meta.roll_w   = 0;
    sol.meta.yaw_w    = 0;
    sol.meta.pitch_w  = 0;
    sol.meta.y_r      = 0;
    sol.meta.z_r      = 0;
    sol.meta.roll_r   = 0;
    sol.meta.xcp_r    = 0;
    sol.meta.ycp_r    = 0;
    sol.meta.zcp_r    = 0;
    sol.meta.deltcp_r = 0;
    sol.meta.xcp_w    = 0;
    sol.meta.ycp_w    = 0;
    sol.meta.zcp_w    = 0;
    sol.meta.npatch   = 0;
    sol.meta.ipatch   = 0;
    sol.meta.rnom_whl = 460;
    sol.meta.rnom_rol = 800;
    if (fmtmat>=1401)
     sol.meta.tim      = tmp(irow,1);
     sol.meta.xcp_r    = tmp(irow,2);
     sol.meta.ycp_r    = tmp(irow,3);
     irow = irow + 1;
    end
   end

   % extract the grid discretisation parameters from the next row

   sol.mx  = tmp(irow, 1); mx = sol.mx;
   sol.my  = tmp(irow, 2); my = sol.my;
   sol.xl  = tmp(irow, 3);
   sol.yl  = tmp(irow, 4);
   sol.dx  = tmp(irow, 5);
   sol.dy  = tmp(irow, 6);
   sol.kincns.chi = tmp(irow, 7);
   sol.kincns.dq  = tmp(irow, 8);
   if (fmtmat>=1901)
    sol.config  = tmp(irow, 9);
    sol.d_digit = tmp(irow, 10);
   elseif (sol.meta.yrw<0)
    sol.config = 0; sol.d_digit = 0;
   else
    sol.config = 1; sol.d_digit = 0;
   end
   if (fmtmat>=2310)
    sol.meta.ynom_whl = tmp(irow,11);
   else
    sol.meta.ynom_whl = 0;
   end

   % x_offset, y_offset: shift, real world coordinates to be used for the
   %              origin of the contact coordinate system

   sol.x_offset = [];
   sol.y_offset = [];

   % construct the grid coordinates

   sol.x = sol.xl + ([1:mx]-0.5)*sol.dx;
   sol.y = sol.yl + ([1:my]-0.5)*sol.dy;

   % extract the material constants from the next row

   if (fmtmat<1220)
    disp('WARNING: the mat-file is created with an old version of CONTACT.');
    disp('         The material constants are not filled in.');
    sol.mater = struct('m_digit',0,'gg',[1 1],'poiss',[1,1]);
    sol.kincns.t_digit = 3;
   else
    irow = irow + 1;
    sol.kincns.t_digit = tmp(irow,1);
    sol.mater.m_digit  = tmp(irow,2);
    sol.mater.gg       = tmp(irow,3:4);
    sol.mater.poiss    = tmp(irow,5:6); % note: poiss(2) n.a. when m_digit==2 or 3
    if (sol.mater.m_digit==1)
     sol.mater.fg      = tmp(irow,7:8);
     sol.mater.tc      = tmp(irow,9:10);
    elseif (sol.mater.m_digit==2 | sol.mater.m_digit==3)
     sol.mater.poiss(2)= sol.mater.poiss(1); % note: poiss(2) n.a. when m_digit==2 or 3
     sol.mater.flx     = tmp(irow,6:8);
     sol.mater.k0_mf   = tmp(irow,9);
     sol.mater.alfamf  = tmp(irow,10);
     sol.mater.betamf  = tmp(irow,11);
    elseif (sol.mater.m_digit==4)
     sol.mater.gg3     = tmp(irow,7);
     sol.mater.laythk  = tmp(irow,8);
     sol.mater.tau_c0  = tmp(irow,9);
     sol.mater.k_tau   = tmp(irow,10);
    end

    g1 = sol.mater.gg(1); nu1 = sol.mater.poiss(1);
    g2 = sol.mater.gg(2); nu2 = sol.mater.poiss(2);
    sol.mater.ga = 2.0 / (1.0/g1 + 1.0/g2);
    sol.mater.si = 0.5 * sol.mater.ga * (nu1/g1 + nu2/g2);
    sol.mater.ak = 0.25 * sol.mater.ga * ((1-2*nu1)/g1 - (1-2*nu2)/g2);

   end
   if (sol.mater.m_digit~=4)
    has_plast = 0;
   else
    has_plast = (sol.mater.tau_c0>1e-10 & sol.mater.tau_c0<1e10);
   end

   % extract the friction law parameters from the next row

   irow = irow + 1;
   frclaw = tmp(irow, 1); sol.fric.frclaw = frclaw;
   if (frclaw<0 | frclaw>6)
    error(sprintf('The value for the L-digit (%d) in %s is incorrect. Maybe the file is generated with an older CONTACT version?',frclaw,fname));
    return
   end
   if (fmtmat<1220)
    idum=0;
    sol.h_digit = 0;
   else
    idum=1; % introduced dummy in 2nd column in version 12.2
    sol.h_digit = tmp(irow,2);
   end
   has_heat   = (sol.h_digit>=1);
   sol.kincns.veloc  = tmp(irow, 2+idum);
   if (frclaw == 0)
    fstat  = tmp(irow, 3+idum); sol.fric.fstat  = fstat;
    fkin   = tmp(irow, 4+idum); sol.fric.fkin   = fkin;
   elseif (frclaw == 2)
    fkin   = tmp(irow, 3+idum); sol.fric.fkin   = fkin;
    flin1  = tmp(irow, 4+idum); sol.fric.flin1  = flin1;
    sabsh1 = tmp(irow, 5+idum); sol.fric.sabsh1 = sabsh1;
    flin2  = tmp(irow, 6+idum); sol.fric.flin1  = flin1;
    sabsh2 = tmp(irow, 7+idum); sol.fric.sabsh2 = sabsh2;
    fstat  = fkin + flin1 + flin2; sol.fric.fstat = fstat;
   elseif (frclaw == 3)
    fkin   = tmp(irow, 3+idum); sol.fric.fkin   = fkin;
    frat1  = tmp(irow, 4+idum); sol.fric.frat1  = frat1;
    sabsh1 = tmp(irow, 5+idum); sol.fric.sabsh1 = sabsh1;
    frat2  = tmp(irow, 6+idum); sol.fric.frat1  = frat1;
    sabsh2 = tmp(irow, 7+idum); sol.fric.sabsh2 = sabsh2;
    fstat  = fkin + frat1 + frat2; sol.fric.fstat = fstat;
   elseif (frclaw == 4)
    fkin   = tmp(irow, 3+idum); sol.fric.fkin   = fkin;
    fexp1  = tmp(irow, 4+idum); sol.fric.fexp1  = fexp1;
    sabsh1 = tmp(irow, 5+idum); sol.fric.sabsh1 = sabsh1;
    fexp2  = tmp(irow, 6+idum); sol.fric.fexp2  = fexp2;
    sabsh2 = tmp(irow, 7+idum); sol.fric.sabsh2 = sabsh2;
    fstat  = fkin + fexp1 + fexp2; sol.fric.fstat = fstat;
   elseif (frclaw == 6)
    fref   = tmp(irow, 3+idum); sol.fric.fref   = fref;
    tref   = tmp(irow, 4+idum); sol.fric.tref   = tref;
    dfheat = tmp(irow, 5+idum); sol.fric.dfheat = dfheat;
    dtheat = tmp(irow, 6+idum); sol.fric.dtheat = dtheat;
    fstat  = fref;              sol.fric.fstat  = fstat;
    fkin   = fref+dfheat;       sol.fric.fkin   = fkin;
   end
   if (frclaw >=2 & frclaw <=4)
    sol.fric.memdst = tmp(irow, 8+idum);
    sol.fric.mem_s0 = tmp(irow, 9+idum);
   elseif (frclaw==6)
    sol.fric.memdst = tmp(irow, 7+idum);
    sol.fric.mem_s0 = tmp(irow, 8+idum);
   end

   % handle case of empty contact area

   exter_el = 0;
   adhes_el = 1;
   slip_el  = 2;
   plast_el = 3;

   irow = irow + 1;
   iend = size(tmp,1);
   ncon = iend - irow + 1;

   if (ncon<=0)
    disp('ERROR: mat-file contains no elements in the contact area');
    sol.eldiv  =  ones(my,mx)*exter_el;
    sol.h      = zeros(my,mx);
    sol.mu     = zeros(my,mx);
    sol.pn     = zeros(my,mx);
    sol.px     = zeros(my,mx);
    sol.py     = zeros(my,mx);
    sol.un     = zeros(my,mx);
    sol.ux     = zeros(my,mx);
    sol.uy     = zeros(my,mx);
    sol.srel   = zeros(my,mx);
    sol.shft   = zeros(my,mx);
    sol.trcbnd = zeros(my,mx);
    return
   end

   % get the numbers of the elements that are inside the contact area

   ir = tmp(irow:end,1);

   % extract the data for all active elements from the table

   sol.eldiv =  ones(mx*my,1)*exter_el;
   sol.h     =  ones(mx*my,1)*2*max(tmp(irow:iend,3));
   sol.mu    = zeros(mx*my,1);
   sol.pn    = zeros(mx*my,1);
   sol.px    = zeros(mx*my,1);
   sol.py    = zeros(mx*my,1);
   sol.un    = zeros(mx*my,1);
   sol.ux    = zeros(mx*my,1);
   sol.uy    = zeros(mx*my,1);
   sol.srel  = zeros(mx*my,1);
   if (has_plast)
    sol.taucrt = zeros(mx*my,1);
    sol.uplsx  = zeros(mx*my,1);
    sol.uplsy  = zeros(mx*my,1);
   end
   if (has_heat)
    sol.temp1 = zeros(mx*my,1);
    sol.temp2 = zeros(mx*my,1);
   end

   if (fmtmat<1220)
    idum=0;
   else
    idum=1; % introduced mu in 4th column in version 12.2
   end

   sol.eldiv(ir) = tmp(irow:iend, 2);
   sol.h(ir)     = tmp(irow:iend, 3);
   if (fmtmat>=1220)
    sol.mu(ir)    = tmp(irow:iend, 4);
   end
   sol.pn(ir)    = tmp(irow:iend,3+idum+1);
   sol.px(ir)    = tmp(irow:iend,3+idum+2);
   sol.py(ir)    = tmp(irow:iend,3+idum+3);
   sol.un(ir)    = tmp(irow:iend,3+idum+4);
   sol.ux(ir)    = tmp(irow:iend,3+idum+5);
   sol.uy(ir)    = tmp(irow:iend,3+idum+6);
   sol.srel(ir)  = tmp(irow:iend,3+idum+7);
   iofs = 3+idum+7;
   if (has_plast)
    sol.taucrt(ir) = tmp(irow:iend,iofs+1);
    sol.uplsx(ir)  = tmp(irow:iend,iofs+2);
    sol.uplsy(ir)  = tmp(irow:iend,iofs+3);
    iofs = iofs + 3;
   end
   if (has_heat)
    sol.temp1(ir) = tmp(irow:iend,iofs+1);
    sol.temp2(ir) = tmp(irow:iend,iofs+2);
    iofs = iofs + 2;
   end

   sol.eldiv = reshape(sol.eldiv,mx,my)';
   sol.h     = reshape(sol.h    ,mx,my)';
   sol.mu    = reshape(sol.mu   ,mx,my)';
   sol.pn    = reshape(sol.pn   ,mx,my)';
   sol.px    = reshape(sol.px   ,mx,my)';
   sol.py    = reshape(sol.py   ,mx,my)';
   sol.un    = reshape(sol.un   ,mx,my)';
   sol.ux    = reshape(sol.ux   ,mx,my)';
   sol.uy    = reshape(sol.uy   ,mx,my)';
   sol.srel  = reshape(sol.srel ,mx,my)';
   if (has_plast)
    sol.taucrt = reshape(sol.taucrt,mx,my)';
    sol.uplsx  = reshape(sol.uplsx ,mx,my)';
    sol.uplsy  = reshape(sol.uplsy ,mx,my)';
   end
   if (has_heat)
    sol.temp1  = reshape(sol.temp1 ,mx,my)';
    sol.temp2  = reshape(sol.temp2 ,mx,my)';
   end

   % in rolling, compute the shift from the slip velocity
   % in shifts, fill the shft-array by copying (dq==1)

   if (sol.kincns.t_digit>=2)
    sol.shft = sol.srel * sol.kincns.dq;
   else
    sol.shft = sol.srel;
   end

   % compute the traction bound

   if (~has_plast)
    sol.trcbnd =     sol.mu .* sol.pn;
   else
    sol.trcbnd = min(sol.mu .* sol.pn, sol.taucrt);
   end

  end % fill_struct

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function [ blk1, blk2, blk3, blk4, blk5, blk6, blk7, blk8 ] = loadstrs( expnam, icase, ipatch )
   %
   % [ blk1, blk2, .. ] = loadstrs( expnam, [icase], [ipatch] )
   %
   % Read a single <expnam>.subs-file and return structures with the contents.
   %
   %   expnam = the experiment name used, i.e. the input-filename without
   %            extension .inp.
   %   icase  = the case-number to be read
   %   ipatch = contact patch number, in case of module 1, wheel/rail contact
   %   blk1   = structure with results of the computation for the first
   %               block (grid) of subsurface stress points
   %   blk2   = structure with results of the computation for the second
   %               block (grid) of subsurface stress points
   %   blk3, .. = consecutive blocks
   %

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % initialize optional arguments

   % #loadstrs -2
   if (nargin<2 | isempty(icase))
    icase = 1;
   end
   if (nargin<3 | isempty(ipatch))
    ipatch = 1;
   end

   % locate and read the file

   fname1=[expnam,'.subs'];
   if (icase<=9999)
    fname2 = [expnam, sprintf('.%04d.subs', icase)];
    fname3 = [expnam, sprintf('.%04d%c.subs', icase, 96+ipatch)];
   else
    fname2 = [expnam, sprintf('.%06d.subs', icase)];
    fname3 = [expnam, sprintf('.%06d%c.subs', icase, 96+ipatch)];
   end

   if (exist(fname1,'file'))
    all_data = load(fname1, '-ascii');
   elseif (exist(fname2,'file'))
    all_data = load(fname2, '-ascii');
   elseif (exist(fname3,'file'))
    all_data = load(fname3, '-ascii');
   else
    disp(['ERROR: cannot find file ',fname1,' nor ',fname2,' or ',fname3])
    sol = [];
    return
   end

   % decompose the file contents into the constitutive blocks

   if (size(all_data,2)==14)
    full_tens= 1;
   elseif (size(all_data,2)==4)
    full_tens=-1;
   else
    full_tens= 0;
   end

   nrows = size(all_data,1);
   iblk  = 0;      % shorthand for sol.nblocks

   irow  = 1;
   while (irow <= nrows)
    iblk = iblk + 1;
    sol = [];    % sol == the current block being processed

    % first line of the block: nx, ny, nz
    nx = all_data(irow,1); sol.nx = nx;
    ny = all_data(irow,2); sol.ny = ny;
    nz = all_data(irow,3); sol.nz = nz;
    % disp(sprintf('block %d has %d points (%d x %d x %d)', iblk, nx*ny*nz, nx, ny, nz));

    % the data lines for the block:
    rng = irow + [1:nx*ny*nz];
    blk_data = all_data(rng,:);

    % sort data such that x runs fastest, then y, then z
    [tx, iperm] = sort(blk_data(:,1));
    blk_data = blk_data(iperm,:);
    [ty, iperm] = sort(blk_data(:,2));
    blk_data = blk_data(iperm,:);
    [tz, iperm] = sort(blk_data(:,3));
    blk_data = blk_data(iperm,:);

    % the data-lines are given with iz running fastest, ix slowest.
    sol.npoints = nx * ny * nz;
    sol.x       = reshape(blk_data(:,1), nx, ny, nz);
    sol.y       = reshape(blk_data(:,2), nx, ny, nz);
    sol.z       = reshape(blk_data(:,3), nx, ny, nz);
    if (full_tens==-1)
     sol.sigvm   = reshape(blk_data(:,4), nx, ny, nz);
    else
     sol.ux      = reshape(blk_data(:,4), nx, ny, nz);
     sol.uy      = reshape(blk_data(:,5), nx, ny, nz);
     sol.uz      = reshape(blk_data(:,6), nx, ny, nz);
     sol.sighyd  = reshape(blk_data(:,7), nx, ny, nz);
     sol.sigvm   = reshape(blk_data(:,8), nx, ny, nz);
    end

    if (full_tens==1)
     sol.sigxx = reshape(blk_data(:, 9), nx, ny, nz);
     sol.sigxy = reshape(blk_data(:,10), nx, ny, nz);
     sol.sigxz = reshape(blk_data(:,11), nx, ny, nz);
     sol.sigyy = reshape(blk_data(:,12), nx, ny, nz);
     sol.sigyz = reshape(blk_data(:,13), nx, ny, nz);
     sol.sigzz = reshape(blk_data(:,14), nx, ny, nz);

     if (1==1)
      % 3D calculation of principal stresses
      princ = zeros(nx*ny*nz,3);
      sigma = blk_data(:,9:14);
      for k=1:nx*ny*nz
       princ(k,:) = sort(eig([sigma(k,1), sigma(k,2), sigma(k,3);
        sigma(k,2), sigma(k,4), sigma(k,5);
        sigma(k,3), sigma(k,5), sigma(k,6)]),'descend');
      end
      sol.sigma1 = reshape(princ(:,1), nx, ny, nz);
      sol.sigma2 = reshape(princ(:,2), nx, ny, nz);
      sol.sigma3 = reshape(princ(:,3), nx, ny, nz);
      sol.sigtr = sol.sigma1 - sol.sigma3;
     else
      % 2D principal stress calculation for Oxz: assuming plane strain
      princ = zeros(nx*ny*nz,2);
      sigma = blk_data(:,9:14);
      for k=1:nx*ny*nz
       princ(k,:) = eig([sigma(k,1), sigma(k,3);
        sigma(k,3), sigma(k,6)]);
      end
      sol.sigma1 = reshape(princ(:,1), nx, ny, nz);
      sol.sigma2 = reshape(princ(:,2), nx, ny, nz);
      sol.sigma3 = 0.28*(sol.sigma1+sol.sigma2);
      sol.sigtr = sol.sigma1 - sol.sigma2;
     end
    end

    sol.x       = squeeze(sol.x(:,1,1));
    sol.y       = squeeze(sol.y(1,:,1)); sol.y = sol.y';
    sol.z       = squeeze(sol.z(1,1,:));

    % advance row-number to first line of next block
    irow = irow + 1 + nx*ny*nz;
    % disp(sprintf('next block (%d) starts at row %d', iblk+1, irow));

    % copy the results to the appropriate output argument
    if (iblk==1 | iblk<=nargout)
     eval(sprintf('blk%d = sol;',iblk));
    else
     % ignore struct sol if no more output arguments are available
     % write a warning upon reading/processing of the whole file
     if (irow>=nrows)
      disp(sprintf('WARNING: this .subs-file contains data for %d blocks of points.',iblk));
      disp(sprintf('         Only %d output arguments are provided, the remaining blocks are ignored.',max(1,nargout)));
     end
    end
   end

   if (iblk<nargout)
    disp(sprintf('WARNING: this .subs-file contains data for %d blocks of points.',iblk));
    disp(sprintf('         %d output arguments are provided, the last %d are not filled in.', nargout, nargout-iblk));
    for jblk= iblk+1 : nargout
     eval(sprintf('blk%d = [];',jblk));
    end
   end
  end

  %% #make

  function [ spl2d ] = make_2dspline( ui, vj, xij, yij, zij, mask_j, use_approx, use_insert, use_cylindr, idebug, show_fig, fig_ofs )

   % [ spl2d ] = make_2dspline( ui, vj, xij, yij, zij, [mask_j], [use_approx], [use_insert],
   %                                                            [use_cylindr], [idebug], [show_fig], [fig_ofs] )
   % can be used with ui==xi, xij==[]
   %
   % compute parametric 2D spline {tui, tvj, cij_x, cij_y, cij_z} (B-form) for data {xij,yij,zij},
   %   ui      = parameter (data) positions in longitudinal (u-) direction
   %   vj      = parameter (scaled data) positions in lateral (v-) direction
   %   mask_j  = mask for inactive points (0) esp. for 'sections' with shorter slices
   %   use_approx  = select approximating spline (1, default) or interpolating spline (0) for x-direction
   %   use_insert  = insert knots at not-a-knot positions (1, default) or not (0)
   %   use_cylindr = define spline on u \in [-pi,pi), using wrap-around in evaluation

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % #make_2dspline
   k      = 4;
   nmeasu = length(ui);
   nmeasv = length(vj);

   if (nargin<6 | isempty(mask_j))
    mask_j = ones(nmeasu,nmeasv);
   end
   if (nargin<7 | isempty(use_approx))
    use_approx = 1;
   end
   if (nargin<8 | isempty(use_insert))
    use_insert = 1;   % option for pure tensor splines: do not(0)/do(1) fill in not-a-knot positions
   end
   if (nargin<9 | isempty(use_cylindr))
    use_cylindr = 0;
   end
   if (nargin<10 | isempty(idebug))
    idebug = 0;
   end
   if (nargin<11 | isempty(show_fig))
    show_fig = [];
   end
   if (nargin<12 | isempty(fig_ofs))
    fig_ofs = 0;
   end

   if (~all(mask_j))
    use_insert = 1;   % filling-in not-a-knots is required when a non-trivial mask is given
   end

   if (nmeasu<k | nmeasv<k)
    disp(sprintf('Error: spline needs >= k=%d data positions in u (nmeasu=%d) and v (nmeasv=%d)', ...
     k, nmeasu, nmeasv));
    return;
   end

   if (idebug>=2)
    disp(' ');
    if (use_approx)
     disp(sprintf('make_2dspline: nmeasu = %d, nmeasv = %d, approximating spline', nmeasu, nmeasv));
    else
     disp(sprintf('make_2dspline: nmeasu = %d, nmeasv = %d, interpolating spline', nmeasu, nmeasv));
    end
   end

   has_xij = (~isempty(xij));
   has_yij = (~isempty(yij));
   has_zij = (~isempty(zij));

   if (size(ui,1)==1), ui = ui'; end         % make column vector ui
   if (size(vj,2)==1), vj = vj'; end         % make row vector vj
   if (has_xij & size(xij,1)~=nmeasu), xij = xij'; end % make xij (nmeasu, nmeasv)
   if (has_yij & size(yij,1)~=nmeasu), yij = yij'; end % make yij (nmeasu, nmeasv)
   if (has_zij & size(zij,1)~=nmeasu), zij = zij'; end % make zij (nmeasu, nmeasv)

   if ((has_xij & any(size(xij)~=[nmeasu nmeasv])) | ...
     (has_yij & any(size(yij)~=[nmeasu nmeasv])) | ...
     (has_zij & any(size(zij)~=[nmeasu nmeasv])))

    disp('Arrays ui, vj, xij, yij, zij have incompatible sizes');
    disp( [size(ui); size(vj); size(xij); size(yij); size(zij)] );
    return;
   end

   % developed originally with use_insert=1, n_ins=2; skip knot-insertion when use_insert=0.

   if (use_insert)
    n_ins = 2;
   else
    n_ins = 0;
   end
   n_hlf = n_ins / 2;

   % determine master knot-vector tui_full with knots at all measurement points + extension at start/end

   use_repl = 0;
   tui_full = cntc.make_knot_vector_atmeas( ui, 'u', 1, use_repl, idebug );

   % determine master knot-vector tvj_full with knots at all measurement points + extension at start/end

   use_repl = 0;
   tvj_full = cntc.make_knot_vector_atmeas( vj, 'v', 1, use_repl, idebug );

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Phase 1: build B-splines per slice in lateral direction v -> (y,z)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % consider 'sections' of slices i0:i1 with same mask per slice

   ci_x = zeros(nmeasu, nmeasv+n_ins);
   ci_y = zeros(nmeasu, nmeasv+n_ins);
   ci_z = zeros(nmeasu, nmeasv+n_ins);

   i1 = 0;
   while(i1 < nmeasu)

    % select slices i0:i1 with the same mask per slice

    i0 = i1 + 1;
    i1 = i0;
    while (i1<nmeasu & all(mask_j(i1+1,:)==mask_j(i0,:)))
     i1 = i1 + 1;
    end
    ni  = i1 - i0 + 1;

    % process slices i0:i1 in one go with ni right hand sides

    j0  = find(mask_j(i0,:), 1, 'first');
    j1  = find(mask_j(i0,:), 1, 'last');
    nj  = j1 - j0 + 1;

    if (idebug>=1)
     disp(sprintf('compute  ci_[xyz] for ix=[%3d,%3d], active jv=[%3d,%3d]', i0,i1, j0,j1));
    end

    % select local knot-vector tvj from tvj_full, skipping not-a-knot positions

    jof = [ 1+[-3:0], 3:nj-2, nj+[0:3] ];  % jof = numbering within selection
    jk  = j0 + jof  + k-2;                 % jk  = converted to overall numbering
    jk1 = j0 + 2    + k-2;
    jk2 = j0 + nj-1 + k-2;
    tvj = tvj_full(jk);

    % disp('knots tvj:')
    % disp(tvj)

    % determine collocation matrix for v-direction: evaluate each B-spline at each measurement location

    [~, ~, ~, Bmat] = cntc.eval_bspline_basisfnc( tvj, vj(j0:j1) );

    % determine 1d spline coefficients ci_x, ci_y, ci_z, solving B [ci_x, ci_y, ci_z] = [xij, yij, zij]

    cnd = condest(Bmat);
    if (cnd>1e10)
     disp(sprintf('make_2dspline: Bmat(v) has size %3d x %3d, condition %6.2e', size(Bmat), cnd));
    end

    % ci_y: [ nmeasu, nmeasv ], Bmat: [ nmeasv, nmeasv ], yij: [ nmeasu, nmeasv ]

    % disp('inp_c_y:'); disp(yij(i0:i1,j0:j1))

    if (has_xij), c_x = Bmat \ xij(i0:i1,j0:j1)'; end
    if (has_yij), c_y = Bmat \ yij(i0:i1,j0:j1)'; end
    if (has_zij), c_z = Bmat \ zij(i0:i1,j0:j1)'; end

    % disp('out_c_y:'); disp(c_y)

    if (use_insert)

     % insert knots at not-a-knot positions

     if (has_xij)
      [ t1, c_x ] = cntc.insert_knot( tvj, tvj_full(jk1), c_x, k, idebug );
      [ t2, c_x ] = cntc.insert_knot( t1, tvj_full(jk2), c_x, k, idebug );
     end
     if (has_yij)
      [ t1, c_y ] = cntc.insert_knot( tvj, tvj_full(jk1), c_y, k, idebug );
      [ t2, c_y ] = cntc.insert_knot( t1, tvj_full(jk2), c_y, k, idebug );
     end
     if (has_zij)
      [ t1, c_z ] = cntc.insert_knot( tvj, tvj_full(jk1), c_z, k, idebug );
      [ t2, c_z ] = cntc.insert_knot( t1, tvj_full(jk2), c_z, k, idebug );
     end

    end

    % store results

    if (has_xij), ci_x(i0:i1, j0:j1+n_ins) = c_x'; end
    if (has_yij), ci_y(i0:i1, j0:j1+n_ins) = c_y'; end
    if (has_zij), ci_z(i0:i1, j0:j1+n_ins) = c_z'; end

   end

   if (any(show_fig==4) & has_zij)
    figure(4+fig_ofs); clf;
    surf(ci_z');
    view([90 -90]);
    shading flat;
    colorbar;
    xlabel('u');
    ylabel('v');
   end

   if (any(show_fig==5) & has_zij)
    figure(5+fig_ofs); clf;
    plot(vj, ci_z(5:6,2:end-1), '-o');
    xlabel('v');
    grid on;
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Phase 2: build B-splines for longitudinal direction x -> (cy,cz)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % consider 'regions' of points j0:j1 with same mask per point

   cij_x    = zeros(nmeasu+n_ins, nmeasv+n_ins);
   cij_y    = zeros(nmeasu+n_ins, nmeasv+n_ins);
   cij_z    = zeros(nmeasu+n_ins, nmeasv+n_ins);
   spl_mask = zeros(nmeasu+n_ins, nmeasv+n_ins);

   j1 = 0;
   while(j1 < nmeasv)

    % select 'interpolation paths' j0:j1 with the same mask per point

    j0 = j1 + 1;
    j1 = j0;
    while (j1<nmeasv & all(mask_j(:,j1+1)==mask_j(:,j0)))
     j1 = j1 + 1;
    end

    % determine corresponding columns jt0:jt1 in arrays ci_x, ci_y, ci_z

    if (~use_insert)
     jt0 = j0;
     jt1 = j1;
    else
     jt0 = j0 + n_hlf;
     jt1 = j1 + n_hlf;

     if (j0<=1 | nnz(mask_j(:,j0))>nnz(mask_j(:,j0-1)))
      jt0 = jt0 - 1;     % add guard band
     else
      jt0 = jt0 + 1;     % move 1st column to neighbour
     end
     if (j1>=nmeasv | nnz(mask_j(:,j1))>nnz(mask_j(:,j1+1)))
      jt1 = jt1 + 1;     % add guard band
     else
      jt1 = jt1 - 1;     % move last column to neighbour
     end
    end
    njt = jt1 - jt0 + 1;

    % determine active i0:i1 and corresponding rows it0:it1 in cij_x, cij_y, cij_z

    i0   = find(mask_j(:,j0), 1, 'first');
    i1   = find(mask_j(:,j0), 1, 'last');
    ni   = i1 - i0 + 1;
    it0  = i0;
    it1  = i1 + 2;

    if (idebug>=1)
     disp(sprintf('compute cij_[xyz] for jt=[%3d,%3d], active ix=[%3d,%3d]', jt0,jt1, i0,i1));
    end

    if (use_approx)

     % approximating spline -- interpret (cx, cy, cz) as control points for 2D spline

     c_x = zeros(ni+n_ins,njt);
     c_y = zeros(ni+n_ins,njt);
     c_z = zeros(ni+n_ins,njt);

     if (~use_insert)

      if (has_xij), c_x(1:ni,:) = ci_x(i0:i1,jt0:jt1); end
      if (has_yij), c_y(1:ni,:) = ci_y(i0:i1,jt0:jt1); end
      if (has_zij), c_z(1:ni,:) = ci_z(i0:i1,jt0:jt1); end

     else

      % select interior 'control points' from ci_x, ci_y, ci_z

      if (has_xij), c_x(2:ni+1,:) = ci_x(i0:i1,jt0:jt1); end
      if (has_yij), c_y(2:ni+1,:) = ci_y(i0:i1,jt0:jt1); end
      if (has_zij), c_z(2:ni+1,:) = ci_z(i0:i1,jt0:jt1); end

      % add fantom control points at rows 1 and ni+2

      dt_ext = tui_full(3+i0)   - tui_full(3+i0-1);
      dt_int = tui_full(3+i0+1) - tui_full(3+i0);
      c_x(1,:)  = c_x(2,:) - dt_ext/dt_int * (c_x(3,:) - c_x(2,:));
      c_y(1,:)  = c_y(2,:) - dt_ext/dt_int * (c_y(3,:) - c_y(2,:));
      c_z(1,:)  = c_z(2,:) - dt_ext/dt_int * (c_z(3,:) - c_z(2,:));
      if (idebug>=1)
       disp(sprintf('start: t_ext = %5.2f, t_bnd = %5.2f, t_int = %5.2f', tui_full(3+i0+[-1:1])));
       disp(sprintf('       c_ext = %5.2f *c_bnd + %5.2f *c_int', [1,0]+[1,-1]*dt_ext/dt_int))
      end

      dt_int = tui_full(3+i1)   - tui_full(3+i1-1);
      dt_ext = tui_full(3+i1+1) - tui_full(3+i1);
      c_x(ni+2,:)  = c_x(ni+1,:) + dt_ext/dt_int * (c_x(ni+1,:) - c_x(ni,:));
      c_y(ni+2,:)  = c_y(ni+1,:) + dt_ext/dt_int * (c_y(ni+1,:) - c_y(ni,:));
      c_z(ni+2,:)  = c_z(ni+1,:) + dt_ext/dt_int * (c_z(ni+1,:) - c_z(ni,:));
      if (idebug>=1)
       disp(sprintf('end:   t_int = %5.2f, t_bnd = %5.2f, t_ext = %5.2f', tui_full(3+i1+[-1:1])));
       disp(sprintf('       c_ext = %5.2f *c_bnd + %5.2f *c_int', [1,0]+[1,-1]*dt_ext/dt_int))
      end
     end

    else

     % interpolating spline -- take (cx, cy, cz) as data points for 2D spline

     % select local knot-vector tui from tui_full

     iof = [ 1+[-3:0], 3:ni-2, ni+[0:3] ]; % iof = numbering within selection
     ik  = i0 + iof  + k-2;               % ik  = converted to overall numbering
     ik1 = i0 + 2    + k-2;
     ik2 = i0 + ni-1 + k-2;
     tui = tui_full(ik);
     % disp('knots tui_full:')
     % disp(tui_full')
     % disp('knots tui:')
     % disp(tui')

     % determine collocation matrix for u-direction: evaluate each B-spline at each measurement location

     [~, ~, ~, Bmat] = cntc.eval_bspline_basisfnc( tui, ui(i0:i1) );

     % determine 2d spline coefficients cij_[xyz], solving B [cij_x, cij_y, cij_z] = [ci_x, ci_y, ci_z]

     cnd = condest(Bmat);
     if (cnd>1e10)
      disp(sprintf('make_2dspline: Bmat(u) has size %3d x %3d, condition %6.2e', size(Bmat), cnd));
     end

     % cij_y: [ nmeasu, nmeasv+2 ], Bmat: [ nmeasu, nmeasu ], ci_y: [ nmeasu, nmeasv+2 ]

     if (has_xij), c_x = Bmat \ ci_x(i0:i1,jt0:jt1); end
     if (has_yij), c_y = Bmat \ ci_y(i0:i1,jt0:jt1); end
     if (has_zij), c_z = Bmat \ ci_z(i0:i1,jt0:jt1); end

     if (use_insert)

      % insert knots at not-a-knot positions

      if (has_xij)
       [ t1, c_x ] = cntc.insert_knot( tui, tui_full(ik1), c_x, k, idebug );
       [ t2, c_x ] = cntc.insert_knot( t1, tui_full(ik2), c_x, k, idebug );
      end
      if (has_yij)
       [ t1, c_y ] = cntc.insert_knot( tui, tui_full(ik1), c_y, k, idebug );
       [ t2, c_y ] = cntc.insert_knot( t1, tui_full(ik2), c_y, k, idebug );
      end
      if (has_zij)
       [ t1, c_z ] = cntc.insert_knot( tui, tui_full(ik1), c_z, k, idebug );
       [ t2, c_z ] = cntc.insert_knot( t1, tui_full(ik2), c_z, k, idebug );
      end
     end

    end % use_approx

    % store results

    if (has_xij), cij_x(i0:i1+n_ins, jt0:jt1) = c_x; end
    if (has_yij), cij_y(i0:i1+n_ins, jt0:jt1) = c_y; end
    if (has_zij), cij_z(i0:i1+n_ins, jt0:jt1) = c_z; end
    spl_mask(i0:i1+n_ins, jt0:jt1) = 1;
   end

   spl2d = struct('ui',ui, 'vj',vj, 'use_cylindr',use_cylindr, 'xij',xij, 'yij',yij, 'zij',zij, ...
    'mask_j',spl_mask, 'tui',tui_full, 'tvj',tvj_full);

   if (~use_insert)
    spl2d.tui = spl2d.tui([1:k, k+2:end-k-1, end-k+1:end]);
    spl2d.tvj = spl2d.tvj([1:k, k+2:end-k-1, end-k+1:end]);
   end

   spl2d.tui_grev = conv(spl2d.tui(2:end-1), [1 1 1]/3, 'valid');
   spl2d.tvj_grev = conv(spl2d.tvj(2:end-1), [1 1 1]/3, 'valid');

   if (has_xij), spl2d.cij_x = cij_x; end
   if (has_yij), spl2d.cij_y = cij_y; end
   if (has_zij), spl2d.cij_z = cij_z; end

   % temporarily, copy intermediate coefficients to spl2d

   if (has_xij), spl2d.ci_x = ci_x; end
   if (has_yij), spl2d.ci_y = ci_y; end
   if (has_zij), spl2d.ci_z = ci_z; end

   if (idebug>=1)

    % check correctness of interpolating spline at meas. positions

    [ ~, x_out, y_out, z_out ] = cntc.eval_2dspline( spl2d, ui, vj, [], idebug );

    disp(sprintf('Max diff(spline y - yij) = %3.1e, max(spline z - zij) = %3.1e', ...
     max(max(abs(y_out-yij))), max(max(abs(z_out-zij))) ))
   end

  end % make_2dspline

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ tj ] = make_knot_vector_atmeas( si, namcoor, use_insert, use_repl, idebug )
   % [ tj ] = make_knot_vector_atmeas( si, namcoor, use_insert, use_repl, idebug )
   % use_insert = 1 == knots at all measurements, 0 == skip 2nd and n-1th for not-a-knot b.c.
   % use_repl   = 1 == repeat 1st/last points, 0 == repeat 1st/last intervals
   % #make_knot_vector_atmeas

   nmeas  = length(si);
   nspl   = nmeas;      % number of B-splines equal to number of measured values
   if (use_insert)
    nknot = nspl + 6; % cubic spline, order k=4: add 3 knots at both ends
   else
    nknot = nspl + 4; % not-a-knot b.c.: remove 2nd and n-1th meas. points
   end

   if (size(si,1)>1)
    tj     = zeros(nknot,1);
   else
    tj     = zeros(1,nknot);
   end

   dt0    = si(2) - si(1);       % compute 1st and last intervals
   dt1    = si(end) - si(end-1);

   if (use_insert & use_repl)

    tj(1:4)           = si(1);                % replicate s(1) and s(end) 4 times for order k=4,
    tj(5:nknot-4)     = si(2:nmeas-1);        % use all measurement points
    tj(nknot-3:nknot) = si(nmeas);

   elseif (use_insert)

    tj(1:4)           = si(1)    +[-3:0]*dt0; % replicate 1st and last intervals 3 times for order k=4,
    tj(5:nknot-4)     = si(2:nmeas-1);        % use all measurement points
    tj(nknot-3:nknot) = si(nmeas)+[ 0:3]*dt1;

   elseif (use_repl)

    tj(1:4)           = si(1);                % replicate s(1) and s(end) 4 times for order k=4,
    tj(5:nknot-4)     = si(3:nmeas-2);        % skip s(2) and s(end-1) for not-a-knot boundaries,
    tj(nknot-3:nknot) = si(nmeas);            %    retain internal s(i) for interpolating spline

   else

    tj(1:4)           = si(1)    +[-3:0]*dt0; % replicate first & last intervals 3 times for order k=4
    tj(5:nknot-4)     = si(3:nmeas-2);        % skip s(2) and s(end-1) for not-a-knot boundaries
    tj(nknot-3:nknot) = si(nmeas)+[ 0:3]*dt1;

   end

   if (idebug>=3)
    disp(sprintf(' knot vector t%c has %d knots, nmeas= %d:', namcoor, nknot, nmeas));
    disp(tj')
    %  disp(tj(end-5:end)')
   end

  end % make_knot_vector_atmeas

  function [ s ] = make_arclength( y, z )

   % [ s ] = make_arclength( y, z )
   %
   % compute arc-length parameterization { s_i } of curve { (y_i, z_i) }:
   % cumulative sum of distances ds between successive points

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % #make_arclength
   dy = diff(y);
   dz = diff(z);
   ds = sqrt(dy.^2 + dz.^2);
   s  = [0; cumsum(ds) ];

   % interpolate [y, s] to y=0, shift s := s - s0

   is_monotonic = (all(diff(y)>0) | all(diff(y)<0));
   ix_zeros     = find(abs(y)<1e-20);
   if (isempty(ix_zeros))
    ix_zeros  = find(y(1:end-1).*y(2:end)<0);
   end

   if (is_monotonic & min(y)<0 & max(y)>0)
    sref = interp1(y, s, 0);
   elseif (isscalar(ix_zeros))
    if (ix_zeros<length(y))
     ix   = ix_zeros + [ 0,1];
    else
     ix   = ix_zeros + [-1,0];
    end
    sref = interp1(y(ix), s(ix), 0);
   else
    sref = 0;
   end
   s  = s - sref;

  end % make_arclength


  function [ spl ] = make_ppform( spl, D1, D2, D3 )

   % [ spl ] = make_ppform( spl, [D1], [D2], [D3] )
   %
   % determine PP-form: value of spline & derivatives at start of each segment

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % compute D1, D2, D3 if not provided

   % #make_ppform
   if (nargin<4 | isempty(D1) | isempty(D2) | isempty(D3))

    % average step sizes over 1/2/3 adjacent intervals

    nknot  = length(spl.tj);
    dtj_k1 =  spl.tj(2:end) - spl.tj(1:nknot-1);
    dtj_k2 = (spl.tj(3:end) - spl.tj(1:nknot-2)) / 2;
    dtj_k3 = (spl.tj(4:end) - spl.tj(1:nknot-3)) / 3;

    % inverse average step sizes, zero at zero step size

    dtj_inv1 = zeros(size(dtj_k1)); ix  = find(dtj_k1>0); dtj_inv1(ix) = 1./ dtj_k1(ix);
    dtj_inv2 = zeros(size(dtj_k2)); ix  = find(dtj_k2>0); dtj_inv2(ix) = 1./ dtj_k2(ix);
    dtj_inv3 = zeros(size(dtj_k3)); ix  = find(dtj_k3>0); dtj_inv3(ix) = 1./ dtj_k3(ix);

    % backward difference matrices using average step lengths 1/2/3

    D1   = spdiags(dtj_inv1, 0, nknot-1, nknot-2) - spdiags(dtj_inv1(2:end), -1, nknot-1, nknot-2);
    D2   = spdiags(dtj_inv2, 0, nknot-2, nknot-3) - spdiags(dtj_inv2(2:end), -1, nknot-2, nknot-3);
    D3   = spdiags(dtj_inv3, 0, nknot-3, nknot-4) - spdiags(dtj_inv3(2:end), -1, nknot-3, nknot-4);
   end

   % determine value of each basis-function at the start of each segment ksi_i

   s = max(spl.s, spl.tj(4));      % avoid warning for extended knots
   s = min(    s, spl.tj(end-3));  %  --> constant extrapolation outside basic interval
   [B1, B2, B3, B4] = cntc.eval_bspline_basisfnc( spl.tj, s );

   % determine PP-form: value of spline & derivatives at start of each segment

   spl.ay0 = B4 * spl.cy;
   spl.ay1 = B3 * D3 * spl.cy;
   spl.ay2 = B2 * D2 * D3 * spl.cy / 2;
   spl.ay3 = B1 * D1 * D2 * D3 * spl.cy / 6;
   spl.az0 = B4 * spl.cz;
   spl.az1 = B3 * D3 * spl.cz;
   spl.az2 = B2 * D2 * D3 * spl.cz / 2;
   spl.az3 = B1 * D1 * D2 * D3 * spl.cz / 6;

  end % make_ppform


  function [ x0, y0, z0, th_rfl ] = make_reflec( spl2d, vec_a, vec_v, ui, vj )

   % [ x0, y0, z0, th_rfl ] = make_reflec( spl2d, vec_a, vec_v, ui, vj  )
   %
   % compute reflection angles theta [deg] for light source direction vec_a and view direction vec_v
   % sampling the spline at u-parameters ui and v-parameters vj
   %
   % vec_a, vec_v can be [x,y,z] or [az,el], with az,el in degrees

   % #make_reflec
   rad = pi / 180;

   if (length(vec_a)==2)
    a_az = vec_a(1);
    a_el = vec_a(2);
    vec_a = [ cos(a_az*rad)*cos(a_el*rad), sin(a_az*rad)*cos(a_el*rad), sin(a_el*rad) ]';
   end
   if (length(vec_v)==2)
    v_az = vec_v(1);
    v_el = vec_v(2);
    vec_v = [ cos(v_az*rad)*cos(v_el*rad), sin(v_az*rad)*cos(v_el*rad), sin(v_el*rad) ]';
   end
   if (size(vec_a,2)>1)
    vec_a = vec_a';
   end
   if (size(vec_v,2)>1)
    vec_v = vec_v';
   end
   vec_a = vec_a / norm(vec_a);
   vec_v = vec_v / norm(vec_v);

   % evaluate spline surface at positions ui x vj

   idebug = 0;
   [~, x0, y0, z0] = cntc.eval_2dspline(spl2d, ui, vj, [], idebug);

   % construct first normalized tangent vector at each ui x vj

   [dx1, dy1, dz1] = cntc.eval_2dspline_deriv(spl2d, ui, vj, 1);
   nrm = sqrt(dx1.^2 + dy1.^2 + dz1.^2);
   t1(:,:,1) = dx1./nrm; t1(:,:,2) = dy1./nrm; t1(:,:,3) = dz1./nrm;

   % construct second normalized tangent vector at each ui x vj

   [dx2, dy2, dz2] = cntc.eval_2dspline_deriv(spl2d, ui, vj, 2);
   nrm = sqrt(dx2.^2 + dy2.^2 + dz2.^2);
   t2(:,:,1) = dx2./nrm; t2(:,:,2) = dy2./nrm; t2(:,:,3) = dz2./nrm;

   % construct normalized normal vector at each ui x vj

   nvec = cross(t1, t2);
   nrm = sqrt(sum(nvec.^2, 3));
   nvec(:,:,1) = nvec(:,:,1) ./ nrm;
   nvec(:,:,2) = nvec(:,:,2) ./ nrm;
   nvec(:,:,3) = nvec(:,:,3) ./ nrm;

   % compute projection of v onto plane P

   v_proj = vec_v - (vec_a'*vec_v)*vec_a; v_proj = v_proj / norm(v_proj);
   v_orth = cross(vec_a, v_proj);

   % compute theta at surface evaluation positions

   nu = length(ui);
   nv = length(vj);
   th_rfl = zeros(nu, nv);

   if (0==1)

    for iu = 1 : nu
     for jv = 1 : nv
      n_i   = squeeze(nvec(iu,jv,:));
      v_t   = vec_v - (n_i'*vec_v)*n_i;
      vec_r = 2 * v_t - vec_v;

      r_proj = vec_r - (vec_a'*vec_r)*vec_a; r_proj = r_proj / norm(r_proj);
      r_p   = (-r_proj'*v_proj);
      r_o   = (-r_proj'*v_orth);
      th_rfl(iu,jv) = atan2(r_o, r_p) / rad;
     end
    end

   else

    one  = ones(nu,nv);

    % n_i   = squeeze(nvec(iu,jv,:));
    % v_t   = vec_v - (n_i'*vec_v) * n_i;
    % vec_r = 2 * v_t - vec_v;

    dot_nv = zeros(nu,nv);
    for k = 1 : 3
     dot_nv = dot_nv + vec_v(k) * nvec(:,:,k);
    end

    vec_r  = zeros(nu,nv,3);
    dot_ar = zeros(nu,nv);
    for k = 1 : 3
     vec_r(:,:,k) = vec_v(k) - 2 * dot_nv .* nvec(:,:,k);
     dot_ar = dot_ar + vec_a(k) * vec_r(:,:,k);
    end

    % r_proj = vec_r - (vec_a'*vec_r)*vec_a;
    % r_p   = (-r_proj'*v_proj);
    % r_o   = (-r_proj'*v_orth);

    r_proj = zeros(nu,nv,3);
    for k = 1 : 3
     r_proj(:,:,k) = vec_r(:,:,k) - dot_ar * vec_a(k);
    end

    dot_rvp = zeros(nu,nv);
    dot_rvo = zeros(nu,nv);
    for k = 1 : 3
     dot_rvp = dot_rvp - r_proj(:,:,k) * v_proj(k);
     dot_rvo = dot_rvp - r_proj(:,:,k) * v_orth(k);
    end

    th_rfl = atan2(dot_rvo, dot_rvp) / rad;

   end

  end % make_reflec

  function [ spl ] = make_spline( s, y, z, lambda, wgt, ikinks, iaccel, use_bspline, ds_bspl, use_deriv, use_repl, idebug )

   % [ spl ] = make_spline( [s], y, [z], lambda, [wgt], [ikinks], [iaccel], [use_bspline], ...
   %                                                           [ds_bspl], [use_deriv], [use_repl], [idebug] )
   %
   % compute parametric smoothing spline {s, a0y-a3y, a0z-a3z} (pp-form) for data {y(s),z(s)},
   %   s         = arc-length, will be computed if not provided
   %   lambda    = weight for smoothness, default 0 (no smoothing)
   %   wgt       = weight for data-points, default all ones.
   %   ikinks    = discontinuous 1st derivatives: sharp corners/kinks at list of indices [ i1, i2, ... ]
   %   iaccel    = discontinuous 2nd derivatives: radius jumps, acceleration points [i1, i2, ...]
   %   ds_bspl   = target step size for B-spline method
   %   use_deriv = 2 or 3 for penalty on 2nd or 3rd derivative
   %   use_repl  = 1=replicate knots at s(1) and s(end) or 0=extend knots, replicating intervals ds(1), ds(end)

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % #make_spline
   if (nargin>=3 & isempty(y) & ~isempty(z))
    y      = ones(size(z));
   end
   if (nargin<3 | isempty(z))
    z      = ones(size(y));
   end
   if (nargin<4 | isempty(lambda))
    lambda = 0;  % default: no smoothing
   end
   if (nargin<5 | isempty(wgt))
    wgt = ones(size(y));
   end
   if (nargin<6)
    ikinks = []; % default: no kinks/sharp corners
   end
   if (nargin<7)
    iaccel = []; % default: no acceleration points
   end
   if (nargin<8 | isempty(use_bspline))
    use_bspline = 1;  % default: new B-spline method
   end
   if (nargin<9 | isempty(ds_bspl))
    ds_bspl = 2;
   end
   if (nargin<10 | isempty(use_deriv))
    use_deriv = 3;
    if (~use_bspline), use_deriv = 2; end
   end
   if (nargin<11 | isempty(use_repl))
    use_repl = 1;
   end
   if (nargin<12 | isempty(idebug))
    idebug = 0;
   end

   if (~use_bspline & ~isempty(iaccel))
    disp('Accelerations are ignored in the PP-form method.');
   end
   if (~use_bspline & use_deriv~=2)
    disp('The PP-spline uses use_deriv=2.');
    use_deriv = 2;
   end

   npnt = length(y);
   if (size(s,1)==1), s = s'; end % make column vector
   if (size(y,1)==1), y = y'; end
   if (size(z,1)==1), z = z'; end
   if (size(wgt,1)==1), wgt = wgt'; end
   if (size(ikinks,1) > size(ikinks,2)) % make row vector
    ikinks = ikinks';
   end
   if (size(iaccel,1) > size(iaccel,2)) % make row vector
    iaccel = iaccel';
   end

   ix = find(iaccel<=0 | iaccel>npnt);
   if (~isempty(ix))
    disp(sprintf('ERROR: %d accelerations out of range [1,%d]', length(ix), npnt));
    disp(iaccel')
    return;
   end

   if (isempty(s))
    s = cntc.make_arclength( y, z );
   end
   if (any(size(s)~=size(y)) | any(size(z)~=size(y)) | any(size(wgt)~=size(y)))
    disp('Arrays s, y, z, wgt have incompatible sizes');
    disp( [size(s); size(y); size(z); size(wgt)] );
    return;
   end

   if (use_bspline)
    spl = cntc.make_bspline( s, y, z, lambda, wgt, ikinks, iaccel, ds_bspl, use_deriv, use_repl, idebug );
   else
    spl = cntc.make_pp_spline( s, y, z, lambda, wgt, ikinks, idebug );
   end

  end % make_spline

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ tjx, a_out ] = insert_knot( tj, tnew, a_in, k, idebug )
   % [ tjx, a_out ] = insert_knot( tj, tnew, a_in, k, [idebug] )
   %
   % insert one new knot tnew into knot vector tj and compute B-spline coefficients a_out

   % #insert_knot -2

   if (nargin<5 | isempty(idebug))
    idebug = 0;
   end

   if (size(tj,1)>1), tj = tj'; end % using row vector tj

   nknot = length(tj);
   nspl  = size(a_in, 1);
   ncol  = size(a_in, 2);
   if (nknot ~= nspl+k)
    disp(sprintf('ERROR: size of tj (%d) does not match with a_in (%d, %d)', nknot, nspl, k));
    return;
   end

   % determine position jnew for tnew, form extended knot vector

   jnew = find( tj > tnew, 1, 'first');

   tjx  = [tj(1:jnew-1) , tnew, tj(jnew:end) ];

   a_out            = zeros(nspl+1,ncol);
   for j = 1 : jnew - k
    a_out(j,:) = a_in(j,:);
   end
   for j = jnew-k+1 : jnew-1
    fj         = (tjx(jnew) - tj(j)) / (tj(j+3) - tj(j));
    a_out(j,:) = (1 - fj) * a_in(j-1,:) + fj * a_in(j,:);
   end
   for j = jnew : nspl+1
    a_out(j,:) = a_in(j-1,:);
   end

  end % insert_knot

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ spl ] = make_bspline( s_in, y_in, z_in, lambda, wgt, ikinks, iaccel, ds_bspl, use_deriv, use_repl, idebug )

   % [ spl ] = make_bspline( s_in, y_in, z_in, lambda, wgt, ikinks, iaccel, ds_bspl, ...
   %                                                               use_deriv, use_repl, [idebug] )
   %
   % compute parametric cubic smoothing spline (y(s),z(s)) using least squares B-spline method with
   % penalty on 2nd or 3rd derivative
   % #make_bspline

   if (nargin<11 | isempty(idebug))
    idebug = 1;
   end
   k  = 4;  % spline order 4, cubic splines

   % check size of ds_bspl: too small value spoils the condition of matrix M

   ds_min = 0.1 * lambda^(1/(2*use_deriv));
   if (ds_bspl < ds_min)
    if (idebug>=1)
     disp(sprintf('Input ds_bspl = %6.1e too small, using ds_min = %6.3f', ds_bspl, ds_min));
    end
    ds_bspl = ds_min;
   end

   % determine appropriate knot-vector tj
   %  - simple algorithm: just apply ds_bspl as much as possible, needs some smoothing
   %  - advanced algorithm: no more than 1 knot per meas.segment, no smaller than ds_bspl, fewer at kinks, accel and boundaries

   use_simple  = 1;      % default
   lmb_min     = 1e-6;
   if (lambda<lmb_min), use_simple = 0; end

   [tj, keep_reduc] = cntc.make_knot_vector( s_in, ikinks, iaccel, ds_bspl, use_simple, use_repl, idebug );
   nknot = length(tj);

   % determine collocation matrix: evaluate each B-spline at each measurement location

   nmeas = length(s_in);
   [~, ~, ~, Bmat] = cntc.eval_bspline_basisfnc( tj, s_in );
   if (idebug>=5 & nmeas<20), disp(full(Bmat)); end

   % weight matrix W for input data

   Wmat = spdiags(wgt, 0, nmeas, nmeas);
   % disp(Wmat)

   % average step sizes over 1/2/3 adjacent intervals

   dtj_k1 =  tj(2:end) - tj(1:nknot-1);
   dtj_k2 = (tj(3:end) - tj(1:nknot-2)) / 2;
   dtj_k3 = (tj(4:end) - tj(1:nknot-3)) / 3;

   % inverse average step sizes, zero at zero step size

   dtj_inv1 = zeros(size(dtj_k1)); ix  = find(dtj_k1>0); dtj_inv1(ix) = 1./ dtj_k1(ix);
   dtj_inv2 = zeros(size(dtj_k2)); ix  = find(dtj_k2>0); dtj_inv2(ix) = 1./ dtj_k2(ix);
   dtj_inv3 = zeros(size(dtj_k3)); ix  = find(dtj_k3>0); dtj_inv3(ix) = 1./ dtj_k3(ix);

   % backward difference matrices using average step lengths 1/2/3

   D1   = spdiags(dtj_inv1, 0, nknot-1, nknot-2) - spdiags(dtj_inv1(2:end), -1, nknot-1, nknot-2);
   D2   = spdiags(dtj_inv2, 0, nknot-2, nknot-3) - spdiags(dtj_inv2(2:end), -1, nknot-2, nknot-3);
   D3   = spdiags(dtj_inv3, 0, nknot-3, nknot-4) - spdiags(dtj_inv3(2:end), -1, nknot-3, nknot-4);

   % penalty matrix D^T * C * D for penalty on 2nd or 3rd derivative

   if (use_deriv==2)
    Dmat = D2(3:nknot-k,2:nknot-k) * D3(2:nknot-k,1:nknot-k);
    Cmat =   2/3 * spdiags( dtj_k2(1:nknot-k),  0, nknot-k, nknot-k ) + ...
     + 1/6 * spdiags( dtj_k1(2:nknot-k), -1, nknot-k, nknot-k ) + ...
     + 1/6 * spdiags( dtj_k1(1:nknot-k),  1, nknot-k, nknot-k );
    Cmat = Cmat(3:nknot-k,3:nknot-k);
   else
    Dmat = D1(4:nknot-k,3:nknot-k) * D2(3:nknot-k,2:nknot-k) * D3(2:nknot-k,1:nknot-k);
    Cmat = spdiags( dtj_k1, 0, nknot-1, nknot-1 );
    Cmat = Cmat(4:nknot-k,4:nknot-k);
   end

   % Pmat = Dmat' * Cmat * Dmat;
   % dt = dtj_k1(1);
   % disp(full(dt^5*Pmat(1:5,1:6)))
   % disp(flipud(fliplr( full(dt^5*Pmat(end-4:end,end-5:end))) ))

   if (idebug>=3)
    for i = 1:min(10,nknot-1)
     tmp = full( Dmat(i,max(1,i-3):min(i,nknot-4) ));
     disp(sprintf('i=%3d: d(i,:) = %10.3e %10.3e %10.3e %10.3e', i, tmp))
    end
   end

   % combine basisfunctions at segments that use a reduced order

   if (any(~keep_reduc))
    if (idebug>=1)
     disp(sprintf('Total %d knots, %d knots removed, %d after reduction', nknot, nnz(~keep_reduc), ...
      nnz(keep_reduc)))
    end

    Bmat = cntc.reduce_matrix(Bmat, tj, keep_reduc, idebug);
    Dmat = cntc.reduce_matrix(Dmat, tj, keep_reduc);
   end

   % determine spline coefficients, solving the least squares optimization equation

   if (idebug>=2)
    if (0==1)
     disp('Bt*W*B:');
     tmp = Bmat'*Wmat*Bmat;
    else
     disp('B:');
     tmp = Bmat;
    end
    if (size(tmp,1)>20)
     disp(full(tmp(1:10, 1:10)))
     disp(full(tmp(end-9:end, end-9:end)))
    else
     disp(full(tmp))
    end
    % figure(6); spy(tmp);

    if (lambda>1e-12)
     for i = 1: min(10,nknot-1)
      tmp = full( Dmat(i,max(1,i-5):min(i-2,nknot-4) ));
      disp(sprintf('i=%3d: dr(i,:) = %10.3e %10.3e %10.3e %10.3e', i, tmp));
     end

     tmp = Dmat'*Cmat*Dmat;
     if (size(tmp,1)>20)
      disp('Dt*C*D (top-left):');
      disp(sprintf('%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n', ...
       full(tmp(1:10, 1:10))))
      disp('Dt*C*D (bottom-right):');
      disp(sprintf('%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n', ...
       full(tmp(end-9:end, end-9:end))))
     else
      disp('Dt*C*D:');
      disp(full(tmp))
     end

     % disp(log10(diag(tmp)))

     disp('[Bt*W*B + lmb * Dt*C*D, rhs]:')
     tmp = [Bmat'*Wmat*Bmat + lambda * Dmat'*Cmat*Dmat, Bmat'*Wmat*[y_in,z_in]];
     if (size(tmp,1)>10)
      disp(full(tmp(1:10, [1:8,end-1:end])))
     else
      disp(full(tmp));
     end
    end

   end

   Mat = (Bmat'*Wmat*Bmat + lambda * Dmat'*Cmat*Dmat);
   cnd = condest(Mat);
   if (cnd>1e10)
    disp(sprintf('make_bspline: matrix condition %6.2e, spline may be unstable', cnd));
   end

   coef_y = Mat \ (Bmat' * Wmat * y_in);
   coef_z = Mat \ (Bmat' * Wmat * z_in);

   if (idebug>=3)
    disp('B-spline coefficients:')
    if (length(coef_y)>20)
     disp([ [1:10,length(coef_y)+[-9:0]]', coef_y([1:10,end-9:end]), coef_z([1:10,end-9:end])])
    else
     disp([ [1:length(coef_y)]', coef_y, coef_z])
    end
   end

   coef_y = cntc.expand_vector( coef_y, tj, keep_reduc );
   coef_z = cntc.expand_vector( coef_z, tj, keep_reduc );

   if (idebug>=3 & ~all(keep_reduc))
    disp('Expanded coefficients:')
    disp([ [1:10,length(coef_y)+[-9:0]]', coef_y([1:10,end-9:end]), coef_z([1:10,end-9:end])])
   end

   ksi = unique(tj);
   spl = struct('s',ksi, 'tj',tj, 'cy',coef_y, 'cz',coef_z);
   spl = cntc.make_ppform( spl, D1, D2, D3 );

  end % make_bspline

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ is_ok ] = check_kink_accel( si, ikinks, iaccel )

   % [ is_ok ] = check_kink_accel( si, ikinks, iaccel )
   %
   % check requirements of B-spline method on kinks and accelerations
   % #check_kink_accel -2

   is_ok  = 0;

   nmeas  = length(si);
   nkink  = length(ikinks);
   naccel = length(iaccel);
   ikinks = sort(ikinks);
   iaccel = sort(iaccel);
   if (size(ikinks,1)>size(ikinks,2)), ikinks = ikinks'; end
   if (size(iaccel,1)>size(iaccel,2)), iaccel = iaccel'; end

   % perform checks on ikinks and iaccel:
   %  - all kinks and accelerations must be different

   if (length(unique([ikinks iaccel])) < length([ikinks iaccel]))
    disp('ERROR: check_kink_accel: ikinks and iaccel must all be different.')
    disp(ikinks);
    disp(iaccel);
    is_ok = -1;
    return;
   end

   %  - all kinks must lie in interior of measurement range

   if (any(ikinks<=1 | ikinks>=nmeas))
    disp(sprintf('ERROR: check_kink_accel: ikinks must be >= 2 and <= nmeas = %d', nmeas))
    disp(ikinks);
    is_ok = -1;
    return;
   end

   %  - all accelerations must lie in interior of measurement range, away from boundaries

   if (any(iaccel<=3 | iaccel>=nmeas-2))
    disp(sprintf('ERROR: check_kink_accel: iaccel must be >= 4 and <= nmeas-3 = %d', nmeas-3))
    disp(iaccel);
    is_ok = -1;
    return;
   end

   %  - all kinks must have 0 or >=2 points between them, accelerations and kinks need >= 2 points between them

   iall  = [  1, nmeas,  ikinks,          iaccel           ];
   itype = [  0,    0 ,  1*ones(1,nkink), 2*ones(1,naccel) ];
   [iall, iperm] = sort(iall);
   itype = itype(:,iperm);
   idist = diff(iall);

   if (any(idist<=2 & (itype(1:end-1)==2 | itype(2:end)==2) ))
    disp('ERROR: check_kink_accel: accelerations must have at least 2 points between them and');
    disp('                         kinks, boundaries and other accelerations');
    disp([iall; itype])
    is_ok = -1;
    return
   end
   if (any(idist==2 & itype(1:end-1)<=1 & itype(2:end)==1))
    disp(sprintf('ERROR: check_kink_accel: kinks must be adjacent or have >=2 points between them'));
    disp(ikinks);
    is_ok = -1;
    return;
   end

   is_ok = 1;

  end % check_kink_accel

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ tj, keep_reduc ] = make_knot_vector( si, ikinks, iaccel, ds_bspl, use_simple,use_repl, idebug )

   % [ tj, keep_reduc ] = make_knot_vector( s_meas, ikinks, iaccel, ds_bspl, use_simple, ...
   %                                                                               use_repl, [idebug] )
   %
   % determine knot-vector with repeated knots at kinks and accelerations
   % reduce knots in intervals between double kinks or kinks adjacent to boundaries
   %  - ds_bspl       = target step size
   %  - use_simple    = true:  apply ds_bspl as much as possible
   %                    false: at most 1 knot per data point
   %  - use_repl      = knots at boundaries: replicated [1,1,1,1]*tj(1) or extended tj(1)+[-3,-2,-1,0]*dt1
   %  - keep_reduc    = flag indicating that tj is/isnt kept in reduced knots t^*

   % #make_knot_vector -1
   if (nargin<6 | isempty(idebug))
    idebug = 1;
   end
   if (idebug>=3)
    make_plot = 6;
   else
    make_plot = 0;
   end

   nmeas  = length(si);
   ikinks = sort(ikinks);
   iaccel = sort(iaccel);
   if (size(ikinks,1)>size(ikinks,2)), ikinks = ikinks'; end
   if (size(iaccel,1)>size(iaccel,2)), iaccel = iaccel'; end

   is_ok = cntc.check_kink_accel( si, ikinks, iaccel );
   if (is_ok<=0)
    disp('ERROR with kinks and/or accelerations')
    return;
   end

   if (use_simple)
    [ tj ] = cntc.make_knot_vector_simple(si, ikinks, iaccel, ds_bspl, use_repl, idebug);
   else
    [ tj ] = cntc.make_knot_vector_advanced(si, ikinks, iaccel, ds_bspl, use_repl, idebug);
   end
   nknot = length(tj);

   % filter knots for segments with reduced order

   iall   = [1, ikinks, nmeas];
   iclose = find(diff(iall)==1);
   isegm  = iall(iclose);

   keep_reduc = ones(size(tj));
   for i = 1 : length(isegm)
    sr   = si(isegm(i));
    jsta = find(tj>sr, 1, 'first') - 4;
    keep_reduc(jsta+[1,2]) = 0;
   end

   % plot knot multiplicities & reduced knots

   if (make_plot)
    ksi  = unique(tj);
    nksi = length(ksi);
    mult_j = zeros(nksi,1);
    mult_r = zeros(nksi,1);
    for iksi = 1 : nksi
     mult_j(iksi) = nnz(tj==ksi(iksi));
     mult_r(iksi) = nnz(tj==ksi(iksi) & keep_reduc);
    end

    figure(make_plot); clf; hold on;
    plot(si, zeros(size(si)), '-o');

    for i = 1 : length(si)
     text(si(i), -0.03, num2str(i), 'horizontalalignment','center', 'verticalalignment','top')
    end

    for m = 1 : 4
     ik = find(mult_j>=m);
     plot(ksi(ik), 0.1*m*ones(size(ik)), '*', 'color',cntc.matlab_color(2));
     ik = find(mult_r>=m);
     plot(ksi(ik), 0.1*m*ones(size(ik)), 'o', 'color',cntc.matlab_color(4));
    end
    v = axis; v(3:4) = [-0.2 1]; axis(v);
    grid on;
    legend('measurement s_i', 'knots t_j', 'reduced knots t^*_j', 'location','northwest');
    xlabel('s / t [mm]');
   end

   if (any(isnan(tj)))
    disp(sprintf('knot vector has %d NaN-values at positions %d, %d, %d, %d, %d, %d, %d, %d', ...
     nnz(isnan(tj)), find(isnan(tj))));
   end
   if (idebug>=3)
    disp(' tj,  keep_reduc:')
    disp([ tj , keep_reduc ]')
    %  disp(tj(end-5:end)')
   end

  end % make_knot_vector

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ tj ] = make_knot_vector_simple( si, ikinks, iaccel, ds_bspl, use_repl, idebug )

   % [ tj ] = make_knot_vector_simple( s_meas, ikinks, iaccel, ds_bspl, use_repl, [idebug] )
   %
   % simple version that just applies ds_bspl as much as possible -- needs lambda>0 to avoid singular systems

   % #make_knot_vector_simple
   nmeas  = length(si);
   nkink  = length(ikinks);
   naccel = length(iaccel);

   % define sections between kinks and acceleration points

   iranges = sort([ 1 , ikinks, iaccel, nmeas ]);
   nsec    = length(iranges) - 1;

   tj   = [];
   tiny = 1e-12;

   % loop over sections

   for isec = 1 : nsec

    % take next section

    sec     = [iranges(isec), iranges(isec+1)];

    % determine actual step dt \approx ds_bspl that fits integer multiple times

    len_sec = si(sec(2)) - si(sec(1));
    if (diff(sec)<=1)
     nintv   = 1;      % double kink or kink next to boundary
    else
     nintv   = max(1, round(len_sec / ds_bspl));
    end
    dt      = len_sec / nintv;

    % add knots for this section including start/end-points

    tj   = [tj; si(sec(1)) + [0:nintv]' * dt ];

    % accelerations will be repeated once, end of isec & start of isec+1
    % kinks will be repeated once more

    if (any(ikinks==sec(2)))
     tj = [ tj ; tj(end) ];
     if (idebug>=1)
      disp(sprintf('kink i = %3d: repeat tj = %7.3f', sec(2), tj(end)));
     end
    end
   end

   if (use_repl)
    % repeat first and last knots three times for cubic spline
    tj   = [ [1;1;1]*tj(1); tj; [1;1;1]*tj(end) ];
   else
    % extend by repeating first and last steps three times for cubic spline
    dt0  = tj(2) - tj(1);
    dt1  = tj(end) - tj(end-1);
    tj   = [ tj(1)+[-3:0]'*dt0; tj(2:end-1); tj(end)+[0:3]'*dt1 ];
   end

  end % make_knot_vector_simple

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ tj ] = make_knot_vector_advanced( si, ikinks, iaccel, ds_bspl, use_repl, idebug )

   % [ tj ] = make_knot_vector_advanced( s_meas, ikinks, iaccel, ds_bspl, use_repl, [idebug] )
   %
   % advanced knot-vector with steps >= ds_bspl, using at most one knot per data-interval,
   %          with measurements skipped for not-a-knot conditions at boundaries, kinks and accelerations

   % define sections between kinks and acceleration points
   % #make_knot_vector_advanced

   nmeas   = length(si);
   iranges = sort([ 1 , ikinks, iaccel, nmeas ]);
   nsec    = length(iranges) - 1;

   tj   = [];
   tiny = 1e-12;

   % loop over sections

   for isec = 1 : nsec

    % take next section

    sec     = [iranges(isec), iranges(isec+1)];

    % determine whether 2nd measurement and before-last measurements should be skipped

    start_kink = (any(ikinks==sec(1)));
    end_kink   = (any(ikinks==sec(2)));
    end_accel  = (any(iaccel==sec(2)));

    skip_2nd   = (isec==1 | start_kink);
    skip_2last = (isec==nsec | end_kink | end_accel);

    % determine step dt \approx ds_bspl that fits integer multiple times
    %    this will be the default for steps in this section

    len_sec = si(sec(2)) - si(sec(1));
    nintv   = round(len_sec / ds_bspl);
    dt      = len_sec / nintv;

    if (idebug>=1)
     disp(sprintf('section %d: ip=[%3d,%3d], s=[%8.3f,%8.3f], dt=%7.4f, nintv=%4d', ...
      isec, sec, si(sec(1)), si(sec(2)), dt, nintv));
     if (isec==1),    disp('           first section - skip 2nd measurement'); end
     if (start_kink), disp('           kink at start - skip 2nd measurement'); end
     if (isec==nsec), disp('           last section  - skip before-last measurement'); end
     if (end_kink),   disp('           kink at end   - skip before-last measurement'); end
     if (end_accel),  disp('           accel at end  - skip before-last measurement'); end
    end

    % add knots for this section

    tj   = [tj; si(sec(1))];

    iprf = sec(1);
    while(iprf < sec(2))
     % disp(sprintf('iprf=%2d, tj(end)=%6.3f, si+dt=%6.3f', iprf, tj(end), si(iprf)+dt));

     % take step dt

     tj_new = tj(end) + dt;

     % enlarge if needed to go to next meas.segment

     tj_new = max(tj_new, si(iprf+1));

     % enlarge if needed to go to third meas.segment

     if (skip_2nd & iprf==sec(1) & sec(2)-sec(1)>=2)
      tj_new = max(tj_new, si(iprf+2));
     end

     % if we are inside the very last two meas.segments, jump to end of section

     if (skip_2last & iprf>=sec(2)-2)
      tj_new = si(sec(2));
     end

     % if step takes us inside the last meas.segment, jump to end of section

     if (skip_2last & tj_new >= si(sec(2)-1))
      tj_new = si(sec(2));
     end

     % if step takes us inside before-last meas.segment, go to beginning of this segment

     if (skip_2last & sec(2)-sec(1)>=2 & tj_new >= si(sec(2)-2) & tj_new <= si(sec(2)-1))
      tj_new = si(sec(2)-2);
     end

     % never go beyond end of section

     tj_new = min(tj_new, si(sec(2)));

     % add to list of knots

     tj = [tj; tj_new];

     % increment iprf to the last measurement that has been covered

     iprf = find(si<=tj(end)+tiny, 1, 'last');
    end

    % accelerations will be repeated once, end of isec / start of isec+1
    % kinks will be repeated once more

    if (end_kink)
     tj = [ tj ; tj(end) ];
     if (idebug>=1)
      disp(sprintf('kink i = %3d: repeat tj = %7.3f', sec(2), tj(end)));
     end
    end
   end

   if (use_repl)
    % repeat first and last knots three times for cubic spline
    tj   = [ [1;1;1]*tj(1); tj; [1;1;1]*tj(end) ];
   else
    % extend by repeating first and last steps three times for cubic spline
    dt0  = tj(2) - tj(1);
    dt1  = tj(end) - tj(end-1);
    tj   = [ tj(1)+[-3:0]'*dt0; tj(2:end-1); tj(end)+[0:3]'*dt1 ];
   end

  end % make_knot_vector_advanced

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ Rmat ] = reduce_matrix(Mat, tj, keep_reduc, idebug)

   % [ Rmat ] = reduce_matrix(Mat, tj, keep_reduc)
   %
   % Combine columns for segments using a reduced order, e.g. linear segments between double kinks
   %   tj         == knot vector (1:nknot)
   %   keep_reduc == flags for knots to be kept in reduced knot vector tr = t^* (1:nreduc)
   % #reduce_matrix -2

   if (nargin<4 | isempty(idebug))
    idebug = 0;
   end

   nknot  = length(tj);

   jsta_reduc = find(keep_reduc(1:end-1) & ~keep_reduc(2:end));
   if (idebug>=2 & ~isempty(jsta_reduc))
    disp('Using reduced order at knots jsta =');
    disp(jsta_reduc');
   end

   % select indices j_orig of knots that are kept in reduced knot vector

   j_orig = find(keep_reduc);
   tr     = tj(j_orig);

   % new number for knots that are kept

   j_reduc = cumsum(keep_reduc) .* keep_reduc;
   nreduc  = length(tr);

   if (size(Mat,2)~=nknot-4)
    disp('Incorrect size Mat')
    disp([size(Mat,2), nknot-4])
    return
   end

   % combine columns of Mat(.,nknot-4) into Rmat(.,nreduc-4)

   %  - copy columns that are kept as a whole

   Rmat  = Mat(:,j_orig(1:nreduc-4));

   %  - add columns that are deleted

   for j0 = jsta_reduc
    Rmat(:,j_reduc(j0)  ) = Rmat(:,j_reduc(j0)  ) + 2/3 * Mat(:,j0+1) + 1/3 * Mat(:,j0+2);
    Rmat(:,j_reduc(j0+3)) = Rmat(:,j_reduc(j0+3)) + 1/3 * Mat(:,j0+1) + 2/3 * Mat(:,j0+2);
   end


  end % reduce_matrix

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ a ] = expand_vector( ra, tj, keep_reduc, idebug )

   % [ a ] = expand_vector( ra, tj, keep_reduc )
   %    ra = coefficients for reduced basisfunctions a^*
   %    tj = full knot vector t_j (1:nknot)
   %    keep_reduc = flags for knots kept in reduced knot vector tr = t^* (1:nreduc)
   % #expand_vector -2

   if (nargin<4 | isempty(idebug))
    idebug = 0;
   end

   nknot   = length(tj);
   nreduc  = nnz(keep_reduc);

   % locate knots jsta_reduc at start of segment with reduced order

   jsta_reduc = find(keep_reduc(1:end-1) & ~keep_reduc(2:end));
   if (idebug>=2 & ~isempty(jsta_reduc))
    disp('Using reduced order at knots jsta =');
    disp(jsta_reduc');
   end

   % select indices j_orig of knots that are kept in reduced knot vector

   j_orig = find(keep_reduc);

   % new number for knots that are kept

   j_reduc = cumsum(keep_reduc) .* keep_reduc;

   % expand coefficients ra(1:nreduc-4) to coefficients a(1:nknot-4)

   a = zeros(nknot-4,1);

   for j = 1 : nknot-4
    if     (~keep_reduc(j) & any(jsta_reduc==j-1))

     a(j) = 2/3 * ra(j_reduc(j-1)) + 1/3 * ra(j_reduc(j+2));

    elseif (~keep_reduc(j) & any(jsta_reduc==j-2))

     a(j) = 1/3 * ra(j_reduc(j-2)) + 2/3 * ra(j_reduc(j+1));

    elseif (keep_reduc(j))

     a(j) = ra(j_reduc(j));

    end
   end

  end % expand_vector

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ spl ] = make_pp_spline( s, y, z, lambda, wgt, ikinks, idebug )

   % [ spl ] = make_pp_spline( s, y, z, lambda, wgt, ikinks, idebug )
   %
   % compute parametric cubic smoothing spline (y(s),z(s)) using PP-form method

   % #make_pp_spline -2

   if (nargin<7 | isempty(idebug))
    idebug = 0;
   end
   use_knot = 1; % using free ends (0) or not-a-knot boundaries (1).

   % arguments have been checked: s, y, z, wgt are all [npnt,1]

   % combine measurements that lie too close together

   ds_min = 1e-3 * max(diff(s));
   [s, y, z, wgt] = cntc.combine_close_meas(s, y, z, wgt, ds_min);
   % disp(sprintf('Ratio hmax / hmin = %3.1e', max(diff(s)) / min(diff(s)) ));

   npnt = length(s);

   % process separate sections between break-points

   nsec    = length(ikinks) + 1;
   iranges = [ 1 , ikinks , npnt ];
   ay0 = []; ay1 = []; ay2 = []; ay3 = [];
   az0 = []; az1 = []; az2 = []; az3 = [];

   for isec = 1 : nsec
    rg = [iranges(isec) : iranges(isec+1)];
    if (idebug>=1)
     disp(sprintf('PPspline: make section %2d, rg=[%4d:%4d], s0=%10.6f', isec, rg(1), rg(end), s(rg(1)) ));
    end

    [a0, a1, a2, a3] = cntc.make_pp_spline_section( s(rg), y(rg), lambda, wgt(rg), use_knot, idebug );
    if (isec<nsec)       % a0,a2 have npnt values, a1, a3 have nseg values
     ay0 = [ay0; a0(1:end-1)]; ay1 = [ay1; a1];
     ay2 = [ay2; a2(1:end-1)]; ay3 = [ay3; a3];
    else
     ay0 = [ay0; a0]; ay1 = [ay1; a1];
     ay2 = [ay2; a2]; ay3 = [ay3; a3];
    end

    [a0, a1, a2, a3] = cntc.make_pp_spline_section( s(rg), z(rg), lambda, wgt(rg), use_knot, idebug );
    if (isec<nsec)
     az0 = [az0; a0(1:end-1)]; az1 = [az1; a1];
     az2 = [az2; a2(1:end-1)]; az3 = [az3; a3];
    else
     az0 = [az0; a0]; az1 = [az1; a1];
     az2 = [az2; a2]; az3 = [az3; a3];
    end
   end

   % add a1, a3 at end-point to make all arrays equal size
   % v (s)  = a0 + a1     * s +   a2     * s^2 +   a3     * s^3
   % v'(s1) =      a1(s0)     + 2*a2(s0) * ds  + 3*a3(s0) * ds^2

   npnt = length(s);
   ds   = s(npnt) - s(npnt-1);
   ay1(npnt) = ay1(npnt-1) + 2 * ds * ay2(npnt-1) + 3 * ds^2 * ay3(npnt-1);
   az1(npnt) = az1(npnt-1) + 2 * ds * az2(npnt-1) + 3 * ds^2 * az3(npnt-1);
   ay3(npnt) = ay3(npnt-1);
   az3(npnt) = az3(npnt-1);

   % create structure with final result

   spl = struct('y',y, 'z',z, 's',s, 'lambda',lambda, 'wgt',wgt, 'ikinks', ikinks, 'iaccel', [], ...
    'ay0',ay0, 'ay1',ay1, 'ay2',ay2, 'ay3',ay3, ...
    'az0',az0, 'az1',az1, 'az2',az2, 'az3',az3);

  end % make_pp_spline

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [a0, a1, a2, a3] = make_pp_spline_section( s, y, lambda, wgt, use_knot, idebug )

   % compute smoothing spline {a0--a3} for data {s,y} with
   % mis-fit weight lambda and data weights wgt
   % #make_pp_spline_section

   if (nargin<6 | isempty(idebug))
    idebug = 0;
   end

   % data   has        npnt   points    numbered ipnt = 1..npnt
   % spline has nseg = npnt-1 segments  numbered iseg = 1..nseg

   npnt = length(s);
   nseg = npnt - 1;
   % disp(sprintf('%d points, %d intervals, %d equations',npnt,nseg,npnt));

   if (size(s,2)>size(s,1)), s = s'; end
   if (size(y,2)>size(y,1)), y = y'; end
   a3 = zeros(nseg,1);
   a2 = zeros(npnt,1);
   a1 = zeros(nseg,1);
   a0 = zeros(npnt,1);

   % compute inverse weight matrix Sigma
   % Note: should this be squared?

   Sigma = spdiags(1./wgt, 0, npnt, npnt);

   % set relative importance of smoothness of the result

   kappa = 2/3 * lambda;

   if (idebug>=5)
    disp(sprintf('kappa = %6.3f, Sigma =', kappa));
    disp(full(Sigma));
   end

   % compute distances between points

   h   = diff(s);

   if (idebug>=5)
    disp('h =');
    disp(h');
   end

   % build R- and Q-matrices: npnt equations, including boundary conditions
   %  R * a2 = Q^T * a0

   R   = sparse(npnt,npnt);
   Q   = sparse(npnt,npnt);

   if (~use_knot | nseg<=1)  % free boundary at start
    R(1, 1) = 2*h(1);
    Q(1, 1) = 0;
   else            % not-a-knot boundary
    R(1, 1) = -h(2);
    R(1, 2) =  h(1) + h(2);
    R(1, 3) = -h(1);
    Q(1, 1) =  0;
    Q(2, 1) =  0;
    Q(3, 1) =  0;
   end

   % equations for interior points

   for i = 2 : npnt-1
    R(i,i-1) =    h(i-1);
    R(i,i)   = 2*(h(i-1) + h(i));
    R(i,i+1) =             h(i);

    Q(i-1,i) =             3 / h(i-1);
    Q(i  ,i) = -3 / h(i) - 3 / h(i-1);
    Q(i+1,i) =  3 / h(i);
   end

   if (~use_knot | nseg<=1)  % free boundary at end
    R(npnt,npnt) = 1;
    Q(npnt,npnt) = 0;
   else
    R(npnt,npnt-2) = -h(nseg);
    R(npnt,npnt-1) =  h(nseg-1) + h(nseg);
    R(npnt,npnt  ) = -h(nseg-1);
    Q(npnt-2,npnt) = 0;
    Q(npnt-1,npnt) = 0;
    Q(npnt  ,npnt) = 0;
   end

   % solve system for a2 and a0

   A   = kappa * Q' * Sigma * Q + R;
   rhs = Q' * y;

   if (idebug>=5)
    disp('R =');
    disp(full(R));
    disp('Qt =');
    disp(full(Q'));
    disp('kappa * Q'' * Sigma * Q =');
    disp(full(kappa*Q'*Sigma*Q));
    disp('A =');
    disp(full(A));
    disp('rhs =');
    disp(rhs');
   end

   % scale free boundaries at start & end
   % note: Q ~ 1/h, wgt ~ 1/h, therefore A ~ 1/h^3

   mx_D         = max(abs(diag(A)));
   A(1,1)       = mx_D * A(1,1);
   A(npnt,npnt) = mx_D * A(npnt,npnt);

   c   = condest(A);       % condition c grows with n^4 * (hmax/hmin)^2 * f(lambda)
   if (idebug>=2)
    disp(sprintf('nmeas = %d, lambda = %4.0f: cond(A) = %3.1e', npnt, lambda, c));
   end
   if (c>1e10)
    disp(sprintf('Warning: cond(A) = %3.1e',c));
   end
   if (isnan(c))
    disp('h(1:6):')
    disp(h(1:6))
    disp('R(1;5,1:6):')
    disp(R(1:5,1:6))
    disp('Q(1;5,1:6):')
    disp(Q(1:5,1:6))
   end

   if (1==1)

    % direct implementation:

    % imid = round(npnt/2);
    % tmp = Q' * Q;
    % disp( R(imid,imid-3:imid+3) )

    a2   = A \ rhs;

   else

    if (~use_knot)
     % eliminate boundary values b_1 = b_np = 0 from equations 2 and npnt-1 for LDLt decomposition
     A(2,1) = 0;
     A(npnt-1,end) = 0;

     % use L-D-Lt factorization

     [L,D,P] = ldl(A);    % Pt * A * P = L * D * Lt
     if (any(P~=eye(npnt)))
      disp('LDLt: Non-trivial permutation');
     end

     tmp1 = L \ rhs;
     tmp2 = D \ tmp1;
     a2   = L' \ tmp2;

    else

     [L,U,P] = lu(A);    % P * A = L * U
     if (any(P~=eye(npnt)))
      disp('LU: Non-trivial permutation');
     end

     % A * a2 = rhs -->  P * A * a2 = P * rhs --> L * U * a2 = P * rhs
     % U * a2 = tmp -->  L * tmp = P * rhs
     tmp = L \ (P * rhs);
     a2  = U \ tmp;

    end

   end

   a0 = y - kappa * Sigma * Q * a2;

   % spline parameters

   a3(1:nseg) = diff(a2(1:npnt)) ./ (3*h(1:nseg));
   a1(1:nseg) = diff(a0(1:npnt)) ./ h(1:nseg) ...
    - a3(1:nseg).*h(1:nseg).^2 - a2(1:nseg).*h(1:nseg);

   if (idebug>=5)
    disp('a0 =');
    disp(a0');
    disp('a1 =');
    disp(a1');
    disp('a2 =');
    disp(a2');
    disp('a3 =');
    disp(a3');
   end

   if (0==1 & is_z)
    ix = 465;
    grad_L1 = -inv(Sigma) * (y - a0);
    mat = Q;
    % mat(ix, ix+[-5:3])
    grad_L2 =  kappa * Q * inv(R) * Q' * a0;
    % [y(ix), a0(ix), grad_L1(ix), -grad_L2(ix)]

    figure(5); clf; hold on;
    % spy(abs(mat)>1e5);
    l = plot(grad_L1+grad_L2, '-*'); % set(l([2,4]),'linestyle','--');
   end

  end % make_pp_spline_section

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ s, y, z, wgt ] = combine_close_meas( s, y, z, wgt, ds_min )

   % [ s, y, z, wgt ] = combine_close_meas( s, y, z, lambda, wgt, ikinks, idebug )
   %
   % combine measurements that lie too close together
   % #combine_close_meas -2

   iter = 0;
   comb = ones(size(s));
   ds   = diff(s);
   ix   = find(ds < ds_min);

   while(~isempty(ix))
    iter = iter + 1;
    disp(sprintf('It %d: merging %d measurements with neighbouring ones', iter, length(ix)));

    % replace s, y, z by weighted averages
    s(ix)   = (comb(ix) .* s(ix) + comb(ix+1) .* s(ix+1)) ./ (comb(ix) + comb(ix+1));
    y(ix)   = (comb(ix) .* y(ix) + comb(ix+1) .* y(ix+1)) ./ (comb(ix) + comb(ix+1));
    z(ix)   = (comb(ix) .* z(ix) + comb(ix+1) .* z(ix+1)) ./ (comb(ix) + comb(ix+1));
    s(ix+1) = []; y(ix+1) = []; z(ix+1) = [];

    % replace wgt by sum of weights
    wgt(ix) = (wgt(ix) + wgt(ix+1));
    wgt(ix+1) = [];

    % compute ds anew, iterate
    ds = diff(s);
    ix = find(ds < ds_min);
   end % while

  end % combine_close_meas

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% #File

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [npatch, ws_pos, tot_forc, cp_pos, cp_creep, cp_force] = parse_out1(fname, ire_out, idebug)
   %
   % [ sol ] = parse_out1(fname, [ire_out], [idebug])
   %    or
   % [ npatch, ws_pos, tot_forc, cp_pos, cp_creep, cp_force, cp_subs ] = ...
   %                                                               parse_out1(fname, [ire_out], [idebug])
   %
   % Reads CONTACT output for w/r contact cases (module 1) for result element 'ire_out' from out-file fname
   % and returns the overall values from it:
   %   npatch   = number of contact patches per case
   %   ws_pos   = (input) wheel-set positions, orientations and velocities
   %              struct with x, y, z, roll, yaw, pitch, vx, ... vpitch
   %   tot_forc = total forces in track coordinates & wheelset coordinates
   %              struct with fx_tr, ... fz_ws, mx_tr, ... mz_ws
   %   cp_pos   = contact reference positions
   %              struct with xtr, ytr, ztr, delttr, yr, zr, xw, yw, zw, ncon, nadh, nslip, nplast, prev_icp,
   %              using npatch rows
   %   cp_creep = list of creepages
   %              struct with pen, cksi, ceta, cphi, veloc - npatch rows
   %   cp_force = forces and moments in contact reference coordinates
   %              struct with fn, fx, fs, mn, el.en, fricw, pmax, temp - npatch rows
   %   cp_subs  = maximum subsurface stresses and locations
   %              [sighyd, sigvm, sigtr, sigma1, sigma3, sigxx, sigyy, sigzz]
   %

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
   % #parse_out1 -2

   if (nargin<2 | isempty(ire_out))
    ire_out = -1; % result element number as used in CONTACT library for which output is wanted, -1 == all
   end
   if (nargin<3 | isempty(idebug))
    idebug=0; % display input when idebug>=2, sections found when idebug>=5
   end

   % open file

   if (isempty(strfind(fname,'.out')) & isempty(strfind(fname,'.ref_out')))
    fname=[deblank(fname),'.out'];
   end
   if (~exist(fname,'file'))
    disp(sprintf('\nERROR: cannot find file %s\n',fname))
    npatch=[]; ws_pos=[]; tot_forc=[]; cp_pos=[]; cp_creep=[]; cp_force=[];
    return
   end

   % list the special lines that precede the output data

   headers = strvcat('Case', ...
    'WHEEL-RAIL CONTACT', ...
    'WHEEL-ROLLER CONTACT', ...
    'WHEEL-SET POSITION', ...
    'FLEXIBLE WHEEL-SET DEVIATIONS', ...
    'RAIL IRREGULARITY', ...
    'MASSLESS RAIL DEFLECTION', ...
    'TOTAL FORCES AND MOMENTS', ...
    'AVERAGE CONTACT POSITION', ...
    'FX(TR)      FY(TR)', ...
    'ELAST.EN.   FRIC', ...
    'DATA FOR CONTACT PATCH', ...
    'CONTACT REFERENCE LOCATION', ...
    'KINEMATIC CONSTANTS', ...
    'TOTAL FORCES, TORSIONAL', ...
    'CONTACT STATISTICS', ...
    'ABSMAX SIGHYD');

   subs_keys = strvcat('SIGHYD', 'SIGVM', 'SIGTR', 'SIGMA1', 'SIGMA3', 'SIGXX', 'SIGYY', 'SIGZZ');
   nsubs_keys = size(subs_keys,1);

   % read contents of file into memory

   f = cntc.read_file(fname, idebug);
   f = cntc.find_all_headers(f, headers, idebug);

   % read first/second lines for checking that it is a .out-file.

   iline = 0;
   [s, iline] = cntc.read_line(f, iline, idebug);
   if (~isempty(strfind(s, 'refresh the license')))
    [s, iline] = cntc.read_line(f, iline, idebug);
   end
   if (~isempty(strfind(s, 'w/r contact on result element')))
    % CONTACT library: case starts at line 1
    iline = 0;
   else
    [s, iline] = cntc.read_line(f, iline, idebug);
    if (isempty(strfind(s, 'detailed investigation of 3D')))
     disp('ERROR: this file doesn''t look like a CONTACT .out-file.');
     return;
    end
   end

   % parse the contents of the file

   icase  = 0;
   max_patch = 10;
   ire    = []; jcase  = []; tim    = []; npatch  = [];
   x_ws   = []; y_ws   = []; z_ws   = []; roll    = []; yaw    = []; pitch  = [];
   vx_ws  = []; vy_ws  = []; vz_ws  = []; vroll   = []; vyaw   = []; vpitch = [];
   dx_whl = []; dy_whl = []; dz_whl = []; drollw  = []; dyaww  = []; dpitchw = [];
   vxwhl  = []; vywhl  = []; vzwhl  = []; vrollw  = []; vyaww  = []; vpitchw = [];
   dyrail = []; dzrail = []; drrail = []; vyrail  = []; vzrail = []; vrrail = [];
   kyrail = []; dydefl = []; fydefl = []; kzrail  = []; dzdefl = []; fzdefl = [];
   fx_tr  = []; fy_tr  = []; fz_tr  = []; fx_ws   = []; fy_ws  = []; fz_ws  = [];
   xcp_tr = []; ycp_tr = []; zcp_tr = []; delt_tr = []; ycp_r  = []; zcp_r  = [];
   xcp_w  = []; ycp_w  = []; zcp_w  = [];
   pen    = []; cksi   = []; ceta   = []; cphi    = []; veloc  = [];
   fn_loc = []; fx_loc = []; fs_loc = []; mn_loc  = []; elen   = []; fric   = [];
   temp1  = []; temp2  = []; pmax   = [];
   ncon   = []; nadh   = []; nslip  = []; nplast  = []; prevcp = [];
   subs_max = [];

   while (iline<f.nline)

    % skip lines until one of the headers is found

    [ s, iline, ihdr ] = cntc.next_header(f, iline, idebug);
    if (ihdr<0), break; end

    % Process the consecutive section

    if (~isempty(strfind(s, 'Case')))

     icase = icase + 1;
     if (idebug>=5)
      disp(sprintf('...Case %d starts at line %d',icase,iline));
     end

     % Case     1
     % Case  29246 for w/r contact on result element  2, t=      0.08706221

     tmp = sscanf(s,' Case %d for w/r contact on result element %d, t= %f');
     jcase(icase) = tmp(1);
     if (~isempty(strfind(s, 'result element')))
      ire(icase) = tmp(2);
     end
     if (~isempty(strfind(s, ', t=')))
      tim(icase) = tmp(3);
     end

     % initialize the patch number for module 3

     ipatch = 1;

     % initialize output-values

     npatch(icase)  = 0;
     x_ws(icase)    = NaN; y_ws(icase)    = NaN; z_ws(icase)    = NaN;
     roll(icase)    = NaN; yaw(icase)     = NaN; pitch(icase)   = NaN;
     vx_ws(icase)   = NaN; vy_ws(icase)   = NaN; vz_ws(icase)   = NaN;
     vroll(icase)   = NaN; vyaw(icase)    = NaN; vpitch(icase)  = NaN;
     dx_whl(icase)  = NaN; dy_whl(icase)  = NaN; dz_whl(icase)  = NaN;
     drollw(icase)  = NaN; dyaww(icase)   = NaN; dpitchw(icase) = NaN;
     vxwhl(icase)   = NaN; vywhl(icase)   = NaN; vzwhl(icase)   = NaN;
     vrollw(icase)  = NaN; vyaww(icase)   = NaN; vpitchw(icase) = NaN;
     fx_tr(icase)   = NaN; fy_tr(icase)   = NaN; fz_tr(icase)   = NaN;
     fx_ws(icase)   = NaN; fy_ws(icase)   = NaN; fz_ws(icase)   = NaN;
     mx_tr(icase)   = NaN; my_tr(icase)   = NaN; mz_tr(icase)   = NaN;
     mx_ws(icase)   = NaN; my_ws(icase)   = NaN; mz_ws(icase)   = NaN;

     al = [1:max_patch];
     xcp_tr(al,icase) = NaN; ycp_tr(al,icase) = NaN; zcp_tr(al,icase) = NaN;
     delt_tr(al,icase)= NaN; ycp_r(al,icase)  = NaN; zcp_r(al,icase)  = NaN;
     xcp_w(al,icase)  = NaN; ycp_w(al,icase)  = NaN; zcp_w(al,icase)  = NaN;
     pen(al,icase)    = NaN; veloc(al,icase)  = NaN;
     cksi(al,icase)   = NaN; ceta(al,icase)   = NaN; cphi(al,icase)   = NaN;
     fn_loc(al,icase) = NaN; fx_loc(al,icase) = NaN; fs_loc(al,icase) = NaN;
     mn_loc(al,icase) = NaN; elen(al,icase)   = NaN; fric(al,icase)   = NaN;
     temp1(al,icase)  = NaN; temp2(al,icase)  = NaN; pmax(al,icase)   = NaN;
     ncon(al,icase)   = NaN; nadh(al,icase)   = NaN; nslip(al,icase)  = NaN;
     nplast(al,icase) = NaN; prevcp(al,al,icase) = 0;

     if (~isempty(subs_max))
      subs_max(1:4, al, icase, 1:nsubs_keys) = NaN;
     end
    end

    if (~isempty(strfind(s, 'WHEEL-RAIL CONTACT')) | ~isempty(strfind(s, 'WHEEL-ROLLER CONTACT')))

     % WHEEL-RAIL CONTACT, LEFT WHEEL,  1 CONTACT PATCH

     if (idebug>=5), disp(sprintf('...Found number of patches at line %d',iline)); end
     npatch(icase) = sscanf(s, 'WHEEL-R%*s CONTACT, %*s WHEEL, %d');
    end

    if (~isempty(strfind(s, 'WHEEL-SET POSITION')))

     % WHEEL-SET POSITION AND VELOCITY
     %     X_WS        Y_WS        Z_WS        ROLL        YAW         PITCH
     %     0.000       0.000      0.2154       0.000       0.000       0.000
     %
     %     VX_WS       VY_WS       VZ_WS       VROLL       VYAW        VPITCH
     %     40.00       1.000E-02   0.000       0.000       0.000      -86.93

     if (idebug>=5)
      disp(sprintf('...Found wheel-set position and velocity at line %d',iline));
     end
     [s, iline] = cntc.read_line(f, iline, idebug);
     [s, iline] = cntc.read_line(f, iline, idebug);
     tmp = sscanf(s, '%f %f %f %f %f %f');
     x_ws(icase)   = tmp(1); y_ws(icase)   = tmp(2); z_ws(icase)   = tmp(3);
     roll(icase)   = tmp(4); yaw(icase)    = tmp(5); pitch(icase)  = tmp(6);

     [s, iline] = cntc.read_line(f, iline, idebug);
     [s, iline] = cntc.read_line(f, iline, idebug);
     [s, iline] = cntc.read_line(f, iline, idebug);
     tmp = sscanf(s, '%f %f %f %f %f %f');
     vx_ws(icase)  = tmp(1); vy_ws(icase)  = tmp(2); vz_ws(icase)  = tmp(3);
     vroll(icase)  = tmp(4); vyaw(icase)   = tmp(5); vpitch(icase) = tmp(6);
    end

    if (~isempty(strfind(s, 'FLEXIBLE WHEEL-SET DEVIATIONS')))

     % FLEXIBLE WHEEL-SET DEVIATIONS
     %     DXWHL       DYWHL       DZWHL       DROLLW      DYAWW       DPITCHW
     %         0.000       0.000       0.000  -4.000E-02       0.000       0.000
     %     VXWHL       VYWHL       VZWHL       VROLLW      VYAWW       VPITCHW
     %         0.000       0.000       0.000       0.000       0.000       0.000

     if (idebug>=5)
      disp(sprintf('...Found flexible wheel-set deviations at line %d',iline));
     end
     [s, iline] = cntc.read_line(f, iline, idebug);
     [s, iline] = cntc.read_line(f, iline, idebug);
     tmp = sscanf(s, '%f %f %f %f %f %f');
     dx_whl(icase) = tmp(1); dy_whl(icase) = tmp(2); dz_whl(icase)  = tmp(3);
     drollw(icase) = tmp(4); dyaww(icase)  = tmp(5); dpitchw(icase) = tmp(6);

     [s, iline] = cntc.read_line(f, iline, idebug);
     [s, iline] = cntc.read_line(f, iline, idebug);
     tmp = sscanf(s, '%f %f %f %f %f %f');
     vxwhl(icase)  = tmp(1); vywhl(icase)  = tmp(2); vzwhl(icase)   = tmp(3);
     vrollw(icase) = tmp(4); vyaww(icase)  = tmp(5); vpitchw(icase) = tmp(6);
    end

    if (~isempty(strfind(s, 'RAIL IRREGULARITY')))

     % RAIL IRREGULARITY
     %     DYRAIL      DZRAIL      DROLLR      VYRAIL      VZRAIL      VROLLR
     %    -2.773E-02      0.2268       0.000      -213.1       0.000       0.000

     if (idebug>=5)
      disp(sprintf('...Found rail irregularities at line %d',iline));
     end
     [s, iline] = cntc.read_line(f, iline, idebug);
     [s, iline] = cntc.read_line(f, iline, idebug);
     tmp = sscanf(s, '%f %f %f %f %f %f');
     dyrail(icase) = tmp(1); dzrail(icase) = tmp(2); drrail(icase) = tmp(3);
     vyrail(icase) = tmp(4); vzrail(icase) = tmp(5); vrrail(icase) = tmp(6);
    end

    if (~isempty(strfind(s, 'MASSLESS RAIL DEFLECTION')))

     % MASSLESS RAIL DEFLECTION
     %     KY_RAIL     DY_DEFL     FY_RAIL     KZ_RAIL     DZ_DEFL     FZ_RAIL
     %      100.0000     -1.0000   1.162E+05      0.0000      0.0000   1.230E+05

     if (idebug>=5)
      disp(sprintf('...Found massless rail deflection at line %d',iline));
     end
     [s, iline] = cntc.read_line(f, iline, idebug);
     [s, iline] = cntc.read_line(f, iline, idebug);
     tmp = sscanf(s, '%f %f %f %f %f %f');
     kyrail(icase) = tmp(1); dydefl(icase) = tmp(2); fydefl(icase) = tmp(3);
     kzrail(icase) = tmp(4); dzdefl(icase) = tmp(5); fzdefl(icase) = tmp(6);
    end

    if (~isempty(strfind(s, 'TOTAL FORCES AND MOMENTS')) | ...
      ~isempty(strfind(s, 'FX(TR)      FY(TR)')) )

     % TOTAL FORCES AND MOMENTS ON RAIL | ROLLER
     %     FX(TR)      FY(TR)      FZ(TR)      FX(WS)      FY(WS)      FZ(WS)
     %       70.09      -40.92   1.000E+04       70.09      -40.92   1.000E+04
     %
     %     MX(TR)      MY(TR)      MZ(TR)      MX(WS)      MY(WS)      MZ(WS)
     %   5.170E+07   1.092E+06  -1.812E+07   5.544E+07   1.492E+06  -1.944E+07

     if (idebug>=5)
      disp(sprintf('...Found global forces and moments at line %d',iline));
     end
     % support for compact .out-files without 'TOTAL FORCES...'
     if (isempty(strfind(s, 'FX(TR)      FY(TR)')))
      [s, iline] = cntc.read_line(f, iline, idebug);
     end
     [s, iline] = cntc.read_line(f, iline, idebug);
     tmp = sscanf(s, '%f %f %f %f %f %f');
     fx_tr(icase) = tmp(1); fy_tr(icase) = tmp(2); fz_tr(icase) = tmp(3);
     fx_ws(icase) = tmp(4); fy_ws(icase) = tmp(5); fz_ws(icase) = tmp(6);

     % check for MX(TR) on following lines
     [s, iline] = cntc.read_line(f, iline, idebug);
     [s, iline] = cntc.read_line(f, iline, idebug);
     if (~isempty(strfind(s, 'MX(TR)      MY(TR)')))
      [s, iline] = cntc.read_line(f, iline, idebug);
      tmp = sscanf(s, '%f %f %f %f %f %f');
      mx_tr(icase) = tmp(1); my_tr(icase) = tmp(2); mz_tr(icase) = tmp(3);
      mx_ws(icase) = tmp(4); my_ws(icase) = tmp(5); mz_ws(icase) = tmp(6);
     else
      iline = iline - 2;
      [s, iline] = cntc.read_line(f, iline, idebug);
     end
    end

    if (~isempty(strfind(s, 'DATA FOR CONTACT PATCH')))

     % DATA FOR CONTACT PATCH  1 (NEW PATCH)
     % DATA FOR CONTACT PATCH  1 (PREV  1,  2)

     if (idebug>=5)
      disp(sprintf('...Found contact patch number at line %d',iline));
     end
     ipatch = sscanf(strtrim(s), '----- DATA FOR CONTACT PATCH %d');

     if (~isempty(strfind(s, 'NEW PATCH')))
      prevcp(ipatch,:,icase) = 0;
     else
      [tmp, nval] = sscanf(strtrim(s), '----- DATA FOR CONTACT PATCH %*d (PREV %d, %d, %d, %d)');
      prevcp(ipatch,1:nval,icase) = tmp;
     end

     % support for compact .out-files without 'WHEEL-RAIL CONTACT...'
     if (ipatch>npatch(icase))
      npatch(icase) = ipatch;
     end
    end

    if (~isempty(strfind(s, 'CONTACT REFERENCE LOCATION')))

     % CONTACT REFERENCE LOCATION
     %     XCP(TR)     YCP(TR)     ZCP(TR)    DELTCP(TR)   YCP(R)      ZCP(R)
     %      0.0000   -751.0869      0.1489      0.0312     -9.4445      0.1490
     %
     %     XCP(W)      YCP(W)      ZCP(W)
     %      0.0000      1.0869     -0.0491

     if (idebug>=5)
      disp(sprintf('...Found contact point location at line %d',iline));
     end
     [s, iline] = cntc.read_line(f, iline, idebug);
     [s, iline] = cntc.read_line(f, iline, idebug);
     tmp = sscanf(s, '%f %f %f %f %f %f');
     xcp_tr(ipatch,icase) = tmp(1); ycp_tr(ipatch,icase)  = tmp(2);
     zcp_tr(ipatch,icase) = tmp(3); delt_tr(ipatch,icase) = tmp(4);
     ycp_r(ipatch,icase)  = tmp(5); zcp_r(ipatch,icase)   = tmp(6);

     [s, iline] = cntc.read_line(f, iline, idebug);
     [s, iline] = cntc.read_line(f, iline, idebug);
     [s, iline] = cntc.read_line(f, iline, idebug);
     tmp = sscanf(s, '%f %f %f');
     xcp_w(ipatch,icase) = tmp(1); ycp_w(ipatch,icase) = tmp(2);
     zcp_w(ipatch,icase) = tmp(3);
    end

    if (~isempty(strfind(s, 'KINEMATIC CONSTANTS')))

     % KINEMATIC CONSTANTS
     %     CHI         DQ          VELOC       CKSI        CETA        CPHI
     %     0.000       1.000           1.000   7.282E-04   7.938E-04  -1.736E-04

     % get the creepages of a contact patch

     if (idebug>=5)
      disp(sprintf('...Found kinematic constants at line %d',iline));
     end
     [s1, iline] = cntc.read_line(f, iline, idebug);
     [s2, iline] = cntc.read_line(f, iline, idebug);
     tmp = sscanf(s2, '%f %f %f %f %f %f');
     veloc(ipatch,icase)  = tmp(3);
     cksi(ipatch,icase)   = tmp(4);
     ceta(ipatch,icase)   = tmp(5);
     cphi(ipatch,icase)   = tmp(6);
    end

    if (~isempty(strfind(s, 'TOTAL FORCES, TORSIONAL')) | ...
      ~isempty(strfind(s, 'ELAST.EN.   FRIC')) )

     % T=1, F=2:
     % TOTAL FORCES, TORSIONAL MOMENT, ELASTIC ENERGY, FRICTIONAL POWER
     %     FN          FX          FS          MN         ELAST.EN.  FRIC.POWER
     %     25.00      0.1885       0.000       0.000       4.736       0.000
     %     FN/G       SHIFT X      SHIFT Y    APPROACH       PMAX
     %     25.00      -2.708E-03  -1.106E-15  6.492E-03    12.34

     % get the creepages or forces

     if (idebug>=5)
      disp(sprintf('...Found contact patch forces at line %d',iline));
     end
     % support for compact .out-files without 'TOTAL FORCES...'
     if (isempty(strfind(s, 'ELAST.EN.   FRIC')))
      [s , iline] = cntc.read_line(f, iline, idebug);
     end
     [s1, iline] = cntc.read_line(f, iline, idebug);
     tmp = sscanf(s1, '%f %f %f %f %f %f');
     fn_loc(ipatch,icase)    = tmp(1);
     fx_loc(ipatch,icase)    = tmp(2); % storing absolute force
     fs_loc(ipatch,icase)    = tmp(3); % storing absolute force
     mn_loc(ipatch,icase)    = tmp(4);
     elen(ipatch,icase)      = tmp(5);
     fric(ipatch,icase)      = tmp(6);

     [s2, iline] = cntc.read_line(f, iline, idebug);
     [s3, iline] = cntc.read_line(f, iline, idebug);
     tmp = sscanf(s3, '%f %f %f %f %f %f');
     if (~isempty([strfind(s2, 'CREEP X'), strfind(s2, 'SHIFT X')]))
      cksi(ipatch,icase)   = tmp(2);
      % else
      %    fx_loc(ipatch,icase) = tmp(2); % storing relative force
     end
     if (~isempty([strfind(s2, 'CREEP Y'), strfind(s2, 'SHIFT Y')]))
      ceta(ipatch,icase)   = tmp(3);
      % else
      %    fs_loc(ipatch,icase) = tmp(3); % storing relative force
     end
     pen(ipatch,icase) = tmp(4);
     if (~isempty(strfind(s2, 'PMAX')))
      pmax(ipatch,icase)  = tmp(5);              % >= v22.1
     end
     if     (~isempty(strfind(s2, 'MAX(T1)')))     % <= v21.1
      temp1(ipatch,icase) = tmp(5);
      temp2(ipatch,icase) = tmp(6);
     elseif (~isempty(strfind(s2, 'MAX(T1,T2)')))  % >= v22.1
      temp1(ipatch,icase) = tmp(6);
     end
    end % TOTAL FORCES, TORSIONAL

    if (~isempty(strfind(s, 'CONTACT STATISTICS')))

     % CONTACT STATISTICS
     % N: NUMBER OF ELEMENTS IN REGION,  I: NUMBER OF ITERATIONS,
     % POT: POTENTIAL CONTACT AREA,  CON: CONTACT,  ADH: ADHESION
     %    NPOT   NCON   NADH  NSLIP  INORM  ITANG
     %   3731   1286    919    367      1      1
     % or
     %    NPOT   NCON   NADH  NSLIP  NPLAST INORM  ITANG
     %      30     26     14      8      4      1      1

     % get the number of elements in contact/adhesion/slip

     if (idebug>=5)
      disp(sprintf('...Found contact statistics at line %d',iline));
     end
     [s , iline] = cntc.read_line(f, iline, idebug);
     [s , iline] = cntc.read_line(f, iline, idebug);
     [s , iline] = cntc.read_line(f, iline, idebug);
     [s1, iline] = cntc.read_line(f, iline, idebug);
     tmp = sscanf(s1, '%f %f %f %f %f %f');
     ncon(ipatch,icase)    = tmp(2);
     nadh(ipatch,icase)    = tmp(3);
     nslip(ipatch,icase)   = tmp(4);
     if (~isempty(strfind(s, 'NPLAST')))
      nplast(ipatch,icase) = tmp(5);
     end
    end % STATISTICS

    if (~isempty(strfind(s, 'ABSMAX SIGHYD')))

     if (isempty(subs_max))
      subs_max(1:4, 1:max_patch, 1:icase, 1:nsubs_keys) = NaN;
     end

     % ABSMAX SIGHYD =    -0.853 AT (X,Y,Z) = (   0.000,   0.000,   0.000)
     % MAX     SIGVM =     0.665 AT (X,Y,Z) = (   0.000,   0.000,   0.400)
     % MAX     SIGTR =     0.665 AT (X,Y,Z) = (   0.000,   0.000,   0.400)
     % MAX    SIGMA1 =     0.008 AT (X,Y,Z) = (   0.000,   0.000,   1.250)
     % MIN    SIGMA3 =    -1.000 AT (X,Y,Z) = (   0.000,   0.000,   0.000)
     % ABSMAX  SIGXX =    -0.780 AT (X,Y,Z) = (   0.000,   0.000,   0.000)
     % ABSMAX  SIGYY =    -0.780 AT (X,Y,Z) = (   0.000,   0.000,   0.000)
     % ABSMAX  SIGZZ =    -1.000 AT (X,Y,Z) = (   0.000,   0.000,   0.000)

     % process all lines containing a maximum value

     while(iline<=f.nline & ~isempty(strfind(s, 'AT (X,Y,Z)')))

      % get the key number
      key = strtrim(s(8:14));
      ikey = 1;
      while(ikey<nsubs_keys & ~strcmp(key, deblank(subs_keys(ikey,:))))
       ikey = ikey + 1;
      end
      if (~strcmp(key, deblank(subs_keys(ikey,:))))
       disp(sprintf('Internal error: unknown key "%s".', key));
       break
      end

      % get the values
      tmp = sscanf(s(17:end),' %f AT (X,Y,Z) = ( %f, %f, %f)');

      % store in subs_max
      subs_max(1:4, ipatch, icase, ikey) = tmp;

      % read next line from file
      [ s, iline ] = cntc.read_line(f, iline, idebug);
     end

    end % ABSMAX, SIGHYD

   end % while (iline<f.nline)

   % select cases for result element ire_out, if provided

   if (ire_out>0 & ~isempty(ire))
    ncase    = length(jcase);
    i_select = find(ire == ire_out);
    disp(sprintf('Total #cases = %d, selected %d cases for ire_out = %d', ncase, length(i_select), ire_out));

    npatch = npatch(i_select);
    ire    = ire   (i_select); jcase  = jcase (i_select); tim     = tim    (i_select);
    x_ws   = x_ws  (i_select); y_ws   = y_ws  (i_select); z_ws    = z_ws   (i_select);
    roll   = roll  (i_select); yaw    = yaw   (i_select); pitch   = pitch  (i_select);
    vx_ws  = vx_ws (i_select); vy_ws  = vy_ws (i_select); vz_ws   = vz_ws  (i_select);
    vroll  = vroll (i_select); vyaw   = vyaw  (i_select); vpitch  = vpitch (i_select);
    dx_whl = dx_whl(i_select); dy_whl = dy_whl(i_select); dz_whl  = dz_whl (i_select);
    drollw = drollw(i_select); dyaww  = dyaww (i_select); dpitchw = dpitchw(i_select);
    vxwhl  = vxwhl (i_select); vywhl  = vywhl (i_select); vzwhl   = vzwhl  (i_select);
    vrollw = vrollw(i_select); vyaww  = vyaww (i_select); vpitchw = vpitchw(i_select);
    fx_tr  = fx_tr (i_select); fy_tr  = fy_tr (i_select); fz_tr   = fz_tr  (i_select);
    fx_ws  = fx_ws (i_select); fy_ws  = fy_ws (i_select); fz_ws   = fz_ws  (i_select);
    mx_tr  = mx_tr (i_select); my_tr  = my_tr (i_select); mz_tr   = mz_tr  (i_select);
    mx_ws  = mx_ws (i_select); my_ws  = my_ws (i_select); mz_ws   = mz_ws  (i_select);

    xcp_tr  = xcp_tr (:,i_select); ycp_tr  = ycp_tr(:,i_select); zcp_tr  = zcp_tr(:,i_select);
    delt_tr = delt_tr(:,i_select); ycp_r   = ycp_r (:,i_select); zcp_r   = zcp_r (:,i_select);
    xcp_w   = xcp_w  (:,i_select); ycp_w   = ycp_w (:,i_select); zcp_w   = zcp_w (:,i_select);
    pen     = pen    (:,i_select); veloc   = veloc (:,i_select);
    cksi    = cksi   (:,i_select); ceta    = ceta  (:,i_select); cphi    = cphi  (:,i_select);
    fn_loc  = fn_loc (:,i_select); fx_loc  = fx_loc(:,i_select); fs_loc  = fs_loc(:,i_select);
    mn_loc  = mn_loc (:,i_select); elen    = elen  (:,i_select); fric    = fric  (:,i_select);
    temp1   = temp1  (:,i_select); temp2   = temp2 (:,i_select); pmax    = pmax  (:,i_select);
    ncon    = ncon   (:,i_select); nadh    = nadh  (:,i_select); nslip   = nslip (:,i_select);
    nplast  = nplast (:,i_select); prevcp  = prevcp(:,:,i_select);
    if (~isempty(subs_max))
     subs_max = subs_max(:,:,i_select,:);
    end
   end

   % group arrays into structures

   ws_pos   = struct('x',     x_ws,   'y',     y_ws,   'z',     z_ws,   ...
    'roll',  roll,   'yaw',   yaw,    'pitch', pitch,  ...
    'vx',    vx_ws,  'vy',    vy_ws,  'vz',    vz_ws,  ...
    'vroll', vroll,  'vyaw',  vyaw,   'vpitch',vpitch);
   if (any(~isnan(dx_whl)))
    ws_pos.dx_whl = dx_whl; ws_pos.dy_whl = dy_whl; ws_pos.dz_whl  = dz_whl;
    ws_pos.drollw = drollw; ws_pos.dyaww  = dyaww;  ws_pos.dpitchw = dpitchw;
    ws_pos.vxwhl  = vxwhl;  ws_pos.vywhl  = vywhl;  ws_pos.vzwhl   = vzwhl;
    ws_pos.vrollw = vrollw; ws_pos.vyaww  = vyaww;  ws_pos.vpitchw = vpitchw;
   end
   if (any(~isnan(dyrail)))
    ws_pos.dyrail = dyrail; ws_pos.dzrail = dzrail;  ws_pos.drrail  = drrail;
    ws_pos.vyrail = vyrail; ws_pos.vzrail = vzrail;  ws_pos.vrrail  = vrrail;
   end
   if (any(~isnan(kyrail)))
    ws_pos.kyrail = kyrail; ws_pos.dydefl = dydefl;  ws_pos.fydefl  = fydefl;
    ws_pos.kzrail = kzrail; ws_pos.dzdefl = dzdefl;  ws_pos.fzdefl  = fzdefl;
   end
   tot_forc = struct('fx_tr', fx_tr,  'fy_tr', fy_tr,  'fz_tr', fz_tr,  ...
    'fx_ws', fx_ws,  'fy_ws', fy_ws,  'fz_ws', fz_ws,  ...
    'mx_tr', mx_tr,  'my_tr', my_tr,  'mz_tr', mz_tr,  ...
    'mx_ws', mx_ws,  'my_ws', my_ws,  'mz_ws', mz_ws);

   max_patch = max(npatch);
   al        = [1:max_patch];
   cp_pos    = struct('xcp_tr', xcp_tr(al,:),  'ycp_tr',  ycp_tr(al,:),  ...
    'zcp_tr', zcp_tr(al,:),  'delt_tr', delt_tr(al,:), ...
    'ycp_r',  ycp_r(al,:),   'zcp_r',   zcp_r(al,:),   ...
    'xcp_w',  xcp_w(al,:),   'ycp_w',   ycp_w(al,:),   ...
    'zcp_w',  zcp_w(al,:),   'ncon',    ncon(al,:),    ...
    'nadh',   nadh(al,:),    'nslip',   nslip(al,:),   ...
    'nplast', nplast(al,:),  'prev_icp', prevcp(al,al,:));
   cp_creep  = struct('pen',    pen(al,:),     'cksi',    cksi(al,:),    ...
    'ceta',   ceta(al,:),    'cphi',    cphi(al,:),    ...
    'veloc',  veloc(al,:));
   cp_force  = struct('fn',     fn_loc(al,:),  'fx',      fx_loc(al,:),  ...
    'fs',     fs_loc(al,:),  'mn',      mn_loc(al,:),  ...
    'elen',   elen(al,:),    'fric',    fric(al,:));
   if (~all(isnan(pmax(al,:))))
    cp_force.pmax  = pmax(al,:);
   end
   if (~all(isnan(temp1(al,:))) & all(isnan(temp2(al,:))))
    cp_force.temp  = temp1(al,:);  % >= v22.1
   elseif (~all(isnan(temp1(al,:))))
    cp_force.temp1 = temp1(al,:);  % <= v21.1
    cp_force.temp2 = temp2(al,:);
   end

   if (isempty(subs_max))
    cp_subs   = [];
   else
    cp_subs   = struct('sighyd', squeeze(subs_max(:,al,:,1)), ...
     'sigvm',  squeeze(subs_max(:,al,:,2)), ...
     'sigtr',  squeeze(subs_max(:,al,:,3)), ...
     'sigma1', squeeze(subs_max(:,al,:,4)), ...
     'sigma3', squeeze(subs_max(:,al,:,5)), ...
     'sigxx',  squeeze(subs_max(:,al,:,6)), ...
     'sigyy',  squeeze(subs_max(:,al,:,7)), ...
     'sigzz',  squeeze(subs_max(:,al,:,8)));
   end

   % convert separate structs into one overall struct, if requested

   if (nargout<=1)
    sol = struct( 'ws_pos',ws_pos, 'tot_forc',tot_forc, 'npatch',  npatch, ...
     'cp_pos',cp_pos, 'cp_creep',cp_creep, 'cp_force',cp_force,...
     'cp_subs',cp_subs);
    if (~isempty(ire)), sol.ire = ire; end
    if (~isempty(tim)), sol.tim = tim; end
    sol.icase = jcase;
    npatch = sol;
   end

  end % parse_out1

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ myfile ] = read_file(fname, idebug)
   % read_file: Helper function for reading entire file into memory
   % [ file_contents ] = read_line(fname, idebug)
   % #read_file -2

   myfile.fname = fname;
   myfile.nline = 0;
   myfile.contents = {};

   f = fopen(fname);
   while(~feof(f))
    myfile.nline = myfile.nline + 1;
    myfile.contents{myfile.nline} = fgets(f);
   end
   fclose(f);

  end % read_file

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [ f ] = find_all_headers(f, headers, idebug)
   % find_all_headers: Helper function for locating headers in file
   % [ f ] = find_all_headers(f, headers, idebug)
   % #find_all_headers -2

   f.headers = zeros(f.nline,1);

   num_hdr = size(headers,1);
   for ihdr = 1 : num_hdr
    ix = strfind(f.contents, deblank(headers(ihdr,:)));
    for iline = 1 : f.nline
     if (~isempty(ix{iline}))
      f.headers(iline) = ihdr;
     end
    end
   end

  end % find_all_headers

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ s, iline ] = read_line(f, iline, idebug)

   % read_line: Helper function for getting next line of input
   % [ s, iline ] = read_line(f, iline, idebug)
   % #read_line-2

   iline = iline + 1;
   s = f.contents{iline};

   if (idebug>=2)
    % debug: print line without trailing \n (?)
    if (length(s)>=2 & double(s(end-1))==13)
     fprintf('line %4d: %s\n', iline, s(1:end-2) )
    elseif (length(s)>=1 & double(s(end))==10)
     fprintf('line %4d: %s\n', iline, s(1:end-1) )
    else
     fprintf('line %4d: %s\n', iline, s)
    end
   end

  end % read_line


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ s, iline, ihdr ] = next_header(f, iline, idebug)

   % next_header: Helper function for locating the next section
   % [ s, iline, ihdr ] = next_header(f, iline, idebug)
   % #next_header -2

   iprev = iline;
   inext = iline+1;
   while(inext<f.nline & f.headers(inext)<=0)
    inext = inext + 1;
   end

   if (f.headers(inext)>0)
    if (idebug>=6)
     disp(sprintf('...Next_header: inext = %d, fast-forword to iline = %d',inext,iline+inext));
    end
    iline = inext;
    s     = strtrim(f.contents{iline});
    ihdr  = f.headers(iline);
   else
    if (idebug>=6)
     disp(sprintf('...Next_header: no more headers'));
    end
    iline = f.nline;
    s     = [];
    ihdr  = -1;
   end

   if (idebug>=2)
    for jline = iprev+1 : iline
     sloc = f.contents{jline};
     % debug: print line without trailing \n
     if (length(sloc)>=2 & double(sloc(end-1))==13)
      disp(sprintf('line %4d: %s', jline, sloc(1:end-2) ))
     elseif (length(sloc)>=1 & double(sloc(end))==10)
      disp(sprintf('line %4d: %s', jline, sloc(1:end-1) ))
     else
      disp(sprintf('line %4d: %s', jline, sloc))
     end
    end
   end

  end % next_header

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [creep, force, subs_max] = parse_out3(fname, use_struct, idebug)
   %
   % [ sol ] = parse_out3(fname, use_struct, [idebug])
   %    or
   % [creep, force, subs_max] = parse_out3(fname, [use_struct], [idebug])
   %
   % Reads CONTACT output for basic contact cases (module 3) from out-file fname and
   % returns the overall values from it:
   %   creep    = struct with creepages    [pen, cksi, ceta, cphi, ncon, nadh, nslip; nplast]
   %   force    = struct with total forces [fn, fx, fy, mz, el.en, fricw, pmax, temp1, temp2]
   %   subs_max = struct with maximum subsurface stresses and locations
   %                                       [sighyd, sigvm, sigtr, sigma1, sigma3, sigxx, sigyy, sigzz]
   % for all cases in the file.
   %
   % if use_struct=0, the outputs are returned in arrays instead of structs.

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
   % #parse_out3 -2

   if (nargin<2 | isempty(use_struct))
    use_struct=1; % return structs with values
   end
   if (nargin<3 | isempty(idebug))
    idebug=1; % display input when idebug>=2
   end

   % open file and read second line for checking that it is a .out-file.

   if (isempty(strfind(fname,'.out')) & isempty(strfind(fname,'.ref_out')))
    fname=[deblank(fname),'.out'];
   end
   if (~exist(fname,'file'))
    disp(['ERROR: cannot find file ',fname])
    creep=[]; force=[];
    return
   end
   f=fopen(fname);

   iline=0;
   [s, iline] = cntc.read_line(f, iline, idebug);
   [s, iline] = cntc.read_line(f, iline, idebug);
   if (~strfind(s, 'detailed investigation of 3D'))
    disp('ERROR: this file doesn''t look like a CONTACT .out-file.');
    return;
   end

   % list the special lines that precede the output data

   headers = strvcat('Case', ...
    'KINEMATIC CONSTANTS', ...
    'TOTAL FORCES, TORSIONAL', ...
    'ELAST.EN.   FRIC.', ...
    'CONTACT STATISTICS', ...
    'ABSMAX SIGHYD');

   subs_keys = strvcat('SIGHYD', 'SIGVM', 'SIGTR', 'SIGMA1', 'SIGMA3', 'SIGXX', 'SIGYY', 'SIGZZ');
   nsubs_keys = size(subs_keys,1);

   % parse the contents of the file

   pen = []; cksi = []; ceta = []; cphi = [];
   fn = []; fx = []; fy = []; mz = []; elen = []; fric = []; pmax = []; temp1 = []; temp2 = [];
   ncon = []; nadh = []; nslip = []; nplast = [];
   subs_max = [];

   while(~feof(f))

    % skip lines until one of the headers is found (case-sensitive)

    [ s, iline, ihdr ] = cntc.find_header(f, iline, headers, 1, idebug);
    if (ihdr<0), break; end

    % Process the consecutive section

    if (~isempty(strfind(s, 'Case')))

     % Case     1

     icase = sscanf(s,' Case %d');
     if (idebug>=5), disp(sprintf('...Case %d starts at line %d',icase,iline)); end

     % initialize output-values for this case

     pen(icase)=NaN;  cksi(icase)=NaN; ceta(icase)=NaN;  cphi(icase)=NaN;
     fn(icase)=NaN;   fx(icase)=NaN;   fy(icase)=NaN;    mz(icase)=NaN;
     elen(icase)=NaN; fric(icase)=NaN; pmax(icase)=NaN;  temp1(icase)=NaN;  temp2(icase)=NaN;
     ncon(icase)=NaN; nadh(icase)=NaN; nslip(icase)=NaN; nplast(icase)=NaN;
     subs_max(icase, 1:4, 1:nsubs_keys)=NaN;
    end

    if (~isempty(strfind(s, 'KINEMATIC CONSTANTS')))

     %  KINEMATIC CONSTANTS, SHIFT T=1, F=1:
     %      DT          VELOC     FX/FSTAT/FN   CETA        CPHI
     %      0.001       1.000     -0.6570       0.000       0.000
     %  KINEMATIC CONSTANTS, ROLLING T=3, F=2:
     %       CHI          DQ        VELOC     FX/FSTAT/FN  FY/FSTAT/FN  CPHI
     %      0.000       1.000       1.000      0.2525E-01   0.000       0.000
     % get the creepages or forces

     if (idebug>=9), disp(sprintf('...Found kinematic constants at line %d',iline)); end
     [s1, iline] = cntc.read_line(f, iline, idebug);
     [s2, iline] = cntc.read_line(f, iline, idebug);
     if (strfind(s1,'CHI'))
      tmp = sscanf(s2, '%*f %*f %*f %f %f %f');
     else
      tmp = sscanf(s2, '%*f %*f %f %f %f');
     end
     if (strfind(s1, 'CKSI'))
      cksi(icase) = tmp(1);
     else
      fx(icase)   = tmp(1);
     end
     if (strfind(s1, 'CETA'))
      ceta(icase) = tmp(2);
     else
      fy(icase)   = tmp(2);
     end
     cphi(icase) = tmp(3);
    end

    if (~isempty(strfind(s, 'TOTAL FORCES, TORSIONAL')))
     if (idebug>=9), disp(sprintf('...Found total forces at line %d',iline)); end
     % just skip one line, next line should contain 'ELAST.EN.   FRIC.POWER')
     %                                           or 'ELAST.EN.   FRIC.WORK')
     [s , iline] = cntc.read_line(f, iline, idebug);
    end

    if (~isempty(strfind(s, 'ELAST.EN.  FRIC.')) | ~isempty(strfind(s, 'ELAST.EN.   FRIC.')))
     % T=1, F=2:
     % TOTAL FORCES, TORSIONAL MOMENT, ELASTIC ENERGY, FRICTIONAL POWER
     %     FN          FX          FY          MZ         ELAST.EN.  FRIC.POWER
     %     25.00      0.1885       0.000       0.000       4.736       0.000
     %     FN/G       SHIFT X      SHIFT Y    APPROACH    MAX(T1)     MAX(T1)
     %     25.00      -2.708E-03  -1.106E-15  6.492E-03    637.0       637.0

     % get the creepages or forces

     if (idebug>=9), disp(sprintf('...Found contact patch forces at line %d',iline)); end
     [s1, iline] = cntc.read_line(f, iline, idebug);
     tmp = sscanf(s1, '%f %f %f %f %f %f');
     fn(icase)   = tmp(1);
     mz(icase)   = tmp(4);
     elen(icase) = tmp(5);
     fric(icase) = tmp(6);

     [s2, iline] = cntc.read_line(f, iline, idebug);
     [s3, iline] = cntc.read_line(f, iline, idebug);
     tmp = sscanf(s3, '%f %f %f %f %f %f');
     if (~isempty([strfind(s2, 'CREEP X'), strfind(s2, 'SHIFT X')]))
      cksi(icase)  = tmp(2);
     else
      fx(icase)    = tmp(2);
     end
     if (~isempty([strfind(s2, 'CREEP Y'), strfind(s2, 'SHIFT Y')]))
      ceta(icase)  = tmp(3);
     else
      fy(icase)    = tmp(3);
     end
     if (~isempty(strfind(s2, 'APPROACH')))
      pen(icase)   = tmp(4);
     end
     if (~isempty(strfind(s2, 'PMAX')))        % >= v22.1
      pmax(icase)  = tmp(5);
     end
     if (~isempty(strfind(s2, 'MAX(T1,T2)')))  % >= v22.1
      temp1(icase) = tmp(6);
     elseif (~isempty(strfind(s2, 'MAX(T1)'))) % <= v21.1
      temp1(icase) = tmp(5);
      temp2(icase) = tmp(6);
     end
    end

    if (~isempty(strfind(s, 'CONTACT STATISTICS')))

     % CONTACT STATISTICS
     % N: NUMBER OF ELEMENTS IN REGION,  I: NUMBER OF ITERATIONS,
     % POT: POTENTIAL CONTACT AREA,  CON: CONTACT,  ADH: ADHESION
     %    NPOT   NCON   NADH  NSLIP  INORM  ITANG
     %   3731   1286    919    367      1      1
     % or
     %    NPOT   NCON   NADH  NSLIP  NPLAST INORM  ITANG
     %      30     26     14      8      4      1      1

     % get the number of elements in contact/adhesion/slip

     if (idebug>=5)
      disp(sprintf('...Found contact statistics at line %d',iline));
     end
     [s , iline] = cntc.read_line(f, iline, idebug); % N:
     [s , iline] = cntc.read_line(f, iline, idebug); % POT:
     if (isempty(strfind(s, 'CON:'))), % support older format
      [s , iline] = cntc.read_line(f, iline, idebug);
     end
     [s , iline] = cntc.read_line(f, iline, idebug); % NPOT
     [s1, iline] = cntc.read_line(f, iline, idebug);
     tmp = sscanf(s1, '%f %f %f %f %f %f');
     ncon(icase)  = tmp(2);
     nadh(icase)  = tmp(3);
     nslip(icase) = tmp(4);
     if (~isempty(strfind(s, 'NPLAST')))
      nplast(icase) = tmp(5);
     end
    end

    if (~isempty(strfind(s, 'ABSMAX SIGHYD')))

     % ABSMAX SIGHYD =    -0.853 AT (X,Y,Z) = (   0.000,   0.000,   0.000)
     % MAX     SIGVM =     0.665 AT (X,Y,Z) = (   0.000,   0.000,   0.400)
     % MAX     SIGTR =     0.665 AT (X,Y,Z) = (   0.000,   0.000,   0.400)
     % MAX    SIGMA1 =     0.008 AT (X,Y,Z) = (   0.000,   0.000,   1.250)
     % MIN    SIGMA3 =    -1.000 AT (X,Y,Z) = (   0.000,   0.000,   0.000)
     % ABSMAX  SIGXX =    -0.780 AT (X,Y,Z) = (   0.000,   0.000,   0.000)
     % ABSMAX  SIGYY =    -0.780 AT (X,Y,Z) = (   0.000,   0.000,   0.000)
     % ABSMAX  SIGZZ =    -1.000 AT (X,Y,Z) = (   0.000,   0.000,   0.000)

     % process all lines containing a maximum value

     while(~feof(f) & ~isempty(strfind(s, 'AT (X,Y,Z)')))

      % get the key number
      key = strtrim(s(8:14));
      ikey = 1;
      while(ikey<nsubs_keys & ~strcmp(key, deblank(subs_keys(ikey,:))))
       ikey = ikey + 1;
      end
      if (~strcmp(key, deblank(subs_keys(ikey,:))))
       disp(sprintf('Internal error: unknown key "%s".', key));
       break
      end

      % get the values
      tmp = sscanf(s(17:end),' %f AT (X,Y,Z) = ( %f, %f, %f)');

      % store in subs_max
      subs_max(icase, 1:4, ikey) = tmp;

      % read next line from file
      [ s, iline ] = cntc.read_line(f, iline, idebug);
     end
    end

   end
   fclose(f);

   if (use_struct)
    sol = struct;
    sol.creep = struct('pen',pen, 'cksi',cksi, 'ceta',ceta, 'cphi',cphi, ...
     'ncon',ncon, 'nadh',nadh, 'nslip',nslip, 'nplast',nplast);
    sol.force = struct('fn',fn, 'fx',fx, 'fy',fy, 'mz',mz, 'elen',elen, ...
     'pmax',pmax, 'temp1',temp1, 'temp2',temp2);

    sol.subs = struct();
    if (~isempty(subs_keys))
     for ikey = 1 : nsubs_keys
      key = deblank(lower(subs_keys(ikey,:)));
      sol.subs = setfield(sol.subs, key, subs_max(:,:,ikey));
     end
    end
    creep = sol;
   else
    creep = [pen; cksi; ceta; cphi; ncon; nadh; nslip; nplast]';
    force = [fn; fx; fy; mz; elen; fric; pmax; temp1; temp2]';
   end

  end % parse_out3

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ s, iline, ihdr ] = find_header(f, iline, headers, case_sens, idebug)

   % find_header: Helper function for locating the next section
   % [ s, iline, ihdr ] = find_header(f, iline, headers, case_sens, idebug)
   % case_sens = 1: case-sensitive matching; 0: case-insensitive
   % #find_header -2

   num_hdr = size(headers,1);
   if (~case_sens)
    headers = lower(headers);
   end
   s       = [];
   is_key  = 0;

   while(~is_key & ~feof(f))

    % read next line from input

    [s, iline] = cntc.read_line(f, iline, idebug);
    if (iline==-39), idebug=100; end

    % convert to lower-case when doing case-insensitive matching

    if (case_sens)
     scmp = s;
    else
     scmp = lower(s);
    end

    % check if it contains one of the headers

    for iheader = 1: num_hdr
     if (~isempty(strfind(scmp, deblank(headers(iheader,:)))))
      is_key = iheader;
      break
     end
    end
   end

   s = strtrim( s );
   if (feof(f))
    ihdr = -1;
   else
    ihdr = is_key;
   end

  end % find_header

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% #show

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function out=sol2curve()
   %% #sol2curve Transform to SDT format -2

   LI=cntc.call;
   C1=struct('X',{{}},'Xlab',{{'comp','igy','igx',{'ire';'iwe';'icp'},'iTime'}},'Y',[]);
   C1.X{5}=(1:length(LI.sol))';C1.X{4}=[1 1 1];

   sol=LI.sol{1} ;
   x=cellfun(@(x)x.x,LI.sol,'uni',0)';
   i1=unique(round(horzcat(x{:})/sol.dx));C1.X{3}=i1(:);
   x=cellfun(@(x)x.y,LI.sol,'uni',0)';
   i1=unique(round(horzcat(x{:})/sol.dy));C1.X{2}=i1(:);

   %% fields at gauss points
   st=fieldnames(sol);i1=[length(sol.y) length(sol.x)];
   for j1=1:length(st)
    if ~isequal(size(sol.(st{j1})),i1);continue;end
    j3=size(C1.Y,1)+1;C1.X{1}(j3,1)=st(j1);
    for j2=1:length(LI.sol)
     [~,i2]=ismember(round(LI.sol{j2}.x/sol.dx),C1.X{3});
     [~,i3]=ismember(round(LI.sol{j2}.y/sol.dy),C1.X{2});
     C1.Y(j3,i3,i2,1,j2)=LI.sol{j2}.(st{j1});
    end
    sol=rmfield(sol,st{j1});
   end
   C1.PlotInfo=ii_plp('PlotInfo2D -type "contour"',C1);
   C1.DimPos=[2 3 5 1 4];C1.name='field';
   C1=sdsetprop(C1,'PlotInfo','ua.axProp',{'yscale','linear'});
   out.Cfield=C1;

   %% profiles
   C1=struct('X',{{[],{'s';'y';'z'},(1:length(LI.sol))'}}, ...
    'Xlab',{{'i','comp','iTime'}},'name','prr');
   r2=cellfun(@(x)x.prr.ProfileS,LI.sol,'uni',0);
   C1.Y=horzcat(r2{:});C1.Y=reshape(C1.Y,size(C1.Y,1),1,size(C1.Y,2));
   r2=cellfun(@(x)x.prr.ProfileY,LI.sol,'uni',0);C1.Y(:,2,:)=horzcat(r2{:});
   r2=cellfun(@(x)x.prr.ProfileZ,LI.sol,'uni',0);C1.Y(:,3,:)=horzcat(r2{:});
   C1.X{1}=(1:size(C1.Y,1))';

   C1.PlotInfo=ii_plp('PlotInfo2D -type "contour"',C1);
   C1.DimPos=[1 3 2];
   C1=sdsetprop(C1,'PlotInfo','ua.axProp',{'yscale','linear'});


   out.prr=C1;

   C1=struct('X',{{[],{'s';'y';'z'},(1:length(LI.sol))'}}, ...
    'Xlab',{{'i','comp','iTime'}},'name','prw');
   r2=cellfun(@(x)x.prw.ProfileS,LI.sol,'uni',0);
   C1.Y=horzcat(r2{:});C1.Y=reshape(C1.Y,size(C1.Y,1),1,size(C1.Y,2));
   r2=cellfun(@(x)x.prw.ProfileY,LI.sol,'uni',0);C1.Y(:,2,:)=horzcat(r2{:});
   r2=cellfun(@(x)x.prw.ProfileZ,LI.sol,'uni',0);C1.Y(:,3,:)=horzcat(r2{:});
   C1.X{1}=(1:size(C1.Y,1))';
   out.prw=C1;

   sol=sdtm.rmfield(sol,'x','y','prr','prw');
   out.sol=sol;

  end


  function [ ] = show_profiles(sol, opt)

   % [ ] = show_profiles(sol, opt)
   %
   % create 3d surf plot of rail and/or wheel profiles
   %   sol     - structure with CONTACT results as returned by loadcase
   %   opt     - structure with options as defined by plot3d.
   %

   % set the x range for plotting
   % #show_profiles

   if (isempty(opt.xrange))
    if (any(strcmp(opt.rw_surfc,{'prr','both'})))
     xmax = max( abs(sol.meta.xcp_r) + 0.40 * sol.mx * sol.dx, 1.20 * sol.mx * sol.dx);
     if (isfield(sol,'slcs')), xmax = 3 * xmax; end
     if (isfield(sol,'slcw')), xmax = sol.meta.rnom_whl; end
    else
     xmax = max( abs(sol.meta.xcp_w) + 0.40 * sol.mx * sol.dx, 1.20 * sol.mx * sol.dx);
    end
    opt.xrange = [-xmax, xmax ];
   end

   % set target step-sizes for plotting

   if (isempty(opt.xysteps))
    xlen = opt.xrange(2) - opt.xrange(1);  ylen = 0;
    if (any(strcmp(opt.rw_surfc,{'prr','both'})))
     ylen = max(ylen, max(sol.prr.ProfileY) - min(sol.prr.ProfileY));
    end
    if (any(strcmp(opt.rw_surfc,{'prw','both'})))
     ylen = max(ylen, max(sol.prw.ProfileY) - min(sol.prw.ProfileY));
    end
    opt.xysteps = max(xlen/4, ylen/15);
   end
   if (isscalar(opt.xysteps)); opt.xysteps = [1 1] * opt.xysteps;end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % if requested, show the rail profile

   want_rail = 1;
   if (any(strcmp(opt.rw_surfc,{'prr','both'})) & strcmp(opt.typplot,'rw_rear'))

    % 2d rear view of y-z-plane

    % form profile curve at x_cp

    xval = sol.meta.xcp_r;
    [ xcurv, ycurv, zcurv ] = cntc.make_3d_surface( sol, opt, want_rail, [xval, xval], [], 1, ...
     [], [], []);

    plot(ycurv, zcurv, 'color',cntc.matlab_color(1));
    set(gca,'ydir','reverse');
    axis equal;
    hold on

   elseif (any(strcmp(opt.rw_surfc,{'prr','both'})) & strcmp(opt.typplot,'rw_side'))

    % 2d side view of x-z-plane

    % s-value with y(ip) closest to y_cp

    [~,ip] = min( abs(sol.prr.ProfileY - sol.meta.ycp_r) );
    sval = sol.prr.ProfileS(ip);

    % form profile curve at given x, s

    totlen = opt.xrange(2) - opt.xrange(1);
    [ xcurv, ycurv, zcurv ] = cntc.make_3d_surface( sol, opt, want_rail, opt.xrange, [], opt.xysteps(1)/10, ...
     [sval, sval], [], 1);

    plot(xcurv, zcurv, 'color',cntc.matlab_color(1));
    set(gca,'ydir','reverse');
    axis equal;
    hold on

   elseif (any(strcmp(opt.rw_surfc,{'prr','both'})))

    % 3d view of rail surface

    % form profile surface at given x, all s - fine sampling

    % disp('...full surface')
    [ xsurf, ysurf, zsurf ] = cntc.make_3d_surface( sol, opt, want_rail, opt.xrange, [], opt.xysteps(1)/10, ...
     [], [], []);

    % plot 90% transparent surface (alpha=0.1); color-value?

    csurf = -1 * ones(size(xsurf));
    l=surf(xsurf, ysurf, zsurf, csurf, 'EdgeAlpha',0.1, 'FaceAlpha',0.1);
    set(l,'tag','prr');

    % set view appropriate for rails with z positive downwards

    set(gca,'xdir','reverse', 'zdir','reverse');
    axis equal;
    hold on

    % determine transverse curves along the surface at fixed steps in x

    % disp('...transverse curves')
    [ xcurv, ycurv, zcurv ] = cntc.make_3d_surface( sol, opt, want_rail, opt.xrange, [], opt.xysteps(1), ...
     [], [], []);

    % plot transverse curves along the surface

    l = plot3(xcurv', ycurv', zcurv', 'color',[.5 .5 .5], 'linewidth',0.5);

    [~,imin] = min(abs(xcurv(:,1)));
    set(l(imin), 'color', [.3 .3 .3], 'linewidth',1);

    % determine longitudinal curves along the surface at fixed steps in s

    % disp('...longitudinal curves')
    [ xcurv, ycurv, zcurv ] = cntc.make_3d_surface( sol, opt, want_rail, opt.xrange, [], opt.xysteps(1)/10, ...
     [], [], opt.xysteps(2));

    % plot curves along surface in longitudinal direction

    clear l ip;
    l = plot3(xcurv, ycurv, zcurv, 'color',[.5 .5 .5], 'linewidth',0.5);
    imid = ceil(length(l)/2);
    set(l(imid), 'color', [.3 .3 .3], 'linewidth',1);

   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % if requested, show the wheel profile

   want_rail = 0;
   if (any(strcmp(opt.rw_surfc,{'prw','both'})) & strcmp(opt.typplot,'rw_rear'))

    % 2d rear view of y-z-plane
    if (sol.meta.npatch>1)
     xval = 0;      % use principal profile; xcp_w can be much different between patches
    else
     xval = sol.meta.xcp_w;
    end
    [ xcurv, ycurv, zcurv ] = cntc.make_3d_surface( sol, opt, want_rail, [xval, xval], [], 1, ...
     [], [], []);

    plot(ycurv, zcurv, 'color',cntc.matlab_color(3));
    set(gca,'ydir','reverse');
    axis equal;
    hold on

   elseif (any(strcmp(opt.rw_surfc,{'prw','both'})) & strcmp(opt.typplot,'rw_side'))

    % 2d side view of x-z-plane

    % s-value with y(ip) closest to y_cp

    [~,ip] = min( abs(sol.prw.ProfileY - sol.meta.ycp_w) );
    sval = sol.prw.ProfileS(ip);

    % form wheel profile curve at given x, s

    [ xcurv, ycurv, zcurv ] = cntc.make_3d_surface( sol, opt, want_rail, opt.xrange, [], opt.xysteps(1)/10, ...
     [sval, sval], [], 1);

    plot(xcurv, zcurv, 'color',cntc.matlab_color(3));
    set(gca,'ydir','reverse');
    axis equal;
    hold on

   elseif (any(strcmp(opt.rw_surfc,{'prw','both'})))

    % 3d view of wheel surface

    % form profile surface at given x, all s - fine sampling

    [xsurf,ysurf,zsurf]= cntc.make_3d_surface(sol,opt,want_rail,opt.xrange,[],opt.xysteps(1)/10, ...
     [], [], []);

    % plot 90% transparent surface (alpha=0.1); color-value?

    csurf = -1 * ones(size(xsurf));
    l=surf(xsurf, ysurf, zsurf, csurf, 'EdgeAlpha',0.1, 'FaceAlpha',0.1);
    set(l,'tag','prw');

    % set view appropriate for rails with z positive downwards

    set(gca,'xdir','reverse', 'zdir','reverse');
    axis equal;
    hold on

    % determine transverse curves along the surface at fixed steps in x

    [ xcurv, ycurv, zcurv ] = cntc.make_3d_surface( sol, opt, want_rail, opt.xrange, [], opt.xysteps(1), ...
     [], [], []);

    % plot transverse curves along the surface

    l = plot3(xcurv', ycurv', zcurv', 'color',[.5 .5 .5], 'linewidth',0.5);

    [~,imin] = min(abs(xcurv(:,1)));
    set(l(imin), 'color', [.3 .3 .3], 'linewidth',1);

    % determine curves along the surface at fixed steps in s

    [ xcurv, ycurv, zcurv ] = cntc.make_3d_surface( sol, opt, want_rail, opt.xrange, [], opt.xysteps(1)/10, ...
     [], [], opt.xysteps(2));

    % plot curves along surface in longitudinal direction

    clear l ip;
    l = plot3(xcurv, ycurv, zcurv, 'color',[.5 .5 .5], 'linewidth',0.5);
    imid = ceil(length(l)/2);
    set(l(imid), 'color', [.3 .3 .3], 'linewidth',1);

   end

  end % show_profiles

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ xsurf, ysurf, zsurf ] = make_3d_surface( sol, opt, want_rail, xrange, dx_true, dx_appx,srange, ds_true, ds_appx, idebug)

   % create 2D arrays (x,y,z) for a prismatic rail (prr), variable rail (slcs), roller (prr, nom_radius),
   %    round wheel (prw, nom_radius), or out-of-round wheel (slcw)
   %  - the profile will be plotted for xrange = [xmin, xmax] in track coordinates
   %  - fixed steps of size dx_true may be requested, or n-time fitting steps of approximately dx_appx
   %  - the whole profile will be used if srange = [smin, smax] is left empty
   %  - the profile will be sampled at sj-values in the profile close to uniform sj = [smin: ds: smax]
   %  - all profile points will be used if both ds_true and ds_appx are left empty
   %  - for variable profiles, x-positions and arc-lengths s are replaced by the (u,v)-parametrization
   % #make_3d_surface -2

   if (nargin< 5), dx_true = []; end
   if (nargin< 6), dx_appx = []; end
   if (nargin< 7), srange  = []; end
   if (nargin< 8), ds_true = []; end
   if (nargin< 9), ds_appx = []; end
   if (nargin<10 | isempty(idebug)), idebug  =  0; end

   if (want_rail)
    prf = sol.prr;
    nom_radius = sol.meta.rnom_rol;
   else
    prf = sol.prw;
    nom_radius = sol.meta.rnom_whl;
   end
   has_slcs   = ( want_rail & isfield(sol,'slcs'));
   has_slcw   = (~want_rail & isfield(sol,'slcw'));

   %% determine longitudinal positions xi (track coords) for evaluation of surface

   if (~isempty(dx_true))

    % if 'final' dx_true is prescribed,
    %    set xi = { i * dx_true } for appropriate i

    i0 = ceil( xrange(1) / dx_true );
    i1 = floor( xrange(2) / dx_true );
    xi = [i0 : i1] * dx_true;

   elseif (~isempty(dx_appx))

    % if 'target' dx_appx is given,
    %    divide [xmin,xmax] in odd #intervals of size dx_true close to dx_appx
    nintv   = floor( (xrange(2) - xrange(1)) / dx_appx );
    if (mod(nintv,2)==1)  % prefer odd #lines for symmetric [-xmax, xmax]
     nintv = nintv + 1;
    end
    dx = (xrange(2) - xrange(1)) / max(1,nintv);
    xi = xrange(1) + [0:nintv] * dx;
   else
    disp('Error: either dx_true or dx_appx must be given');return
   end

   % determine lateral positions sj for evaluation of surface

   if (has_slcs);    profile_s = sol.slcs.vj;
   elseif (has_slcw);profile_s = sol.slcw.vj;
   else ;            profile_s = prf.ProfileS;
   end

   if (isempty(srange))
    srange = [ profile_s(1), profile_s(end) ];
   end

   if (~isempty(ds_true))

    % if 'final' ds_true is prescribed,
    %    set si = { j * ds_true } for appropriate j
    j0 = ceil( (srange(1)-profile_s(1)) / ds_true );
    j1 = floor( (srange(2)-profile_s(1)) / ds_true );
    sj = [j0 : j1] * ds_true;

   elseif (~isempty(ds_appx))

    % if 'target' ds_appx is given,
    %    divide [smin,smax] in odd #intervals of size ds_true close to ds_appx

    if (has_slcs)
     len = 1.3 * max(max(sol.slcs.ysurf)) - min(min(sol.slcs.ysurf));
     nintv = floor ( len / ds_appx );
    elseif (has_slcw)
     len = 1.3 * max(max(sol.slcw.ysurf)) - min(min(sol.slcw.ysurf));
     nintv = floor ( len / ds_appx );
    else
     nintv = floor( (srange(2) - srange(1)) / ds_appx );
    end
    if (mod(nintv,2)==1)
     nintv = nintv + 1;  % prefer odd #lines for symmetric [-smax, smax]
    end
    ds = (srange(2) - srange(1)) / max(1,nintv);
    sj = srange(1) + [0:nintv] * ds;

   else
    % if neither is given, use all profile points in range srange
    sj = profile_s(profile_s>=srange(1) & profile_s<=srange(2) );
   end

   % determine indices js closest to each sj

   js = zeros(size(sj));
   for j = 1 : length(sj)
    [~,js(j)] = min(abs(profile_s - sj(j)));
   end

   % form profile surface at given xi and js

   nx    = length(xi);
   ns    = length(js);

   if (want_rail & ~has_slcs & abs(nom_radius)<1e-3)

    % prismatic rail

    xsurf = xi' * ones(1,ns);
    ysurf = ones(nx,1) * prf.ProfileY(js)';
    zsurf = ones(nx,1) * prf.ProfileZ(js)';

   elseif (want_rail & cntc.myisfield(sol, 'slcs.spl2d'))

    % variable rail profile, using 2D spline algorithm

    % form profile surface at given ui (==x_fc) and vj

    sf_min = min(sol.slcs.spl2d.ui);
    sf_max = max(sol.slcs.spl2d.ui);
    ui     = max(sf_min, min(sf_max, sol.meta.s_ws + xi));
    vj     = sj;

    [ ~, xsurf, ysurf, zsurf ] = cntc.eval_2dspline( sol.slcs.spl2d, ui, vj );

   elseif (want_rail & has_slcs)

    xsurf = xi' * ones(1,ns);
    ysurf = []; zsurf = [];
    for ix = 1 : nx
     x_cur = sol.meta.s_ws + xi(ix);
     tmp_prr = cntc.get_profile_slice(sol.slcs, x_cur);
     ysurf = [ysurf; tmp_prr.ProfileY(js)'];
     zsurf = [zsurf; tmp_prr.ProfileZ(js)'];
    end

   elseif (want_rail)

    % roller surface: circle x^2 + (rnom - z)^2 = (rnom - z(0,y))^2

    xsurf = xi' * ones(1,ns);
    ysurf = ones(nx,1) * prf.ProfileY(js)';
    r_y   = nom_radius - prf.ProfileZ(js)';
    zsurf = nom_radius - sqrt( max(0, (ones(nx,1)*r_y).^2 - xsurf.^2) );

   elseif (~want_rail & ~has_slcw)

    % round wheel surface: circle x^2 + (rnom + z)^2 = (rnom + z(0,y))^2
    % #Prw.Rev origin -2
    xsurf = xi' * ones(1,ns);
    ysurf = ones(nx,1) * prf.ProfileY(js)';
    r_y   =  nom_radius + prf.ProfileZ(js)';
    zsurf = -nom_radius + sqrt( max(0, (ones(nx,1)*r_y).^2 - xsurf.^2) );

   elseif (~want_rail & cntc.myisfield(sol, 'slcw.spl2d'))

    % variable wheel profile, using 2D spline algorithm

    % form profile surface at given ui (==th_w) and vj

    th_min = min(sol.slcw.spl2d.ui);
    th_max = max(sol.slcw.spl2d.ui);
    th_w   = -sol.meta.th_ws + mean(xi)/sol.meta.rnom_whl;
    ui     = max(th_min, min(th_max, -sol.meta.th_ws + xi/sol.meta.rnom_whl));
    vj     = sj;

    [ ~, thsurf, ysurf, drsurf ] = cntc.eval_2dspline( sol.slcw.spl2d, ui, vj );
    [xsurf, ysurf, zsurf] = cntc.wcyl_to_wprof_coords(sol, thsurf, ysurf, drsurf);

   else

    disp('make_3d_surface: ERROR: case not supported.')
    disp([want_rail, has_slcs, has_slcw, cntc.myisfield(sol,'slcw.spl2d'), exist('eval_2dspline')])
    return;

   end

   % profiles are given for right-side w/r combination
   % mirror y for left-side w/r combination

   if (cntc.is_left_side(sol))
    ysurf = -ysurf;
   end

   % transform to track coordinates if needed

   if (want_rail & strcmp(opt.rw_surfc,'both'))

    %xyzsurf=cat(3,xsurf,ysur,zsurf)
    [xsurf, ysurf, zsurf] = cntc.rail_to_track_coords(sol, xsurf, ysurf, zsurf);

   elseif (strcmp(opt.rw_surfc,'both'))

    [xsurf, ysurf, zsurf] = cntc.wheel_to_track_coords(sol, xsurf, ysurf, zsurf);

   end

  end % make_3d_surface

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ ] = show_scalar_field(sol, field, opt)
   %
   % [ ] = show_scalar_field(sol, field, opt)
   %
   % generic routine for plotting a grid with cell-centered values
   %   sol     - structure with CONTACT results as returned by loadcase
   %   field   - the actual array to be plotted, or a field name within sol
   %   opt     - structure with options as defined by plot3d.
   %

   % delegate work to separate functions depending on type of plot
   % #show_scalar_field

   if (strcmp(opt.typplot, 'rw_rear'))
    cntc.show_rw_rear_view(sol, field, opt);
    if (isfield(sol, 'subs'))
     cntc.show_subs_field(sol, opt);
    end
   elseif (strcmp(opt.typplot, 'rw_side'))
    cntc.show_rw_side_view(sol, field, opt);
   else
    cntc.show_3d_field(sol, field, opt);
   end

  end % show_scalar_field

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ ] = show_subs_field(sol, opt)
   %
   % [ ] = show_subs_field(sol, field, opt)
   %
   % add plot of subsurface stresses to rw_rear plot
   % #show_subs_field

   if (0==1)
    % use slice closest to xc=0
    [~,ix] = min(abs(sol.subs.x));
    fld = squeeze( sol.subs.sigvm(ix,:,:) );
   else
    fld = squeeze( max(abs(sol.subs.sigvm), [], 1) );
   end

   % plot the potential contact area

   [YC, ZC] = meshgrid(sol.subs.y, sol.subs.z);
   if (strcmp(opt.rw_surfc,'prr'))
    [ ~, ypln, zpln ] = cntc.to_rail_coords(sol, [], YC, ZC);
   elseif (strcmp(opt.rw_surfc,'prw'))
    [ ~, ypln, zpln ] = cntc.to_wheel_coords(sol, [], YC, ZC);
   elseif (strcmp(opt.rw_surfc,'both'))
    [ ~, ypln, zpln ] = cntc.to_track_coords(sol, [], YC, ZC);
   end

   [C, h] = contourf( ypln, zpln, fld' );
   set(h, 'linewidth', 0.01);

  end % show_subs_field

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ ] = show_3d_field(sol, field, opt)

   % determine version information on the Matlab graphics system
   % #show_3d_field

   % get field to be plotted

   if (ischar(field)); eval(['tmp = sol.',field,';']);
   else; tmp = field;
   end

   % set element numbers ix, iy to be plotted: interval of [1:mx(my)]

   [ix_plot, iy_plot] = cntc.plot_ranges(sol, opt);

   % select the range to be plotted

   tmp    = tmp(iy_plot,ix_plot);
   eldiv  = sol.eldiv(iy_plot,ix_plot);
   xcornr = sol.xl + [ix_plot(1)-1, ix_plot]*sol.dx;
   ycornr = sol.yl + [iy_plot(1)-1, iy_plot]*sol.dy;
   if (~isempty(sol.x_offset)), xcornr = xcornr + sol.x_offset; end
   if (~isempty(sol.y_offset)), ycornr = ycornr + sol.y_offset; end
   xcornr = opt.xstretch * xcornr;

   % fill exterior area with default values

   if (~isempty(opt.exterval))
    ii = find(eldiv==0);
    tmp(ii) = opt.exterval;
   end

   if (isempty(opt.exterval) | ~isnan(opt.exterval))

    % expand one row/column for shading flat for cell-centered values

    tmp = [tmp, tmp(:,end)];
    tmp = [tmp; tmp(end,:)];

   elseif (exist('repelem','builtin'))

    % with flat shading, Matlab doesn't show faces where one of the corners
    % has value NaN --> duplicate values, shift-duplicate corner coordinates

    tmp = repelem(tmp,2,2);
    xcornr = [xcornr(1), repelem(xcornr(2:end-1),2), xcornr(end)];
    ycornr = [ycornr(1), repelem(ycornr(2:end-1),2), ycornr(end)];

   else

    % fall-back for Matlab versions < R2015a

    % expand one row/column for shading flat for cell-centered values

    tmp = [tmp, tmp(:,end)];
    tmp = [tmp; tmp(end,:)];

    %  perform dilatation:
    %  at exterior elements that lie to the right or above of interior
    %  elements, the tmp-value must be non-NaN

    sz=size(tmp);
    % rows 1:end-1 ==> rows 2:end
    [i,j]=find(eldiv(1:end-1,:)~=0 & eldiv(2:end,:)==0);
    ii0=sub2ind(sz,i,j); ii1=sub2ind(sz,i+1,j  ); tmp(ii1)=tmp(ii0);

    % row end ==> row end+1
    [i,j]=find(eldiv(end,:)~=0); i = i + size(eldiv,1) - 1;
    ii0=sub2ind(sz,i,j); ii1=sub2ind(sz,i+1,j  ); tmp(ii1)=tmp(ii0);
    ii1=sub2ind(sz,i+1,j+1); tmp(ii1)=tmp(ii0);
    % columns 1:end-1 ==> columns 2:end
    [i,j]=find(eldiv(:,1:end-1)~=0 & eldiv(:,2:end)==0);
    ii0=sub2ind(sz,i,j); ii1=sub2ind(sz,i  ,j+1); tmp(ii1)=tmp(ii0);

    % column end ==> column end+1
    [i,j]=find(eldiv(:,end)~=0); j = j + size(eldiv,2) - 1;
    ii0=sub2ind(sz,i,j); ii1=sub2ind(sz,i+1,j  ); tmp(ii1)=tmp(ii0);
    ii1=sub2ind(sz,i+1,j+1); tmp(ii1)=tmp(ii0);

    % diagonally iend,jend -> iend+1,jend+1
    [i,j]=find(eldiv(1:end-1,1:end-1)~=0 & eldiv(2:end,2:end)==0);
    ii0=sub2ind(sz,i,j); ii1=sub2ind(sz,i+1,j+1); tmp(ii1)=tmp(ii0);
   end

   % mark lower-left corner

   if (~strcmp(field,'ptarg'))
    minval   = min(min(tmp));
    tmp(1,1) = -0.2*max(max(abs(tmp)));
   end

   % expand x,y-coordinates to 2d arrays

   nx = length(xcornr);
   ny = length(ycornr);
   xcornr = ones(ny,1) * xcornr;
   ycornr = ycornr' * ones(1,nx);

   % set z-coordinates for plotting
   %% #XXXgae Position -2
   if (strcmp(opt.rw_surfc,'none'))
    zcornr = tmp;
   elseif (strcmp(opt.rw_surfc,'prr'))
    [ xcornr, ycornr, zcornr ] =cntc.to_rail_coords(sol, xcornr, ycornr);
    ii = find(isnan(tmp)); zcornr(ii) = NaN;
   elseif (strcmp(opt.rw_surfc,'prw'))
    [ xcornr, ycornr, zcornr ] =cntc.to_wheel_coords(sol, xcornr, ycornr);
    ii = find(isnan(tmp)); zcornr(ii) = NaN;
   elseif (strcmp(opt.rw_surfc,'both'))
    [ xcornr, ycornr, zcornr ] =cntc.to_track_coords(sol, xcornr, ycornr);
    ii = find(isnan(tmp)); zcornr(ii) = NaN;
   end

   % make plot, adjust colormap, z-range
   prop={'tag',opt.field};
   if (strcmp(opt.colormap,'none') | strcmp(opt.colormap,'black'))
    mesh( xcornr, ycornr, zcornr, tmp ,prop{:});
    view(opt.view);
    colormap([0 0 0]);
   else
    if (strcmp(opt.typplot,'contourf'))
     contourf( xcornr, ycornr, tmp, 'ShowText','on'  ,prop{:});
    else
     % surf( xcornr, ycornr, tmp );
     surf( xcornr, ycornr, zcornr, tmp  ,prop{:});
    end
    view(opt.view);
    if (strcmp(opt.rw_surfc,'none'))
     shading faceted;
    else
     shading flat;
    end
    colormap(opt.colormap);
    if (~isempty(opt.zrange))
     % use prescribed clim
     set(gca,'clim',opt.zrange);
    else
     % avoid spoiling of clim by special value at point (1,1):
     cl=get(gca,'clim');
     if (minval<cl(2))      % ensure that clim is increasing
      cl(1)=max(minval,cl(1));
      set(gca,'clim',cl);
     end
    end
   end

   % adjust axis

   nw_ax = [ min(min(xcornr)), max(max(xcornr)), min(min(ycornr)), max(max(ycornr)) ];
   if (opt.addplot==1)  % allow growing of plot region only
    pv_ax = axis;
    nw_ax = [ min(pv_ax(1),nw_ax(1)), max(pv_ax(2),nw_ax(2)) ...
     min(pv_ax(3),nw_ax(3)), max(pv_ax(4),nw_ax(4)) ];
   end

   if (strcmp(opt.rw_surfc,'none'))
    axis( nw_ax );
   end

   % optionally override automatic ticks, set ticks at cell-centers

   if (0==1)
    ixstep = max(1, round( (ix_plot(end)-ix_plot(1)+1) / 5 ));
    iystep = max(1, round( (iy_plot(end)-iy_plot(1)+1) / 5 ));
    x_tick = sol.xl + (0.5+[ix_plot(1):ixstep:ix_plot(end)])*sol.dx;
    y_tick = sol.yl + (0.5+[iy_plot(1):iystep:iy_plot(end)])*sol.dy;
    set(gca,'xtick', x_tick * opt.xstretch, ...
     'xticklabel', num2str(x_tick', '%5.3f'));
    set(gca,'ytick', y_tick, 'yticklabel', num2str(y_tick', '%5.3f'));
   elseif (opt.xstretch~=1)
    %     set tick-labels to unstretched values
    xt = get(gca,'xtick');
    xv = xt / opt.xstretch;
    %     #digits before decimal point, negative: zeros after decimal point -3
    lg = ceil(log10(max(abs(xv))));
    if (lg>=2)
     fmt = '%2.0f';
    elseif (lg==1 | lg==0)
     fmt = '%3.1f';
    elseif (lg==-1)
     fmt = '%4.2f';
    else
     fmt = '%5.3f';
    end
    set(gca, 'xtick', xt', 'xticklabel', num2str(xv', fmt));
   end

   % set grid, labels

   grid on
   cntc.set_xyzlabels(sol, opt);

   % draw colorbar, set label, if not done so before
   % note: no colorbar on mesh plot

   if (opt.addplot<=0 & ...
     ~(strcmp(opt.colormap,'none') | strcmp(opt.colormap,'black')))

    h = colorbar;
    if (strcmp(opt.field,'pn'))
     ylabel(h, 'Pressure [N/mm^2]');
    elseif (strcmp(opt.field,'ptabs') | ...
      strcmp(opt.field,'ptabs+vec') | ...
      strcmp(opt.field,'px') | ...
      strcmp(opt.field,'py'))
     ylabel(h, 'Traction [N/mm^2]');
    elseif (strcmp(opt.field,'un') | strcmp(opt.field,'ux') | ...
      strcmp(opt.field,'uy') | strcmp(opt.field,'h'))
     ylabel(h, 'Displacement [mm]');
    elseif (strcmp(opt.field,'ptarg'))
     ylabel(h, 'Direction [deg]');
    elseif (strcmp(opt.field,'shft') | (sol.kincns.t_digit<=1 & ...
      (strcmp(opt.field,'sx')  | strcmp(opt.field,'sy'))))
     ylabel(h, 'Shift distance [mm]');
    elseif (strcmp(opt.field,'uplsx') | strcmp(opt.field,'uplsy') | strcmp(opt.field,'upls+vec'))
     ylabel(h, 'Plastic deformation [mm]');
    elseif (strcmp(opt.field,'taucrt'))
     ylabel(h, 'Tangential yield stress [N/mm^2]');
    elseif (strcmp(opt.field,'srel') | ...
      (sol.kincns.t_digit>=2 & (strcmp(opt.field,'sx')  | strcmp(opt.field,'sy'))) )
     ylabel(h, 'Relative slip velocity [-]');
    elseif (strcmp(opt.field,'sabs'))
     ylabel(h, 'Absolute slip velocity [mm/s]');
    elseif (strcmp(opt.field,'fricdens'))
     ylabel(h, 'Frictional power density [W/mm^2]');
    elseif (strcmp(opt.field,'temp1') | strcmp(opt.field,'temp2'))
     ylabel(h, 'Temperature (increase) [{}^\circ{}C]');
     % provisional support for user-added fields:
    elseif (strcmp(opt.field,'wx') | strcmp(opt.field,'wy'))
     ylabel(h, 'Rigid slip velocity [-]');
    elseif (strcmp(opt.field,'kyield'))
     ylabel(h, 'Yield strength K [N/mm^2]');
    elseif (strcmp(opt.field,'rcf'))
     ylabel(h, 'RCF index [-]');
    end
   end

  end % show_3d_field

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ ] = show_vec_field(sol, vx, vy, opt)
   %
   % [ ] = show_vec_field(sol, vx, vy, opt)
   %
   % generic routine for a quiver plot for cell-centered values
   %   sol     - structure with CONTACT results as returned by loadcase
   %   vx, vy  - the vector data to be plotted
   %   opt     - structure with options as defined by plot3d.
   %
   % #show_vec_field

   if (any(strcmp(opt.typplot,{'rw_rear','rw_side'})))
    return
   end

   % determine whether there are just vectors in the plot or magnitudes as well

   has_mgn = ~isempty(strfind(opt.field,'+vec'));

   % add to existing plot if magnitudes are already there

   if (has_mgn)
    hold on;
   end

   % set the vector color and line width, using defaults depending on
   % field to be plotted

   if (has_mgn)
    col='k';
    wid=2;
   else
    col='b';
    wid=1;
   end
   if (~isempty(opt.veccolor))
    col=opt.veccolor;
   end
   if (~isempty(opt.vecwidth))
    wid=opt.vecwidth;
   end

   % draw vectors in maximally about numvecx x numvecy points

   [ix_plot, iy_plot] = cntc.plot_ranges(sol, opt);
   ixstep = max(1, round( (ix_plot(end)-ix_plot(1)+1) / opt.numvecx ));
   iystep = max(1, round( (iy_plot(end)-iy_plot(1)+1) / opt.numvecy ));
   ix_vec = [ix_plot(1) : ixstep : ix_plot(end)];
   iy_vec = [iy_plot(1) : iystep : iy_plot(end)];
   x  = ones(size(iy_vec'))*sol.x(ix_vec);
   y  = sol.y(iy_vec)'*ones(size(ix_vec));
   if (~isempty(sol.x_offset)), x = x + sol.x_offset; end
   if (~isempty(sol.y_offset)), y = y + sol.y_offset; end
   x = x * opt.xstretch;

   % select components for actual vectors

   vx = vx(iy_vec,ix_vec);
   vy = vy(iy_vec,ix_vec);

   % eliminate vectors for points outside the contact area

   ix = find(sol.eldiv(iy_vec,ix_vec)==0); x(ix)=NaN;

   % set z-coordinates for plotting

   if (strcmp(opt.rw_surfc,'none'))
    zval = 1e6 * sign(opt.view(2)); z = zval * ones(size(x));
    vz = zeros(size(vx));
   elseif (strcmp(opt.rw_surfc,'prr'))
    [ x, y, z, deltar ] =cntc.to_rail_coords(sol, x, y);
    vz = sin(deltar) * vy;
    vy = cos(deltar) * vy;
   elseif (strcmp(opt.rw_surfc,'prw'))
    [ x, y, z, deltaw ] =cntc.to_wheel_coords(sol, x, y);
    vz = sin(deltaw) * vy;
    vy = cos(deltaw) * vy;
   elseif (strcmp(opt.rw_surfc,'both'))
    [ x, y, z, deltar ] =cntc.to_track_coords(sol, x, y);
    vz = sin(deltar) * vy;
    vy = cos(deltar) * vy;
   end

   % draw the actual vectors using quiver

   if (isempty(opt.vecscale) | opt.vecscale<=0)
    h = quiver3(x, y, z, vx, vy, vz, col);
   else
    % manual scaling of vectors, needed for addplot option:
    scl = opt.vecscale;
    h = quiver3(x, y, z, scl*vx, scl*vy, scl*vz, 0, col);
   end
   set(h,'linewidth',wid);

   % set huge z-data (+ or -) for the vectors when there are magnitudes as well

   % if (has_mgn)
   %    set(h,'zdata',sign(opt.view(2))*1e5*ones(size(get(h,'vdata'))));
   %    set(h,'wdata',    zeros(size(get(h,'vdata'))));
   % end

   % change axis range in 3D plots

   if (~cntc.is_2d_view(gca, opt.view) & strcmp(opt.rw_surfc,'none'))
    v=axis; axis([v -1e-6 1e-6]); axis equal;
    set(gca,'ztick',[]);
   end
   view(opt.view);

   % set labels

   cntc.set_xyzlabels(sol, opt);

  end % show_vec_field

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ ] = show_eldiv_contour(sol, opt)
   %
   % [ ] = show_eldiv_contour(sol, opt)
   %
   % routine for plotting outlines of contact area if opt.addeldiv>=1, adhesion/plasticity area (>=2)
   %   sol     - structure with CONTACT results as returned by loadcase
   %   opt     - structure with options as defined by plot3d.
   %
   % #show_eldiv_contour

   if (opt.addeldiv<=0)
    return;
   end
   if (sol.mx<=1 | sol.my<=1)
    if (0==1)
     disp('WARNING: option "eldiv_contour" is not available for 2D grids');
    end
    return
   end

   % set the line width, using default depending on field to be plotted
   if (~isempty(opt.eldivwid))
    wid = opt.eldivwid;
   elseif (strcmp(opt.field,'eldiv_contour'))
    wid = 1;
   else
    wid = 2;
   end

   [ix_plot, iy_plot] = cntc.plot_ranges(sol, opt);
   xcentr = sol.x(ix_plot); ycentr = sol.y(iy_plot);
   if (~isempty(sol.x_offset)), xcentr = xcentr + sol.x_offset; end
   if (~isempty(sol.y_offset)), ycentr = ycentr + sol.y_offset; end
   xcentr = xcentr * opt.xstretch;

   eldiv = cntc.derived_data(sol, opt, 'eldiv');
   eldiv = mod(eldiv,10); % diffcase: difference is stored in tens-digit
   eldiv = eldiv(iy_plot,ix_plot);

   % expand selection 'eldiv' and x/ycentr when there are elements in
   % contact in first/last rows/columns
   if (any(eldiv(1,:)))
    eldiv = [zeros(1,size(eldiv,2)) ; eldiv];
    ycentr = [ ycentr(1)-sol.dy, ycentr ];
   end
   if (any(eldiv(end,:)))
    eldiv = [eldiv ; zeros(1,size(eldiv,2))];
    ycentr = [ ycentr, ycentr(end)+sol.dy ];
   end
   if (any(eldiv(:,1)))
    eldiv = [zeros(size(eldiv,1),1) , eldiv];
    xcentr = [ xcentr(1)-sol.dx, xcentr ];
   end
   if (any(eldiv(:,end)))
    eldiv = [eldiv , zeros(size(eldiv,1),1)];
    xcentr = [ xcentr, xcentr(end)+sol.dx ];
   end

   zval = 1e5;
   l1 = []; l2 = []; l3 = [];

   % plot double line around adhesion area
   if (opt.addeldiv>=2)
    mask = 1.0*(eldiv==1);
    col  = opt.eldivcol(1,:);
    if (nnz(mask)>0 & ~strcmp(opt.field,'pn'))
     ContourLvl = [0.45, 0.55]';
     c = contourc(xcentr, ycentr, mask, ContourLvl);
     nx = size(c,2); ix = 1;
     while(ix < nx);
      npnt = c(2,ix); c(:,ix) = NaN; ix = ix + 1 + npnt;
     end

     % set xyz-coordinates for plotting

     if (strcmp(opt.rw_surfc,'none'))
      xcrv = c(1,:); ycrv = c(2,:); zcrv = zval*[1;-1]*ones(1,nx);
     elseif (strcmp(opt.rw_surfc,'prr'))
      [ xcrv, ycrv, zcrv ] =cntc.to_rail_coords(sol, c(1,:), c(2,:));
     elseif (strcmp(opt.rw_surfc,'prw'))
      [ xcrv, ycrv, zcrv ] =cntc.to_wheel_coords(sol, c(1,:), c(2,:));
     elseif (strcmp(opt.rw_surfc,'both'))
      [ xcrv, ycrv, zcrv ] =cntc.to_track_coords(sol, c(1,:), c(2,:));
     end

     l1 = plot3(xcrv, ycrv, zcrv, 'color',col, 'linewidth',wid, 'Tag','Adhes');
    end
   end

   % plot double line around plasticity area
   if (opt.addeldiv>=2)
    mask = 1.0*(eldiv==3);
    col  = opt.eldivcol(3,:);
    if (nnz(mask)>0 & ~strcmp(opt.field,'pn'))
     ContourLvl = [0.45, 0.55]';
     c = contourc(xcentr, ycentr, mask, ContourLvl);
     nx = size(c,2); ix = 1;
     while(ix < nx);
      npnt = c(2,ix); c(:,ix) = NaN; ix = ix + 1 + npnt;
     end

     % set xyz-coordinates for plotting

     if (strcmp(opt.rw_surfc,'none'))
      xcrv = c(1,:); ycrv = c(2,:); zcrv = zval*[1;-1]*ones(1,nx);
     elseif (strcmp(opt.rw_surfc,'prr'))
      [ xcrv, ycrv, zcrv ] =cntc.to_rail_coords(sol, c(1,:), c(2,:));
     elseif (strcmp(opt.rw_surfc,'prw'))
      [ xcrv, ycrv, zcrv ] =cntc.to_wheel_coords(sol, c(1,:), c(2,:));
     elseif (strcmp(opt.rw_surfc,'both'))
      [ xcrv, ycrv, zcrv ] =cntc.to_track_coords(sol, c(1,:), c(2,:));
     end

     l3 = plot3(xcrv, ycrv, zcrv, 'color',col, 'linewidth',wid, 'Tag','Plast');
    end
   end

   % plot single line around contact area
   if (opt.addeldiv>=1)
    mask = 1.0*(eldiv>=1);
    col = opt.eldivcol(2,:);
    if (nnz(mask)>0)
     ContourLvl = [0.5, 9.99]';
     c = contourc(xcentr, ycentr, mask, ContourLvl);
     nx = size(c,2); ix = 1;
     while(ix < nx);
      npnt = c(2,ix); c(:,ix) = NaN; ix = ix + 1 + npnt;
     end

     % set xyz-coordinates for plotting

     if (strcmp(opt.rw_surfc,'none'))
      xcrv = c(1,:); ycrv = c(2,:); zcrv = zval*[1;-1]*ones(1,nx);
     elseif (strcmp(opt.rw_surfc,'prr'))
      [ xcrv, ycrv, zcrv ] =cntc.to_rail_coords(sol, c(1,:), c(2,:));
     elseif (strcmp(opt.rw_surfc,'prw'))
      [ xcrv, ycrv, zcrv ] =cntc.to_wheel_coords(sol, c(1,:), c(2,:));
     elseif (strcmp(opt.rw_surfc,'both'))
      [ xcrv, ycrv, zcrv ] =cntc.to_track_coords(sol, c(1,:), c(2,:));
     end

     l2 = plot3(xcrv, ycrv, zcrv, 'color',col, 'linewidth',wid, 'Tag','Contact');
    end
   end

   % adapt axis in order to fit the plotted range of the contact area
   if (strcmp(opt.rw_surfc,'none'))
    xcornr = sol.xl + [ ix_plot(1)-1, ix_plot ]*sol.dx;
    ycornr = sol.yl + [ iy_plot(1)-1, iy_plot ]*sol.dy;
    if (~isempty(sol.x_offset)), xcornr = xcornr + sol.x_offset; end
    if (~isempty(sol.y_offset)), ycornr = ycornr + sol.y_offset; end
    nw_ax = [ min(xcornr), max(xcornr), min(ycornr), max(ycornr) ];
    if (opt.addplot==1)  % allow growing of plot region only
     pv_ax = axis;
     nw_ax = [ min(pv_ax(1),nw_ax(1)), max(pv_ax(2),nw_ax(2)) ...
      min(pv_ax(3),nw_ax(3)), max(pv_ax(4),nw_ax(4)) ];
    end
    axis( nw_ax );

    % in 3D plots, concentrate on plane Oxy using tiny z-range
    %    except when the eldiv is added to another 3D plot
    if (~cntc.is_2d_view(gca, opt.view) & strcmp(opt.field,'eldiv_contour'))
     v=axis; axis([v -1e-6 1e-6]); axis equal;
     set(gca,'ztick',[]);
    end
   end

   % set view and labels
   view(opt.view);
   cntc.set_xyzlabels(sol, opt);

  end % show_eldiv_contour

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function show_rw_rear_view(sol, field, opt)

   % show_rw_rear_view: plot w/r profiles, contact plane & contact results

   % plot the potential contact area
   % #show_rw_rear_view

   yc = sol.y; zc = zeros(1,sol.my);
   if (strcmp(opt.rw_surfc,'prr'))
    [ ~, ypln, zpln ] =cntc.to_rail_coords(sol, [], yc, zc);
   elseif (strcmp(opt.rw_surfc,'prw'))
    [ ~, ypln, zpln ] =cntc.to_wheel_coords(sol, [], yc, zc);
   elseif (strcmp(opt.rw_surfc,'both'))
    [ ~, ypln, zpln ] =cntc.to_track_coords(sol, [], yc, zc);
   end

   plot(ypln, zpln, '.-', 'markersize',12, 'color',cntc.matlab_color(5));

   % get field to be plotted, replace values in exterior elements
   if (ischar(field))
    eval(['val = sol.',field,';']);
   else
    val = field;
   end
   if (~isempty(opt.exterval))
    eldiv = cntc.derived_data(sol, opt, 'eldiv');
    adjac = eldiv;
    adjac(1:end-1,:) = max(adjac(1:end-1,:), eldiv(2:end  ,:));
    adjac(2:end  ,:) = max(adjac(2:end  ,:), eldiv(1:end-1,:));
    adjac(:,1:end-1) = max(adjac(:,1:end-1), eldiv(:,2:end  ));
    adjac(:,2:end  ) = max(adjac(:,2:end  ), eldiv(:,1:end-1));
    ii = find(adjac==0);
    val(ii) = opt.exterval;
   end

   % convert field to 2d view, aggregating values over columns
   % TODO: provide max, absmax, maxabs, sum
   val   = max(val, [], 2)';
   scale = cntc.get_vecscale( sol, opt, val );

   yc = sol.y; zc = -scale*val;
   if (strcmp(opt.rw_surfc,'prr'))
    [ ~, yval, zval ] =cntc.to_rail_coords(sol, [], yc, zc);
   elseif (strcmp(opt.rw_surfc,'prw'))
    [ ~, yval, zval ] =cntc.to_wheel_coords(sol, [], yc, zc);
   elseif (strcmp(opt.rw_surfc,'both'))
    [ ~, yval, zval ] =cntc.to_track_coords(sol, [], yc, zc);
   end

   col = cntc.matlab_color(6);
   if (~isempty(opt.veccolor))
    col = opt.veccolor;
   end
   plot( yval, zval, 'linewidth',1, 'color',col );
   plot( [ypln; yval], [zpln; zval], 'linewidth',1, 'color',col );

   % plot contact reference point, origin of contact local coordinate system

   axlen = 0.1 * (sol.y(end) - sol.y(1));
   yc = sol.meta.spinyo+[0,0,1]*axlen; zc = [1,0,0]*axlen;
   if (strcmp(opt.rw_surfc,'prr'))
    [ ~, yval, zval, delta ] =cntc.to_rail_coords(sol, [], yc, zc);
   elseif (strcmp(opt.rw_surfc,'prw'))
    [ ~, yval, zval, delta ] =cntc.to_wheel_coords(sol, [], yc, zc);
   elseif (strcmp(opt.rw_surfc,'both'))
    [ ~, yval, zval, delta ] =cntc.to_track_coords(sol, [], yc, zc);
   end

   if (exist('plot_cm','file'))
    l  = cntc.plot_cm([yval(2),zval(2)], 0.3, axlen, delta*180/pi, 'k');
    set(l, 'Tag','CM');
   else
    l1 = plot(yval, zval, 'k');
    l2 = plot(yval(2), zval(2), 'k.', 'markersize',18);
    set([l1,l2], 'Tag','CM');
   end

   cntc.set_yzlabels(sol, opt);
   grid on;

  end % show_rw_rear_view

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function show_rw_side_view(sol, field, opt)

   % show_rw_side_view: plot w/r profiles, contact plane & contact results

   % plot the potential contact area
   % #show_rw_side_view

   xc = sol.x; zc = zeros(1,sol.mx);
   if (strcmp(opt.rw_surfc,'prr'))
    [ xpln, ~, zpln ] =cntc.to_rail_coords(sol, xc, [], zc);
   elseif (strcmp(opt.rw_surfc,'prw'))
    [ xpln, ~, zpln ] =cntc.to_wheel_coords(sol, xc, [], zc);
   elseif (strcmp(opt.rw_surfc,'both'))
    [ xpln, ~, zpln ] =cntc.to_track_coords(sol, xc, [], zc);
   end

   plot(xpln, zpln, '.-', 'markersize',12, 'color',cntc.matlab_color(5));

   % get field to be plotted, replace values in exterior elements
   if (ischar(field))
    eval(['val = sol.',field,';']);
   else
    val = field;
   end
   if (~isempty(opt.exterval))
    eldiv = cntc.derived_data(sol, opt, 'eldiv');
    adjac = eldiv;
    adjac(1:end-1,:) = max(adjac(1:end-1,:), eldiv(2:end  ,:));
    adjac(2:end  ,:) = max(adjac(2:end  ,:), eldiv(1:end-1,:));
    adjac(:,1:end-1) = max(adjac(:,1:end-1), eldiv(:,2:end  ));
    adjac(:,2:end  ) = max(adjac(:,2:end  ), eldiv(:,1:end-1));
    ii = find(adjac==0);
    val(ii) = opt.exterval;
   end

   % convert field to 2d view, aggregating values over rows
   % TODO: provide max, absmax, maxabs, sum

   val   = max(val, [], 1);
   scale = cntc.get_vecscale( sol, opt, val );

   xc = sol.x; zc = -scale*val;
   if (strcmp(opt.rw_surfc,'prr'))
    [ xval, ~, zval ] =cntc.to_rail_coords(sol, xc, [], zc);
   elseif (strcmp(opt.rw_surfc,'prw'))
    [ xval, ~, zval ] =cntc.to_wheel_coords(sol, xc, [], zc);
   elseif (strcmp(opt.rw_surfc,'both'))
    [ xval, ~, zval ] =cntc.to_track_coords(sol, xc, [], zc);
   end

   col = cntc.matlab_color(6);
   if (~isempty(opt.veccolor))
    col = opt.veccolor;
   end

   plot( xval, zval, 'linewidth',1, 'color',col );
   plot( [xpln; xval], [zpln; zval], 'linewidth',1, 'color',col );

   % plot contact reference point, origin of contact local coordinate system

   v = axis;
   axlen = 0.03 * (v(2) - v(1));
   xc = sol.meta.spinxo+[0,0,1]*axlen; zc = [1,0,0]*axlen;
   if (strcmp(opt.rw_surfc,'prr'))
    [ xval, ~, zval ] =cntc.to_rail_coords(sol, xc, [], zc);
   elseif (strcmp(opt.rw_surfc,'prw'))
    [ xval, ~, zval ] =cntc.to_wheel_coords(sol, xc, [], zc);
   elseif (strcmp(opt.rw_surfc,'both'))
    [ xval, ~, zval ] =cntc.to_track_coords(sol, xc, [], zc);
   end

   if (exist('plot_cm','file'))
    roty = 0; % hack for 'grade' normal [nx,0,nz] instead of [0,0,1]
    if (isfield(sol.meta, 'roty')), roty = sol.meta.roty*180/pi; end
    l = cntc.plot_cm([xval(2),zval(2)], 0.15*axlen, axlen, -roty, 'k');
    set(l, 'Tag','CM');
   else
    l1 = plot(xval, zval, 'k');
    l2 = plot(xval(2), zval(2), 'k.', 'markersize',18);
    set([l1,l2], 'Tag','CM');
   end

   cntc.set_xyzlabels(sol, opt);
   grid on;

  end % show_rw_side_view

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ is2d ] = is_2d_view(ax, view);

   % return 1 if the current plot uses 2D view and 0 in case of a 3D view
   % 2D view: myopt.view = [ *, +/-90 ];  axis returns vector of 4 elements.
   % #is_2d_view -2

   is2d = (abs(90 - abs(view(2))) <= 1e-3);

  end % is_2d_view

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ ] = set_yzlabels(sol, opt)
   % #set_yzlabels -2

   if (strcmp(opt.rw_surfc,'prr'))
    xlabel ('y_r [mm]');
    ylabel ('z_r [mm]');
   elseif (strcmp(opt.rw_surfc,'prw'))
    xlabel ('y_w [mm]');
    ylabel ('z_w [mm]');
   elseif (strcmp(opt.rw_surfc,'both'))
    xlabel ('y_{tr} [mm]');
    ylabel ('z_{tr} [mm]');
   end

  end % set_yzlabels

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ ] = set_xyzlabels(sol, opt)
   % #set_xyzlabels -2

   if (strcmp(opt.typplot,'rw_rear'))

    if (strcmp(opt.rw_surfc,'prr'))
     xlabel ('y_r [mm]');
     ylabel ('z_r [mm]');
    elseif (strcmp(opt.rw_surfc,'prw'))
     xlabel ('y_w [mm]');
     ylabel ('z_w [mm]');
    elseif (strcmp(opt.rw_surfc,'both'))
     xlabel ('y_{tr} [mm]');
     ylabel ('z_{tr} [mm]');
    end

   elseif (strcmp(opt.typplot,'rw_side'))

    if (strcmp(opt.rw_surfc,'prr'))
     xlabel ('x_r [mm]');
     ylabel ('z_r [mm]');
    elseif (strcmp(opt.rw_surfc,'prw'))
     xlabel ('x_w [mm]');
     ylabel ('z_w [mm]');
    elseif (strcmp(opt.rw_surfc,'both'))
     xlabel ('x_{tr} [mm]');
     ylabel ('z_{tr} [mm]');
    end

   else

    if (strcmp(opt.rw_surfc,'prr'))
     xlabel ('x_r [mm]');
     ylabel ('y_r [mm]');
    elseif (strcmp(opt.rw_surfc,'prw'))
     xlabel ('x_w [mm]');
     ylabel ('y_w [mm]');
    elseif (strcmp(opt.rw_surfc,'both'))
     xlabel ('x_{tr} [mm]');
     ylabel ('y_{tr} [mm]');
    else
     xlabel ('x_c [mm]');
     if (sol.d_digit==4)
      ylabel ('s_c [mm]');
     else
      ylabel ('y_c [mm]');
     end
     %     else
     %        xlabel ('X-coordinate [mm]');
     %        ylabel ('Y-coordinate [mm]');
    end
   end

  end % set_xyzlabels

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ix_plot, iy_plot] = plot_ranges(sol, opt);

   % set element numbers ix, iy to be plotted: interval of [1:mx(my)]
   % #plot_ranges -2

   if (isempty(opt.ixrange) | ...
     (ischar(opt.ixrange) & strcmp(opt.ixrange,'all')) )
    ix_plot = [1, sol.mx];
   elseif (ischar(opt.ixrange) & strcmp(opt.ixrange,'tight'))
    num_inter = sum(sol.eldiv,1);
    ix_plot = find(num_inter);
    ix_plot = [ix_plot(1)-1 : ix_plot(end)+1 ];
   elseif (ischar(opt.ixrange) & strcmp(opt.ixrange,'auto'))
    num_inter = sum(sol.eldiv,1);
    ix_plot = find(num_inter);
    if (isempty(ix_plot))
     ix_plot = [1, sol.mx];
    else
     ix_plot = [ix_plot(1)-3 : ix_plot(end)+3 ];
    end
   else
    ix_plot = opt.ixrange;
   end
   if (isempty(opt.iyrange) | ...
     (ischar(opt.iyrange) & strcmp(opt.iyrange,'all')) )
    iy_plot = [1, sol.my];
   elseif (ischar(opt.iyrange) & strcmp(opt.iyrange,'tight'))
    num_inter = sum(sol.eldiv,2);
    iy_plot = find(num_inter);
    iy_plot = [iy_plot(1)-1 : iy_plot(end)+1 ];
   elseif (ischar(opt.iyrange) & strcmp(opt.iyrange,'auto'))
    num_inter = sum(sol.eldiv,2);
    iy_plot = find(num_inter);
    if (isempty(iy_plot))
     iy_plot = [1, sol.my];
    else
     iy_plot = [iy_plot(1)-3 : iy_plot(end)+3 ];
    end
   else
    iy_plot = opt.iyrange;
   end
   ix_plot = max(1,min(ix_plot)) : min(sol.mx,max(ix_plot));
   iy_plot = max(1,min(iy_plot)) : min(sol.my,max(iy_plot));

  end % plot_ranges

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ f1, f2, f3, f4, f5, f6 ] =derived_data(sol, opt, n1, n2, n3, n4, n5, n6)
   % #derived_data -2

   px = sol.px; py = sol.py; eldiv = sol.eldiv;

   nfield = nargin - 2;
   for ifield = 1 : nfield
    % get name of requested derived field, by copying n1, n2, .. or n6
    eval(sprintf('nam = n%d;',ifield));

    % compute the derived field
    if (strcmp(nam,'eldiv'))
     field = eldiv;
    elseif (strcmp(nam,'px'))
     field = px;
    elseif (strcmp(nam,'py'))
     field = py;
    elseif (strcmp(nam,'pt'))
     field = sqrt(px.^2 + py.^2);
    elseif (strcmp(nam,'ptarg'))
     field = 180/pi*atan2(py, px);   % range [-180,180]
    elseif (strcmp(nam,'pxdir'))
     field = sol.px ./ max(1d-9,sqrt(sol.px.^2 + sol.py.^2));
    elseif (strcmp(nam,'sabs'))
     field = sol.kincns.veloc * sol.srel;
    else
     disp(sprintf('Error: unknown derived field %s.',nam));
    end

    % store derived field, by copying to f1, f2, .. or f6
    eval(sprintf('f%d = field;',ifield));
   end

  end % derived_data

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ scale ] = get_vecscale( sol, opt, val )
   % #get_vecscale -2

   if (~isempty(opt.vecscale))
    scale = opt.vecscale;
   else
    scale = 10 / max(abs(val));
   end

  end % get_vecscale

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% #Coordinates
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [ xtr, ytr, ztr, delta_tr ] = to_track_coords(sol, xc, yc, zc);

   % transform contact [xc,yc/s]-coordinates to track [xtr,ytr,ztr] coordinates
   % #to_track_coords -2

   if (nargin<3), yc = []; end
   if (nargin<4), zc = []; end
   sz = max( [ size(xc) ; size(yc) ; size(zc) ] );
   if (isempty(xc)), xc = sol.meta.spinxo*ones(sz); end
   if (isempty(yc)), yc = sol.meta.spinyo*ones(sz); end
   if (isempty(zc)), zc = zeros(sz); end

   % reshape inputs to row vectors

   m = size(xc,1); n = size(xc,2);
   xc = reshape(xc, 1, m*n);
   yc = reshape(yc, 1, m*n);
   zc = reshape(zc, 1, m*n);
   coords = [xc - sol.meta.spinxo; yc - sol.meta.spinyo; zc];

   % form transformation matrices and vectors

   Rref = cntc.rotx(sol.meta.deltcp_r*180/pi);
   if (isfield(sol.meta, 'roty')) % hack for 'grade' normal [nx,0,nz] instead of [0,0,1]
    Rref = Rref * cntc.roty(sol.meta.roty*180/pi);
   end
   oref = [sol.meta.xcp_r; sol.meta.ycp_r; sol.meta.zcp_r];

   R_r  = cntc.rotx(sol.meta.roll_r*180/pi);
   o_r  = [0; sol.meta.y_r; sol.meta.z_r];

   if (sol.d_digit~=4) %if planar approach

    % change coordinates from contact to rail to track coordinates

    coords = o_r*ones(1,m*n) + R_r * oref*ones(1,m*n) + R_r * Rref * coords;

   else

    % mirror ProfileY for left-side w/r pairs

    prf_s = sol.prr.ProfileS; prf_y = sol.prr.ProfileY; prf_z = sol.prr.ProfileZ;
    if (cntc.is_left_side(sol))
     prf_s = -flipud(prf_s); prf_y = -flipud(prf_y); prf_z =  flipud(prf_z);
    end

    % find oref in the rail profile

    dst = sqrt( (prf_y-sol.meta.ycp_r).^2 + (prf_z-sol.meta.zcp_r).^2 );
    [~,ix] = min(dst);
    rg = [max(1,ix-5) : min(length(prf_y),ix+5)];

    sref = interp1(prf_y(rg), prf_s(rg), sol.meta.ycp_r);
    if (isnan(sref))
     disp(sprintf('ERROR: prr seems different from prr used in simulation (needs mirroring?).'))
     sref = prf_s(ix);
    end

    % find grid points in the rail profile

    sc = yc;
    if (any(sol.kincns.t_digit==[1 2]))
     sr =        sc;
    else
     sr = sref + sc;
    end
    yr_at_s = interp1(prf_s, prf_y, sr);
    zr_at_s = interp1(prf_s, prf_z, sr);

    % determine n-vectors on the rail profile

    ds = 0.1;
    dy = (interp1(prf_s, prf_y, sr+ds) - yr_at_s) / ds;
    dz = (interp1(prf_s, prf_z, sr+ds) - zr_at_s) / ds;
    nvec = [zeros(1,m*n); -dz; dy] ./ ([1;1;1] * sqrt(dy.^2+dz.^2));

    % hack: slight raise to put the data above the surfaces
    % zc = zc - 0.3;

    coords = [xc+oref(1); yr_at_s; zr_at_s] + nvec .* ([1;1;1] * zc);

    % change coordinates from rail to track coordinates

    coords = o_r*ones(1,m*n) + R_r * coords;

   end

   % extract components & reshape to the original array sizes

   xtr = reshape(coords(1,:), m, n);
   ytr = reshape(coords(2,:), m, n);
   ztr = reshape(coords(3,:), m, n);
   % TODO: correction for \Delta\phi_rr
   % delta_tr = sol.meta.deltcp_r + sol.meta.roll_r;
   delta_tr = sol.meta.deltcp_r + sol.meta.roll_r;

  end % cntc.to_track_coords

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ xr, yr, zr, deltar]=to_rail_coords(sol, xc, yc, zc);

   % transform contact [xc,yc/s]-coordinates to rail [xr,yr,zr] coordinates
   % #to_rail_coords -2

   if (nargin<3), yc = []; end
   if (nargin<4), zc = []; end
   sz = max( [ size(xc) ; size(yc) ; size(zc) ] );
   if (isempty(xc)), xc = sol.meta.spinxo*ones(sz); end
   if (isempty(yc)), yc = sol.meta.spinyo*ones(sz); end
   if (isempty(zc)), zc = zeros(sz); end

   % reshape inputs to row vectors

   m = size(xc,1); n = size(xc,2);
   xc = reshape(xc, 1, m*n);
   yc = reshape(yc, 1, m*n);
   zc = reshape(zc, 1, m*n);
   coords = [xc - sol.meta.spinxo; yc - sol.meta.spinyo; zc];

   % form transformation matrices and vectors

   Rref = cntc.rotx(sol.meta.deltcp_r*180/pi);
   if (isfield(sol.meta, 'roty')) % hack for 'grade' normal [nx,0,nz] instead of [0,0,1]
    Rref = Rref * cntc.roty(sol.meta.roty*180/pi);
   end

   % for variable rails, x_r-coordinates are defined by the slices-file rather than aligned with track x_tr

   has_slcs = isfield(sol,'slcs');
   if (has_slcs)
    oref = [sol.meta.s_ws+sol.meta.xcp_r; sol.meta.ycp_r; sol.meta.zcp_r];
   else
    oref = [              sol.meta.xcp_r; sol.meta.ycp_r; sol.meta.zcp_r];
   end

   if (sol.d_digit~=4)

    % change coordinates from contact to rail coordinates

    coords = oref*ones(1,m*n) + Rref * coords;

   else

    if (has_slcs), disp('Rail coords not supported for conformal on variable profile'); end

    % mirror ProfileY for left-side w/r pairs

    prf_s = sol.prr.ProfileS; prf_y = sol.prr.ProfileY; prf_z = sol.prr.ProfileZ;
    if (cntc.is_left_side(sol))
     prf_s = -flipud(prf_s); prf_y = -flipud(prf_y); prf_z =  flipud(prf_z);
    end

    % find oref in the rail profile

    dst = sqrt( (prf_y-sol.meta.ycp_r).^2 + (prf_z-sol.meta.zcp_r).^2 );
    [~,ix] = min(dst);
    rg = [max(1,ix-5), min(length(prf_y),ix+5)];

    sref = interp1(prf_y(rg), prf_s(rg), sol.meta.ycp_r);

    % find grid points in the rail profile

    sc = yc;
    if (any(sol.kincns.t_digit==[1 2]))
     sr =        sc;
    else
     sr = sref + sc;
    end
    yr_at_s = interp1(prf_s, prf_y, sr);
    zr_at_s = interp1(prf_s, prf_z, sr);

    % determine n-vectors on the rail profile

    ds = 0.1;
    dy = (interp1(prf_s, prf_y, sr+ds) - yr_at_s) / ds;
    dz = (interp1(prf_s, prf_z, sr+ds) - zr_at_s) / ds;
    nvec = [zeros(1,m*n); -dz; dy] ./ ([1;1;1] * sqrt(dy.^2+dz.^2));

    coords = [xc; yr_at_s; zr_at_s] + nvec .* ([1;1;1] * zc);

   end

   % extract components & reshape to the original array sizes

   xr = reshape(coords(1,:), m, n);
   yr = reshape(coords(2,:), m, n);
   zr = reshape(coords(3,:), m, n);
   deltar = sol.meta.deltcp_r;

  end % cntc.to_rail_coords

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ xr, yr, zr, deltar]=to_rail_coords_test(sol, A);

   % transform contact [xc,yc/s]-coordinates to rail [xr,yr,zr] coordinates
   % #to_rail_coords_test -2
   % xxxgae coord

   if (nargin<3), yc = []; end
   if (nargin<4), zc = []; end
   sz = max( [ size(xc) ; size(yc) ; size(zc) ] );
   if (isempty(xc)), xc = sol.meta.spinxo*ones(sz); end
   if (isempty(yc)), yc = sol.meta.spinyo*ones(sz); end
   if (isempty(zc)), zc = zeros(sz); end

   % reshape inputs to row vectors

   m = size(xc,1); n = size(xc,2);
   xc = reshape(xc, 1, m*n);
   yc = reshape(yc, 1, m*n);
   zc = reshape(zc, 1, m*n);
   coords = [xc - sol.meta.spinxo; yc - sol.meta.spinyo; zc];

   % form transformation matrices and vectors

   Rref = cntc.rotx(sol.meta.deltcp_r*180/pi);
   if (isfield(sol.meta, 'roty')) % hack for 'grade' normal [nx,0,nz] instead of [0,0,1]
    Rref = Rref * cntc.roty(sol.meta.roty*180/pi);
   end

   % for variable rails, x_r-coordinates are defined by the slices-file rather than aligned with track x_tr

   has_slcs = isfield(sol,'slcs');
   if (has_slcs)
    oref = [sol.meta.s_ws+sol.meta.xcp_r; sol.meta.ycp_r; sol.meta.zcp_r];
   else
    oref = [              sol.meta.xcp_r; sol.meta.ycp_r; sol.meta.zcp_r];
   end

   if (sol.d_digit~=4)

    % change coordinates from contact to rail coordinates

    coords = oref*ones(1,m*n) + Rref * coords;

   else

    if (has_slcs), disp('Rail coords not supported for conformal on variable profile'); end

    % mirror ProfileY for left-side w/r pairs

    prf_s = sol.prr.ProfileS; prf_y = sol.prr.ProfileY; prf_z = sol.prr.ProfileZ;
    if (cntc.is_left_side(sol))
     prf_s = -flipud(prf_s); prf_y = -flipud(prf_y); prf_z =  flipud(prf_z);
    end

    % find oref in the rail profile

    dst = sqrt( (prf_y-sol.meta.ycp_r).^2 + (prf_z-sol.meta.zcp_r).^2 );
    [~,ix] = min(dst);
    rg = [max(1,ix-5), min(length(prf_y),ix+5)];

    sref = interp1(prf_y(rg), prf_s(rg), sol.meta.ycp_r);

    % find grid points in the rail profile

    sc = yc;
    if (any(sol.kincns.t_digit==[1 2]))
     sr =        sc;
    else
     sr = sref + sc;
    end
    yr_at_s = interp1(prf_s, prf_y, sr);
    zr_at_s = interp1(prf_s, prf_z, sr);

    % determine n-vectors on the rail profile

    ds = 0.1;
    dy = (interp1(prf_s, prf_y, sr+ds) - yr_at_s) / ds;
    dz = (interp1(prf_s, prf_z, sr+ds) - zr_at_s) / ds;
    nvec = [zeros(1,m*n); -dz; dy] ./ ([1;1;1] * sqrt(dy.^2+dz.^2));

    coords = [xc; yr_at_s; zr_at_s] + nvec .* ([1;1;1] * zc);

   end

   % extract components & reshape to the original array sizes

   xr = reshape(coords(1,:), m, n);
   yr = reshape(coords(2,:), m, n);
   zr = reshape(coords(3,:), m, n);
   deltar = sol.meta.deltcp_r;

  end % cntc.to_rail_coords_test


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ xw, yw, zw, deltaw]=to_wheel_coords(sol, xc, yc, zc);

   % transform contact [xc,yc/s]-coordinates to wheel [xw,yw,zw] coordinates
   % #to_wheel_coords -2

   if (nargin<3), yc = []; end
   if (nargin<4), zc = []; end
   sz = max( [ size(xc) ; size(yc) ; size(zc) ] );
   if (isempty(xc)), xc = sol.meta.spinxo*ones(sz); end
   if (isempty(yc)), yc = sol.meta.spinyo*ones(sz); end
   if (isempty(zc)), zc = zeros(sz); end

   % reshape inputs to row vectors

   m = size(xc,1); n = size(xc,2);
   xc = reshape(xc, 1, m*n);
   yc = reshape(yc, 1, m*n);
   zc = reshape(zc, 1, m*n);
   coords = [xc - sol.meta.spinxo; yc - sol.meta.spinyo; zc];

   % form transformation matrices and vectors

   Rref = cntc.rotx(sol.meta.deltcp_r*180/pi);
   oref = [sol.meta.xcp_r; sol.meta.ycp_r; sol.meta.zcp_r];

   R_r  = cntc.rotx(sol.meta.roll_r*180/pi);
   o_r  = [0; sol.meta.y_r; sol.meta.z_r];

   R_w  = cntc.rotx(sol.meta.roll_w*180/pi) * cntc.rotz(sol.meta.yaw_w*180/pi);
   o_w  = [sol.meta.x_w; sol.meta.y_w; sol.meta.z_w];

   if (sol.d_digit~=4)

    % change coordinates from contact to rail to track to wheel coordinates

    xtr = o_r*ones(1,m*n) + R_r * oref*ones(1,m*n) + R_r * Rref * coords;
    coords = R_w' * ( xtr - o_w*ones(1,m*n) );

   else

    % mirror ProfileY for left-side w/r pairs

    prf_s = sol.prw.ProfileS; prf_y = sol.prw.ProfileY; prf_z = sol.prw.ProfileZ;
    if (cntc.is_left_side(sol))
     prf_s = -flipud(prf_s); prf_y = -flipud(prf_y); prf_z =  flipud(prf_z);
    end

    % change oref from rail to track to wheel coordinates

    wref = R_w' * ( o_r + R_r * oref - o_w );

    % find wref.y in the wheel profile

    dst = sqrt( (prf_y-wref(2)).^2 + (prf_z-wref(3)).^2 );
    [~,ix] = min(dst);
    rg = [max(1,ix-5), min(length(prf_y),ix+5)];

    sref = interp1(prf_y(rg), prf_s(rg), wref(2));

    % find grid points in the wheel profile

    sc = yc;
    if (any(sol.kincns.t_digit==[1 2]))
     sw =       -sc;
    else
     sw = sref - sc;    % sw increasing to the left
    end
    yw_at_s = interp1(prf_s, prf_y, sw);
    zw_at_s = interp1(prf_s, prf_z, sw);

    % determine n-vectors on the wheel profile

    ds = 0.1;
    dy = (interp1(prf_s, prf_y, sw-ds) - yw_at_s) / ds;
    dz = (interp1(prf_s, prf_z, sw-ds) - zw_at_s) / ds;
    nvec = [zeros(1,m*n); -dz; dy] ./ ([1;1;1] * sqrt(dy.^2+dz.^2));

    coords = [xc; yw_at_s; zw_at_s] + nvec .* ([1;1;1] * zc);
   end

   % extract coordinates, add wheel curvature

   fac_r = 1 / (2*sol.meta.rnom_whl);
   xw = coords(1,:);
   yw = coords(2,:);
   zw = coords(3,:) - fac_r * xw.^2;

   % reshape to the original array sizes

   xw = reshape(xw, m, n);
   yw = reshape(yw, m, n);
   zw = reshape(zw, m, n);

   deltaw = sol.meta.deltcp_r + sol.meta.roll_r - sol.meta.roll_w;

  end % cntc.to_wheel_coords

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ xtr, ytr, ztr ] = rail_to_track_coords(sol, xr, yr, zr);

   % transform rail [xr,yr,zr]-coordinates to track [xtr,ytr,ztr] coordinates
   % #rail_to_track_coords -2

   % form transformation matrices and vectors
   R_r  = cntc.rotx(sol.meta.roll_r*180/pi);
   % for variable rails, x_r-coordinates are defined by the slices-file rather than aligned with track x_tr

   if isfield(sol,'slcs');o_r  = [-sol.meta.s_ws; sol.meta.y_r; sol.meta.z_r];
   else;                  o_r  = [          0   ; sol.meta.y_r; sol.meta.z_r];
   end

   if nargin==1; xtr=struct('o',o_r,'R',R_r,'name','c_tr_r');end
   % reshape inputs to row vectors

   m = size(xr,1); n = size(xr,2);
   xr = reshape(xr, 1, m*n);
   yr = reshape(yr, 1, m*n);
   zr = reshape(zr, 1, m*n);
   coords = [xr; yr; zr];

   coords = o_r * ones(1,m*n) + R_r * coords;

   % extract & reshape to the original array sizes

   xtr = reshape(coords(1,:), m, n);
   ytr = reshape(coords(2,:), m, n);
   ztr = reshape(coords(3,:), m, n);

  end % rail_to_track_coords

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ xw, yw, zw ] = wcyl_to_wprof_coords(sol, thw, yw, drw);

   % transform cylindrical wheel [thw,yw,drw]-coordinates to wheel profile [xw,yw,zw] coordinates
   % #wcyl_to_wprof_coords -2

   % reshape inputs to row vectors

   m = size(thw,1); n = size(thw,2);
   thw = reshape(thw, 1, m*n);
   yw  = reshape(yw,  1, m*n);
   drw = reshape(drw, 1, m*n);

   % add nominal radius rnom to offset dr

   rw  = sol.meta.rnom_whl + drw;

   % convert cylindrical to cartesian wheel profile coordinates -- lowest point at -theta_ws

   xw  =                      rw .* sin(thw+sol.meta.th_ws);
   zw  = -sol.meta.rnom_whl + rw .* cos(thw+sol.meta.th_ws);

   % reshape to the original array sizes

   xw = reshape(xw, m, n);
   yw = reshape(yw, m, n);
   zw = reshape(zw, m, n);

  end % wcyl_to_wprof_coords

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ xtr, ytr, ztr ] = wheel_to_track_coords(sol, xw, yw, zw);

   % transform wheel [xw,yw,zw]-coordinates to track [xtr,ytr,ztr] coordinates
   % #wheel_to_track_coords -2

   % reshape inputs to row vectors

   m = size(xw,1); n = size(xw,2);
   xw = reshape(xw, 1, m*n);
   yw = reshape(yw, 1, m*n);
   zw = reshape(zw, 1, m*n);
   coords = [xw; yw; zw];

   % form transformation matrices and vectors from rw/ws to tr coordinates

   R_w_tr = cntc.rotx(sol.meta.roll_w*180/pi) * cntc.rotz(sol.meta.yaw_w*180/pi);
   o_w_tr = [sol.meta.x_w; sol.meta.y_w; sol.meta.z_w];

   % change coordinates from wheel to track coordinates

   coords = o_w_tr * ones(1,m*n) + R_w_tr * coords;

   % extract & reshape to the original array sizes

   xtr = reshape(coords(1,:), m, n);
   ytr = reshape(coords(2,:), m, n);
   ztr = reshape(coords(3,:), m, n);

  end % wheel_to_track_coords

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ xr, yr, zr ] = wheel_to_rail_coords(sol, xw, yw, zw);

   % transform wheel [xw,yw,zw]-coordinates to rail [xr,yr,zr] coordinates
   % #wheel_to_rail_coords -2

   % reshape inputs to row vectors

   m = size(xw,1); n = size(xw,2);
   xw = reshape(xw, 1, m*n);
   yw = reshape(yw, 1, m*n);
   zw = reshape(zw, 1, m*n);
   coords = [xw; yw; zw];

   % form transformation matrices and vectors

   R_r  = cntc.rotx(sol.meta.roll_r*180/pi);
   o_r  = [0; sol.meta.y_r; sol.meta.z_r];

   R_w  = cntc.rotx(sol.meta.roll_w*180/pi) * cntc.rotz(sol.meta.yaw_w*180/pi);
   o_w  = [sol.meta.x_w; sol.meta.y_w; sol.meta.z_w];

   % change coordinates from wheel to track to rail coordinates

   coords = R_r' * (o_w-o_r)*ones(1,m*n) + R_r' * R_w * coords;

   % extract & reshape to the original array sizes

   xr = reshape(coords(1,:), m, n);
   yr = reshape(coords(2,:), m, n);
   zr = reshape(coords(3,:), m, n);

  end % wheel_to_rail_coords

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ rot ] = rotx(roll_deg)
   % #rotx -2

   sn = sin(roll_deg*pi/180);
   cs = cos(roll_deg*pi/180);
   rot = [ 1,  0,   0;
    0, cs, -sn;
    0, sn,  cs];

  end % rotx

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ rot ] = roty(pitch_deg)
   % #roty -2

   sn = sin(pitch_deg*pi/180);
   cs = cos(pitch_deg*pi/180);
   rot = [ cs,  0,  sn;
    0,  1,   0;
    -sn,  0,  cs];

  end % roty

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ rot ] = rotz(yaw_deg)
   % #rotz -2

   sn = sin(yaw_deg*pi/180);
   cs = cos(yaw_deg*pi/180);
   rot = [ cs, -sn, 0;
    sn,  cs, 0;
    0,   0, 1];

  end % rotz

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ is_left ] = is_left_side(sol)
   % #is_left_side -2

   if (length(sol)>1)
    is_left = (sol(1).config==0 | sol(1).config==4);
   else
    is_left = (sol.config==0 | sol.config==4);
   end

  end % is_left_side

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% #Plot

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ opts2 ] = plot2d(sol, opt2)

   %
   % [ opts2 ] = plot2d([sol], [opt2])
   %
   % Make a plot of tangential tractions Px or Py along a slice with constant x
   % or constant y.
   %
   %  sol   == structure with CONTACT results as returned by loadcase.
   %
   %  opt2, opts2 == (in/output) structures with options for plot2d.
   %           A valid structure opts2 is returned when plot2d is called with
   %           no input arguments.
   %
   %   - orie  = flag for orientation of selected slices,
   %               1 = grid column(s) with constant x, 2 = grid row(s), const y.
   %                   Default: 2.
   %   - pxory = direction of tractions to plot, 'x' or 'y' for Px resp. Py.
   %             Default: 'x'.
   %   - xslc  = x-coordinate(s) for a column-wise plot (orie=1). Note: the
   %             grid column(s) closest to xslc is/are used. Default: 0.0.
   %   - yslc  = y-coordinate(s) for a row-wise plot (orie=2). Note: the grid
   %             row(s) closest to yslc is/are used. Default: 0.0.
   %   - negpn = flag for vertical range to be plotted.
   %              1 = show positive traction bound  "fstat*pn"
   %              0 = show positive and negative traction bounds
   %             -1 = show negative traction bound "-fstat*pn"
   %   - facpt = multiplication factor for negating Px/Py.
   %              1 = show original Px/Py.
   %             -1 = show -Px (pxory=='x') or -Py (pxory=='y').
   %   - xlim  = range of x-coordinates plotted (orie=2). Default: [x0-dx:x1+dx].
   %   - ylim  = range of y-coordinates plotted (orie=1). Default: [y0-dy:y1+dy].
   %   - plim  = vertical range of plot. Defaults: negpn=1: [0-eps:pmax+eps],
   %             negpn=0: [-pmax-eps:pmax+eps], negpn=-1: [-pmax-eps:0+eps].
   %   - pn_linestyle = Matlab linestyle for traction bound (default '--')
   %   - pt_linestyle = Matlab linestyle for tangential tractions (default '-o')
   %

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % construct a local struct "myopt" in which default values are filled in for
   % all available options.
   % #plot2d

   myopt = struct( ...
    'orie',    2, ...
    'pxory', 'x', ...
    'xslc',   0., ...
    'yslc',   0., ...
    'negpn',   1, ...
    'facpt',   1, ...
    'xlim',   [], ...
    'ylim',   [], ...
    'plim',   [], ...
    'pn_linestyle', '--', ...
    'pt_linestyle', '-o'  ...
    );

   % If the user has not supplied any arguments,
   %    return default options as first and only output argument

   if (nargin<1)
    opts2 = myopt;
    return
   end

   % If the user has not supplied an opt2-struct, use the default

   if (nargin<2 | isempty(opt2))
    opt2 = myopt;
   end

   % Check whether user-supplied opt2-struct contains unknown options

   useropts = fieldnames(opt2);
   ierror = 0;
   for i = 1:length(useropts)
    if (~isfield(myopt, useropts{i}))
     ierror = ierror + 1;
     if (isempty(sol) | nargout>=1)
      disp(sprintf('Removing unknown option "%s"',useropts{i}));
      opt2=rmfield(opt2,useropts{i});
     else
      disp(sprintf('Unknown option "%s" will be ignored',useropts{i}));
      if (ierror==1)
       disp(sprintf('You may use "opt2=plot2d([],opt2);" to remove the option',useropts{i}));
      end
     end
    end
   end

   % Overwrite all values in "myopt" with user-supplied values

   myopts = fieldnames(myopt);
   for i = 1:length(myopts)
    if (isfield(opt2, myopts{i}) & ~isempty(cntc.mygetfield(opt2,myopts{i})))
     myopt = setfield(myopt, myopts{i}, cntc.mygetfield(opt2,myopts{i}));
    end
   end

   % Return options myopt: completed with defaults, removed unknown fields

   if (isempty(sol))
    opts2 = myopt;
    return
   end

   % Add defaults for missing components of sol - for convenience of C.library

   if (~isfield(sol, 'x'))
    sol.x = sol.xl + ([1:sol.mx]-0.5)*sol.dx;
   end
   if (~isfield(sol, 'y'))
    sol.y = sol.yl + ([1:sol.my]-0.5)*sol.dy;
   end
   if (~isfield(sol,'fric')), sol.fric.fstat = 0.3; end
   if (~isfield(sol, 'mu'))
    sol.mu = sol.fric.fstat * ones(size(sol.pn));
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % start actual processing to produce the requested plot
   %

   %------------------------------------------------------------------------
   % determine indices islc corresponding to requested x or y-value(s)

   islc = [];

   if (myopt.orie == 1)          % plot tractions along grid column with const x

    if (~isfield(myopt,'xslc') | isempty(myopt.xslc))
     disp('ERROR: when plotting along grid columns, opt2.xslc may not be empty');
     return
    end
    for j = 1:length(myopt.xslc)
     [xmin,ix] = min(abs(sol.x-myopt.xslc(j)));
     if (abs(sol.x(ix)-myopt.xslc(j))>0.001)
      disp(sprintf('Using grid column ix=%d with x=%f, nearest to x=%f', ...
       ix, sol.x(ix), myopt.xslc(j)));
     end
     islc(j) = ix;
    end

   elseif (myopt.orie == 2)      % plot tractions along grid row with const y

    if (~isfield(myopt,'yslc') | isempty(myopt.yslc))
     disp('ERROR: when plotting along grid rows, opt2.yslc may not be empty');
     return
    end
    for j = 1:length(myopt.yslc)
     [ymin,iy] = min(abs(sol.y-myopt.yslc(j)));
     if (abs(sol.y(iy)-myopt.yslc(j))>0.001)
      disp(sprintf('Using grid row iy=%d with y=%f, nearest to y=%f', ...
       iy, sol.y(iy), myopt.yslc(j)));
     end
     islc(j) = iy;
    end

   else

    disp(['ERROR: invalid value for orientation: ',num2str(myopt.orie)]);
    return

   end

   %------------------------------------------------------------------------
   % get information for horizontal axis for appropriate orientation

   if (myopt.orie == 1)          % plot tractions along grid column with const x

    %  determine range of indices ia0/1 with non-zero pn (px/py)

    %  the following is not supported in Matlab 5.3:
    %  ia0 = find( any(sol.pn(:,islc)>0, 2), 1, 'first')-1;
    ia0 = find( any(sol.pn(:,islc)>0, 2) ) - 1;
    ia1 = find( any(sol.pn(:,islc)>0, 2) ) + 1;
    if (isempty(ia0)), ia0=1; ia1=sol.my; end;
    ia0 = max(1, ia0(1));
    ia1 = min(sol.my, ia1(end));
    [ia0 ia1]

    %  determine plot limits for horizontal axis

    if (~isfield(myopt,'ylim') | isempty(myopt.ylim))
     axlim = [ sol.yl-sol.dy, sol.yl+(sol.my+1)*sol.dy ];
    else
     axlim = myopt.ylim;
    end

    axcoor = sol.y;
    axname = 'Y-coordinate [mm]';

   elseif (myopt.orie == 2)      % plot tractions along row y=const

    %  determine range of indices ia with non-zero pn (px/py)

    ia0 = find( any(sol.pn(islc,:)>0, 1) ) - 1;
    ia1 = find( any(sol.pn(islc,:)>0, 1) ) + 1;
    if (isempty(ia0)), ia0=1; ia1=sol.mx; end;
    ia0 = max(1, ia0(1));
    ia1 = min(sol.mx, ia1(end));

    %  determine plot limits for horizontal axis

    if (~isfield(myopt,'xlim') | isempty(myopt.xlim))
     axlim = [ sol.xl-sol.dx, sol.xl+(sol.mx+1)*sol.dx ];
    else
     axlim = myopt.xlim;
    end

    axcoor = sol.x;
    axname = 'X-coordinate [mm]';

   else
    disp(['ERROR: invalid value for orientation: ',num2str(myopt.orie)]);
    return
   end

   %------------------------------------------------------------------------
   % get tangential tractions of appropriate direction (px='x', py='y').

   if (strcmp(myopt.pxory,'x'))          % plot tangential tractions px
    pt = sol.px;
    pname = 'P_x';
   elseif (strcmp(myopt.pxory,'y'))      % plot tangential tractions py
    pt = sol.py;
    pname = 'P_y';
   else
    disp(['ERROR: invalid value for p-direction: ',myopt.pxory]);
    return
   end

   % get selected slices
   if (myopt.orie == 1)          % plot columns x=const, ix=islc
    trcbnd = sol.mu(ia0:ia1, islc) .* sol.pn(ia0:ia1, islc);
    pt     =     pt(ia0:ia1, islc);
   else                        % plot rows y=const, iy=islc
    trcbnd = sol.mu(islc, ia0:ia1)' .* sol.pn(islc, ia0:ia1)';
    pt     =     pt(islc, ia0:ia1)';
   end

   %------------------------------------------------------------------------
   % determine appropriate vertical-range plim for plot:

   if (isfield(myopt,'plim') & ~isempty(myopt.plim))

    %  vertical plot limits provided by user:

    pmin = myopt.plim(1);
    pmax = myopt.plim(2);

   else

    %  derive vertical plot limits from solution.

    pmax = max(max(trcbnd));

    %  find nice round value for pstep such that pmax ~= 10*pstep
    pdec = 1e-9;
    pfac = [1, 2, 2.5, 5];
    found = 0;
    while (~found)
     for k = 1 : length(pfac)
      if (pmax <= 12*pfac(k)*pdec)
       found = 1;
       pstep = pfac(k)*pdec;
       break;
      end
     end
     if (~found)
      pdec = pdec * 10;
     end
    end

    %  set limits according to expected sign of tangential tractions
    if (myopt.negpn>0)
     pmax = ceil(pmax / pstep)*pstep;
     pmin = -0.5*pstep;
    elseif (myopt.negpn==0)
     pmax = ceil(pmax / pstep)*pstep;
     pmin = -pmax;
    else
     pmin = -ceil(pmax / pstep)*pstep;
     pmax = 0.5*pstep;
    end
   end

   plim = [ pmin, pmax ];

   %------------------------------------------------------------------------
   % layout configuration parameters:

   xtickheight   = 0.01*(plim(2)-plim(1));
   xticklabelsep = 1*xtickheight;
   ytickwidth    = 0.007*(axlim(2)-axlim(1));
   yticklabelsep = 2*ytickwidth;

   %------------------------------------------------------------------------
   % plot normal traction, tangential traction

   clf
   hold on

   % plot traction bound, positive and/or negated

   if (myopt.negpn>=0)
    plot(axcoor(ia0:ia1),  trcbnd, myopt.pn_linestyle);
   end
   if (myopt.negpn<=0)
    plot(axcoor(ia0:ia1), -trcbnd, myopt.pn_linestyle);
   end

   % plot tangential tractions

   plot(axcoor(ia0:ia1), myopt.facpt*pt, myopt.pt_linestyle);

   % use scaling derived from solution or provided through arguments:

   axis([axlim plim]);

   % disable Matlab's axis

   xtick = get(gca,'xtick');
   ytick = get(gca,'ytick');
   axis off

   % plot own horizontal axis with ticks and tick-labels

   plot(axlim,[0 0],'k-');
   for x = xtick,
    if (abs(x)<1e-6), x=0; end
    plot([x x],[0 1]*xtickheight,'k');
    text(x, -xticklabelsep, num2str(x), 'verticalalignment','top', ...
     'horizontalalignment','center');
   end

   % plot vertical axis with ticks and tick-labels

   plot([0 0],[pmin pmax]*1.1,'k-');
   for y = ytick,
    if (abs(y)>1e-6)
     plot([-1 1]*ytickwidth,[y y],'k-');
     text(-yticklabelsep, y, num2str(y), 'horizontalalignment','right');
    end
   end

   % set labels on plot

   text(axlim(2), -6*xticklabelsep, axname, ...
    'verticalalignment','top', 'horizontalalignment','right');
   text(yticklabelsep, 1.1*ytick(end)-0.1*ytick(end-1), ...
    ['Tractions ',pname,', f_{stat}*P_n [N/mm^2]']);

   if (nargout>=1)
    opts2 = myopt;
   end
  end

  function PlotEvt(obj,evt)
   %% #PlotEvt cntc.PlotEvt(gcf,evt);
   gf=ancestor(obj,'figure');
   sdt=getappdata(gf,'sdt');
   sdt.cntc=evt;
   cntc.scroll(gf,struct);
  end

  function scroll(obj,evt)
   %% #scroll / key channel changes -2

   gf=ancestor(obj,'figure'); sdt=getappdata(gf,'sdt');
   LI=cntc.call;
   if ~isfield(sdt,'Scroll')
    %% prepare interactivity (defines the sdt appdata)
    iimouse('interacturn',gf, ...
     {'Key.leftarrow{cntc.scroll,"previous"}';
     'Normal+Scroll.@ga{cntc.scroll,"change step"}'
     'Key.rightarrow{cntc.scroll,"previous"}'})
    sdt=getappdata(gf,'sdt');
   end
   if isfield(sdt,'cntc'); ev1=sdt.cntc;
   else;ev1=struct('ch',1);
   end

   if isfield(evt,'Key')
    switch evt.Key
     case {'rightarrow'}
      ev1.ch=min(ev1.ch+1,length(LI.sol));
     case {'leftarrow'}
      ev1.ch=max(1,ev1.ch-1);

    end
   elseif isfield(evt,'VerticalScrollCount')
    ev1.ch=max(min(ev1.ch+evt.VerticalScrollCount,length(LI.sol)),1);
   end

   if isfield(ev1,'ch')
    sdt.cntc=ev1;
    sol=LI.sol{ev1.ch};ev1.opt.ch=ev1.ch;
    if isfield(ev1,'doSol')
     evt=ev1;
     delete(findobj(gf,evt.prop{1:2})); % do not overlay multiple
     evt.doSol(sol);  title(eval(evt.LabFcn))
    else
     feval(ev1.do,sol,ev1.opt)
    end
   end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ opts3 ] = plot3d(sol, opt3, prr, prw, subs)

   %
   % [ opts3 ] = plot3d([sol], [opt3], [prr], [prw], [subs])
   %
   % Make a plot of different quantities in the contact area.
   %
   %  sol    == structure with CONTACT results as returned by loadcase.
   %            may be empty; may be array of sol structures (#patches>1).
   %
   %  opt3, opts3 == structure with options for plot3d.
   %           A valid structure opts3 is returned when plot3d is called with
   %           no input arguments.
   %
   %  prr, prw == rail and wheel profiles (structs) or profile filenames.
   %           prr may be a variable profile in slcs-format, prw may be in slcw-format.
   %
   %  subs   == structure with CONTACT results as returned by loadstrs.
   %            may be empty; may be array of subs structures (#patches>1, 1 block/patch).
   %
   %  opt3.field:
   %   'eldiv', 'eldiv_spy', 'eldiv_contour': different ways of presenting the element division;
   %   'h':        the undeformed distance of the two bodies;
   %   'mu', 'taucrt':   the actual local coefficient of friction and the critical shear stress
   %               for the tangential plasticity model;
   %   'pn', 'px', 'py': the normal and tangential tractions acting on body (1),
   %               i.e. the rail in module 1;
   %   'ptabs', 'ptarg': the magnitude and direction of the tangential tractions;
   %   'ptvec':    a vector-plot of the tangential tractions;
   %   'ptabs+vec': show magnitude and direction of tangential tractions in one plot;
   %   'un', 'ux', 'uy': the displacement differences in normal and tangential directions;
   %   'uplsx', 'uplsy': the components of the accumulated plastic deformation (if M=4);
   %   'utabs+vec', 'upls+vec': elastic and plastic displacements in tangential directions,
   %               magnitude and orientation;
   %   'sx', 'sy': the components of the local shift (T=1) or relative micro-slip velocity (T=2,3);
   %   'shft', 'sabs', 'srel': the magnitude of the shift or absolute or relative slip velocity;
   %   'shft+vec', 'sabs+vec', 'srel+vec': show in one plot the magnitude and direction of the
   %               local shift resp. the absolute or relative slip velocity;
   %   'fricdens': the frictional power density;
   %   'temp1', 'temp2': surface temperatures of bodies 1 and 2 (if H>=1).
   %
   %  Other options:
   %    rw_surfc: 'none' for flat view, 'prr' for rail coordinates, 'prw' for wheel surface view,
   %               or 'both' for track coordinates;
   %    exterval: the value to be plotted at points of the exterior area;
   %    typplot:  type of plot: 'surf' (default), 'contourf', 'rw_rear' or 'rw_side';
   %    view:     the view direction, e.g. [30 30], or 'rail' for [90 -90], 'rdyn' for [180 90];
   %    xrange:   the range of the x-axis to be displayed in 3D plots
   %    xysteps:  one or two step sizes for sampling surfaces along the x- and y-axes
   %    zrange:   the range of the z-axis to be displayed;
   %    ixrange:  the selection of the elements [1,mx] in x-direction to be displayed in the plot;
   %    iyrange:  the selection of the elements [1,my] in y-direction to be displayed in the plot;
   %    numvecx:  the maximum number of vectors to be displayed in x-direction;
   %    numvecy:  same as numvecx for y-direction;
   %    xstretch: stretching of x-axis esp. for vectors on elongated contact patches (b>>a)
   %    vecscale: manual scaling factor ([mm] per [N/mm2]) for vectors on vector-plots.
   %               used also for quantities in 'rw_rear' and 'rw_side' views
   %    veccolor: color specification for vectors on vector-plots. Default: 'b' for ptvec, 'k' for ptabs+vec;
   %               used also for quantitites in 'rw_rear' and 'rw_side' views, default matlab_color(6)
   %    vecwidth: line-width for vectors on vector-plots;
   %    addeldiv: show contours of element divisions on top of results. <=0: no contours; >=1: show contact
   %              area, >=2: show adhesion area. Available in 2D plots, e.g. view [0 90]
   %    eldivcol: color specification for element division, e.g. ['b';'m';'g']
   %              or [0 0 0.56; 0.5 0 0; 0.8 0.25 0.42]
   %    eldivwid: line-width for contours for element division;
   %    colormap: this changes the colormap for 3D plots;
   %    addplot:  clear current subplot (-1), clear figure (0) or do not clear (1) the figure before plotting.
   %

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % determine version information on the Matlab graphics system
   % #plot3d

   old_graphics = verLessThan('matlab', '8.4.0'); % 8.4 == R2014b
   new_graphics = ~old_graphics;

   % construct a local struct "myopt" in which default values are filled in for
   % all available options.

   myopt = struct( ...
    'field',    'default', ...
    'rw_surfc', 'none', ...
    'exterval',  NaN, ...
    'typplot',  'surf', ...
    'view',     'default', ...
    'xrange',    [], ...
    'xysteps',   [], ...
    'zrange',    [], ...
    'ixrange',  'auto', ...
    'iyrange',  'auto', ...
    'numvecx',   15, ...
    'numvecy',   15, ...
    'xstretch',  1, ...
    'vecscale',  [], ...
    'veccolor',  [], ...
    'vecwidth',  [], ...
    'addeldiv',  [], ...
    'eldivcol',  [0 0 0.56; 0.5 0 0; 0.8 0.25 0.42], ...
    'eldivwid',  [], ...
    'colormap', 'parula', ...
    'addplot',    0,  ... % clear plot (0) or add to exist.plot (1, experimental)
    'ch',1 ... % current channel
    );
   if (old_graphics); myopt.colormap = 'jet'; end

   % If the user has not supplied any arguments,
   %    return default options as first and only output argument

   if (nargin<1); opts3 = myopt;return
   elseif isempty(sol);
    %% #plot3d.ch automated multiple plots
    LI=cntc.call;sol=LI.sol;
    if ischar(opt3)
     % cntc.plot3d([],'{typpplot=surf,field=px,rw_surfc=prw,ch4,gf10}')
     [~,opt3]=sdtm.urnPar(opt3,'{}{ch%g,gf%g}');
     opt3=sdth.sfield('addmissing',opt3,myopt);
     r3=sdtm.findTok(opt3.Other,'([^=]*)=(.*)');
     for j3=1:size(r3,1);opt3.(r3{j3}{1})=r3{j3}{2};end
    end
    if isfield(opt3,'ch')
     if ~isfield(opt3,'gf');opt3.gf=0;
     else;opt3.gf=opt3.gf(1)+opt3.ch;
     end % Figure Offset
     opt=sdtm.rmfield(opt3,'gf','Other');
     for j1=1:length(opt3.ch)
      gf=figure(opt3.gf(j1));opt.ch=opt3.ch(j1);
      cntc.plot3d(LI.sol{opt.ch},opt)
      sdt.cntc=struct('do',@cntc.plot3d,'ch',opt3.ch(j1),'opt',opt);
     end
     return
    end
   elseif (~isstruct(sol) | ~isfield(sol,'pn'))
    disp('ERROR(plot3d): argument sol should be a structure as obtained from loadcase.');
    return;
   end

   % If the user has not supplied an opt3-struct, use the default

   if (nargin<2 | isempty(opt3));opt3 = myopt;
   elseif (~isstruct(opt3))
    disp('ERROR(plot3d): argument opt3 should be a structure as obtained from plot3d.');
    return;
   end

   if (nargin<3)
    if (isfield(sol,'prr'));prr = sol.prr;
    else; prr = [];
    end
   end
   if (nargin<4)
    if (isfield(sol,'prw'));prw = sol.prw;else;prw = [];end
   end
   if (nargin<5);subs = [];end

   % Check whether user-supplied opt3-struct contains unknown options

   useropts = fieldnames(opt3);
   ierror = 0;
   for i = 1:length(useropts)
    if (~isfield(myopt, useropts{i}))
     ierror = ierror + 1;
     disp(sprintf('Unknown option "%s" will be ignored',useropts{i}));
     if (ierror==1)
      disp(sprintf('You may use "opt3=rmfield(opt3,''%s'');" to remove the option',useropts{i}));
     end
    end
   end

   % Overwrite all values in "myopt" with user-supplied values

   myopts = fieldnames(myopt);
   for i = 1:length(myopts)
    if (isfield(opt3, myopts{i}))
     r1= cntc.mygetfield(opt3,myopts{i});
     myopt.(myopts{i})=r1;
    end
   end

   % Set default field

   if (strcmp(myopt.field, 'default'))
    if (isempty(sol) | sol(1).kincns.t_digit<=0); myopt.field = 'pn';
    else;     myopt.field = 'ptabs+vec';
    end
   end

   % Check whether requested field is valid

   ok_fields=strvcat('eldiv', 'eldiv_spy', 'eldiv_contour', 'h', 'mu', ...
    'pn', 'px', 'py', 'ptabs', 'ptarg', 'ptvec', 'ptabs+vec', ...
    'un', 'ux', 'uy', 'utabs+vec', ...
    'taucrt', 'uplsx', 'uplsy', 'upls+vec', ...
    'sx', 'sy', 'shft', 'sabs', 'srel', ...
    'shft_vec', 'sabs_vec', 'srel_vec', ...
    'shft+vec', 'sabs+vec', 'srel+vec', 'fricdens', ...
    'temp1', 'temp2');
   found=0;
   for istr=1:size(ok_fields,1)
    if (strcmp(myopt.field, deblank(ok_fields(istr,:))))
     found=1;
    end
   end
   user_added_fld = 0;
   if (~found & isfield(sol, myopt.field))
    % disp(sprintf('Plotting user-added field "%s"', myopt.field));
    user_added_fld = 1;
    found = 1;
   end
   if (~found)
    disp(sprintf('Unknown field requested="%s"; Available options:',myopt.field));
    disp(ok_fields)
    return;
   end

   % Return if no solution is given (handy for w/r contact: call plot3d for fixed #patches)

   if (isempty(sol))
    % disp('empty sol, returning')
    return;
   end

   % prepare the rail profile, if provided

   if (ischar(prr))
    is_wheel = 0; mirror_y = 0; % note: Miniprof profiles usually need to be mirrored
    prr = cntc.read_profile(prr, is_wheel, mirror_y);
   end
   if (~isempty(prr) & ~isfield(prr,'nslc'))
    if (isempty(prr.ProfileY) | isempty(prr.ProfileZ))
     disp(sprintf('ERROR: rail profile does not contain any data.'));
     prr = [];
    elseif (prr.ProfileY(end) - prr.ProfileY(1) < 0)
     disp(sprintf('ERROR: rail profile Y-coordinates must be ascending.'));
     return
    end
   end

   % prepare the wheel profile, if provided

   if (ischar(prw))
    is_wheel = 1; mirror_y = 0;
    prw = cntc.read_profile(prw, is_wheel, mirror_y);
   end
   if (~isempty(prw) & ~isfield(prw,'nslc'))
    if (isempty(prw.ProfileY) | isempty(prw.ProfileZ))
     disp(sprintf('ERROR: wheel profile does not contain any data.'));
     prw = [];
    end
    if (~isempty(prw) & prw.ProfileY(end) - prw.ProfileY(1) > 0)
     disp(sprintf('ERROR: wheel profile Y-coordinates must be descending.'));
     return
    end
   end

   % Check whether requested rail or wheel surface is available

   if (~any(strcmp(myopt.rw_surfc,{'none','prr','prw','both'})))
    disp(sprintf('Unknown option rw_surfc="%s" will be ignored',myopt.rw_surfc));
    disp(sprintf('Available values are "none" (default), "prr", "prw", and "both".'));
    myopt.rw_surfc = 'none';
   end

   if (isempty(prw) & any(strcmp(myopt.rw_surfc,{'prw','both'})))
    disp(sprintf('No wheel profile provided; ignoring rw_surfc="%s". ', myopt.rw_surfc));
    if (strcmp(myopt.rw_surfc, 'both'))
     myopt.rw_surfc = 'prr';
    else
     myopt.rw_surfc = 'none';
    end
   end
   if (isempty(prr) & any(strcmp(myopt.rw_surfc,{'prr','both'})))
    disp(sprintf('No rail profile provided; ignoring rw_surfc="%s".', myopt.rw_surfc));
    if (strcmp(myopt.rw_surfc, 'both') & ~isempty(prw))
     myopt.rw_surfc = 'prw';
    else
     myopt.rw_surfc = 'none';
    end
   end

   if (~strcmp(myopt.rw_surfc,'none') & any(strcmp(myopt.field, {'eldiv_spy'})))
    disp(sprintf('Option "%s" for rw_surfc is not available for field "%s"', myopt.rw_surfc, myopt.field));
    myopt.rw_surfc = 'none';
   end

   % Check whether requested typplot is valid, set default when needed

   if (~any(strcmp(myopt.typplot,{'surf','contourf','rw_rear','rw_side'})))
    %    typplot:  type of plot: 'surf' (default), 'contourf', 'rw_rear' or 'rw_side';
    disp(sprintf('Unknown value "%s" for option typplot will be ignored.',myopt.typplot));
    disp(sprintf('Available values are "surf" (default), "contourf", "rw_rear", and "rw_side".'));
    myopt.typplot = 'surf';
   end
   if (any(strcmp(myopt.typplot,{'rw_rear','rw_side'})) & ~any(strcmp(myopt.rw_surfc,{'prr','prw','both'})))
    disp(sprintf('Option rw_surfc="%s" therefore option typplot="%s" will be ignored', myopt.rw_surfc, ...
     myopt.typplot));
    myopt.typplot = 'surf';
   end

   % Check whether requested view is valid, set default when needed

   if (isempty(myopt.view)), myopt.view = 'default'; end
   if (ischar(myopt.view))
    if (strcmp(myopt.view,'rail') & strcmp(myopt.typplot,'surf') & ~strcmp(myopt.rw_surfc,'none'))
     % 3d surface plot should not use 2d view
     myopt.view = 'default';
    end
    if (strcmp(myopt.view,'rail'))
     % wheel-rail view: 2D, rolling direction == plot y-direction
     myopt.view=[90 -90];
    elseif (strcmp(myopt.view,'rdyn'))
     % RecurDyn view: 2D, x- and y-axes reversed
     myopt.view=[180 90];
    elseif (strcmp(myopt.view,'default') & ~strcmp(myopt.rw_surfc,'none'))
     % default: 3D wheel-rail profile view

     if cntc.is_left_side(sol)
      myopt.view=[115 30];
     else
      myopt.view=[65 30];
     end
    else
     % default: 2D or 3D view depending on field to be plotted
     if (~strcmp(myopt.view,'default'))
      disp(sprintf('Unknown option view="%s" will be ignored',myopt.view));
      opt.view = 'default';
     end
     if (strcmp(myopt.field,'h')        || ...
       strcmp(myopt.field,'fricdens') );
      myopt.view=[-37.5 30];
     else
      myopt.view=[0 90]; % traditional view, rolling to the right
      % myopt.view=[90 -90]; % rail view, rolling up in the picture
     end
    end
   end

   if (isempty(myopt.colormap))
    myopt.colormap = 'default'; % matlab built-in
   end

   % fill in default for addeldiv (eldiv_contour)

   if (isempty(myopt.addeldiv))
    if (strcmp(myopt.field,'eldiv') | strcmp(myopt.field,'eldiv_spy'))
     myopt.addeldiv = 0;
    else
     myopt.addeldiv = 2;
    end
   end
   if (strcmp(myopt.field,'eldiv_contour'))
    myopt.addeldiv = max(1, myopt.addeldiv);
   end
   if (any(strcmp(myopt.typplot,{'rw_rear','rw_side'})))
    myopt.addeldiv = 0;
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Plot arrays of sol structures using recursion (esp. w/r contacts with multiple patches)

   if (isstruct(sol) & length(sol)>1)
    if (strcmp(myopt.rw_surfc,'none') | (isempty(prr) & isempty(prw)))
     disp(sprintf('ERROR: plotting %d patches at once requires using rail or wheel profiles', length(sol)));
     return
    end

    for isol = 1 : length(sol)
     pmax(isol) = max(max(sol(isol).pn));
    end
    myopt.vecscale = cntc.get_vecscale( sol, myopt, pmax ); % fix vecscale -- same for all patches

    [~,jsol] = sort(pmax, 'descend');    % process highest->lowest -- axes based on first one
    for i = 1 : length(sol)
     cntc.plot3d( sol(jsol(i)), myopt, prr, prw );
     myopt.addplot = 1;
    end
    return
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % Add profiles to the sol struct

   if (any(strcmp(myopt.rw_surfc, {'prr','both'})))
    if (isfield(prr,'nslc')) % variable profile
     sol.slcs = prr;
     sol.prr = cntc.get_profile_slice(sol.slcs, sol.meta.s_ws+sol.meta.xcp_r);
    else
     sol.prr = prr;
    end
   end
   if (any(strcmp(myopt.rw_surfc, {'prw','both'})))
    if (isfield(prw,'nslc')) % variable profile
     sol.slcw = prw;
     sol.prw = cntc.get_profile_slice(sol.slcw, -sol.meta.th_ws+sol.meta.xcp_w/sol.meta.rnom_whl,1);
    else
     sol.prw = prw;
    end
   end
   clear prr prw;

   % Add subsurface data to the sol struct

   if (~isempty(subs))
    sol.subs = subs;
   end

   % Add defaults for missing components of sol - for convenience of C.library

   if (~isfield(sol, 'x_offset')), sol.x_offset = []; end
   if (~isfield(sol, 'y_offset')), sol.y_offset = []; end
   if (~cntc.myisfield(sol, 'kincns.t_digit')), sol.kincns.t_digit = 3; end
   if (~isfield(sol, 'x'))
    sol.x = sol.xl + ([1:sol.mx]-0.5)*sol.dx;
   end
   if (~isfield(sol, 'y'))
    sol.y = sol.yl + ([1:sol.my]-0.5)*sol.dy;
   end

   % Check sizes of main solution arrays - for convenience of C. library

   if (any(size(sol.pn)~=[sol.my sol.mx]))
    disp(sprintf(['ERROR: incorrect size for array pn (%d,%d), ', ...
     'expecting (%d,%d)'], size(sol.pn), sol.my, sol.mx));
    return
   end
   if (any(size(sol.eldiv)~=[sol.my sol.mx]))
    disp(sprintf(['ERROR: incorrect size for array eldiv (%d,%d), ', ...
     'expecting (%d,%d)'], size(sol.eldiv), sol.my, sol.mx));
    return
   end

   % Prefix for figure title, dependent on type of problem

   pfxtit = '';

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % start actual processing to produce the requested plot
   %

   if (myopt.addplot==0 | myopt.addplot<-1)
    clf;  % restart figure
   elseif (myopt.addplot==-1 | myopt.addplot==2)
    cla;  % restart subplot?
   end

   % start by plotting the rail and/or wheel surfaces

   show_slcs = (isfield(sol,'slcs') & any(strcmp(myopt.rw_surfc,{'prr','both'})));
   show_slcw = (isfield(sol,'slcw') & any(strcmp(myopt.rw_surfc,{'prw','both'})));

   %% Show and plot the profile xxxgae P
   if (myopt.addplot<=0 & ~strcmp(myopt.rw_surfc,'none'))  % create new plot
    cntc.show_profiles(sol, myopt);
   elseif ( any(strcmp(myopt.typplot,{'rw_side','rw_rear'})) & (show_slcs | show_slcw) )
    cntc.show_profiles(sol, myopt);   % add new cross-section to existing plot
   end
   hold on;

   % set element numbers ix, iy to be plotted: interval of [1:mx(my)]
   [ix_plot, iy_plot] = cntc.plot_ranges(sol, myopt);

   % make the plot, depending on the field that is requested
   evt=struct('ch',myopt.ch,'LabFcn','sprintf(''ch %i, %s'',evt.ch,evt.title)');
   if isfield(myopt,'ch');stCh=sprintf('ch%i ',myopt.ch);else;stCh='';end

   evt.doSol=@(sol)cntc.show_scalar_field(sol,myopt.field, myopt);
   evt.prop={'tag',myopt.field};
   switch myopt.field
    case 'h';  evt.title='Undeformed distance H';
    case 'mu'; evt.title='Actual local friction coefficient \mu(x)';
    case 'pn'; evt.title='Normal pressure P_n';
    case 'px'
     evt.doSol=@(sol)cntc.show_scalar_field(sol, ...
      cntc.derived_data(sol, myopt, 'px'), myopt);
     evt.title=[pfxtit, 'Tangential traction P_x'];
    case 'py'
     evt.doSol=@(sol)cntc.show_scalar_field(sol, ...
      cntc.derived_data(sol, myopt, 'py'), myopt);
     evt.title=[pfxtit, 'Tangential traction P_y'];
   end
   if isfield(evt,'title')
    evt.gf=gcf; cntc.PlotEvt(gcf,evt);
   end
   if (strcmp(myopt.field,'ptabs') | strcmp(myopt.field,'ptabs+vec'))

    pt = cntc.derived_data(sol, myopt, 'pt');
    cntc.show_scalar_field(sol, pt, myopt);
    title ([pfxtit, 'Magnitude of tangential tractions |P_t|']);

   end
   if (strcmp(myopt.field,'ptvec') | strcmp(myopt.field,'ptabs+vec'))

    [px, py] = cntc.derived_data(sol, myopt, 'px', 'py');
    cntc.show_vec_field(sol, px, py, myopt);
    title ([pfxtit, 'Tangential tractions p_t']);

   end
   if (strcmp(myopt.field,'ptarg'))

    [eldiv, ptarg] = cntc.derived_data(sol, myopt, 'eldiv', 'ptarg');
    % ptarg: range -180:180 deg
    % shift lower bound (-180) to "arg_low", to avoid wrap-around in the plot
    msk=(eldiv>0);
    cnt(1) = nnz(msk.*(ptarg> 135 | ptarg<-135));
    cnt(2) = nnz(msk.*(ptarg>-135 & ptarg<-45));
    cnt(3) = nnz(msk.*(ptarg> -45 & ptarg< 45));
    cnt(4) = nnz(msk.*(ptarg>  45 & ptarg<135));
    [mn,ix]=min(cnt); arg_low=-270+90*ix;
    ix=find(ptarg<arg_low); ptarg(ix)=ptarg(ix)+360;
    ix=find(ptarg>arg_low+360); ptarg(ix)=ptarg(ix)-360;

    cntc.show_scalar_field(sol, ptarg, myopt);
    if (isempty(myopt.zrange))
     ix = find(eldiv>=1);
     rng=[ min(min(ptarg(ix))), max(max(ptarg(ix))) ];
     rng=rng + 0.05*[-1 1]*(rng(2)-rng(1));
     set(gca,'clim', rng);
    end
    title ([pfxtit, 'Direction of tangential tractions arg(P_t)']);

   end
   if (strcmp(myopt.field,'un'))

    cntc.show_scalar_field(sol, 'un', myopt);
    title ('Normal displacement difference U_n');

   end
   if (strcmp(myopt.field,'ux'))

    cntc.show_scalar_field(sol, 'ux', myopt);
    title ([pfxtit, 'Tangential displacement difference U_x']);

   end
   if (strcmp(myopt.field,'uy'))

    cntc.show_scalar_field(sol, 'uy', myopt);
    title ([pfxtit, 'Tangential displacement difference U_y']);

   end
   if (strcmp(myopt.field,'utabs+vec'))

    ut = sqrt(sol.ux.^2 + sol.uy.^2);
    cntc.show_scalar_field(sol, ut, myopt);
    cntc.show_vec_field(sol, sol.ux, sol.uy, myopt);
    title ('Tangential displacements u_t');

   end
   if (strcmp(myopt.field,'uplsx'))

    cntc.show_scalar_field(sol, 'uplsx', myopt);
    title ('Plastic deformation u_{pl,x}');

   end
   if (strcmp(myopt.field,'uplsy'))

    cntc.show_scalar_field(sol, 'uplsy', myopt);
    title ('Plastic deformation u_{pl,y}');

   end
   if (strcmp(myopt.field,'upls+vec'))

    upls = sqrt(sol.uplsx.^2 + sol.uplsy.^2);
    cntc.show_scalar_field(sol, upls, myopt);
    cntc.show_vec_field(sol, sol.uplsx, sol.uplsy, myopt);
    title ('Plastic deformation u_{pl}');

   end
   if (strcmp(myopt.field,'taucrt'))

    cntc.show_scalar_field(sol, 'taucrt', myopt);
    title ('Tangential yield stress \tau_{crt}');

   end

   % compute sx, sy when needed

   if (strcmp(myopt.field,'sx')       | strcmp(myopt.field,'sy') | ...
     strcmp(myopt.field,'shft_vec') | strcmp(myopt.field,'shft+vec') | ...
     strcmp(myopt.field,'sabs_vec') | strcmp(myopt.field,'sabs+vec') | ...
     strcmp(myopt.field,'srel_vec') | strcmp(myopt.field,'srel+vec'))

    if (isfield(sol,'sx') & isfield(sol,'sy'))
     sx = sol.sx;
     sy = sol.sy;
    else
     pxdir = sol.px ./ max(1d-9,sqrt(sol.px.^2 + sol.py.^2));
     pydir = sol.py ./ max(1d-9,sqrt(sol.px.^2 + sol.py.^2));
     if (sol.kincns.t_digit>=2)
      sx = -sol.srel .* pxdir;
      sy = -sol.srel .* pydir;
     else
      sx = -sol.shft .* pxdir;
      sy = -sol.shft .* pydir;
      if (strfind(myopt.field,'sabs'))
       sx = sx * sol.kincns.veloc;
       sy = sy * sol.kincns.veloc;
      end
     end
    end
   end

   if (strcmp(myopt.field,'sx'))

    cntc.show_scalar_field(sol, sx, myopt);
    if (sol.kincns.t_digit>=2)
     title ([pfxtit, 'Relative slip component s_x']);
    else
     title ([pfxtit, 'Shift component S_x']);
    end
   end
   if (strcmp(myopt.field,'sy'))

    cntc.show_scalar_field(sol, sy, myopt);
    if (sol.kincns.t_digit>=2)
     title ([pfxtit, 'Relative slip component s_y']);
    else
     title ([pfxtit, 'Shift component S_y']);
    end

   end
   if (strcmp(myopt.field,'shft') | strcmp(myopt.field,'shft+vec'))

    cntc.show_scalar_field(sol, 'shft', myopt);
    title ([pfxtit, 'Shift distance']);

   end
   if (strcmp(myopt.field,'sabs') | strcmp(myopt.field,'sabs+vec'))

    sabs = cntc.derived_data(sol, myopt, 'sabs');
    cntc.show_scalar_field(sol, sabs, myopt);
    title ([pfxtit, 'Absolute slip velocity']);

   end
   if (strcmp(myopt.field,'srel') | strcmp(myopt.field,'srel+vec'))

    cntc.show_scalar_field(sol, 'srel', myopt);
    title ([pfxtit, 'Relative slip velocity']);

   end
   if (strcmp(myopt.field,'shft_vec') | strcmp(myopt.field,'shft+vec') | ...
     strcmp(myopt.field,'sabs_vec') | strcmp(myopt.field,'sabs+vec') | ...
     strcmp(myopt.field,'srel_vec') | strcmp(myopt.field,'srel+vec'))

    cntc.show_vec_field(sol, sx, sy, myopt)

    if (strcmp(myopt.field,'shft_vec') | strcmp(myopt.field,'shft+vec'))
     title ([pfxtit, 'Tangential shift S_\tau']);
    elseif (strcmp(myopt.field,'sabs_vec') | strcmp(myopt.field,'sabs+vec'))
     title ([pfxtit, 'Absolute slip velocity s_a']);
    elseif (strcmp(myopt.field,'srel_vec') | strcmp(myopt.field,'srel+vec'))
     title ([pfxtit, 'Relative slip velocity s_r']);
    end

   end
   if (strcmp(myopt.field,'fricdens'))

    % unit [W / mm2] = ([mm/s]/1000) * [-]  * [ N/mm2],   [W] = [J/s] = [N.m/s]
    % note: srel is dimensionless, using V.dt=dq=1 in shifts (T=1)
    fricdens = sol.kincns.veloc/1000 * sol.srel .* sqrt(sol.px.^2 + sol.py.^2);
    cntc.show_scalar_field(sol, fricdens, myopt);
    title ([pfxtit, 'Frictional power density']);

   end
   if (strcmp(myopt.field,'temp1'))

    cntc.show_scalar_field(sol, 'temp1', myopt);
    title ('Surface temperature T^{(1)}');

   end
   if (strcmp(myopt.field,'temp2'))

    cntc.show_scalar_field(sol, 'temp2', myopt);
    title ('Surface temperature T^{(2)}');

   end
   if (user_added_fld & isfield(sol, myopt.field))

    cntc.show_scalar_field(sol, myopt.field, myopt);
    title (sprintf('User-added field "%s"', myopt.field));

   end
   if (strcmp(myopt.field,'eldiv'))

    myopt.zrange = [-0.1 2.1];
    eldiv = cntc.derived_data(sol, myopt, 'eldiv');
    ix = find(eldiv==3); eldiv(ix) = 1.5; % set Plast(3) between Adhes(1) and Slip(2)
    cntc.show_scalar_field(sol, eldiv, myopt);
    h = get(gcf,'children'); ix = strmatch('Colorbar', get(h,'tag'));
    if (isscalar(ix))
     jx = find(sol.eldiv==3);
     if (~isempty(jx))
      set(h(ix),'ytick',[0,1,1.5,2],'yticklabel',['exter';'stick';'plast';'slip ']);
     else
      set(h(ix),'ytick',[0,1,2],'yticklabel',['exter';'stick';'slip ']);
     end
     if (isnan(myopt.exterval))
      set(h(ix),'ylim',[0.9 2.1]);
     end
    end
    title ([pfxtit, 'The element divisions']);

   end
   if (strcmp(myopt.field,'eldiv_spy'))

    eldiv = cntc.derived_data(sol, myopt, 'eldiv');
    mask = (eldiv==1); spy(mask, 'b.', 6); % ones in the adhesion area
    mask = (eldiv==2); spy(mask, 'r.', 6); % ones in the slip area
    mask = (eldiv==3); spy(mask, 'g.', 6); % ones in the platicity area
    title ([pfxtit, 'The element divisions (red=slip, blue=adhesion)']);
    axis xy
    ix = get(gca,'xtick')'; iy = get(gca,'ytick')';
    xcentr = sol.xl+(ix+0.5)*sol.dx;
    ycentr = sol.yl+(iy+0.5)*sol.dy;
    if (~isempty(sol.x_offset)), xcentr = xcentr + sol.x_offset; end
    if (~isempty(sol.y_offset)), ycentr = ycentr + sol.y_offset; end
    set(gca,'xticklabel', num2str(xcentr, '%5.3f'));
    set(gca,'yticklabel', num2str(ycentr, '%5.3f'));
    cntc.set_xyzlabels(sol, myopt);
    if (~cntc.is_2d_view(gca, myopt.view))
     v=axis; axis([v -1e-6 1e-6]); axis equal;
     set(gca,'ztick',[]);
    end
    view(myopt.view);
   end
   if ( strcmp(myopt.field,'eldiv_contour') | myopt.addeldiv )
    % eldiv_contour can be used in 2D view & on rail/wheel surface
    if (cntc.is_2d_view(gca,myopt.view) | ~strcmp(myopt.rw_surfc,'none') )
     cntc.show_eldiv_contour(sol, myopt);
    end
   end
   if (any(strcmp(myopt.rw_surfc,{'prr','both'})) & ~any(strcmp(myopt.typplot,{'rw_rear','rw_side'})))
    h = findobj(gcf, 'Tag','prr');
    clim = get(gca,'clim');
    set(h,'cdata', clim(1)*ones(size(get(h,'xdata'))));
   end
   if (any(strcmp(myopt.rw_surfc,{'prw','both'})) & ~any(strcmp(myopt.typplot,{'rw_rear','rw_side'})))
    h = findobj(gcf, 'Tag','prw');
    clim = get(gca,'clim');
    set(h,'cdata', clim(2)*ones(size(get(h,'xdata'))));
   end

  end % plot3d


  function [ myopt ] = plot_2dspline( slcs, opt )

   % [ opt ] = plot_2dspline( slcs, opt )
   %
   % create 3d surface plot of variable rail profile
   %   slcs    - structure with variable profile as returned by read_slices
   %   opt     - plot configuration; called without options, a struct is returned with default settings
   %
   % plot options:
   %   typplot:     type of plot: 'surf' (default), 'refl', 'topol'
   %   urange:      range of surface parameter u for longitudinal direction
   %   vrange:      range of surface parameter v for lateral direction
   %   xysteps:     (approx) distance between drawn lines on surface
   %   show_slc:    list of slice numbers that are highlighted
   %   slc_color:   color-spec for slices that are highlighted
   %   show_feat:   list of feature numbers that are highlighted
   %   feat_color:  color-spec for features that are highlighted
   %   view:        matlab view direction, e.g. [50 25] (azimuth, elevation), can be 'rail' or 'default'
   %   refl_avec:   light source direction for reflection plot, e.g. [0,1,0] or [100, 30] (az, el)
   %   refl_dth:    interval spacing d theta for reflection plot, e.g. 30 deg (default)
   %   zoom:        data aspect ratio, e.g. [1000 1 1] (rail) or [0.1 1 1] (wheel, theta)
   %   addplot:     clear (0) or do not clear (1) the figure before plotting.

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % #plot_2dspline
   if (nargin<1 | isempty(slcs))
    slcs = struct();
   end
   if (nargin<2 | isempty(opt))
    opt = struct();
   end

   % construct a local struct "myopt" in which default values are filled in for all available options.

   myopt = struct( ...
    'typplot',    'surf',  ...
    'urange',     [],  ...
    'vrange',     [],  ...
    'xysteps',    [],  ...
    'show_slc',   [],  ...
    'slc_color',  'r', ...
    'show_feat',  [],  ...
    'feat_color', 'r', ...
    'view',       'default', ...
    'refl_avec',  [0,1,0], ...
    'refl_dth',   30, ...
    'zoom',       [], ...
    'addplot',    0   ...
    );

   if (isempty(fields(slcs)))   % called with no arguments: return default options
    return;
   end
   if (~isfield(slcs, 'u'))
    disp('Error: 2nd argument does not look like a slcs-structure, no u-positions')
    return;
   end
   if (~isfield(slcs, 'spl2d'))
    disp('Error: slcs does not have 2d spline data');
    return;
   end

   % Check whether user-supplied opt3-struct contains unknown options

   useropts = fieldnames(opt);
   ierror = 0;
   for i = 1:length(useropts)
    if (~isfield(myopt, useropts{i}))
     ierror = ierror + 1;
     disp(sprintf('Unknown option "%s" will be ignored',useropts{i}));
     if (ierror==1)
      disp(sprintf('You may use "opt=rmfield(opt,''%s'');" to remove the option',useropts{i}));
     end
    end
   end

   % Overwrite all values in "myopt" with user-supplied values

   myopts = fieldnames(myopt);
   for i = 1:length(myopts)
    if (isfield(opt, myopts{i}))
     myopt = setfield(myopt, myopts{i}, cntc.mygetfield(opt,myopts{i}));
    end
   end

   if (~any(size(myopt.slc_color,2)==[1 3]))
    disp('ERROR: slc_color should be column vector with 1 or 3 entries per row');
    myopt.slc_color = 'r';
   end
   if (~any(size(myopt.feat_color,2)==[1 3]))
    disp('ERROR: feat_color should be column vector with 1 or 3 entries per row');
    myopt.feat_color = 'r';
   end

   if (strcmp(myopt.view,'rail'))
    % wheel-rail view: 2D, rolling direction == plot y-direction
    myopt.view=[90 -90];
   elseif (strcmp(myopt.view,'default'))
    if (strcmp(myopt.typplot,'topol'))
     myopt.view=[90 -90];
    else
     myopt.view=[110 20];
    end
   end

   % start plotting

   if (~myopt.addplot)
    clf;
   end

   % set target step-sizes for plotting

   if (isempty(myopt.urange))
    myopt.urange = [ min(slcs.spl2d.ui), max(slcs.spl2d.ui) ];
   end
   if (isempty(myopt.vrange))
    myopt.vrange = [ min(slcs.spl2d.vj), max(slcs.spl2d.vj) ];
   end
   if (isempty(myopt.xysteps))
    xlen = max(max(slcs.xsurf)) - min(min(slcs.xsurf));
    ylen = max(slcs.spl2d.vj) - min(slcs.spl2d.vj);
    myopt.xysteps = [xlen/4, ylen/10]; % using 4/10 intervals == 5/11 lines
   end
   if (isempty(myopt.zoom))
    if (strcmp(myopt.typplot,'topol'))
     xlen = max(max(slcs.u))  - min(min(slcs.u));
     ylen = max(max(slcs.vj)) - min(min(slcs.vj));
    else
     xlen = max(max(slcs.xsurf)) - min(min(slcs.xsurf));
     ylen = max(max(slcs.ysurf)) - min(min(slcs.ysurf));
    end
    myopt.zoom = [ 0.5*xlen/ylen 1 1 ];       % wheel: data aspect ratio, e.g. [0.05 1 1]
    if (~slcs.is_wheel)
     myopt.zoom(1) = max(1, myopt.zoom(1)); % rail: data aspect ratio could be [1000 1 1]
    end
   end
   if (isscalar(myopt.xysteps))
    myopt.xysteps = [1 1] * myopt.xysteps;
   end

   % form profile surface at given x, all u - fine sampling

   [ xsurf, ysurf, zsurf ] = cntc.make_3d_surface( slcs, myopt, myopt.urange, [], myopt.xysteps(1)/20, ...
    myopt.vrange, [], []);

   % plot 90% transparent surface (alpha=0.1); color-value?

   if (any(strcmp(myopt.typplot,{'surf','topol'})))
    csurf = -1 * ones(size(xsurf));
    l=surf(xsurf, ysurf, zsurf, csurf, 'EdgeAlpha',0.1, 'FaceAlpha',0.1);
   end

   % compute and plot reflection lines

   if (strcmp(myopt.typplot, 'refl'))
    ix = find( slcs.spl2d.tui>= myopt.urange(1) & slcs.spl2d.tui<=myopt.urange(end) );
    ui = cntc.refine_intv( slcs.spl2d.tui(ix), 4 );
    ix = find( slcs.spl2d.tvj>= myopt.vrange(1) & slcs.spl2d.tvj<=myopt.vrange(end) );
    vj = cntc.refine_intv( slcs.spl2d.tvj(ix), 4 );
    [xsurf, ysurf, zsurf, th_rfl] = cntc.make_reflec( slcs.spl2d, myopt.refl_avec, myopt.view, ui, vj);
    th_col = floor(th_rfl/myopt.refl_dth) * myopt.refl_dth;

    % tmp   = parula(12);
    % m_bin = repmat(tmp([3,12],:), n_bin/2, 1);
    n_bin = round(360 / myopt.refl_dth);
    tmp   = [0 0 0; 1 1 .6];
    m_bin = repmat(tmp, n_bin/2, 1);

    l=surf(xsurf, ysurf, zsurf, th_col, 'EdgeColor','none', 'FaceAlpha',0.5);
    colormap(m_bin);
   end

   set(l,'Tag','prr');

   % set view appropriate for rails with z positive downwards

   if (~strcmp(myopt.typplot,'topol'))
    set(gca,'xdir','reverse', 'zdir','reverse');
   end
   axis equal;
   view(myopt.view);
   hold on

   if (strcmp(myopt.typplot,'topol'))
    xlabel('u [-]');
    ylabel('v [-]');
   elseif (slcs.is_wheel)
    xlabel('\theta_{w} [rad]');
    ylabel('y_{w} [mm]');
    zlabel('dr_{w} [mm]');
   else
    xlabel('u_{fc} [mm]');
    ylabel('y_{r} [mm]');
    zlabel('z_{r} [mm]');
   end
   set(gca,'dataaspectratio',myopt.zoom);

   % determine transverse curves along the surface at selected u (x) positions

   [ xcurv, ycurv, zcurv ] = cntc.make_3d_surface( slcs, myopt, myopt.urange, [], myopt.xysteps(1), ...
    [], [], []);
   % plot transverse curves along the surface

   l = plot3(xcurv', ycurv', zcurv', 'color',[.5 .5 .5], 'linewidth',0.5);

   [~,imin] = min(abs(xcurv(:,1)));
   set(l(imin), 'color', [.3 .3 .3], 'linewidth',1);

   % determine longitudinal curves along the surface at selected v (y) positions

   [ xcurv, ycurv, zcurv ] = cntc.make_3d_surface( slcs, myopt, ...
    myopt.urange, [], myopt.xysteps(1)/20, [], [], myopt.xysteps(2));

   % plot curves along surface in longitudinal direction

   clear l ip;
   l = plot3(xcurv, ycurv, zcurv, 'color',[.5 .5 .5], 'linewidth',0.5);
   imid = ceil(length(l)/2);
   set(l(imid), 'color', [.3 .3 .3], 'linewidth',1);

   % plot transverse curves at selected slices

   if (~isempty(myopt.show_slc))
    islc = unique( sort( max(1, min(slcs.nslc, myopt.show_slc)) ));
    nslc = length(islc);
    ncol = size(myopt.slc_color,1);
    jv   = [ max(1, round(myopt.vrange(1))) : min(slcs.npnt,round(myopt.vrange(2))) ];
    nv   = length(jv);
    for i = 1 : nslc
     jslc = islc(i);
     if (slcs.u(jslc)>=myopt.urange(1) & slcs.u(jslc)<=myopt.urange(end))
      icol = mod((i-1), ncol) + 1;
      col  = myopt.slc_color(icol,:);
      if (isscalar(col) & isnumeric(col)), col = cntc.matlab_color(col); end
      if (1==1)
       [ xcurv, ycurv, zcurv ] = cntc.make_3d_surface( slcs, myopt, slcs.u(jslc)*[1 1], [], 1, ...
        myopt.vrange, [], myopt.xysteps(2)/20);
      else
       xcurv = slcs.xsurf(jslc,jv); ycurv = slcs.ysurf(jslc,jv); zcurv = slcs.zsurf(jslc,jv);
      end
      plot3( xcurv, ycurv, zcurv, 'color',col, 'linewidth',1, 'Tag','slice');
     end
    end
   end

   % plot longitudinal curves at selected features

   if (~isempty(myopt.show_feat))
    ifeat = unique( sort( max(1, min(slcs.nfeat, myopt.show_feat)) ));
    nfeat = length(ifeat);
    ncol  = size(myopt.feat_color,1);
    for j = 1 : nfeat
     jp   = slcs.iseg_p( ifeat(j) );
     [ xcurv, ycurv, zcurv ] = cntc.make_3d_surface( slcs, myopt, ...
      myopt.urange, [], myopt.xysteps(1)/20, slcs.vj(jp)*[1 1], [], myopt.xysteps(2));
     jcol = mod((j-1), ncol) + 1;
     col  = myopt.feat_color(jcol,:);
     if (isscalar(col) & isnumeric(col)), col = cntc.matlab_color(col); end
     % plot3( slcs.xsurf(:,jp), slcs.ysurf(:,jp), slcs.zsurf(:,jp), 'color',col );
     plot3( xcurv, ycurv, zcurv, 'color',col, 'linewidth',1, 'Tag','feature');
    end
   end

  end % plot_2dspline



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ v_rfn ] = refine_intv( v_in, k_rfn )

   % refine intervals v_in(1) : v_in(2), v_in(2) : v_in(3), ... into k_rfn>=2 pieces each
   % #refine_intv -2

   n_in  = length( v_in );
   n_rfn = k_rfn * (n_in-1) + 1;
   v_rfn = zeros(n_rfn,1);
   f_rfn = [0 : 1/k_rfn : 1];

   for i_in = 1 : n_in-1
    i_out = k_rfn * (i_in-1) + [1 : k_rfn+1];
    dv_in = v_in(i_in+1) - v_in(i_in);
    v_rfn(i_out) = v_in(i_in) + f_rfn * dv_in;
   end
  end % refine_intv

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% #plot

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ h ] = plot_arrow( pos0, vec, col, scale, width, pln )
   %
   % [ h ] = plot_arrow( pos0, vec, col, scale, width, [pln] )
   %
   % pos0  = position of tail of arrow (may be 3-vector for plot3)
   % vec   = vector from tail to head of arrow (may be a 3-vector for plot3)
   % col   = color-spec for arrow; default='b';
   % scale = scale factor for arrow head; default=1: 25% of arrow size
   % width = width factor for arrow head; default=1: about 30 deg.
   % pln   = binormal direction used in 3d plot to define plane for arrow-head
   %

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
   % #plot_arrow

   if (nargin<3 | isempty(col))
    col = 'b';
   end
   if (nargin<4 | isempty(scale))
    scale = 1;
   end
   if (nargin<5 | isempty(width))
    width = 1;
   end
   if (nargin<6 | isempty(pln))
    pln = [ 0, 0, 1 ];
   end

   % if col is just a single value, it is interpreted as a color index

   if (isnumeric(col) & isscalar(col))
    col = cntc.matlab_color(col);
   end

   % convert arguments into column vectors

   if (size(pos0,2)>1), pos0=pos0'; end
   if (size(vec ,2)>1), vec =vec'; end
   if (size(pln ,2)>1), pln =pln'; end
   pln = pln / norm(pln);

   % use plot3 instead of plot when pos is a 3-vector

   if (length(pos0)~=length(vec))
    disp('plot_arrow: Inconsistent lengths');
    return;
   end
   is_3d = (length(pos0)==3);

   % adapt binormal pln if in direction of vec

   if (is_3d)
    vp = vec'*pln / norm(vec);
    if (abs(vp)>0.99)
     if (abs(pln'*[0;0;1])>0.99)
      pln = [0;1;0];
     else
      pln = [0;0;1];
     end
    end
   end

   % set line to be plotted

   lin  = [pos0, pos0+vec];

   % determine arrow head

   tng  = vec;
   if (is_3d)
    nrm  = cross(tng, pln);
   else
    asp  = get(gca,'dataaspectratio') ./ get(gca,'plotboxaspectratio');
    nrm  = [-tng(2)*asp(1)/asp(2); tng(1)*asp(2)/asp(1)];
   end

   head = [ pos0+vec-0.25*scale*tng+0.125*width*scale*nrm,  ...
    pos0+vec, ...
    pos0+vec-0.25*scale*tng-0.125*width*scale*nrm];

   if (~is_3d)
    h(1) = plot(lin(1,:), lin(2,:), 'color',col);
    hold on;
    h(2) = plot(head(1,:), head(2,:), 'color',col);
   else
    nl = size(lin,2); nh = size(head,2);
    h(1) = plot3(lin(1,:), lin(2,:), lin(3,:), 'color',col);
    hold on;
    h(2) = plot3(head(1,:), head(2,:), head(3,:), 'color',col);
   end

   %stel norm = [ 1, 0 ]
   %     aspc = [ 2, 1 ]
   % ==> transpose = [ 0, 1 ] = half zo lang.
  end

  function [ h ] = plot_axes(O, scl, theta, zpos, axnams, col, mrksiz,siz_o, fontsiz, txtcol, th_txt, fac_txt, no_orie, no_head, left_hnd)

   % [ h ] = plot_axes(O, scl, theta, zpos, axnams, col, mrksiz, ...
   %                            siz_o, fontsiz, txtcol, th_txt, fac_txt, ...
   %                            no_orie, no_head, left_hnd)
   %
   % plot coordinate axes Oxyz when looked upon from either +inf or -inf
   % at one coordinate axis.
   %
   % O        = position of origin on the current axes
   % scl      = length of vectors on the current axes
   % theta    = angle from Matlab x-axis to first coordinate axis [deg]
   % zpos     = height of viewpoint on 3rd axis
   % axnams   = labels to be put on axes x1,x2,x3
   %                axnams = 'xyz', zpos <0: plot Oxy viewed from z = -inf
   %                axnams = 'xyz', zpos>=0: plot Oxy viewed from z =  inf
   %                axnams = 'yzx', zpos <0: plot Oyz viewed from x = -inf
   %                axnams = 'zxy', zpos <0: plot Oxz viewed from y = -inf, theta=-90: z pointing upwards
   %                axnams = 'zyx', zpos>=0: plot Oxz viewed from y =  inf, theta=-90: z pointing downwards
   % col      = color-spec
   % mrksiz   = marker size (dot in circle for origin, default 8)
   % siz_o    = relative size of circle for origin (default 0.12)
   % fontsiz  = font size
   % txtcol   = color-spec for text labels
   % th_txt   = rotation angle [deg] for text labels (default theta)
   % fac_txt  = scale factor(s) for label positioning, can be 3x1 or 3x2
   % no_orie  = option to suppress sense of rotation
   % no_head  = option to suppress arrow heads
   % left_hnd = option to use left-handed sense of rotation
   %
   % output h(1) = marker; h(2:4) = lines center; h(5:8) = arrows; h(9:12) = sense of rotation; h(13:15) = labels

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
   % #plot_axes

   if (nargin<1 | isempty(O))
    O = [0; 0];
   end
   if (size(O,2)>size(O,1))
    O = O';
   end
   if (nargin<2 | isempty(scl))
    scl = 1;
   end
   if (nargin<3 | isempty(theta))
    theta = 0;
   end
   if (nargin<4 | isempty(zpos))
    zpos = 1;
   end
   if (nargin<5 | isempty(axnams))
    axnams = 'xyz';
   end
   if (size(axnams,1)==1)
    axnams = axnams';
   end
   if (nargin<6 | isempty(col))
    col = 'b';
   end
   if (nargin<7 | isempty(mrksiz))
    mrksiz = 8;
   end
   if (nargin<8 | isempty(siz_o))
    siz_o   = 0.12;
   end
   if (nargin<9 | isempty(fontsiz))
    fontsiz = 16;
   end
   if (nargin<10 | isempty(txtcol))
    txtcol = 'k';
   end
   if (nargin<11 | isempty(th_txt))
    th_txt = theta;
   end
   if (nargin<12 | isempty(fac_txt))
    fac_txt = 1;
   end
   if (nargin<13 | isempty(no_orie))
    no_orie = 0;
   end
   if (nargin<14 | isempty(no_head))
    no_head = 0;
   end
   if (nargin<15 | isempty(left_hnd))
    left_hnd = 0;
   end

   % associate ivec with matlab x, jvec with matlab y, kvec out of the screen

   % use left-handed convention if Matlab ydir is reversed

   gca_ydir = get(gca,'ydir');
   if (strcmp(gca_ydir,'reverse'))
    fac_y = -1;
   else
    fac_y =  1;
   end

   % reverse sense of rotation when left_hnd is used

   if (left_hnd)
    fac_y = -fac_y;
   end

   % reverse angles if z-position < 0

   if (zpos >= 0)
    fac_ang =  1;
   else
    fac_ang = -1;
   end
   theta  = fac_y * fac_ang * theta;
   th_txt =         fac_ang * th_txt;

   if (isnumeric(col) & isscalar(col))
    col = cntc.matlab_color(col);
   end
   if (isnumeric(txtcol) & isscalar(txtcol))
    txtcol = cntc.matlab_color(txtcol);
   end
   if (no_head)
    scl_head = 0.001;
   else
    scl_head = 1;
   end

   % expand fac_txt to 3x2 matrix
   if (size(fac_txt,2)==3), fac_txt = fac_txt'; end
   if (size(fac_txt,1)==1), fac_txt = ones(3,1)*fac_txt; end
   if (size(fac_txt,2)==1), fac_txt = fac_txt*ones(1,2); end

   % text labels: pos(1) * ivec + pos(2) * jvec
   pos_txt = scl * fac_txt .* [1.00, 0.25;   0.25, 1.00;   -0.10, -0.20];

   % compute unit vectors ivec, jvec in Matlab x,y-coordinates

   ivec = [ cos(theta*pi/180); sin(theta*pi/180) ];
   jvec = fac_y * fac_ang * [ -ivec(2); ivec(1) ];

   % initialize output-array: handles

   h = zeros(14,1);

   hold on;

   th=[0:10:360]*pi/180;
   sn=sin(th); cs=cos(th);

   % put circle with cross or dot at origin O

   r = siz_o*scl;
   h(4) = plot(O(1)+r*cs, O(2)+r*sn,'-','color',col);

   if (zpos>=0)
    % looking down x3-axis: dot
    h(1) = plot(O(1)+0*scl, O(2)+0*scl, '.','color',col, 'markersize',mrksiz);
   else
    % looking up x3-axis: cross
    ij = (ivec+jvec) / sqrt(2);
    h(2) = plot(O(1)+r*[-1,1]*ij(1), O(2)+r*[-1,1]*ij(2), 'color',col);
    ij = (ivec-jvec) / sqrt(2);
    h(3) = plot(O(1)+r*[-1,1]*ij(1), O(2)+r*[-1,1]*ij(2), 'color',col);
   end

   intrp = 'tex';
   if (strfind(reshape(axnams,1,prod(size(axnams))), '$'))
    intrp = 'latex';
   end

   % plot first axis ivec

   O1 = O + r*ivec;
   h(5:6) = cntc.plot_arrow(O1, (scl-r)*ivec, col, scl_head);

   l1 = O + pos_txt(1,1)*ivec + pos_txt(1,2)*jvec;
   h(13) = text(l1(1), l1(2), deblank(axnams(1,:)), 'rotation', th_txt, ...
    'horizontalalignment','center', 'interpreter', intrp, ...
    'fontsize',fontsiz, 'color', txtcol);

   % plot second axis jvec

   O2 = O + r*jvec;
   h(7:8) = cntc.plot_arrow(O2, (scl-r)*jvec, col, scl_head);

   l2 = O + pos_txt(2,1)*ivec + pos_txt(2,2)*jvec;
   h(14) = text(l2(1), l2(2), deblank(axnams(2,:)), 'rotation', th_txt, ...
    'horizontalalignment','center', 'interpreter', intrp, ...
    'fontsize',fontsiz, 'color', txtcol);

   % plot label for third axis kvec

   l3 = O + pos_txt(3,1)*ivec + pos_txt(3,2)*jvec;
   h(15) = text(l3(1), l3(2), deblank(axnams(3,:)), 'rotation',th_txt, ...
    'horizontalalignment','center', 'interpreter', intrp, ...
    'fontsize',fontsiz, 'color', txtcol);

   % plot optional sense of positive rotation

   if (no_orie<=0)
    Orot = O + 0.3*scl*ivec - 0.3*scl*jvec;
    f = 0.05*scl*[-1,1];
    h( 9) = plot(Orot(1)+ivec(1)*f, Orot(2)+ivec(2)*f, 'color',col);
    h(10) = plot(Orot(1)+jvec(1)*f, Orot(2)+jvec(2)*f, 'color',col);
    if (fac_ang*fac_y>0)
     h(11:12) = cntc.circ_arrow(Orot, 0.15*scl, theta-120, theta+150, col, ...
      0.7*scl_head, 0.8);
    else
     h(11:12) = cntc.circ_arrow(Orot, 0.15*scl, theta+120, theta-150, col, ...
      0.7*scl_head, 0.8);
    end
   end
  end

  function [ h ] = plot_axes3(O, scl, rot, axnams, col, mrksiz, fontsiz, txtcol, th_txt, fac_txt, no_head, left_hnd)

   % [ h ] = plot_axes3(O, scl, rot, axnams, col, mrksiz, ...
   %                            fontsiz, txtcol, th_txt, fac_txt, ...
   %                            no_head, left_hnd)
   %
   % 3d plot of coordinate axes Oxyz rotated by [roll,yaw,pitch].
   %
   % O        = position of origin on the current axes
   % scl      = length of vectors on the current axes
   % rot      = Euler angles [r,y,p] in [deg] or 3x3 rotation matrix
   % axnams   = labels to be put on axes x1,x2,x3
   % col      = color-spec
   % mrksiz   = marker size (dot in circle for origin, default 8)
   % fontsiz  = font size
   % txtcol   = color-spec for text labels
   % th_txt   = rotation angle [deg] for text labels
   % fac_txt  = scale factor(s) for label positioning, can be 3x1 or 3x3
   % no_head  = option to suppress arrow heads
   % left_hnd = option to use left-handed sense of rotation

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
   % #plot_axes3

   if (nargin<1 | isempty(O))
    O = [0; 0];
   end
   if (size(O,2)>size(O,1))
    O = O';
   end
   if (nargin<2 | isempty(scl))
    scl = 1;
   end
   if (nargin<3 | isempty(rot))
    rot = eye(3);
   end
   if (nargin<4 | isempty(axnams))
    axnams = 'xyz';
   end
   if (size(axnams,1)==1)
    axnams = axnams';
   end
   if (nargin<5 | isempty(col))
    col = 'b';
   end
   if (nargin<6 | isempty(mrksiz))
    mrksiz = 8;
   end
   if (nargin<7 | isempty(fontsiz))
    fontsiz = 16;
   end
   if (nargin<8 | isempty(txtcol))
    txtcol = 'k';
   end
   if (nargin<9 | isempty(th_txt))
    th_txt = 0;
   end
   if (nargin<10 | isempty(fac_txt))
    fac_txt = 1;
   end
   if (nargin<11 | isempty(no_head))
    no_head = 0;
   end
   if (nargin<12 | isempty(left_hnd))
    left_hnd = 0;
   end

   % associate ivec with matlab x, jvec with matlab y, kvec out of the screen

   % use left-handed convention if Matlab ydir is reversed

   gca_ydir = get(gca,'ydir');
   if (strcmp(gca_ydir,'reverse'))
    fac_y = -1;
   else
    fac_y =  1;
   end

   % reverse sense of rotation when left_hnd is used

   if (left_hnd)
    fac_y = -fac_y;
   end

   % convert color-spec to [r g b]

   if (isnumeric(col) & isscalar(col))
    col = cntc.matlab_color(col);
   end
   if (isnumeric(txtcol) & isscalar(txtcol))
    txtcol = cntc.matlab_color(txtcol);
   end

   % set size of arrow head

   if (no_head)
    scl_head = 0.001;
   else
    scl_head = 1;
   end

   % expand fac_txt to 3x2 matrix

   if (size(fac_txt,2)==3), fac_txt = fac_txt'; end
   if (size(fac_txt,1)==1), fac_txt = ones(3,1)*fac_txt; end
   if (size(fac_txt,2)==1), fac_txt = fac_txt*ones(1,3); end

   % text labels: pos(1) * ivec + pos(2) * jvec

   pos_txt = scl * fac_txt .* [1.00, 0.25, 0;   0.25, 1.00, 0;   0.20, 0.20, 1.00];

   % compute unit vectors ivec, jvec, kvec in Matlab x,y,z-coordinates

   if (any(size(rot)~=[3 3]))
    roll  = rot(1);
    yaw   = rot(2);
    pitch = rot(3);
    rot   = cntc.rotx(roll) * cntc.rotz(yaw) * cntc.roty(pitch);
   end
   ivec = rot(:,1);
   jvec = rot(:,2);
   kvec = rot(:,3);

   % initialize output-array: handles [dot, line,line,txt, line,line,txt, line,line,txt]

   h = zeros(10,1);

   % put dot at origin O

   h(1) = plot3(O(1), O(2), O(3), '.','color',col, 'markersize',mrksiz);
   hold on;

   intrp = 'tex';
   if (strfind(reshape(axnams,1,prod(size(axnams))), '$'))
    intrp = 'latex';
   end

   % plot first axis ivec

   h(2:3) = cntc.plot_arrow(O, scl*ivec, col, scl_head, [], jvec);

   l1 = O + pos_txt(1,1)*ivec + pos_txt(1,2)*jvec + pos_txt(1,3)*kvec;
   h(4) = text(l1(1), l1(2), l1(3), deblank(axnams(1,:)), 'rotation', th_txt, ...
    'horizontalalignment','center', 'interpreter', intrp, ...
    'fontsize',fontsiz, 'color', txtcol);

   % plot second axis jvec

   h(5:6) = cntc.plot_arrow(O, scl*jvec, col, scl_head, [], ivec);

   l2 = O + pos_txt(2,1)*ivec + pos_txt(2,2)*jvec + pos_txt(2,3)*kvec;
   h(7) = text(l2(1), l2(2), l2(3), deblank(axnams(2,:)), 'rotation', th_txt, ...
    'horizontalalignment','center', 'interpreter', intrp, ...
    'fontsize',fontsiz, 'color', txtcol);

   % plot third axis kvec

   h(8:9) = cntc.plot_arrow(O, scl*kvec, col, scl_head, [], jvec);

   l3 = O + pos_txt(3,1)*ivec + pos_txt(3,2)*jvec + pos_txt(3,3)*kvec;
   h(10) = text(l3(1), l3(2), l3(3), deblank(axnams(3,:)), 'rotation',th_txt, ...
    'horizontalalignment','center', 'interpreter', intrp, ...
    'fontsize',fontsiz, 'color', txtcol);

  end % plot_axes3




  function [ h ] = plot_cm(O, r, veclen, theta, col, no_vec, no_head)

   % [ h ] = plot_cm(O, r, veclen, theta, col, no_vec, no_head)
   %
   % plot a CM marker
   %  O       - position of CM wrt Matlab axes
   %  r       - radius of symbol
   %  veclen  - length of arrows - starting at r
   %  theta   - rotation wrt coordinate system [deg]
   %  col     - color
   %  no_vec  - disable vectors
   %  no_head - disable vector heads

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
   % #plot_cm

   if (nargin<1 | isempty(O))
    O = [0; 0];
   end
   if (size(O,2)>size(O,1))
    O = O';
   end
   if (nargin<2 | isempty(r))
    r = 1;
   end
   if (nargin<3 | isempty(veclen))
    veclen = 2*r;
   end
   if (nargin<4 | isempty(theta))
    theta = 0;
   end
   if (nargin<5 | isempty(col))
    col = 'b';
   end
   if (nargin<6 | isempty(no_vec))
    no_vec = 0;
   end
   if (nargin<7 | isempty(no_head))
    no_head = 0;
   end

   if (isnumeric(col) & isscalar(col))
    col = cntc.matlab_color(col);
   end
   if (no_head)
    scl_head = 0.001;
   else
    scl_head = 1;
   end

   h = zeros(7,1);

   hold on;

   % put circle at origin O

   th = [0:5:360]*pi/180;
   sn = sin(th); cs=cos(th);
   h(1) = plot(O(1)+r*cs, O(2)+r*sn,'-','color',col);

   % put filled quarter circles in 2nd & 4th quadrants

   th = (theta + [90:5:180])*pi/180;
   sn = sin(th); cs=cos(th);
   h(2) = fill(O(1)+[0,r*cs,0], O(2)+[0,r*sn,0], col, 'edgecolor','none');

   th = (theta + [270:5:360])*pi/180;
   sn = sin(th); cs=cos(th);
   h(3) = fill(O(1)+[0,r*cs,0], O(2)+[0,r*sn,0], col, 'edgecolor','none');


   if (~no_vec)
    % plot x-axis

    n = [ cos(theta*pi/180); sin(theta*pi/180)];
    O1 = O + r*n;
    h(4:5) = cntc.plot_arrow(O1, veclen*n, col, scl_head);

    % plot y-axis

    t = [ -n(2) ; n(1) ];
    O2 = O + r*t;
    h(6:7) = cntc.plot_arrow(O2, veclen*t, col, scl_head);
   end
  end

  function [ dst_max, l, n ] = plot_update(prf1, prf2, fac, plot_n, use_angl2, norm_dist, show_angl, ...
    add_plot, force_spline)

   % [ dst_max, l, n ] = plot_update(prf1, prf2, fac, plot_n, use_angl2, norm_dist, show_angl,
   %                                                       add_plot, force_spline)
   %
   % show difference of two profiles at magnification 'fac'
   %
   % fac       == magnification factor, negative/default: automatic setting
   % plot_n    == connect profile and update with normal lines at points [1 : plot_n : end]
   % use_angl2 == use surface angle of second profile
   % norm_dist == plot tangential updates (<0), full updates (0) or normal updates (>0)
   % show_angl == add markers at specified surface inclinations [ show_angl ] (deg)
   % add_plot  == clear current figure (0) or add plot to figure (1)
   % force_spline == assume different s parameterization even if #points is equal

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
   % #plot_update

   if (nargin<3 | isempty(fac))
    fac       = -1;
   end
   if (nargin<4 | isempty(plot_n))
    plot_n    = 5;
   end
   if (nargin<5 | isempty(use_angl2))
    use_angl2 = 0;
   end
   if (nargin<6 | isempty(norm_dist))
    norm_dist = 1;
   end
   if (nargin<7);
    show_angl = [];
   end
   if (nargin<8 | isempty(add_plot))
    add_plot  = 0;
   end
   if (nargin<9 | isempty(force_spline))
    force_spline = 0;
   end

   % get (y,z)-data for the two profiles

   y1_in = prf1.ProfileY;
   z1_in = prf1.ProfileZ;
   n1    = length(y1_in);
   if (min(size(y1_in))>1)
    disp('ERROR in plot_update: y1_in should be an n x 1 vector');
    return;
   end
   if (size(y1_in,1)<size(y1_in,2))
    y1_in = y1_in'; z1_in = z1_in';  % make column vectors
   end

   y2    = prf2.ProfileY;
   z2    = prf2.ProfileZ;
   n2    = length(y2);
   if (min(size(y2))>1)
    disp('ERROR in plot_update: y2 should be an n x 1 vector');
    return;
   end
   if (size(y2,1)<size(y2,2))
    y2 = y2'; z2 = z2';              % make column vectors
   end

   % different #points: resample using spline interpolation

   if (n1~=n2 | force_spline)
    % form spline for profile 1
    lambda = 0; use_bspline = 0; wgt = []; ikinks = []; iaccel = [];
    spl  = cntc.make_spline([], y1_in, z1_in, lambda, wgt, ikinks, iaccel, use_bspline);

    % determine coordinates s1 in spline 1 for points (y2,z2)
    % figure(5); clf; hold on;
    % plot(y1_in, z1_in, '-o');
    % plot(y2(774:776), z2(774:776), '*');
    % set(gca,'ydir','reverse');
    % grid on;
    % axis equal;
    % axis([-0.5 0.5 1.4 2.2]);

    idebug = -2; ifig = -5;
    si   = map_pts_to_spline(spl, y2, z2, idebug, ifig);

    % determine (y1,z1) at s1-positions
    [~, y1_int, z1_int] = eval_spline(spl, si);
   else
    y1_int = y1_in;
    z1_int = z1_in;
   end

   % tangent and inward normal on profile surface at positions used for prf2

   if (use_angl2)
    dy = diff(y2);
    dz = diff(z2);
   else
    dy = diff(y1_int);
    dz = diff(z1_int);
   end
   ds = sqrt(dy.^2 + dz.^2);

   ty = dy ./ ds;
   tz = dz ./ ds;
   ny = -tz;
   nz =  ty;

   % move normals from segments to end-points

   ty = [ty(1); (ty(1:end-1)+ty(2:end))/2; ty(end)];
   tz = [tz(1); (tz(1:end-1)+tz(2:end))/2; tz(end)];
   ny = [ny(1); (ny(1:end-1)+ny(2:end))/2; ny(end)];
   nz = [nz(1); (nz(1:end-1)+nz(2:end))/2; nz(end)];
   alph = atan2(tz, ty);

   % dst_n, dst_t: signed distance in normal/tangential direction between input profiles

   dst_y = (y2 - y1_int);
   dst_z = (z2 - z1_int);
   dst_n = ny .* dst_y + nz .* dst_z;
   dst_t = ty .* dst_y + tz .* dst_z;

   if (norm_dist>0)
    [dst_max, ixm] = max(abs(dst_n));
   elseif (norm_dist==0)
    [dst_max, ixm] = max(dst_y.^2 + dst_z.^2);
   else
    [dst_max, ixm] = max(abs(dst_t));
   end

   % fac<0: automatic scaling factor, fac * dst_max ~=~ 4% to 8% * tot_width

   if (fac<0)
    fac_list  = [ -1e-6 0.00005 0.0001 0.00025 0.0005 0.001 0.0025 0.005 0.01 0.025 0.05 0.1 0.25 0.5 1.0;
     1000   20000  10000    4000   2000  1000    400   200  100    40   20  10    4   2   1];
    tot_width = max(y1_in) - min(y1_in);
    ifac      = find(fac_list(1,:)<25*dst_max/tot_width, 1, 'last');
    fac       = fac_list(2,ifac);
   end

   % plot original profiles 1 + 2

   if (~add_plot)
    clf; hold on;
   end
   set(gca,'ydir','reverse');
   axis equal;
   grid on;

   l(1) = plot(y1_in, z1_in);
   l(2) = plot(y2, z2, '--');

   % plot magnified updates 2-1 and maximum update

   if (norm_dist>0)        % updates in normal direction
    l(3) = plot(y1_int+fac*dst_n.*ny, z1_int+fac*dst_n.*nz);
    plot(y1_int(ixm)+fac*dst_n(ixm)*ny(ixm), z1_int(ixm)+fac*dst_n(ixm)*nz(ixm), '*');
   elseif (norm_dist==0)   % full updates normal + tangential
    l(3) = plot(y1_int+fac*dst_y, z1_int+fac*dst_z);
    plot(y1_int(ixm)+fac*dst_y(ixm), z1_int(ixm)+fac*dst_z(ixm), '*');
   else                    % updates in tangential direction, normal to curve
    l(3) = plot(y1_int+fac*dst_t.*ny, z1_int+fac*dst_t.*nz);
    plot(y1_int(ixm)+fac*dst_t(ixm)*ny(ixm), z1_int(ixm)+fac*dst_t(ixm)*nz(ixm), '*');
   end

   % connect 1 with magnified updates at each profile point i

   n = []; j = 0;
   if (plot_n)
    if (norm_dist>0)
     for i = [1 : plot_n : length(dst_n)]
      if (fac*abs(dst_n(i))>0.01)
       j = j + 1;
       n(j) = plot(y1_int(i)+[0 fac*dst_n(i)*ny(i)], z1_int(i)+[0 fac*dst_n(i)*nz(i)], ...
        'color',cntc.matlab_color(6), 'linewidth',1);
      end
     end
    elseif (norm_dist==0)
     dst = sqrt(dst_y.^2+dst_z.^2);
     for i = [1 : plot_n : length(dst)]
      if (fac*dst(i)>0.01)
       j = j + 1;
       n(j) = plot(y1_int(i)+[0 fac*dst_y(i)], z1_int(i)+[0 fac*dst_z(i)], ...
        'color',cntc.matlab_color(4), 'linewidth',1);
      end
     end
    else
     for i = [1 : plot_n : length(dst_t)]
      if (fac*abs(dst_t(i))>0.01)
       j = j + 1;
       n(j) = plot(y1_int(i)+[0 fac*dst_t(i)*ny(i)], z1_int(i)+[0 fac*dst_t(i)*nz(i)], ...
        'color',cntc.matlab_color(5), 'linewidth',1);
      end
     end
    end
   end

   % pull smoothed profile to foreground

   if (1==1)
    c = get(gca,'children');
    c = [c(end-1); c(1:end-2); c(end)];
    set(gca,'children',c);
   end

   % plot angle indicators

   if (~isempty(show_angl))
    if (size(show_angl,1)>size(show_angl,2)), show_angl = show_angl'; end
    alph_deg = alph * 180/pi;
    for angl = show_angl
     % disp(sprintf('Show angle marker at %4d deg', angl));
     % take profile points at start of segments with zero crossing
     ix = find( (alph_deg(1:end-1)-angl).*(alph_deg(2:end)-angl) <= 0 );
     if (~isempty(ix))
      % disp(sprintf('  found at ix=%4d %4d %4d %4d',ix));
      for j = 1 : length(ix)
       yi = y1_int(ix(j));
       zi = z1_int(ix(j));
       nvec = [-sin(alph(ix(j))), cos(alph(ix(j)))]; % inward normal
       plot(yi+5*nvec(1)*[-1,1], zi+5*nvec(2)*[-1,1], '-','color',cntc.matlab_color(5),'linewidth',1);
       text(yi-5*nvec(1), zi-5*nvec(2), sprintf('$%d^\\circ$',angl), 'interpreter','latex', ...
        'fontsize',10, 'horizontalalignment','center', 'verticalalignment','bottom');
      end
     end
    end
   end

   % print biggest updates on flange

   if (0==1)
    ix = find( (y1_int>5 & y1_int<40) | (y1_int>125 & y1_int<138));
    for j = ix'
     if (abs(dst_n(j))==max(abs(dst_n([j-15:j+15]))))
      dst = fac * dst_n(j) + 1 * sign(dst_n(j));
      text(y1_int(j)+dst*ny(j), z1_int(j)+dst*nz(j), ...
       sprintf('%5.3f',dst_n(j)), 'horizontalalignment','center', 'verticalalignment','middle');
     end
    end
   end

   % show magnification factor

   if (1==1) % northeast / northwest
    text(0.95, 0.93, sprintf('$%d\\times$',round(fac)), 'units','normalized', ...
     'horizontalalignment','right', 'interpreter','latex');
    text(0.05, 0.93, sprintf('max update %5.3f mm',dst_max), 'units','normalized', ...
     'horizontalalignment','left', 'interpreter','latex');
   else % southeast / southwest
    text(0.95, 0.07, sprintf('$%d\\times$',round(fac)), 'units','normalized', ...
     'horizontalalignment','right', 'interpreter','latex');
    text(0.05, 0.07, sprintf('max update %5.3f mm',dst_max), 'units','normalized', ...
     'horizontalalignment','left', 'interpreter','latex');
   end
  end

  function [ opt_out ] = plotstrs(sol, opt)

   %
   % [ opt ] = plotstrs(sol, opt)
   %
   % Make a plot of quantities in the interiors of the contacting bodies.
   %
   %  sol   == structure with subsurface results as returned by loadstrs.
   %
   %  opt   == (input/output) structure with options for plotstrs.
   %           A valid structure opt is returned when plotstrs is called with
   %           no input arguments.
   %
   %   - block = selects a single block of subsurface points from the input
   %   - field = quantity to plot,
   %      'ux', 'uy', 'uz' for displacements,
   %      'sighyd' or 'hydro' for the mean hydrostatic stress sighyd = I_1/3,
   %      'sigvm' or 'mises' for the von Mises stress sigvm = sqrt(3*J_2),
   %      'sigma1','sigma2','sigma3' for the principal stresses,
   %      'sigtr' or 'tresca' for the maximum shear stress sigtr = sigma1-sigma3,
   %      'sigxx', 'sigxy', 'sigxz', 'sigyy', 'sigyz', 'sigzz' for the
   %               components of the stress tensor.
   %     Default: 'mises'.
   %   - dir   = orientation of selected slice, 'y' == plot 'Oxz'-plane (default)
   %   - yslc  = y-coordinate(s) for an Oxz-plot (dir='y'). Default: 0.0.
   %             Note: the grid column(s) closest to yslc is/are used.
   %             Setting 'max' shows the maximum value over all slices.
   %   - xslc, zslc : see yslc.
   %   - addplot = clear/do not clear the figure before plotting;
   %   - typplot = 'contour', 'contourf' (filled), 'surf'. Default: 'contourf'.
   %   - cntrlvl = values at which contours are required. Default: 'auto'.
   %   - clabel  = 'on' or [values]: add labels on contours
   %   - scale   = 'linear' or 'log'. Note: log-scale takes absolute value of
   %               the data, and works best with user-defined contour levels
   %               (cntrlvl). Default: 'linear'.
   %   - colormap = Matlab color-map to be used, e.g. 'cool', 'hot', 'jet', ...

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % determine version information on the Matlab graphics system
   % #plotstrs

   % construct a local struct "myopt" in which default values are filled in for all available options.

   myopt = struct( ...
    'field',    'mises', ...
    'dir',      'y', ...
    'xslc',     0., ...
    'yslc',     0., ...
    'zslc',     0., ...
    'addplot',  0, ... % clear plot (0) or add to exist.plot (1, experimental)
    'typplot',  'contourf', ...
    'cntrlvl',  'auto', ...
    'clabel',   'off', ...
    'scale',    'linear', ...
    'colormap', 'parula' ...
    );
   if (old_graphics)
    myopt.colormap = 'jet';
   end

   % If the user has not supplied any arguments,
   %    return default options as first and only output argument

   if (nargin<1)
    opt_out = myopt;
    return
   end

   % If the user has not supplied an opt-struct, use the default

   if (nargin<2 | isempty(opt))
    opt = myopt;
   end

   % Check whether user-supplied opt-struct contains unknown options

   useropts = fieldnames(opt);
   ierror = 0;
   for i = 1:length(useropts)
    if (~isfield(myopt, useropts{i}))
     ierror = ierror + 1;
     disp(sprintf('Unknown option "%s" will be ignored',useropts{i}));
     if (ierror==1)
      disp(sprintf('You may use "opt=rmfield(opt,''%s'');" to remove the option',useropts{i}));
     end
    end
   end

   % Overwrite all values in "myopt" with user-supplied values

   myopts = fieldnames(myopt);
   for i = 1:length(myopts)
    if (isfield(opt, myopts{i}) & ~isempty(cntc.mygetfield(opt,myopts{i})))
     myopt = setfield(myopt, myopts{i}, cntc.mygetfield(opt,myopts{i}));
    end
   end

   % Check whether requested field is valid

   ok_fields=strvcat('ux', 'uy', 'uz', 'sighyd', 'hydro', 'sigvm', 'mises', ...
    'sigma1', 'sigma2', 'sigma3', 'sigtr', 'tresca', ...
    'sigxx', 'sigxy', 'sigxz', 'sigyy', 'sigyz', 'sigzz');
   found=0;
   for istr=1:size(ok_fields,1)
    if (strcmp(myopt.field, deblank(ok_fields(istr,:))))
     found=1;
    end
   end
   if (~found)
    disp(sprintf('Unknown field requested="%s"; Available options:',myopt.field));
    disp(ok_fields)
    return;
   end

   % Check whether clabel is valid

   is_ok = (isempty(myopt.clabel) | strcmp(myopt.clabel,'auto') | strcmp(myopt.clabel,'on') | ...
    strcmp(myopt.clabel,'off') | isnumeric(myopt.clabel) );
   if (~is_ok)
    disp(sprintf('Unknown value clabel="%s"; Available options: on|off|auto or [values]',myopt.clabel));
    myopt.clabel = 'off';
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % start actual processing to produce the requested plot
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   if (myopt.addplot<=0)
    clf;
   end
   hold on;

   if (strcmp(myopt.field,'ux'))

    cntc.show_2d_slice(sol, 'ux', myopt);
    title (['Displacement U_x in X-direction', get(get(gca,'title'),'string')]);
    hc=colorbar; log_colorbar(hc, opt);
    ylabel(hc, 'Displacement U_x [mm]');

   end

   if (strcmp(myopt.field,'uy'))

    cntc.show_2d_slice(sol, 'uy', myopt);
    title (['Displacement U_y in Y-direction', get(get(gca,'title'),'string')]);
    hc=colorbar; log_colorbar(hc, opt);
    ylabel(hc, 'Displacement U_y [mm]');

   end

   if (strcmp(myopt.field,'uz'))

    cntc.show_2d_slice(sol, 'uz', myopt);
    title (['Displacement U_z in Z-direction', get(get(gca,'title'),'string')]);
    hc=colorbar; log_colorbar(hc, opt);
    ylabel(hc, 'Displacement U_z [mm]');

   end

   if (strcmp(myopt.field,'sighyd') | strcmp(myopt.field,'hydro'))

    cntc.show_2d_slice(sol, 'sighyd', myopt);
    title (['Mean hydrostatic stress SIGHYD', get(get(gca,'title'),'string')]);
    hc=colorbar; log_colorbar(hc, opt);
    ylabel(hc, 'stress invariant \sigma_{hyd} [N/mm^2]');

   end

   if (strcmp(myopt.field,'sigvm') | strcmp(myopt.field,'mises'))

    cntc.show_2d_slice(sol, 'sigvm', myopt);
    title (['Von Mises stress \sigma_{vm}', get(get(gca,'title'),'string')]);
    hc=colorbar; log_colorbar(hc, opt);
    ylabel(hc, 'Von Mises stress \sigma_{vm} [N/mm^2]');

   end

   if (strcmp(myopt.field,'sigxx') | strcmp(myopt.field,'sigxy') | strcmp(myopt.field,'sigxz') | ...
     strcmp(myopt.field,'sigyy') | strcmp(myopt.field,'sigyz') | strcmp(myopt.field,'sigzz') )

    if (~isfield(sol, myopt.field))
     disp(['ERROR: no data available for ',myopt.field,'.']);
    else
     cntc.show_2d_slice(sol, myopt.field, myopt);
     title (['Stress component ',myopt.field,' ',get(get(gca,'title'),'string')]);
     hc=colorbar; log_colorbar(hc, opt);
     ylabel(hc, ['Stress component ',myopt.field,' [N/mm^2]']);
    end
   end

   if (strcmp(myopt.field,'sigma1') | strcmp(myopt.field,'sigma2') | strcmp(myopt.field,'sigma3'))

    if (~isfield(sol, myopt.field))
     disp(['ERROR: no data available for ',myopt.field,'.']);
    else
     cntc.show_2d_slice(sol, myopt.field, myopt);
     title (['Principal stress \sigma_',myopt.field(6),' ',get(get(gca,'title'),'string')]);
     hc=colorbar; log_colorbar(hc, opt);
     ylabel(hc, ['principal stress \sigma_',myopt.field(6),' [N/mm^2]']);
    end
   end

   if (strcmp(myopt.field,'sigtr') | strcmp(myopt.field,'tresca'))

    if (~isfield(sol, 'sigtr'))
     disp(['ERROR: no data available for ',myopt.field,'.']);
    else
     cntc.show_2d_slice(sol, 'sigtr', myopt);
     title (['Maximum shear stress \sigma_{tresca}',get(get(gca,'title'),'string')]);
     hc=colorbar; log_colorbar(hc, opt);
     ylabel(hc, ['maximum shear stress \sigma_{tresca} [N/mm^2]']);
    end
   end

  end

  function [ rgb ] = matlab_color( num )

   % [ rgb ] = matlab_color( num )
   %
   % return the rgb-triplet for color order index num

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
   % #matlab_color( num ) -3

   matlab_colors = [
    0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
    ];
   num_colors = size(matlab_colors,1);

   if (nargin<1 | isempty(num) | ~isnumeric(num))
    num = 1;
   end

   % reshape to column vector
   num = reshape(num, prod(size(num)), 1);

   % clip
   num = mod(num-1, num_colors) + 1;

   % get rows
   rgb = matlab_colors(num, :);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ h ] = circ_arrow( pos_c, r, th0, th1, col, scale, width, rot )
   % [ h ] = circ_arrow( pos_c, r, th0, th1, col, scale, width, [rot] )
   %
   % pos_c = position of arrow center of rotation, 2-vector or 3-vector
   % r     = circle radius;  negative r: reversed rotation direction.
   %         vectorial: [rx,ry], separate radius for x and y-direction
   % th0   = angle of tail of arrow [deg]
   % th1   = angle of head of arrow
   % col   = color-spec for arrow; default='b';
   % scale = scale factor for arrow head; default=1: 25% of arrow size
   % width = width factor for arrow head; default=1: about 30 deg.
   % rot   = for 3d plots: rotation matrix from [x,y,0] to desired orientation
   %         [0,0,1; 1,0,0; 0,1,0] for arc in yz-plane
   %

   % Copyright 2008-2023 by Vtech CMCC.
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % #circ_arrow -2

   if (nargin<5 | isempty(col))
    col = 1;
   end
   if (nargin<6 | isempty(scale))
    scale = 1;
   end
   if (nargin<7 | isempty(width))
    width = 1;
   end
   if (nargin<8 | isempty(rot))
    rot   = eye(3);
   end

   % if r is a scalar, copy to make a 2-vector [rx,ry]
   if (max(size(r))==1)
    r = [r,r];
   end
   rx=r(1); ry=r(2);

   % if col is just a single value, it is interpreted as a color index
   if (isnumeric(col) & isscalar(col))
    col = cntc.matlab_color(col);
   end

   if (size(pos_c,2)>1), pos_c=pos_c'; end % make column vectors

   th0 = th0 * pi/180;             % convert to radians
   th1 = th1 * pi/180;             % convert to radians

   % use plot3 instead of plot when pos is a 3-vector

   is_3d = (length(pos_c)==3);
   if (is_3d)
    pos_c = rot' * pos_c;
   else
    pos_c = [pos_c; 0];
   end

   hold on;

   % compute circular arc with center pos_c, radius [rx,ry,0] and angles theta (th)
   th   = th0 + [0:30]/30 * (th1-th0);
   lin  = pos_c*ones(size(th)) + [rx*cos(th); ry*sin(th); 0*th];

   pos1 = lin(:,end);
   % plot(pos_c(1), pos_c(2), 'b*');

   % compute nice theta-value for getting the arrow head direction
   th_h = th1 - min(1,0.25*scale) * (th1-th0);

   % compute point on the curve for the back of the arrow head
   pos_h = pos_c + [rx*cos(th_h); ry*sin(th_h); 0*th_h];

   % compute tangent vector at th_h == orientation of head;
   % make size of arrowhead proportional to length of curve (th1-th0)
   tng  = [-ry*sin(th_h); rx*cos(th_h); 0] * (th1-th0);

   % compute the normal direction
   nrm  = [-tng(2); tng(1); 0];

   % reverse direction when needed
   if (th0>th1), nrm=-nrm; end

   % compute arrow head points, on the inside closer to the curve
   head = [ pos_h+0.10*width*scale*nrm,  ...
    pos1, ...
    pos_h-0.13*width*scale*nrm];

   if (is_3d)
    lin  = rot * lin;
    head = rot * head;
    h(1) = plot3(lin(1,:), lin(2,:), lin(3,:), 'color',col);
    h(2) = plot3(head(1,:), head(2,:), head(3,:), 'color',col);
   else
    h(1) = plot(lin(1,:), lin(2,:), 'color',col);
    h(2) = plot(head(1,:), head(2,:), 'color',col);
   end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ ] = show_2d_slice(sol, field, opt)
   %
   % [ ] = show_2d_slice(sol, field, opt)
   %
   % generic routine for a 2d contour plot of 3d data
   %   sol     - structure with subsurface results as returned by loadstrs
   %   field   - the actual array to be plotted, or a field name within sol
   %   opt     - structure with options as defined by plotstrs.
   %
   % #show_2d_slice -2

   % get field to be plotted

   nx  = sol.nx;
   ny  = sol.ny;
   nz  = sol.nz;

   if (ischar(field))
    eval(['tmp = sol.',field,';']);
   else
    tmp = field;
   end

   % get the data depending on the viewing direction

   if (strcmp(opt.dir,'x'))
    slc = opt.xslc; c0 = sol.x;
    dir1 = 'y'; n1 = ny; c1 = sol.y;  % to be plotted on 1-x-axis of graph
    dir2 = 'z'; n2 = nz; c2 = sol.z;  % to be plotted on 2-y-axis of graph
   elseif (strcmp(opt.dir,'y'))
    slc = opt.yslc; c0 = sol.y;
    dir1 = 'x'; n1 = nx; c1 = sol.x;
    dir2 = 'z'; n2 = nz; c2 = sol.z;
   elseif (strcmp(opt.dir,'z'))
    slc = opt.zslc; c0 = sol.z;
    dir1 = 'x'; n1 = nx; c1 = sol.x;
    dir2 = 'y'; n2 = ny; c2 = sol.y;
   else
    disp(sprintf('ERROR: illegal value for opt.dir=%s',opt.dir));
    return;
   end

   % check size for requested coordinate direction
   if ((n1==1 | n2==1) & ...
     (strcmp(opt.typplot,'contour') | strcmp(opt.typplot,'contourf')))
    disp(sprintf('ERROR: the block has %d x %d points per %c-slice;', n1, n2, opt.dir));
    disp(sprintf('       %s-plot not available for %c-slices.', opt.typplot, opt.dir));
    return
   end

   % determine requested slice number
   if (isempty(slc))
    disp(sprintf('ERROR: no %c-coordinate provided for %c%c-slice',opt.dir,dir1,dir2));
    return;
   end
   if (ischar(slc))
    if (strcmp(lower(slc),'max'))
     islc = -1;
    else
     disp(sprintf('ERROR: unknown option slc=''%c''',slc));
     return;
    end
   else
    [mn, islc] = min(abs(c0-slc));
    if (abs(c0-slc)>1e-4);
     disp(sprintf('Showing results for %c=%5.3f, nearest to %cslc=%5.3f',...
      opt.dir, c0(islc), opt.dir, slc));
    end
   end

   % get the data for the slice
   if (islc<=0)
    idir = double(lower(opt.dir)) - double('x') + 1;
    dta = squeeze( max(abs(tmp), [], idir ));
   else
    if (strcmp(opt.dir,'x'))
     dta = reshape( tmp(islc,:,:), n1,n2);
    elseif (strcmp(opt.dir,'y'))
     dta = reshape( tmp(:,islc,:), n1,n2);
    elseif (strcmp(opt.dir,'z'))
     dta = reshape( tmp(:,:,islc), n1,n2);
    end
   end

   if (~isempty(opt.cntrlvl) && ~strcmp(opt.cntrlvl,'auto'))
    cntrlvl = opt.cntrlvl;
   end

   % prepare for logarithmic color scale
   if (strcmp(opt.scale,'log'))
    dta     = log10(abs(dta));
    if (~isempty(opt.cntrlvl) && ~strcmp(opt.cntrlvl,'auto'))
     ix      = find(cntrlvl>0);
     cntrlvl = log10(cntrlvl(ix));
    end
   end

   % for surf-plot: add one column and row to the data
   dta = [dta, dta(:,end)];
   dta = [dta; dta(end,:)];

   % avoid constant data
   if (max(max(dta))-min(min(dta))<1e-20)
    dta(1,1) = dta(1,1)+1e-20;
   end

   % interpret c1, c2 as cell-centers, compute corners of the cells for surf-plot
   if (isscalar(c1))
    d1 = (c2(end)-c2(1))/n2;
    cornr1 = [ c1(1)-d1(1)/2; c1(end)+d1(end)/2];
   else
    d1 = c1(2:end) - c1(1:end-1); % step size
    cornr1 = [ c1(1)-d1(1)/2; c1(1:end-1)+d1(1:end)/2; c1(end)+d1(end)/2];
   end
   if (isscalar(c2))
    d2 = (c1(end)-c1(1))/n1;
    cornr2 = [ c2(1)-d2(1)/2; c2(end)+d2(end)/2];
   else
    d2 = c2(2:end) - c2(1:end-1); % step size
    cornr2 = [ c2(1)-d2(1)/2; c2(1:end-1)+d2(1:end)/2; c2(end)+d2(end)/2];
   end
   % do not extend with half a cell-size at z=0:
   if (strcmp(dir2,'z') & abs(c2(1))<1e-5), cornr2(1) = 0; end;
   if (strcmp(dir2,'z') & length(c2)>1 & abs(c2(end))<1e-5), cornr2(end) = 0; end;

   % plot contours or surf plot
   if (strcmp(opt.typplot,'contour'))
    if (isempty(opt.cntrlvl) || strcmp(opt.cntrlvl,'auto'))
     [C, h] = contour(c1, c2, dta(1:end-1,1:end-1)');
    else
     [C, h] = contour(c1, c2, dta(1:end-1,1:end-1)', cntrlvl);
    end
   elseif (strcmp(opt.typplot,'contourf'))
    if (isempty(opt.cntrlvl) || strcmp(opt.cntrlvl,'auto'))
     [C, h] = contourf(c1, c2, dta(1:end-1,1:end-1)');
    else
     [C, h] = contourf(c1, c2, dta(1:end-1,1:end-1)', cntrlvl);
    end
   else
    % size(cornr1), size(cornr2), size(dta)
    l1 = surf(cornr1, cornr2, dta');
    % shading faceted
    % set(l1,'linestyle','--');
    % l2=mesh(c1,c2,10*ones(size(dta)-1));
    % set(l2,'facecolor','none');
    % set(l2,'edgecolor','k');
    view([0 90]);
    axis([min(cornr1), max(cornr1), min(cornr2), max(cornr2)]);
   end

   if (isnumeric(opt.clabel) & (strcmp(opt.typplot,'contour') | strcmp(opt.typplot,'contourf')))
    clabel(C, h, opt.clabel);
   elseif (strcmp(opt.clabel,'on') & (strcmp(opt.typplot,'contour') | strcmp(opt.typplot,'contourf')))
    clabel(C, h);
   end

   if (~isempty(opt.colormap))
    colormap(opt.colormap);
   end

   xlabel(sprintf('%c_c [mm]',dir1));
   ylabel(sprintf('%c_c [mm]',dir2));
   % xlabel(sprintf('%c-coordinate [mm]',dir1));
   % ylabel(sprintf('%c-coordinate [mm]',dir2));
   if (islc>=1)
    title(sprintf(' at %c=%5.3f', opt.dir, c0(islc)));
   else
    title(sprintf(', max in %c', opt.dir));
   end

   % optionally adjust ticks

   % set(gca,'xtick', xcentr, 'xticklabel', num2str(xcentr', '%5.3f'));
   % set(gca,'ytick', ycentr, 'yticklabel', num2str(ycentr', '%5.3f'));

  end % show_2d_slice

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ ] = log_colorbar(hc, opt)
   %
   % [ ] = log_colorbar(hc, opt)
   %
   % change colorbar to logarithmic scale when necessary
   % #log_colorbar -2

   if (strcmp(opt.scale,'log'))
    if (isempty(opt.cntrlvl) || strcmp(opt.cntrlvl,'auto'))
     yt=get(hc,'ytick');
    else
     ix = find(opt.cntrlvl>0);
     yt = log10(opt.cntrlvl(ix));
    end
    set(hc,'ytick',yt,'yticklabel',10.^yt);
   end
  end % log_colorbar

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function [ dif ] = diffcase(sol1, sol2)

   %
   % [ dif ] = diffcase(sol1, sol2)
   %
   % Compute the difference of the results for two calculations sol1 and sol2.
   %
   % dif         - output struct "sol1 - sol2"
   % sol1, sol2  - structs with surface tractions as defined by loadcase.
   %

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % #diffcase -2
   dif = [];

   % Check equality of the structures w.r.t. grids used

   if (sol1.mx ~= sol2.mx || sol1.my ~= sol2.my)
    disp(sprintf('The two structs concern different grids: %dx%d, %dx%d',...
     sol1.mx, sol1.my, sol2.mx, sol2.my)); %#ok<*DSPSP>
    return
   end

   grd1 = [sol1.xl, sol1.yl, sol1.dx, sol1.dy];
   grd2 = [sol2.xl, sol2.yl, sol2.dx, sol2.dy];
   tol  = 0.01 * max(sol1.dx,sol1.dy);
   if (any( abs(grd1-grd2) > tol ))
    disp(sprintf(['The two structs concern different grids:\n',...
     '     grid 1: (xl,yl)=(%6.3f,%6.3f), (dx,dy)=(%6.3f,%6.3f)\n', ...
     '     grid 2: (xl,yl)=(%6.3f,%6.3f), (dx,dy)=(%6.3f,%6.3f)\n'], ...
     grd1, grd2));
    return
   end

   % Copy the administration from sol1

   dif.d_digit = sol1.d_digit;
   dif.mx      = sol1.mx;
   dif.my      = sol1.my;
   dif.xl      = sol1.xl;
   dif.yl      = sol1.yl;
   dif.dx      = sol1.dx;
   dif.dy      = sol1.dy;
   dif.x       = sol1.x;
   dif.y       = sol1.y;

   dif.kincns.t_digit = sol1.kincns.t_digit;
   dif.kincns.chi     = sol1.kincns.chi;
   dif.kincns.dq      = sol1.kincns.dq;
   dif.kincns.veloc   = sol1.kincns.veloc;

   dif.config  = sol1.config;
   dif.h_digit = sol1.h_digit;
   dif.mater   = sol1.mater;
   dif.fric    = sol1.fric;

   dif.x_offset  = sol1.x_offset;
   dif.y_offset  = sol1.y_offset;

   % copy element division from sol1 and encode the difference

   dif.eldiv = sol1.eldiv + 10*(sol1.eldiv-sol2.eldiv);

   % determine differences of solution variables

   dif.h      = sol1.h      - sol2.h;
   sol.mu     = sol1.mu     - sol2.mu;
   dif.pn     = sol1.pn     - sol2.pn;
   dif.px     = sol1.px     - sol2.px;
   dif.py     = sol1.py     - sol2.py;
   dif.un     = sol1.un     - sol2.un;
   dif.ux     = sol1.ux     - sol2.ux;
   dif.uy     = sol1.uy     - sol2.uy;
   dif.srel   = sol1.srel   - sol2.srel;
   dif.shft   = sol1.shft   - sol2.shft;
   dif.trcbnd = sol1.trcbnd - sol2.trcbnd;
   if (isfield(sol1,'taucrt') && isfield(sol2,'taucrt'))
    dif.taucrt = sol1.taucrt - sol2.taucrt;
    dif.uplsx  = sol1.uplsx  - sol2.uplsx;
    dif.uplsy  = sol1.uplsy  - sol2.uplsy;
   end
   if (isfield(sol1,'temp1') && isfield(sol2,'temp1'))
    dif.temp1 = sol1.temp1 - sol2.temp1;
    dif.temp2 = sol1.temp2 - sol2.temp2;
   end

   % determine the difference in magnitude pt

   dif.pt     = sqrt(sol1.px.^2 + sol1.py.^2) - sqrt(sol2.px.^2 + sol2.py.^2);
  end

  function [ dif ] = diffstrs(sol1, sol2)

   %
   % [ dif ] = diffstrs(sol1, sol2)
   %
   % Compute the difference of the results for two subsurface stress
   % calculations sol1 and sol2.
   %
   % dif         - output struct "sol1 - sol2"
   % sol1, sol2  - structs with subsurface stresses as defined by loadstrs.
   %

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % #diffstrs -2
   dif = [];

   % Check equality of the grids w.r.t. subsurface points used

   if (sol1.nx~=sol2.nx)
    disp(sprintf('The two structs concern different subsurface points (nx).'));
    return
   end
   if (sol1.ny~=sol2.ny)
    disp(sprintf('The two structs concern different subsurface points (ny).'));
    return
   end
   if (sol1.nz~=sol2.nz)
    disp(sprintf('The two structs concern different subsurface points (nz).'));
    return
   end

   %  check x,y,z-coordinates of block

   difc = max(abs(sol1.x-sol2.x));
   if (difc>1e-10)
    disp(sprintf('Warning: the x-coordinates differ by at most %f;\n         continuing with coords of first case', difc));
   end
   difc = max(abs(sol1.y-sol2.y));
   if (difc>1e-10)
    disp(sprintf('Warning: the y-coordinates differ by at most %f;\n         continuing with coords of first case', difc));
   end
   difc = max(abs(sol1.z-sol2.z));
   if (difc>1e-10)
    disp(sprintf('Warning: the z-coordinates differ by at most %f;\n         continuing with coords of first case', difc));
   end

   %  copy administration from sol1

   dif.nx      = sol1.nx;
   dif.ny      = sol1.ny;
   dif.nz      = sol1.nz;
   dif.npoints = sol1.npoints;
   dif.x       = sol1.x;
   dif.y       = sol1.y;
   dif.z       = sol1.z;

   %  compute differences of displacements and stress invariants

   dif.ux = sol1.ux - sol2.ux;
   dif.uy = sol1.uy - sol2.uy;
   dif.uz = sol1.uz - sol2.uz;

   dif.sighyd = sol1.sighyd - sol2.sighyd;
   dif.sigvm  = sqrt(sol1.sigvm.^2 - sol2.sigvm.^2);

   if (isfield(sol1,'sigxx') && isfield(sol2,'sigxx'))
    dif.sigxx = sol1.sigxx - sol2.sigxx;
    dif.sigxy = sol1.sigxy - sol2.sigxy;
    dif.sigxz = sol1.sigxz - sol2.sigxz;
    dif.sigyy = sol1.sigyy - sol2.sigyy;
    dif.sigyz = sol1.sigyz - sol2.sigyz;
    dif.sigzz = sol1.sigzz - sol2.sigzz;
   end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% #Read
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ p ] = read_miniprof(fname, idebug, make_plot)

   % [ p ] = read_miniprof(fname, [idebug], [make_plot])
   %
   % Lower-level routine for reading wheel/rail profiles in Miniprof format.
   % Does not do automatic corrections like mirroring or reversing order.
   %
   %   fname     - filename of input file
   %   make_plot - figure number to use for plotting the profile

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
   % #read_miniprof -1

   if (nargin<2 | isempty(idebug))
    idebug = 1;
   end
   if (nargin<3 | isempty(make_plot))
    make_plot = 0;
   end

   keep_keyword_values = 0;

   known_fields = struct( ...
    'AlarmFailureCount',             'int', ...
    'AlarmWarningCount',             'int', ...
    'AlignDA',                       'unknown', ...
    'AlignDX',                       'unknown', ...
    'AlignDY',                       'unknown', ...
    'AlignOX',                       'unknown', ...
    'AlignOY',                       'unknown', ...
    'AreaGain',                      'unknown', ...
    'AreaLoss',                      'unknown', ...
    'AreaRef',                       'unknown', ...
    'AxleNo',                        'int', ...
    'BaseOffset',                    'unknown', ...
    'BaseTilt',                      'unknown', ...
    'Battery_Capacity',              'real', ...
    'Battery_Level',                 'real', ...
    'Bogie',                         'int', ...
    'CarNo',                         'int', ...
    'Chainage',                      'str', ...
    'Chars',                         'str', ...
    'ColumnDef',                     'str', ...
    'Comment',                       'str', ...
    'CrownR',                        'real', ...
    'CrownRadius',                   'real', ...
    'Curve',                         'str', ...
    'Date',                          'str', ...
    'DiameterFlange',                'real', ...
    'DiameterFlange_AlarmStatus',    'int', ...
    'DiameterTaperline',             'real', ...
    'DiameterTaperline_AlarmStatus', 'int', ...
    'Direction',                     'str', ...
    'Elevation',                     'real', ...
    'EmplNum',                       'int', ...
    'FileName',                      'str', ...
    'Flag',                          'unknown', ...
    'Flat',                          'unknown', ...
    'Gauge',                         'real', ...
    'Gauge_AlarmStatus',             'int', ...
    'GaugeRodLength',                'real', ...
    'GaugeSensorDistance',           'real', ...
    'GaugeSensorSamples',            'real', ...
    'GaugeSensorStdDev',             'real', ...
    'Grade',                         'real', ...
    'GradeRatio',                    'real', ...
    'Instrument_Connection',         'int', ...
    'Instrument_Firmware',           'str', ...
    'KP',                            'unknown', ...
    'Line',                          'str', ...
    'Location_Altitude',             'real', ...
    'Location_HorizontalAccuracy',   'real', ...
    'Location_Latitude',             'real', ...
    'Location_Longitude',            'real', ...
    'Location_VerticalAccuracy',     'real', ...
    'Mileage',                       'str', ...
    'MPCalDate',                     'str', ...
    'MPCalTime',                     'str', ...
    'MPD1',                          'int', ...
    'MPD2',                          'int', ...
    'MPD3',                          'int', ...
    'MPD4',                          'int', ...
    'MPD5',                          'int', ...
    'MPD6',                          'int', ...
    'MPD7',                          'int', ...
    'MPD8',                          'int', ...
    'MPD9',                          'int', ...
    'MPD10',                         'int', ...
    'MPD11',                         'int', ...
    'MPD12',                         'int', ...
    'MPD13',                         'int', ...
    'MPD14',                         'int', ...
    'MPD15',                         'int', ...
    'MPD16',                         'int', ...
    'MPD17',                         'int', ...
    'MPD18',                         'int', ...
    'MPD19',                         'int', ...
    'MPDC',                          'int', ...
    'MPDN',                          'str', ...
    'MPSerNo',                       'int', ...
    'MPTypeNo',                      'int', ...
    'OriginalHeader',                'unknown', ...
    'Position',                      'str', ...
    'ProfileAlignment',              'unknown', ...
    'ProfileShaping',                'unknown', ...
    'ProgramDate',                   'str', ...
    'ProgramName',                   'str', ...
    'ProgramVer',                    'str', ...
    'qR',                            'real', ...
    'qR_AlarmStatus',                'int', ...
    'Rail',                          'str', ...
    'RailAlignType',                 'int', ...
    'RailAngle',                     'real', ...
    'RailAngle_AlarmStatus',         'int', ...
    'RailGaugePoint',                'str', ...
    'RailRodLength',                 'str', ...
    'RCF',                           'unknown', ...
    'ReferenceProfile',              'str', ...
    'RefPoint1',                     'real', ...
    'RefPoint2',                     'real', ...
    'RefPoint3',                     'real', ...
    'RefPoint4',                     'real', ...
    'RefPoint5',                     'real', ...
    'RefPoint6',                     'real', ...
    'RefPoint7',                     'real', ...
    'RefPoint8',                     'real', ...
    'RefPoint9',                     'real', ...
    'RodLength',                     'real', ...
    'Sd',                            'real', ...
    'Sd_AlarmStatus',                'int', ...
    'SectionNo',                     'unknown', ...
    'Sh',                            'real', ...
    'Sh_AlarmStatus',                'int', ...
    'SuperElevation',                'real', ...
    'SuperElevationHeight',          'real', ...
    'Side',                          'str', ...
    'Stock',                         'str', ...
    'SurfaceCondition',              'unknown', ...
    'Temperature',                   'real', ...
    'Tilt',                          'unknown', ...
    'Time',                          'str', ...
    'Track',                         'str', ...
    'TrackSection',                  'unknown', ...
    'Transformation',                'str', ...
    'TransformState',                'str', ...
    'UOF1',                          'str', ...
    'UOF2',                          'str', ...
    'UOF3',                          'str', ...
    'UOF4',                          'str', ...
    'UserName',                      'str', ...
    'W1',                            'real', ...
    'W1_AlarmLevel',                 'real', ...
    'W1_AlarmStatus',                'int', ...
    'W2',                            'real', ...
    'W2_AlarmStatus',                'int', ...
    'W3',                            'real', ...
    'W3_AlarmStatus',                'int', ...
    'WheelAdjustPoint',              'real', ...
    'WheelDiameterTaperline',        'real', ...
    'WheelID',                       'int', ...
    'Windows_ComputerName',          'str', ...
    'Windows_NetFramework',          'str', ...
    'Windows_ServicePack',           'str', ...
    'Windows_UserName',              'str', ...
    'Windows_Version',               'str', ...
    'Xoffset',                       'real', ...
    'Yoffset',                       'real', ...
    'XYPoints',                      'int'  ...
    );


   % Initialize output structure

   p                  = struct;
   p.ProfileData      = [];
   nrows_data         = 0;
   p.ProfileS         = [];
   p.ProfileY         = [];
   p.ProfileZ         = [];
   p.ProfileAngle     = [];
   p.ProfileCurvature = [];
   p.ProfileKYield    = [];

   if (~exist(fname,'file'))
    disp(['ERROR: cannot find file ',fname]);
    return
   end
   f = fopen(fname,'r');
   iline = 0;
   in_profile = 0;

   FORMAT_UNKNOWN = 0;
   FORMAT_OLD = 1;
   FORMAT_BAN = 2; % .ban and .whl files
   FORMAT_MPT = 3; % .mpt files supporting multiple profiles or profile sections
   file_format = FORMAT_UNKNOWN;

   while(~feof(f))
    % get next line from input
    s  = deblank(fgets(f));
    iline = iline + 1;
    if (idebug>=5 & iline<10)
     disp(sprintf('%4d: %s', iline, s));
    end

    % determine type of Miniprof file
    if (iline<=3 & file_format==FORMAT_UNKNOWN)
     if (strfind(s, 'MINIPROF'))
      file_format = FORMAT_OLD;
      if (idebug>=1), disp('Old style Miniprof file'); end
     elseif (strfind(s, '[FileInfo]'))
      file_format = FORMAT_MPT;
      if (idebug>=1), disp('Multi-part Miniprof file'); end
     end
    end
    if (iline>3 & file_format==FORMAT_UNKNOWN)
     file_format = FORMAT_BAN;
     if (idebug>=1), disp('Regular Miniprof file'); end
    end

    % remove comments, everything after ", !, or %
    ix = strfind(s, '"'); if (~isempty(ix)), s=s(1:ix-1); end
    ix = strfind(s, '!'); if (~isempty(ix)), s=s(1:ix-1); end
    ix = strfind(s, '%'); if (~isempty(ix)), s=s(1:ix-1); end

    if (~isempty(s) & s(1)=='[' & idebug>=1)
     disp(sprintf('Line %5d: start of section %s', iline, strtrim(s)));
    end

    % OLD, BAN-formats: data starts after empty line

    if (any(file_format==[FORMAT_OLD, FORMAT_BAN]) & isempty(strtrim(s)))

     in_profile = 1;

     % MPT-format: next profile section starts after [Profile*Points], add row of NaNs between sections

    elseif (file_format==FORMAT_MPT & strfind(s, 'Points]'))

     in_profile = in_profile + 1;
     if (idebug>=2)
      disp(sprintf('Line %5d: starting profile section %d', iline, in_profile))
     end
     if (in_profile > 1)
      nrows_data = nrows_data + 1;
      p.ProfileData(nrows_data,:) = NaN;
     end
    end

    % look for an equals sign
    ix = strfind(s,'=');

    % if found: interpret keyword=value syntax
    if (~isempty(ix))
     ix = ix(1);
     field_name = deblank(s(1:ix-1));
     field_name = strrep(field_name, '.', '_');
     field_name = strrep(field_name, ' ', '_');
     field_val  = deblank(s(ix+1:end));
     % check if the keyword is known - convert to case-sensitive name
     my_field_name = cntc.check_field(known_fields, field_name);
     if (~keep_keyword_values)
      % do not store file header keyword/value information
     elseif (~isempty(my_field_name))
      % interpret numerical value
      if (strcmp(known_fields.(my_field_name), 'int'))
       p.(my_field_name) = str2num(field_val);
      elseif (strcmp(known_fields.(my_field_name), 'real'))
       p.(my_field_name) = str2num(field_val);
      elseif (strcmp(known_fields.(my_field_name), 'num'))
       p.(my_field_name) = str2num(field_val);
       % store string value
      elseif (strcmp(known_fields.(my_field_name), 'str'))
       p.(my_field_name) = field_val;
       % unknown type: store string value if it's non-empty
       %               empty fields of type 'unknown' will be discarded
      elseif (~isempty(field_val))
       if (idebug>=2)
        disp(sprintf('Unknown field-type %s: value="%s"',field_name,field_val));
       end
       p.(my_field_name) = field_val;
      end
      % unknown keyword: report & store string value
     else
      if (idebug>=2)
       disp(sprintf('Unknown field "%s"',field_name));
      end
      p.(field_name) = field_val;
     end

     % no keyword found - interpret data line
    elseif (~isempty(s))

     ix = strfind(s, ','); if (~isempty(ix)), s(ix)=' '; end
     [tmp, cnt] = sscanf(s, '%f %f %f %f %f %f %f %f');

     if (cnt>0 & (in_profile>0 | file_format~=FORMAT_OLD))
      nrows_data = nrows_data + 1;
      p.ProfileData(nrows_data, 1:cnt) = tmp;
     end

    end % ~isempty(ix), keyword handling
   end % while (~feof)
   fclose(f);

   if (~isfield(p, 'ColumnDef'))
    p.ColumnDef = 'X,Y,*';
    if (idebug>=2)
     disp(sprintf('Assuming columndef "%s"', p.ColumnDef));
    end
   end

   ncol = size(p.ProfileData,2);
   KnownCols = ['X';'Y';'A';'C';'N';'K'];
   for icol = 1 : length(KnownCols)
    ipos = (findstr(KnownCols(icol), p.ColumnDef) + 1) / 2;
    if (ipos)
     if (idebug>=3)
      disp(sprintf('Found column "%s" at ipos = %d', KnownCols(icol), ipos));
     end
     if     (KnownCols(icol)=='X')
      p.ProfileY         = p.ProfileData(:,ipos);
     elseif (KnownCols(icol)=='Y')
      p.ProfileZ         = p.ProfileData(:,ipos);
     elseif (KnownCols(icol)=='A')
      p.ProfileAngle     = p.ProfileData(:,ipos);
     elseif (KnownCols(icol)=='C')
      p.ProfileCurvature = p.ProfileData(:,ipos);
     elseif (KnownCols(icol)=='N')
      % ignored
     elseif (KnownCols(icol)=='K')
      p.ProfileKYield    = p.ProfileData(:,ipos);
     else
      disp(sprintf('Unknown column icol=%d, "%s"', icol, KnownCols(icol)));
     end
    elseif (icol<=2)
     disp(sprintf('Invalid ColumnDef="%s", should specify at least X,Y', p.ColumnDef));
    end
   end

   if (make_plot)
    figure(make_plot); clf; hold on;
    plot(p.ProfileY, p.ProfileZ, '-o');
    grid on;
    xlabel('y_{prf} [mm]'); ylabel('z_{prf} [mm]');
    set(gca,'ydir','reverse');
   end

  end % read_miniprof

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ my_field ] = check_field( known_fields, field_in )
   % #check_field -2

   fn = fieldnames(known_fields);
   chk = strcmpi(field_in, fn);
   i = find(chk);
   if (isempty(i))
    my_field = [];
   else
    my_field = fn{i};
   end

  end % check_field

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ p ] = read_profile(fname, is_wheel, mirror_y, mirror_z, scale_yz, rgt_side, idebug, make_plot)

   % [ p ] = read_profile(fname, [is_wheel], [mirror_y], [mirror_z], [scale_yz], [rgt_side], ...
   %                                                                          [idebug], [make_plot])
   %
   % Main routine for reading wheel/rail profiles.
   % Switches between Slices, SIMPACK, MiniProf and Vampire files using the filename extension.
   % Uses modify_profiles to apply automatic corrections.
   %
   % is_wheel    -  0 for rail, 1 for wheel
   % mirror_y    - -1 or 0 for no, 1 for yes
   % mirror_z    - -1 for no, 0 for automatic, 1 for yes
   % scale_yz    - scale-factor to convert to mm, e.g. 1000. for data in meters
   % rgt_side    - Vampire format: select right side (1, default) or left side (0)
   % make_plot   - <=0 for no, fig.number for yes

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
   % #read_profile

   % parameters
   if nargin<3&&ischar(fname)
    % package SDT options
    DoOpt=['iswheel(#%g#"0 for rail, 1 for wheel")' ...
     'fname(#%s#"Name of the file")' ...
     'mirrory(0#%g#"Do the symmetry w.r.t the y axis")' ...
     'mirrorz(0#%g#"Do the symmetry w.r.t the z axis")' ...
     'sclfac(0#%g#"scale-factor to convert to mm, e.g. 1000. for data in meters")' ...
     'rgt_side(0#%g#"Vampire format: select right side (1, default) or left side (0)")' ...
     'make_plot(0#%g#"<=0 for no, fig.number for yes")'];

    if nargin==0; CAM='';else; CAM=fname;end
    LI=cntc.call;
    LI.ReadProfile=cingui('paramedit -doclean2',DoOpt,{struct,CAM});

    fname=sdth.sfield('addselected',struct,LI.ReadProfile,{'fname'});
    fname=struct2cell(fname);fname=fname{1};
    params={'iswheel','mirrory','mirrorz','sclfac','rgt_side','idebug','make_plot'};
    params=sdth.sfield('addselected',struct,LI.ReadProfile,params);
    params=struct2cell(params);%params=horzcat(params{:});

    p=cntc.read_profile(fname,params{:});
    if nargout==0
     [~,fname,ext]=fileparts(p.Fname);
     C1=struct('X',{{p.ProfileData(:,1)/p.spline.units_len_f, ...
      {'Z','m',struct('DispUnit','mm','coef',p.spline.units_len_f)}}}, ...
      'Xlab', ...
      {{{'Y','m',struct('DispUnit','mm','coef',p.spline.units_len_f)},'val'}}, ...
      'Y',p.ProfileData(:,2)/p.spline.units_len_f,'Ylab',2,'name',[fname ext]);
     cdm.urnVec(C1,'{x,1}{y,1}{gf10}');
    end
    return

   end
   % end SDT packaging

   %% read_profile(fname, is_wheel, mirror_y, mirror_z, scale_yz, rgt_side, idebug, make_plot)
   LI=cntc.call;
   if ~exist(fname,'file'); fname=fullfile(LI.ProjectWd,fname);end
   fname=sdtu.f.safe(fname);

   [~,~,ext] = fileparts(fname);

   if (nargin<2 | isempty(is_wheel))
    is_wheel = (strcmp(lower(ext), '.slcw') | strcmp(lower(ext), '.prw') | ...
     strcmp(lower(ext), '.whe') | strcmp(lower(ext), '.whl')); % default: rail
   end
   if (nargin<3 | isempty(mirror_y))
    mirror_y = 0;
   end
   if (nargin<4 | isempty(mirror_z))
    mirror_z = 0;
   end
   if (nargin<5 | isempty(scale_yz))
    scale_yz = 1;
   end
   if (scale_yz<1e-11)
    disp('ERROR(read_profile): scale_yz must be > 0');
    return;
   end
   if (nargin<6 | isempty(rgt_side))
    rgt_side = 1;
   end
   if (nargin<7 | isempty(idebug))
    idebug = 0;
   end
   if (nargin<8 | isempty(make_plot))
    make_plot = 0;
   end

   % switch to appropriate routine based on filename extension

   is_slices  = 0;
   is_simpack = 0;
   is_vampire = 0;

   if (any(strcmpi(ext, {'.slcs', '.slcw'}) ))

    p = cntc.read_slices(fname, mirror_y, mirror_z, scale_yz, idebug);
    is_slices  = 1;

   elseif (any(strcmpi(ext, {'.prr', '.prw'}) ))

    p = cntc.read_simpack(fname, idebug, make_plot);
    is_simpack = 1;

   elseif (any(strcmpi(ext, {'.rai', '.whe'}) ))

    p = cntc.read_vampire(fname, rgt_side, idebug, make_plot);
    is_vampire = 1;

   elseif (any(strcmpi(ext, {'.ban', '.whl'}) ))

    p = cntc.read_miniprof(fname, idebug, make_plot);

   else

    if (idebug>=1)
     disp(sprintf('Unknown file extension "%s", trying 2-column format.', ext))
    end
    p = cntc.read_miniprof(fname, idebug, make_plot);

   end

   if (isempty(p) | (~is_slices & ~isfield(p, 'ProfileY')) | (~is_slices & isempty(p.ProfileY)))
    if (idebug>=1)
     disp('No profile, returning');
    end
    return
   end

   p.is_wheel = is_wheel;

   % make automatic adjustments: mirroring, swapping

   if (~is_slices)

    y = p.ProfileY; ny = length(y);

    if (is_simpack) %#ok<IFBDUP>
     delete_foldback = 0;
    else
     delete_foldback = 0;
    end

    if (mirror_z==0) % automatic mode
     if (is_simpack)
      mirror_z = -1;
     else
      if (is_wheel)
       % mirror z if lowest z-value occurs in interior (flange)
       [minz,ix0]  = min(p.ProfileZ);
       [maxz,ix1]  = max(p.ProfileZ);
       if (ix0<0.05*ny | ix0>0.95*ny)
        mirror_z = -1;  % lowest z at start or end
       elseif (ix1<0.05*ny | ix1>0.95*ny)
        mirror_z =  1;  % highest z at start or end
       else
        ix_m = round(mean([ix0 ix1]));
        z_m = mean([minz, maxz]);
        mirror_z = (p.ProfileZ(ix_m) > z_m); % mirror when profile z > straight line
       end
      else
       % mirror z if lowest z-value does not occur on the tread of the rail
       [~,ix]   = min(p.ProfileZ);
       mirror_z = (ix<0.1*ny | ix>0.9*ny);
      end
     end
    end

    if (is_wheel)
     % reverse order if first y-value < last y-value, after mirroring
     reverse_order = ( (mirror_y<=0 & y(1)<y(end)) | (mirror_y>=1 & y(1)>y(end)) );
    else
     % reverse order if last y-value < first y-value, after mirroring
     reverse_order = ( (mirror_y<=0 & y(end)<y(1)) | (mirror_y>=1 & y(end)>y(1)) );
    end

    if (scale_yz~=1)
     point_dist_min = 1e-4 / scale_yz;
    elseif (max(abs(p.ProfileY))>1) % data in [mm]
     point_dist_min = 1e-4; % mm
    else
     point_dist_min = 1e-7; % m
    end

    p = cntc.modify_profile(p, is_wheel, mirror_y, mirror_z, reverse_order, point_dist_min, ...
     delete_foldback, scale_yz, make_plot, idebug);

    % check for big gaps in the data

    dy = diff(p.ProfileY);
    dz = diff(p.ProfileZ);
    ds = sqrt(dy.^2 + dz.^2);

    if (idebug>=1 & max(ds) > 20*mean(ds))
     disp(sprintf('Warning: max(ds)=%6.2f >> avg(ds)=%6.2f, missing data?', max(ds), mean(ds)));
    end

    % compute s-coordinate along profile

    p.ProfileS = cntc.make_arclength(p.ProfileY, p.ProfileZ);
   end % not slices

   p.Fname = fname;

  end % read_profile

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ p_out ] = modify_profile(p_in, is_wheel, mirror_y, mirror_z, reverse_order, point_dist_min,delete_foldback, scale_yz, make_plot, idebug)

   % [ p_out ] = modify_profile(p_in, is_wheel, mirror_y, mirror_z, reverse_order, point_dist_min, ...
   %                                     delete_foldback, scale_yz, make_plot, idebug)
   %
   % apply a number of common filters to the input profile p_in
   %   - mirror_y         - change y(i)  --> -y(i)    <=0 means no, >=1 yes
   %   - mirror_z         - change z(i)  --> -z(i)    <=0 means no, >=1 yes
   %   - reverse_order    - change [1:n] --> [n:-1:1]
   %   - point_dist_min   - delete points that lie too close together
   %   - delete_foldback  - identify & remove points where y(i+1)<y(i) (rail) or y(i+1)>y(i) (wheel)
   %   - scale_yz         - scale y(i), z(i) --> scale*y(i), scale*z(i)

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % #modify_profile -2

   if (nargin<2 | isempty(is_wheel))
    is_wheel = 0;
   end
   if (nargin<3 | isempty(mirror_y))
    mirror_y = 0;
   end
   if (nargin<4 | isempty(mirror_z))
    mirror_z = 0;
   end
   if (nargin<5 | isempty(reverse_order))
    reverse_order = 0;
   end
   if (nargin<6 | isempty(point_dist_min))
    point_dist_min = -1;
   end
   if (nargin<7 | isempty(delete_foldback))
    delete_foldback = 0;
   end
   if (nargin<8 | isempty(scale_yz))
    scale_yz = 1;
   end
   if (nargin<9 | isempty(make_plot))
    make_plot = 0;
   end
   if (nargin<10 | isempty(idebug))
    idebug = 1;
   end

   % Initialize output structure

   p_out          = p_in;
   if (~isfield(p_out, 'ProfileAngle') | isempty(p_out.ProfileAngle))
    p_out.ProfileAngle = zeros(size(p_out.ProfileY));
   end
   if (~isfield(p_out, 'ProfileCurvature') | isempty(p_out.ProfileCurvature))
    p_out.ProfileCurvature = zeros(size(p_out.ProfileY));
   end
   if (~isfield(p_out, 'ProfileData'))
    p_out.ProfileData = [ p_out.ProfileY, p_out.ProfileZ, ...
     p_out.ProfileAngle, p_out.ProfileCurvature];
   end
   p_out.OrigData = p_out.ProfileData;

   % Apply modifications to original input data

   % modification 1: mirror horizontally w.r.t. y=0, changing y into -y

   if (mirror_y>=1)
    if (idebug>=1)
     disp('Mirroring y-values from left-side to right-side profile...');
    end
    p_out.ProfileY     = -p_out.ProfileY;
    p_out.ProfileAngle = -p_out.ProfileAngle;
   end

   % modification 2: mirror vertically w.r.t. z=0, changing z into -z

   if (mirror_z>0)
    if (idebug>=1)
     disp('Mirroring z-values to get z positive downwards...');
    end
    p_out.ProfileZ     = -p_out.ProfileZ;
    p_out.ProfileAngle = -p_out.ProfileAngle;
   end

   % modification 3: reverse the order of the data points

   if (reverse_order)
    if (idebug>=1)
     disp('Reversing the order of the points, get inside at right hand...');
    end
    p_out.ProfileY         = flipud(p_out.ProfileY);
    p_out.ProfileZ         = flipud(p_out.ProfileZ);
    p_out.ProfileAngle     = flipud(p_out.ProfileAngle);
    p_out.ProfileCurvature = flipud(p_out.ProfileCurvature);
   end

   % modification 4: remove points that lie too close together

   if (point_dist_min>=1e-6)
    np   = length(p_out.ProfileY);
    dy   = diff(p_out.ProfileY);
    dz   = diff(p_out.ProfileZ);
    dist = sqrt(dy.^2 + dz.^2);
    cum  = [0; cumsum(dist)];

    keep = ones(np,1);
    ik   = 1; % highest number that is kept so far
    for ip = 2 : np
     if (cum(ip)-cum(ik) < point_dist_min)
      keep(ip) = 0;
     else
      ik = ip;
     end
    end
    ix0 = find(keep==0); ix1 = find(keep==1);

    if (idebug>=1 & ~isempty(ix0))
     for i = ix0'
      disp(sprintf('Deleting (%7.2f,%7.2f) (ip=%4d), too close to adjacent', ...
       p_out.ProfileY(i), p_out.ProfileZ(i), i));
     end
    end
    if ((idebug>=1 & ~isempty(ix0)) | idebug>=2)
     disp(sprintf('Point_dist_min=%6.0e: deleted %d points from profile, %d remaining', point_dist_min, ...
      length(ix0), length(ix1)));
    end

    p_out.ProfileY         = p_out.ProfileY(ix1);
    p_out.ProfileZ         = p_out.ProfileZ(ix1);
    p_out.ProfileAngle     = p_out.ProfileAngle(ix1);
    p_out.ProfileCurvature = p_out.ProfileCurvature(ix1);
   end

   % modification 5: remove points where the profile folds back
   %                 (i.e. where the function y --> z(y) is multi-valued)

   if (~delete_foldback)
    p_out.Mask             = ones(size(p_out.ProfileY));
   else
    if (idebug>=2)
     disp('Checking for points where the profile folds back...');
    end

    % copy data to work variables
    points0 = [p_out.ProfileY, p_out.ProfileZ, p_out.ProfileAngle, ...
     p_out.ProfileCurvature];
    np0 = size(points0,1);

    % find "midpoint" of profile: tread of rail, bottom of flange on wheel
    % z is assumed positive downwards: tread == minimum z, flange == maximum z
    if (is_wheel)
     [~,i_mid0] = max(points0(:,2));
    else
     [~,i_mid0] = min(points0(:,2));
    end

    % set threshold "dy" below which y-values are "the same"
    dy_thresh = 1e-6 * (max(points0(:,1))-min(points0(:,1)));

    % delete points where the profile is "vertical", i.e. |y(i+1)-y(i)|<thresh
    % where (i<i_mid0) keep i, where (i>i_mid0) keep i+1.

    dy   = diff(points0(:,1));
    idel = find(abs(dy)<dy_thresh);
    if (isempty(idel))
     points1 = points0; np1 = np0; ikeep1 = [1:np0]; i_mid1=i_mid0;
    else
     j = find(idel<i_mid0); idel(j) = idel(j) + 1;
     keep1 = ones(np0,1); keep1(idel) = 0; ikeep1 = find(keep1);
     points1 = points0(ikeep1,:); np1 = size(points1,1);
     i_mid1 = nnz(keep1(1:i_mid0));
     if (idebug>=1)
      disp(sprintf('Identified %d points with vertical slope', length(idel)));
     end
    end

    % delete points where the profile folds back, i.e. y(i+1) < yrmx(i) (rail)
    %                                              or  y(i+1) > yrmx(i) (wheel)

    keep2 = zeros(np1,1); keep2(i_mid1) = 1;
    y = points1(:,1);
    if (is_wheel), sgn = -1; else, sgn = 1; end

    % rail: check right of mid-point (wheel: left of mid-point)
    y_runmax = sgn*y(i_mid1);
    for i = i_mid1+1 : np1
     if (sgn*y(i) >= y_runmax+dy_thresh)
      keep2(i) = 1;
     end
     y_runmax = max(y_runmax, sgn*y(i));
    end

    % rail: check left of mid-point (wheel: right)
    y_runmin = sgn*y(i_mid1);
    for i = i_mid1-1 : -1 : 1
     if (sgn*y(i) <= y_runmin-dy_thresh)
      keep2(i) = 1;
     end
     y_runmin = min(y_runmin, sgn*y(i));
    end

    idel = find(keep2==0);
    if (isempty(idel))
     points2 = points1; np2 = np1; ikeep2 = [1:np1];
    else
     ikeep2 = find(keep2);
     points2 = points1(ikeep2,:); np2 = size(points2,1);
     if (idebug>=1)
      disp(sprintf('Identified %d points where profile folds back', length(idel)));
     end
    end

    % get mask that combines keep1 and keep2

    msk = zeros(np0,1);
    msk( ikeep1(ikeep2) ) = 1;

    if (delete_foldback<=1)
     p_out.Mask             = msk;
     p_out.ProfileY         = points0(:,1);
     p_out.ProfileZ         = points0(:,2);
     p_out.ProfileAngle     = points0(:,3);
     p_out.ProfileCurvature = points0(:,4);
    else
     p_out.Mask             = ones(np2,1);
     p_out.ProfileY         = points2(:,1);
     p_out.ProfileZ         = points2(:,2);
     p_out.ProfileAngle     = points2(:,3);
     p_out.ProfileCurvature = points2(:,4);
    end
   end

   % modification 6: scale y,z-coordinates

   if (scale_yz~=1)
    if (idebug>=1)
     disp(sprintf('Scaling y,z-values with factor %3.1f to convert to mm...',scale_yz));
    end
    p_out.ProfileY     = scale_yz * p_out.ProfileY;
    p_out.ProfileZ     = scale_yz * p_out.ProfileZ;
    p_out.ProfileCurvature = p_out.ProfileCurvature / scale_yz;
   end

   % combine the columns into array ProfileData

   p_out.ProfileData = [p_out.ProfileY, p_out.ProfileZ, p_out.ProfileAngle, p_out.ProfileCurvature];
   p_out.XYPoints = size(p_out.ProfileData,1);

   % plot the modified profile if requested

   if (make_plot)
    figure(make_plot); clf; hold on;
    plot(p_out.OrigData(:,1), p_out.OrigData(:,2), '-*');
    plot(p_out.ProfileY, p_out.ProfileZ, '-o');
    % axis([-70 70 -5 30]);
    grid on;
    xlabel('y_{prf} [mm]'); ylabel('z_{prf} [mm]');
    set(gca,'ydir','reverse');
    % legend('Original data','Modified data','location','NorthWest');
   end

  end % modify_profile


  function [ p ] = read_simpack(fname, idebug, make_plot)

   % [ p ] = read_simpack(fname, [idebug], [make_plot])
   %
   % Lower-level routine for reading wheel/rail profiles in SIMPACK format.
   % Does not do automatic corrections like mirroring or reversing order.

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
   % #read_simpack -1

   if (nargin<2 | isempty(idebug))
    idebug = 1;
   end
   if (nargin<3 | isempty(make_plot))
    make_plot = 0;
   end

   known_fields_header = struct( ...
    'version',           'int', ...
    'type',              'int',...
    'meas_pos_e',        'real', ...
    'meas_pos_qr',       'real'  ...
    );

   known_fields_spline = struct( ...
    'approx_smooth',     'real', ...
    'file',              'str', ...
    'file_mtime',        'int', ...
    'comment',           'str', ...
    'type',              'int', ...
    'point_dist_min',    'real', ...
    'shift_y',           'real', ...
    'shift_z',           'real', ...
    'rotate',            'real', ...
    'bound_y_min',       'real', ...
    'bound_y_max',       'real', ...
    'bound_z_min',       'real', ...
    'bound_z_max',       'real', ...
    'mirror_y',          'int', ...
    'mirror_z',          'int', ...
    'inversion',         'int', ...
    'units_len',         'str', ...
    'units_ang',         'str', ...
    'units_len_f',       'real', ...
    'units_ang_f',       'real' ...
    );

   % Initialize output structure

   p = struct;
   p.OrigData  = [0 0 0];
   nrows_data  = 0;
   p.ProfileY  = [];
   p.ProfileZ  = [];
   p.SplineWgt = [];

   if (~exist(fname,'file'))
    disp(['ERROR: cannot find file ',fname]);
    return
   end
   f = fopen(fname,'r');
   iline = 0;
   in_header = 0;
   in_spline = 0;
   in_point = 0;

   while(~feof(f))

    % get next line from input

    s  = strtrim(fgets(f));
    iline = iline + 1;
    if (idebug>=3)
     disp(sprintf('%4d: %s',iline,s));
    end

    % remove comments -- TODO: support comments in text strings

    ix = findstr('!', s);
    if (~isempty(ix))
     s = s(1:ix(1)-1);
    end

    % look for section start

    ix = findstr('.begin',s);
    if (~isempty(ix))
     section = strtrim(s(1:ix-1));
     if (strcmp(section,'header'))
      in_header = 1;
     elseif (strcmp(section,'spline'))
      in_spline = 1;
     elseif (strcmp(section,'point'))
      in_point = 1;
     else
      disp(['Unknown section: ',section]);
      return
     end
     continue % go read next line
    end

    % look for section end

    ix = findstr('.end',s);
    if (~isempty(ix))
     section = strtrim(s(1:ix-1));
     if (strcmp(section,'header'))
      in_header = 0;
     elseif (strcmp(section,'spline'))
      in_spline = 0;
     elseif (strcmp(section,'point'))
      in_point = 0;
     else
      disp(['Unknown section: ',section]);
      return
     end
     continue % go read next line
    end

    % look for an equals sign
    % if found: interpret keyword=value syntax

    ix = findstr('=',s);

    if (~isempty(ix))
     ix = ix(1);
     field_name = strtrim(s(1:ix-1));
     % replace dots in keyword by underscores
     jx = findstr('.',field_name);
     field_name(jx) = '_';
     field_val  = strtrim(s(ix+1:end));

     if (in_header)
      section = 'header';
      known_fields = known_fields_header;
     elseif (in_spline)
      section = 'spline';
      known_fields = known_fields_spline;
     else
      disp('ERROR: keyword outside a section');
      return
     end

     % check if the keyword is known
     if (isfield(known_fields, field_name))
      % interpret real numerical value
      if (strcmp(known_fields.(field_name), 'num') | strcmp(known_fields.(field_name), 'real'))
       p.(section).(field_name) = str2num(field_val);
       % interpret integer value
      elseif (strcmp(known_fields.(field_name), 'int'))
       p.(section).(field_name) = round(str2num(field_val));
       % store string value
      elseif (strcmp(known_fields.(field_name), 'str'))
       p.(section).(field_name) = field_val;
       % unknown type: store string value if it's non-empty
      elseif (~isempty(field_val))
       disp(sprintf('Unknown field-type %s: value="%s"',field_name,field_val));
       p.(section).(field_name) = field_val;
      end
      % unknown keyword: report & store string value
     else
      disp(sprintf('Unknown field "%s"',field_name));
      p.(section).(field_name) = field_val;
     end

     % no keyword found - interpret data line

    elseif (in_spline & in_point)

     [tmp, cnt] = sscanf(s, '%f %f %f %f %f');
     if (cnt > 0)
      nrows_data = nrows_data + 1;
      p.OrigData(nrows_data, 1:cnt) = tmp;
     end

    end % ~isempty(ix), keyword handling
   end % while (~feof)
   fclose(f);

   % Apply modifications to original input data

   opts   = p.spline;
   points = p.OrigData;
   np     = size(points,1);

   % modification 1: remove points too close together

   if (isfield(opts,'point_dist_min') & opts.point_dist_min>=1e-12)
    tmp = points;
    inew = 1;
    for ip = 2:np
     if (abs(tmp(inew,1)-points(ip,1))>opts.point_dist_min)
      inew = inew + 1;
      tmp(inew,:) = points(ip,:);
     end
    end
    points = tmp(1:inew,:);
    np     = inew;
   end

   % modifications 2 & 3: shift y and z

   if (isfield(opts,'shift_y'))
    points(:,1) = points(:,1) + opts.shift_y;
   end
   if (isfield(opts,'shift_z'))
    points(:,2) = points(:,2) + opts.shift_z;
   end

   % modification 10: convert angle to radians

   if (isfield(opts,'units_ang_f') & isfield(opts,'rotate'))
    opts.rotate = opts.rotate / opts.units_ang_f;
   end

   % modification 4: rotate profile about profile origin

   if (isfield(opts,'rotate'))
    rot = [cos(opts.rotate), -sin(opts.rotate); sin(opts.rotate), cos(opts.rotate)];
    points(:,1:2) = points(:,1:2) * rot';
   end

   % modification 5: remove points outside [bound_y_min, bound_y_max]

   if (isfield(opts,'bound_y_min') & isfield(opts,'bound_y_max') & opts.bound_y_min<opts.bound_y_max)
    ix = find(points(:,1)<opts.bound_y_min | points(:,1)>opts.bound_y_max);
    if (~isempty(ix))
     points(ix,:) = [];
     np = size(points,1);
    end
   end

   % modification 6: clipping to [bound_z_min, bound_z_max]

   if (isfield(opts,'bound_z_min') & isfield(opts,'bound_z_max') & opts.bound_z_min<opts.bound_z_max)
    points(:,2) = max(opts.bound_z_min, min(opts.bound_z_max, points(:,2)));
   end

   % modification 7: mirror y

   if (isfield(opts,'mirror_y') & opts.mirror_y==1)
    points(:,1) = - points(:,1);
   end

   % modification 8: mirror z

   if (isfield(opts,'mirror_z') & opts.mirror_z==1)
    points(:,2) = - points(:,2);
   end

   % modification 9: invert data order

   if (isfield(opts,'inversion') & opts.inversion==1)
    points = flipud(points);
   end

   % modification 10: conversion of lengths unit

   if (isfield(opts,'units_len_f'))
    points(:,1:2) = 1000 * points(:,1:2) / opts.units_len_f;
   end

   % store end-result in structure

   p.ProfileY  = points(:,1);
   p.ProfileZ  = points(:,2);
   p.SplineWgt = points(:,3);

   if (make_plot)
    figure(make_plot); clf; hold on;
    plot(p.OrigData(:,1), p.OrigData(:,2), '-*');
    plot(p.ProfileY, p.ProfileZ, '-o');
    axis([-70 70 -5 30]);
    grid on;
    xlabel('y_{prf} [mm]'); ylabel('z_{prf} [mm]');
    set(gca,'ydir','reverse');
    legend('Original data','Modified data','location','NorthWest');
   end

  end % read_simpack

  function [ slcs ] = read_slices(fname, mirror_y, mirror_z, scale_yz, idebug)

   % [ slcs ] = read_slices(fname, [mirror_y], [mirror_z], [scale_yz], [idebug])
   %
   % Lower-level routine for reading wheel/rail profiles in slices format.
   % mirror_y    - -1 or 0 for no, 1 for yes
   % mirror_z    - -1 for no, 0 for automatic, 1 for yes
   % scale_yz    - scale-factor to convert profile data to mm, e.g. 1000. for data in meters

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
   % #read_slices -1

   if (nargin<2 | isempty(mirror_y))
    mirror_y = 0;
   end
   if (nargin<3 | isempty(mirror_z))
    mirror_z = 0;
   end
   if (nargin<4 | isempty(scale_yz))
    scale_yz = 1;
   end
   if (nargin<5 | isempty(idebug))
    idebug = 1;
   end
   dsmax_2d = 0.5;

   show_fig = [ ];
   fig_ofs  = 4;
   % figs 1 - 4: resampling of slices
   % figs 5 - 6: slices before/after interruption

   % using '.slcs' for rails and '.slcw' for wheels

   [~,~,ext] = fileparts(fname);
   is_wheel  = strcmp(lower(ext), '.slcw');

   if (is_wheel)           % 1st surface parameter u == theta or u == s along track curve
    nam_udir = 'theta';
   else
    nam_udir = 'U';
   end

   % Initialize output structure

   slcs = struct;
   slcs.slc_file = fname;
   slcs.is_wheel = is_wheel;
   slcs.nslc     = 0;
   slcs.u_offset = 0;
   slcs.u_scale  = 1;
   slcs.u        = [];
   slcs.fnames   = [];

   % read file contents

   if (~exist(fname,'file'))
    disp(['ERROR: cannot find file ',fname]);
    return
   end

   f = fopen(fname,'r'); iline  = 0;
   [slcs, iline, ierror] = cntc.read_slices_header(slcs, nam_udir, f, iline, idebug);
   if (ierror==0)
    [slcs, iline, ierror] = cntc.read_slices_fnames(slcs, nam_udir, f, iline, idebug);
   end
   if (ierror==0)
    [slcs, iline, ierror] = cntc.read_feature_info(slcs, nam_udir, f, iline, idebug);
   end
   fclose(f);
   if (ierror), return; end

   % apply u-offset and scaling

   slcs.u = slcs.u_scale * (slcs.u + slcs.u_offset);

   % read each of the profiles

   [slcs, ierror] = cntc.read_slice_profiles( fname, slcs, mirror_y, mirror_z, scale_yz, idebug );

   % resample all slices to the same number of points

   if (ierror==0)
    [slcs, ierror] = cntc.resample_slices( slcs, dsmax_2d, idebug, show_fig, fig_ofs );
   end

   % add 2d spline representation
   % rails:  parameter (u,v) --> cartesian (x,y,z) with x==u.
   % wheels: parameter (u,v) --> cylindrical (th,y,dr) with th==u, basic interval [-pi,pi) + wrap-around
   % TODO: smoothing in x/theta-direction

   if (ierror==0 & exist('make_2dspline'))
    use_approx = (slcs.u_intpol==2); use_insert = 1; use_cylindr = is_wheel; idebug = 0;
    slcs.spl2d = cntc.make_2dspline(slcs.u, slcs.vj, [], slcs.ysurf, slcs.zsurf, slcs.mask_j, ...
     use_approx, use_insert, use_cylindr, idebug);
   end

  end % read_slices

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ slcs, iline, ierror ] = read_slices_header(slcs, nam_udir, f, iline, idebug)

   % [ slcs, iline ] = read_slices_header(slcs, nam_udir, f, iline, [idebug])
   %
   % Lower-level routine for reading contents of slcs-file
   % #read_slices_header

   ierror    = 0;
   in_header = 1;

   % get next line from input

   while(~feof(f) & in_header<=4)

    s  = strtrim(fgets(f));
    iline = iline + 1;
    if (idebug>=3)
     disp(sprintf('%4d: %s',iline,s));
    end

    % remove comments, everything after %

    ix = strfind(s, '%'); if (~isempty(ix)), s=s(1:ix-1); end

    % process line

    if (~isempty(s))

     if (in_header==1)

      % 1st line: [S,TH]_OFFSET, [S,TH]_SCALE

      [tmp, nval] = sscanf(s, '%f %f');
      if (nval<2)
       disp(sprintf('ERROR: first line should have two values, %s_OFFSET, %s_SCALE, obtained %d',...
        nam_udir, nam_udir, nval));
       ierror = 1;
      else
       slcs.u_offset = tmp(1);
       slcs.u_scale  = tmp(2);
      end
      in_header = in_header + 1;

     elseif (in_header==2)

      % 2nd line: NSLC

      slcs.nslc = sscanf(s, '%f');
      in_header = in_header + 1;

     elseif (in_header==3)

      % 3rd line: NFEAT, NKINK, NACCEL

      [tmp, nval] = sscanf(s, '%f %f %f');
      if (nval<3), tmp(3) = 0; end

      slcs.nfeat  = tmp(1);
      slcs.nkink  = tmp(2);
      slcs.naccel = tmp(3);
      in_header = in_header + 1;

     elseif (in_header==4)

      % 4th line: [U,TH]_INTPOL

      slcs.u_intpol = sscanf(s, '%f');
      in_header = in_header + 1;
      if (~any(slcs.u_intpol==[1 2]))
       disp(sprintf('ERROR: %s_INTPOL should be 1 or 2', nam_udir));
       disp(['input: ', s]);
       ierror = 4;
      end

     end

    end % ~isempty(s)
   end % while

   if (idebug>=-2)
    disp(sprintf('Slices-file "%s": %d slices, %d features, %d parts, %d kinks, %d accel', ...
     slcs.slc_file, slcs.nslc, slcs.nfeat, max(0,slcs.nfeat-1), slcs.nkink, slcs.naccel));
   end

  end % read_slices_header

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ slcs, iline, nerror ] = read_slices_fnames(slcs, nam_udir, f, iline, idebug)

   % [ slcs, iline, nerror ] = read_slices_fnames(slcs, nam_udir, f, iline, [idebug])
   %
   % Lower-level routine for reading filenames from slcs-file
   % #read_slices_fnames

   nslc   = 0; % #slices_read_from file -3
   nerror = 0;

   while(~feof(f) & nslc<slcs.nslc)

    s  = strtrim(fgets(f));
    iline = iline + 1;
    if (idebug>=3)
     disp(sprintf('%4d: %s',iline,s));
    end

    % remove comments, everything after %

    ix = strfind(s, '%'); if (~isempty(ix)), s=s(1:ix-1); end

    if (~isempty(s))

     % NSLC lines: S/TH_SLC, 'R/WFNAME' or S/TH_SLC, "R/WFNAME"

     ix = findstr(s, '''');
     if (length(ix)==0)
      ix = findstr(s, '"');
     end

     if (length(ix)==2)
      nslc = nslc + 1;
      val            = sscanf(s(1:ix(1)-1), '%f');
      if (isscalar(val))
       slcs.u(nslc,1) = val;
      end
      nam            = s(ix(1)+1:ix(2)-1);
      if (~isempty(nam))
       slcs.fnames    = strvcat( slcs.fnames, nam );
      end
     end
     if (length(ix)~=2 | length(val)~=1 | isempty(nam))
      disp(sprintf('ERROR: could not interpret line %d:',iline));
      disp(sprintf('       "%s"',s));
      disp(sprintf('       line should provide %s_SLC, ''RFNAME''', nam_udir));
      nerror = nerror + 1;
     end

    end % ~isempty(s)
   end % while (~feof & islc<nslc)

   if (nslc<slcs.nslc)
    disp(sprintf('ERROR: got %d slices, expected NSLC = %d', nslc, slcs.nslc));
    nerror = nerror + 1;
   end

   % check that u-positions are in strictly increasing order, du >= 0.001 mm or >= 0.001 rad

   tiny_du = 0.001 / slcs.u_scale;
   du = diff(slcs.u);
   if (any(du<tiny_du))
    ix = find(du<tiny_du, 1, 'first');
    disp(sprintf('ERROR: slice %s-positions should be strictly increasing.', nam_udir));
    disp(sprintf('       islc= %d, %d: %s=%8.3f, %8.3f', ix, ix+1, nam_udir, slcs.u(ix), slcs.u(ix+1)));
    nerror = nerror + 1;
   end

   % for wheels, require that |u_end - u_1| < 2*pi

   du = (slcs.u(end) - slcs.u(1)) * slcs.u_scale;
   if (slcs.is_wheel & du>2*pi)
    disp(sprintf('ERROR: slice %s-positions should be at most 2*pi apart.', nam_udir));
    disp(sprintf('       islc= %d, %d: %s=%8.3f, %8.3f', 1, slcs.nslc, nam_udir, slcs.u(1), slcs.u(end)));
    nerror = nerror + 1;
   end

  end % read_slices_fnames

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ slcs, iline, nerror ] = read_feature_info(slcs, nam_udir, f, iline, idebug)

   % [ slcs, iline, nerror ] = read_feature_info(slcs, nam_udir, f, iline, [idebug])
   %
   % Lower-level routine for reading feature information from a slcs-file

   % TODO: add 2 lines: P_KINK, P_ACCEL
   % #read_feature_info

   nerror = 0;

   % NFEAT <= 1: using default feature information

   if (slcs.nfeat<=1)

    slcs.s_feat = ones(slcs.nslc,1) * [0, 1e6];

   else

    % next NSLC lines: u_slc, s_f, f = 1 : nfeat

    nslc_p = 0; % #slices_with_features_read from file -3
    slcs.s_feat = zeros(slcs.nslc, slcs.nfeat);

    while(~feof(f) & nslc_p<slcs.nslc)

     % get next line from input

     s  = strtrim(fgets(f));
     iline = iline + 1;
     if (idebug>=3)
      disp(sprintf('%4d: %s',iline,s));
     end

     % remove comments, everything after %

     ix = strfind(s, '%'); if (~isempty(ix)), s=s(1:ix-1); end

     % process line

     if (~isempty(s))

      nslc_p  = nslc_p + 1;
      is_ok  = 1;

      [tmp, nval] = sscanf(s, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');

      if (nval>=1)
       u_slc = tmp(1);
       i_slc = find(abs(slcs.u-u_slc)<1e-3);

       if (length(i_slc)<1)
        is_ok = 0;
        disp(sprintf('ERROR: could not find slice at %s_slc = %4.2f in list of slices', ...
         nam_udir, u_slc));
        nerror = nerror + 1;
       elseif (length(i_slc)>1)
        is_ok = 0;
        disp(sprintf('ERROR: multiple slices with %s_slc = %4.2f in list of slices', ...
         nam_udir, u_slc));
        if (~isempty(i_slc))
         disp(sprintf('       possible slices i = %d %d %d %d %d', i_slc));
        end
        nerror = nerror + 1;
       elseif (any(slcs.s_feat(i_slc,:)))
        is_ok = 0;
        disp(sprintf('ERROR: repeated definition of features for slice with %s_slc = %4.2f', ...
         nam_udir, u_slc));
        nerror = nerror + 1;
       end
      end

      nbreak = nval - 1;

      if (is_ok & nbreak~=slcs.nfeat)
       disp(sprintf('ERROR: incorrect number of feature positions for slice at %s_slc = %4.2f', ...
        nam_udir, u_slc));
       disp(sprintf('       obtained %d values, expecting n_feat = %d', nbreak, slcs.nfeat));
       is_ok = 0;
       nerror = nerror + 1;
      end

      if (is_ok)
       slcs.s_feat(i_slc, 1:nbreak) = tmp(2:nbreak+1);
      end

     end % ~isempty(s)
    end % while (~feof)

    if (nslc_p<slcs.nslc)
     disp(sprintf('ERROR: got feature information for %d slices, expected NSLC = %d', nslc_p, slcs.nslc));
     nerror = nerror + 1;
    end
   end % nfeat<=1

  end % read_feature_info

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ slcs, ierror ] = read_slice_profiles( fname, slcs, mirror_y, mirror_z, scale_yz, idebug )

   % [ slcs, ierror ] = read_slice_profiles( fname, slcs, mirror_y, mirror_z, scale_yz, idebug )
   %
   % read each of the per-slice profiles
   % #read_slice_profiles

   ierror = 0;

   % determine relative folder where slices-file resides

   filepath = fileparts(fname);

   nslc = slcs.nslc;

   for is = 1 : nslc
    if (idebug>=1 | (idebug>=-1 & is==1) | (idebug>=-1 & is==nslc))
     disp(sprintf('Slice %3d: read file "%s"', is, strtrim(slcs.fnames(is,:)) ));
    end

    slc_file = fullfile(filepath, strtrim(slcs.fnames(is,:)));
    prf = cntc.read_profile(slc_file, slcs.is_wheel, mirror_y, mirror_z, scale_yz, [], idebug-1);

    % in case of an error (missing file), prf will be a struct with empty members

    if (isempty(prf.ProfileY))
     ierror = 1;
    else
     prf.npnt = length(prf.ProfileY);
     prf.ltot = prf.ProfileS(end) - prf.ProfileS(1);
    end

    if (is>1 & ~isempty(prf.ProfileY))

     % set empty for fields present in slcs.prf / missing in new prf

     slcs_fields = fieldnames(slcs.prf(1));
     for i = 1:length(slcs_fields)
      if (~isfield(prf, slcs_fields{i}))
       % disp(sprintf('prf(%d) does not have field "%s", left empty', is, slcs_fields{i}));
       prf = setfield(prf, slcs_fields{i}, []);
      end
     end

     % set empty for fields present in prf / missing in previous slcs.prf

     prf_fields = fieldnames(prf);
     for i = 1:length(prf_fields)
      if (~isfield(slcs.prf(1), prf_fields{i}))
       % disp(sprintf('prf(%d) introduces field "%s", empty for previous slices', is, prf_fields{i}));
       slcs.prf = setfield(slcs.prf, prf_fields{i}, []);
      end
     end
    end % is>1 & prf

    if (~isempty(prf.ProfileY))
     % disp(sprintf('Adding prf(%2d) = "%s"', is, prf.Fname));
     slcs.prf(is,1) = prf;
    end
   end % for is

  end % read_slice_profiles

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ p ] = read_vampire(fname, right_side, idebug, make_plot)

   % [ p ] = read_vampire(fname, [right_side], [idebug], [make_plot])
   %
   % Lower-level routine for reading wheel/rail profiles in Vampire whe/rai format.
   % Does not do automatic corrections like mirroring or reversing order.
   %
   %   fname      - filename of input file
   %   right_side - select right side profile (1, default) or left side (0)
   %   make_plot  - figure number to use for plotting the profile

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
   % #read_vampire(fname, right_side, idebug, make_plot) -1

   if (nargin<2 | isempty(right_side))
    right_side = 1;
   end
   if (nargin<3 | isempty(idebug))
    idebug = 1;
   end
   if (nargin<4 | isempty(make_plot))
    make_plot = 0;
   end

   % Initialize output structure

   p                  = struct;
   p.ProfileData      = [];
   nrows_data         = 0;
   p.ProfileY         = [];
   p.ProfileZ         = [];
   inch               = 25.4; % [mm/inch]

   if (~exist(fname,'file'))
    disp(['ERROR: cannot find file ',fname]);
    return
   end

   f = fopen(fname,'r');
   iline = 0;
   in_profile = 0;

   % read first line, should contain 'WHEELPROFILE' or 'RAILPROFILE'

   [ s, iline ] = cntc.get_next_line(f, iline, idebug);
   if (strcmp(s, 'WHEELPROFILE'))
    p.is_wheel = 1;
   elseif (strcmp(s, 'RAILPROFILE'))
    p.is_wheel = 0;
   else
    disp(sprintf('ERROR: file ''%s'' is not in Vampire format.', fname));
    return;
   end

   % second line: 80-character title

   [ s, iline ] = cntc.get_next_line(f, iline, idebug);
   p.Title = s;

   % third line: tread datum (wheel) or gauge point (rail)

   if (p.is_wheel)

    % third line, wheel: tread datum

    [ s, iline ] = cntc.get_next_line(f, iline, idebug);
    ix = length('TREADDATUM') + 1;
    p.TreadDatum = strtrim(s(ix:end));

    ix = strfind('CUSTOM', p.TreadDatum);
    if (~isempty(ix))
     p.FlangeBackPos = sscanf(p.TreadDatum(ix+6:end), ' %f');
    elseif (strfind('EUROPEAN', p.TreadDatum))
     p.FlangeBackPos = 70;
    elseif (strfind('US1', p.TreadDatum))
     p.FlangeBackPos = (2 + 27/32) * inch;
    elseif (strfind('US2', p.TreadDatum))
     p.FlangeBackPos = (3 + 1/16) * inch;
    else
     disp(['Error: unknown tread datum ',p.TreadDatum]);
    end
    p.FlangeBackPos = -p.FlangeBackPos;     % CONTACT: from datum to flange back

   else

    % third line, rail: gauge point

    [ s, iline ] = cntc.get_next_line(f, iline, idebug);
    ix = length('GAUGEPOINT') + 1;
    p.GaugePoint = strtrim(s(ix:end));

    ix = strfind('CUSTOM', p.GaugePoint);
    if (~isempty(ix))
     p.GaugeHeight = sscanf(p.GaugePoint(ix+6:end), ' %f');
    elseif (strfind('EUROPEAN', p.GaugePoint))
     p.GaugeHeight = 14;
    elseif (strfind('US', p.GaugePoint))
     p.GaugeHeight = 5/8 * inch;
    else
     disp(['Error: unknown gauge point ',p.GaugePoint]);
    end
   end

   % fourth line: flange back distance (wheel) or gauge (rail)

   [ s, iline ] = cntc.get_next_line(f, iline, idebug);

   if (p.is_wheel)
    ix = length('FLANGEBACK') + 1;
    p.FlangeBackDist = sscanf(s(ix:end), '%f');
   else
    ix = length('GAUGE') + 1;
    p.GaugeWidth = sscanf(s(ix:end), '%f');
   end

   while(~feof(f))

    % get next line from input

    [ s, iline ] = cntc.get_next_line(f, iline, idebug);

    % remove comments, everything after ", !, or %
    ix = strfind(s, '"'); if (~isempty(ix)), s=s(1:ix-1); end
    ix = strfind(s, '!'); if (~isempty(ix)), s=s(1:ix-1); end
    ix = strfind(s, '%'); if (~isempty(ix)), s=s(1:ix-1); end

    % get (yl,zl), (yr,zr) positions on left and right profiles

    ix = strfind(s, ','); if (~isempty(ix)), s(ix)=' '; end
    [tmp, cnt] = sscanf(s, '%f %f %f %f %f %f %f %f');

    if (cnt>0)
     nrows_data = nrows_data + 1;
     p.ProfileData(nrows_data, 1:cnt) = tmp;
    end

   end % while (~feof)
   fclose(f);

   % select left or right side profile, convert to right side profile

   if (right_side)
    p.ProfileY =  p.ProfileData(:,3);
    p.ProfileZ =  p.ProfileData(:,4);
   else
    p.ProfileY = -p.ProfileData(:,1);
    p.ProfileZ =  p.ProfileData(:,2);
   end

   % renumber rails from track center to field side

   if (~p.is_wheel & (p.ProfileY(end)<p.ProfileY(1)))
    p.ProfileY = flipud(p.ProfileY);
    p.ProfileZ = flipud(p.ProfileZ);
   end

   % shift the datum to zero

   if (p.is_wheel)
    fb_dist = p.FlangeBackDist;
    fb_pos  = p.FlangeBackPos;
    p.ProfileY = p.ProfileY - fb_dist/2 + fb_pos;
   else
    p.ProfileZ = -p.ProfileZ;    % Vampire rails: z positive upwards

    zmin = min(p.ProfileZ);      % interpolate to get the gauge point
    ix = find(p.ProfileZ-zmin < p.GaugeHeight, 1, 'first');
    if (~isempty(ix) & ix>1)
     ygauge = interp1( p.ProfileZ(ix+[-1,0])-zmin, p.ProfileY(ix+[-1,0]), p.GaugeHeight );
    else
     ygauge = p.ProfileY(1);
    end
    % shift rail to put the gauge point at y_r = 0
    p.ProfileY = p.ProfileY - ygauge;
    p.ProfileZ = p.ProfileZ - zmin;
    p.rail_y0  = ygauge;
    p.rail_z0  = zmin;
   end

   if (make_plot)
    figure(make_plot); clf; hold on;
    plot(p.ProfileY, p.ProfileZ, '-o');
    grid on;
    set(gca,'ydir','reverse');
    if (p.is_wheel)
     xlabel('y_w [mm]'); ylabel('z_w [mm]');
    else
     xlabel('y_r [mm]'); ylabel('z_r [mm]');
    end
    title(p.Title)
   end

  end % read_vampire

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ s, iline ] = get_next_line(f, iline, idebug)
   % #get_next_line -2

   s = deblank(fgets(f));
   iline = iline + 1;
   if (idebug>=5 & iline<10)
    disp(sprintf('%4d: %s', iline, s));
   end

  end % get_next_line

  function [ slcs, ierror ] = resample_slices( slcs, ds_max2d, idebug, show_fig, fig_ofs )

   % [ slcs, ierror ] = resample_slices( slcs, ds_max2d, [idebug], [show_fig], [fig_ofs] )
   %
   % resample slices on the basis of feature information, using the same number of points per part in
   % each slice.
   %
   % in:  slcs.prf    - profiles per slice, with different #points(islc) and arc-length L(islc)
   %      slcs.s_feat - break points per slice [ s0, s1, .., sn ] corresponding to 'geometric features'
   %      ds_max2d    - maximum step size ds after resampling
   %
   % s_feat(islc,:) = [ -1, -1, 5.0, 10.3, 20.0, 999.0, -1, -1 ]
   %         -->  this has 8 breaks, 7 parts
   %         -->  for this slice 'islc', breaks 1, 2, 7, 8 are non-existent
   %         -->  break 3 == 5.0 means that 5 mm will be trimmed off from the start of the profile
   %         -->  break 6 == 999.0 will be changed to s(end), the end-point of the profile

   % variable naming focuses on the use of rail slices. Wheels use (x,y,z) to mean (th,y,dr).

   % method:
   %  - for this slice, len_part(islc,:) = [ NaN, NaN, 5.3, 9.7, s(end)-20.0, NaN, NaN ]
   %  - the maximum max(len_part) is computed over all slices,
   %                         say max_len = [ 5.2, 5.6, 5.8, 9.7, 15.0, 10.0, 9.0 ]
   %  - a target step size ds_max2d is given, say ds_max2d = 1.0
   %  - the number of points per part is set accordingly, e.g. n_p = [ 6, 6, 6, 10, 15, 10, 9 ]
   %  - each slice is resampled using its own step sizes ds, e.g. ds_p = [ NaN, NaN, 5.3/(6-1), ... ]

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
   % #resample_slices -2

   if (nargin<3 | isempty(idebug))
    idebug   = 0;
   end
   if (nargin<4)
    show_fig = [];
   end
   if (nargin<5 | isempty(fig_ofs))
    fig_ofs = 0;
   end
   is_debug = 0; % slice for which information is printed
   if (idebug>=3 & ~isempty(show_fig))
    is_debug = 8;
   end
   ierror = 0;

   s_feat = slcs.s_feat;
   nfeat  = size(s_feat, 2);
   nslc   = size(s_feat, 1);

   % check size of s_feat array

   if (nslc~=slcs.nslc | nfeat<=1)
    disp(sprintf('Incorrect break points s_feat, needs at least 2 s-values for each of %d slices', ...
     slcs.nslc));
    disp(size(s_feat));
    ierror = 1;
    return;
   end

   % check and adapt the feature information

   for is = 1 : nslc

    % - mark non-used breaks (features) for this slice by NaN-value

    ib  = find(s_feat(is,:) < 0);
    s_feat(is,ib) = NaN;

    % - determine active 'breaks' [ib0:ib1] for each slice

    ib0 = find(s_feat(is,:) >= 0, 1, 'first');
    ib1 = find(s_feat(is,:) >= 0, 1, 'last');
    slcs.slc_ib(is,1:2) = [ib0, ib1];

    % - shift s_feat from [0, L] as used in the parts-file to [s_0, s_end] as used in the profile

    s_feat(is,ib0:ib1) = s_feat(is,ib0:ib1) + slcs.prf(is).ProfileS(1);

    % - clip last value after s_end

    s_feat(is,ib1) = min(s_feat(is,ib1), slcs.prf(is).ProfileS(end));

    % - check that s_feat are in strictly increasing order and in range of profile

    if (any(diff(s_feat(is,ib0:ib1))<=0))
     disp(sprintf('Incorrect break points s_feat for slice %d, must be strictly increasing.',is))
     disp([ib0, ib1])
     disp(s_feat(is,:))
     m = min(diff(s_feat(is,ib0:ib1)));
     disp(sprintf('Minimum s_p - s_{p-1} = %3.1e', min(m)));
     ierror = 2;
     return;
    end

    % - determine the length of each part in this slice

    len_part(is,ib0:ib1-1) = diff(s_feat(is,ib0:ib1));

   end % for is: adapt parts-info

   % determine the maximum length per part over all of the slices

   L_part = max( len_part );

   % check scaling

   tot_len = sum(L_part);
   if (tot_len < 5*ds_max2d)
    disp(sprintf('Warning: step ds_max2d = %5.2f may be too large for total length = %5.2f', ...
     ds_max2d, tot_len));
   end

   % determine number of points and step size du for each of the parts

   nseg_p = ceil( L_part / ds_max2d );

   % determine the first point number i_p for each of the parts after resampling

   i_p    = cumsum([1, nseg_p]);

   % set total number of points after resampling

   n_pnt  = i_p(end);

   % set the sampling positions vj

   vj  = [1 : n_pnt];

   if (idebug>=2 | n_pnt<=4)
    for ipart = 1 : nfeat-1
     disp(sprintf('part %d: %3d segments, points [%3d,%3d]', ipart, nseg_p(ipart), i_p(ipart+[0,1])));
    end
   end

   % create output arrays for resampled surface

   slcs.npnt   = n_pnt;
   slcs.iseg_p = i_p;
   slcs.vj     = vj;
   slcs.mask_j =      zeros(slcs.nslc, n_pnt);
   slcs.xsurf  = NaN * ones(slcs.nslc, n_pnt);  % Note: wheel (th,y,dr) stored in (x,y,z)
   slcs.ysurf  = NaN * ones(slcs.nslc, n_pnt);
   slcs.zsurf  = NaN * ones(slcs.nslc, n_pnt);

   % create mask array for use in spline computation

   for is = 1 : slcs.nslc
    ib0 = slcs.slc_ib(is,1);  % active features [ib0:ib1]
    ib1 = slcs.slc_ib(is,2);
    i0 = slcs.iseg_p( ib0 );  % active points [i0:i1]
    i1 = slcs.iseg_p( ib1 );
    slcs.mask_j(is, i0:i1) = 1;
   end

   % loop over slices, resample each slice

   for is = 1 : slcs.nslc

    slc = slcs.prf(is);
    if (idebug>=3 & is==is_debug)
     disp(sprintf(['Slice %d: s=[%5.1f,%5.1f], s_p={%5.1f,',...
      '%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f}'], ...
      is, slcs.prf(is).ProfileS([1,end]), s_feat(is,:)));
    end

    % - determine breaks [ib0:ib1] available for this slice

    ib0 = find(~isnan(s_feat(is,:)), 1, 'first');
    ib1 = find(~isnan(s_feat(is,:)), 1, 'last');

    % - set the positions s_j corresponding to v_j

    sj     = NaN * ones(n_pnt,1);
    for ipart = ib0 : ib1-1

     i0  = i_p(ipart);
     i1  = i_p(ipart+1);
     v_len = vj(i1) - vj(i0);
     s_ofs = s_feat(is,ipart);
     s_len = s_feat(is,ipart+1) - s_feat(is,ipart);
     sj(i0:i1) = s_ofs + s_len * (vj(i0:i1) - vj(i0)) / v_len;

     if (idebug>=3 & is==is_debug)
      disp(sprintf(' part %d: s=[%6.1f,%6.1f], sj={%6.1f,%6.1f ..%6.1f}', ...
       ipart, s_feat(is,ipart+[0:1]), sj([i0,i0+1,i1]) ));
     end
    end

    % resample at selected sj-positions

    i0  = i_p(ib0);
    i1  = i_p(ib1);

    if (exist('make_spline'))

     % make spline for the slice, evaluate at positions sj

     lambda = 0; wgt = []; ikinks = []; iaccel = []; use_bspline = 1; ds_bspline = 0.5;
     if (length(slc.ProfileY)==7), ikinks = [2,3,4,5,6]; end    % hack for V-groove test-case
     spl = cntc.make_spline( slc.ProfileS, slc.ProfileY, slc.ProfileZ, lambda, wgt, ikinks, iaccel, ...
      use_bspline, ds_bspline);
     [~, yj, zj] = eval_spline(spl, sj(i0:i1));
     slcs.xsurf(is,i0:i1) = slcs.u(is);
     slcs.ysurf(is,i0:i1) = yj;
     slcs.zsurf(is,i0:i1) = zj;

    else
     % linear interpolation. Note: s_i(end) can be larger than ProfileS(end) by round-off error

     slcs.xsurf(is,:) = slcs.u(is);
     slcs.ysurf(is,:) = interp1( slc.ProfileS, slc.ProfileY, sj, 'linear', 'extrap' )';
     slcs.zsurf(is,:) = interp1( slc.ProfileS, slc.ProfileZ, sj, 'linear', 'extrap' )';
    end

    nnan = nnz( isnan(slcs.ysurf(is,i0:i1))+isnan(slcs.zsurf(is,i0:i1)) );
    if (nnan>0)
     disp(sprintf('resample_slices: slice %3d has %d nan-values',is, nnan));
     ierror = 3;
    end

    if (idebug>=3 & is==is_debug)
     if (any(show_fig==2))
      figure(2+fig_ofs); clf; hold on;
      plot(slcs.ysurf(is,:));
      plot(i_p, slcs.ysurf(is,i_p), 'o');
      grid on;
      xlabel('point j'); ylabel('s [mm]');
      title(sprintf('Slice %d: resampled y-coordinate',is));
     end

     if (any(show_fig==3))
      figure(3+fig_ofs); clf; hold on;
      plot(slcs.ysurf(is,:), slcs.zsurf(is,:));
      plot(slcs.ysurf(is,i_p), slcs.zsurf(is,i_p), 'o');
      set(gca,'ydir','reverse');
      grid on;
      xlabel('point j'); ylabel('s [mm]');
      title(sprintf('Slice %d: resampled profile',is));
     end
    end
   end

   if (any(show_fig==4))
    figure(4+fig_ofs); clf;
    surf(slcs.xsurf', slcs.ysurf', slcs.zsurf');
    set(gca, 'ydir','reverse', 'zdir','reverse');
    % v = axis; v(1:2) = slcs.u([5,end-4]); axis(v);
    v = get(gca,'dataaspectratio'); v(3) = v(2); set(gca,'dataaspectratio',v);
    shading flat; colorbar;
    if (slcs.is_wheel)
     xlabel('\theta_{w} [rad]');
     ylabel('y_{w} [mm]');
     zlabel('\delta{}r_{w} [mm]');
    else
     xlabel('u_{fc} [mm]');
     ylabel('y_{r} [mm]');
     zlabel('z_{r} [mm]');
    end
   end

  end % resample_slices

  function [ prr ] = get_profile_slice( slcs, u_i, make_plot )
   % #get_profile_slice

   if (nargin<3)
    make_plot = 0;
   end
   i0 = find(slcs.u <  u_i, 1, 'last');
   i1 = find(slcs.u >= u_i, 1, 'first');
   % disp([u_i, i0, i1])
   if (i0 >= slcs.nslc)
    yi = slcs.ysurf(end,:); zi = slcs.zsurf(end,:);
    % disp(sprintf('u_i = %5.1f > u(end) = %5.1f, using slice %d', u_i, slcs.u(end), i1));
   elseif (i1 <= 1)
    yi = slcs.ysurf(1,:); zi = slcs.zsurf(1,:);
    % disp(sprintf('u_i = %5.1f < u(1) = %5.1f, using slice %d', u_i, slcs.u(1), i0));
   else
    fac0 = (slcs.u(i1) - u_i) / (slcs.u(i1) - slcs.u(i0));
    fac1 = (u_i - slcs.u(i0)) / (slcs.u(i1) - slcs.u(i0));
    if (fac0>0.99 & nnz(slcs.mask_j(i0,:))>nnz(slcs.mask_j(i1,:)))
     % disp(sprintf('Using longer slice i0=%d', i0));
     yi = slcs.ysurf(i0,:); zi = slcs.zsurf(i0,:);
    elseif (fac1>0.99 & nnz(slcs.mask_j(i1,:))>nnz(slcs.mask_j(i0,:)))
     % disp(sprintf('Using longer slice i1=%d', i1));
     yi = slcs.ysurf(i1,:); zi = slcs.zsurf(i1,:);
    else
     % disp(sprintf('Using %5.3f * slice %d + %5.3f * slice %d',fac0, i0, fac1, i1));
     yi = fac0 * slcs.ysurf(i0,:) + fac1 * slcs.ysurf(i1,:);
     zi = fac0 * slcs.zsurf(i0,:) + fac1 * slcs.zsurf(i1,:);
    end

    if (make_plot)
     tmpfig = gcf;
     figure(make_plot); clf; hold on;
     plot([slcs.ysurf([i0:i1],:);yi]', [slcs.zsurf([i0:i1],:);zi]');
     set(gca,'ydir','reverse'); grid on; axis equal;
     legend(sprintf('slice %d, u=%6.2f',i0,slcs.u(i0)), sprintf('slice %d, u=%6.2f',i1,slcs.u(i1)), ...
      sprintf('interpolated, u=%6.2f',u_i), 'location','southeast')
     figure(tmpfig);
    end
   end

   % set arc-length s, accounting for NaN's in 'missing parts'
   % s  = slcs.vj; ix = find(isnan(yi)); s(ix) = NaN;

   prr = struct('ProfileY',yi', 'ProfileZ',zi', 'ProfileS',slcs.vj);

  end % get_profile_slice


  function [ p_out, spl ] = smooth_profile( p_in, l_filt, lambda, ds_sampl, ikinks, iaccel, use_wgt, use_bspline, ds_bspl, idebug );

   % [ p_out, spl ] = smooth_profile( p_in, [l_filt], [lambda], [ds_sampl], [ikinks], [iaccel],
   %                                                           [use_wgt], [use_bspline], [ds_bspl], [idebug] );
   %
   % compute smoothed profile using parametric smoothing spline approximation with smoothness
   % parameters l_filt [mm] or lambda.
   %
   % p_in may be a struct ('ProfileY', 'ProfileZ') or array [y_prf, z_prf]
   %
   % re-sample at step-size ds_sampl, if given, or use s-positions of input profile
   %
   % enable detection of kinks (-1, default), disable ([]), or set manually [kink1, kink2,...]
   % enable detection of accelerations (-1, default), disable ([]) or set manually [accel1, accel2, ...]
   % use uniform weighting (use_wgt<=0), by segment length (1, default), increase wgt for curvature (2)
   % using PP-spline (use_bspline=0) or B-spline (use_bspline=1, default) with target step ds_bspl.

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
   % #smooth_profile -2

   if (nargin<2);l_filt = [];end
   if (nargin<3);lambda = [];end
   if (nargin<4);ds_sampl = [];end
   if (nargin<5)
    ikinks = -1;
   end
   if (nargin<6)
    iaccel = [];
   end
   if (nargin<7 | isempty(use_wgt))
    use_wgt = 1;
   end
   if (nargin<8 | isempty(use_bspline))
    use_bspline = 1;
   end
   if (nargin<9 | isempty(ds_bspl))
    ds_bspl = 2;
   end
   if (nargin<10 | isempty(idebug))
    idebug = 0;
   end

   is_wheel = 0;
   if (isstruct(p_in))
    is_wheel = p_in.is_wheel;
   end

   if (use_bspline)
    use_deriv = 3;
   else
    use_deriv = 2;
   end
   use_repl = 1;

   % compute lambda when l_filt is given

   if     ( isempty(l_filt) &  isempty(lambda))
    disp('Either L_filt or lambda must be given.');
    return;
   elseif (~isempty(l_filt) & ~isempty(lambda))
    disp('L_filt and lambda may not both be given.');
    return;
   elseif (~isempty(l_filt))
    if (use_deriv==3)
     lambda = l_filt^6 / (64*pi^6);
    else
     lambda = l_filt^4 / (16*pi^4);
    end
   end

   % obtain y, z and s (optional) from input profile

   if (isstruct(p_in))
    y_in = p_in.ProfileY;
    z_in = p_in.ProfileZ;
    if (isfield(p_in, 'ProfileS'))
     s_in = p_in.ProfileS;
    else
     s_in = [];
    end
    if (size(y_in,2)>1), y_in = y_in'; end
    if (size(z_in,2)>1), z_in = z_in'; end
    if (size(s_in,2)>1), s_in = s_in'; end
   else
    if (size(p_in,2)>size(p_in,1))
     p_in = p_in';
    end
    y_in = p_in(:,1);
    z_in = p_in(:,2);
    s_in = [];
   end

   % compute arc-lengths s if not provided in input profile

   if (isempty(s_in))
    s_in = cntc.make_arclength(y_in, z_in);
   end

   % estimate curvatures needed for adaptive weighting, using moving window [-sw,sw]

   imeth = 0; s_windw = 4; k_hlf = 4;
   curv = estimate_curv(p_in, imeth, [s_windw, k_hlf]);

   % compute weights on data
   %  0: uniform
   %  1: using w_i = ds_i
   %  2: w_i = ds_i * f, with f = f(rho)

   if (use_wgt<=0)
    wgt = ones(size(y_in));
    if (~use_bspline)
     wgt([1,end  ]) = wgt([1,end  ]) * 20;     % PP-spline: increase weight at end points
    end
   else
    ds_in = diff(s_in);
    ds_in = [ds_in(1); (ds_in(1:end-1)+ds_in(2:end))/2; ds_in(end)];
    wgt   = ds_in;
    if (~use_bspline)
     wgt([1,end  ]) = wgt([1,end  ]) * 20;     % PP-spline: increase weight at end points
    end

    if (use_wgt>=2)
     ix = find(abs(curv)>0.1);
     fac = 50; wgt(ix) = fac * wgt(ix);
     ix = find(abs(curv)>0.07 & abs(curv)<=0.1);
     fac = 15; wgt(ix) = fac * wgt(ix);
     ix = find(abs(curv)>0.035 & abs(curv)<=0.07);
     fac = 8; wgt(ix) = fac * wgt(ix);
     ix = find(abs(curv)>0.020 & abs(curv)<=0.035);
     fac = 8; wgt(ix) = fac * wgt(ix);
    end
   end

   % ikinks=-1: autodetect kinks between successive segments

   if (isscalar(ikinks) & ikinks<0)
    dalph_thrs_high = pi/6;
    dalph_thrs_low  = dalph_thrs_high / 5;
    dst_max         = 0.0;
    scale_z         = 1.0;
    is_wheel        = [];
    make_plot       = 0;
    ikinks = detect_kinks( s_in, y_in, z_in, dalph_thrs_high, dalph_thrs_low, dst_max, scale_z, ...
     is_wheel, make_plot);

    % PP-spline: increase weight at kink points
    if (~isempty(ikinks) & ~use_bspline)
     disp('Increasing weight * 100 at break points')
     wgt(ikinks) = 100 * wgt(ikinks);
    end
   end

   % create spline representation

   spl = cntc.make_spline( s_in, y_in, z_in, lambda, wgt, ikinks, iaccel, use_bspline, ds_bspl, ...
    use_deriv, use_repl, idebug );

   % evaluate spline at output positions

   if (~isempty(ds_sampl))
    is_min = ceil(min(spl.s)/ds_sampl);
    is_max = floor(max(spl.s)/ds_sampl);
    s_out = [is_min : is_max]' * ds_sampl;
    % disp(sprintf('sampling at is=[%d:%d], s=[%5.2f,%5.2f]', is_min, is_max, s_out(1), s_out(end)));
   elseif ( use_bspline )
    s_out = s_in;
   else
    s_out = spl.s;
   end

   [ ~, y_out, z_out ] = eval_spline( spl, s_out );

   % create profile structure

   if (isstruct(p_in))
    p_out = p_in;
   else
    p_out = struct();
   end
   p_out.ProfileS = s_out;
   p_out.ProfileY = y_out;
   p_out.ProfileZ = z_out;

   % estimate surface inclination
   % rail:   in [-pi, pi], wheel:  in [0, 2*pi]

   alph = atan2( diff(z_out), diff(y_out) );
   if (is_wheel)
    ix = find(alph<0); alph(ix) = alph(ix) + 2*pi;
   end
   alph = [alph(1); 0.5*(alph(1:end-1)+alph(2:end)); alph(end)];
   p_out.ProfileAngle = alph;

   % reset curvature

   p_out.ProfileCurvature = zeros(size(p_out.ProfileS));

  end % smooth_profile

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% #Solve
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [ u_out, found ] = solve_cubic_eq(ppcoef, useg, yev, idebug, iout, jout, jseg)

   % [ u_out, found ] = solve_cubic_eq(ppcoef, useg, yev, idebug, iout, jout, jseg)
   %
   % solve 3rd degree equation f(u) = y, returning one root in interval [useg(1),useg(2)]
   % ppcoef = [a0, a1, a2, a3],  d = a0-yev, c = a1, b = a2, a = a3
   % #solve_cubic_eq

   if (nargin<4 | isempty(idebug))
    idebug = 0;
   end
   if (nargin<5 | isempty(iout))
    iout   = 0;
   end
   if (nargin<6 | isempty(jout))
    jout   = 0;
   end
   if (nargin<7 | isempty(jseg))
    jseg   = 0;
   end
   if (jseg==-335), idebug = 5; end

   if (1==1)
    [ u_out1, found ] = cntc.solve_cubic_cardano(ppcoef, useg, yev, idebug, iout, jout, jseg);
    [ u_out2, found ] = cntc.solve_cubic_newton(ppcoef, useg, yev, idebug, iout, jout, jseg);
    if (abs(u_out1-u_out2)>1e-6)
     disp(sprintf('  jseg = %d: different solutions u_newt = %8.4f, u_card = %8.4f, diff = %3.1e', ...
      jseg, u_out2, u_out1, abs(u_out2-u_out1)));
    end
    u_out = u_out1;
   else
    [ u_out, found ] = cntc.solve_cubic_cardano(ppcoef, useg, yev, idebug, iout, jout, jseg);
    % [ u_out, found ] = cntc.solve_cubic_newton(ppcoef, useg, yev, idebug, iout, jout, jseg);
   end

  end % solve_cubic_eq

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ u_out, found ] = solve_cubic_newton(ppcoef, useg, yev, idebug, iout, jout, jseg)

   % [ u_out, found ] = #solve_cubic_newton(ppcoef, useg, yev, idebug, iout, jout, jseg)
   %
   % solve 3rd degree equation f(u) = y, returning one root in interval [useg(1),useg(2)]

   % solve f(u) = a3 ul^3 + a2 ul^2 + a1 ul + a0 - yev = 0, with ul the local coordinate in the segment
   % using Newton-Raphson, cannot jump over extremal values, needs appropriate u^0
   % #solve_cubic_newton

   a0 = ppcoef(1);
   a1 = ppcoef(2);
   a2 = ppcoef(3);
   a3 = ppcoef(4);

   if (idebug>=3)
    ul = [-100:0.01:20];
    u4 = useg(2) - useg(1);
    f_u = a3 * ul.^3 +   a2 * ul.^2 + a1 * ul + a0 - yev;

    figure(12); clf; hold on;
    plot(ul, f_u);
    v = axis;
    plot([0 0 NaN u4 u4], [v(3:4) NaN v(3:4)], '--');
    plot(v(1:2), 0*[1 1], '--');
    grid on;
    xlabel(sprintf('u_l in seg %d', jseg))
    ylabel('y_l = f(u_l)')
    title(sprintf('Newton for (i,j)=(%d,%d)',iout,jout));
    l = legend(sprintf('cubic, jseg=%d',jseg));
    % set(l, 'autoupdate','off');
   end

   tiny_a2 = 1e-10;
   tiny_a3 = 1e-10;

   % check for local extremal values:
   % zeros of quadratic f'(ul) = 3a ul^2 + 2b ul + c = 0

   discr = (2*a2)^2 - 4 * (3*a3) * a1;

   if (idebug>=3)
    disp(sprintf('  jseg = %d: a0 = %6.3e, a1 = %6.3e, a2 = %6.3e, a3 = %6.3e, D = %6.3e', ...
     jseg, a0, a1, a2, a3, discr));
   end

   u0 = 0;
   u4 = useg(2) - useg(1);

   if (abs(a2)<tiny_a2 & abs(a3)<tiny_a3)

    % segment function is linear, within roundoff precision; use initial estimate ul = 0

    if (idebug>=5), disp('  case 1'); end
    ul    = u0;
    f_u   = a0 - yev;
    df_du = a1;

   elseif (discr<=0)

    % no local extremal values: use initial estimate ul = 0

    if (idebug>=5), disp('  case 2'); end
    ul    = u0;
    f_u   = a0 - yev;
    df_du = a1;

   elseif (abs(a3)<tiny_a3)

    % segment function is quadratic, within roundoff precision

    if (idebug>=5), disp('  case 3'); end
    ul    = u0;
    f_u   = a0 - yev;
    df_du = a1;

    % check extremal value, should have opposite sign to a2

    u1    = -a1 / (2 * (2*a2));  % extremum
    f1    = a2 * u1^2 + a1 * u1 + a0 - yev;
    if (a2*f1 > 0)
     % disp('case 3: no solution')
     ul = NaN;
    end

   else

    u1 = ( -(2*a2) - sqrt(discr) ) / (2 * (3*a3));
    u2 = ( -(2*a2)               ) / (2 * (3*a3));  % inflection point
    u3 = ( -(2*a2) + sqrt(discr) ) / (2 * (3*a3));
    % swap u1,u3 when u1>u3
    if (u1>u3), tmp=u1; u1=u3; u3=tmp; end

    f0 = a0 - yev;
    f1 = a3 * u1^3 + a2 * u1^2 + a1 * u1 + a0 - yev;
    f2 = a3 * u2^3 + a2 * u2^2 + a1 * u2 + a0 - yev;
    f3 = a3 * u3^3 + a2 * u3^2 + a1 * u3 + a0 - yev;
    f4 = a3 * u4^3 + a2 * u4^2 + a1 * u4 + a0 - yev;

    if (idebug>=3)
     v = axis; dy = (v(4)-v(3)) / 5;
     plot(u0*[1 1], f0+dy*[-1,1],'b--'); text(u0, f0-1.1*dy, 's_0');
     plot(u1*[1 1], f1+dy*[-1,1],'r--'); text(u1, f1-1.1*dy, 's_1');
     plot(u2      , f2          ,'o'  ); text(u2, f2-1.1*dy, 's_2');
     plot(u3*[1 1], f3+dy*[-1,1],'r--'); text(u3, f3-1.1*dy, 's_3');
     plot(u4*[1 1], f4+dy*[-1,1],'b--'); text(u4, f4-1.1*dy, 's_4');
     % disp([a2, discr, f0, f1, f3, f4])
    end

    % there can be one zero before u1, one between u1 and u3, and one after u3
    % select initial estimate: start-point or end-point of segment jseg

    i1_in_u04 = ((u4-u1)*(u1-u0) > 0);
    i2_in_u04 = ((u4-u3)*(u3-u0) > 0);

    if (f1*f3>0 & a3*f1>0)
     % one zero, to the left of u1
     if (u1<=u0) % no solution within [u0,u4]
      if (idebug>=5), disp('  case 4.1.a'); end
      ul = u1+(u1-u0)+(u1-u2); f_u = a3 * ul.^3 +   a2 * ul.^2 + a1 * ul + a0 - yev;
     else
      if (idebug>=5), disp('  case 4.1.b'); end
      ul = u0; f_u = f0;
     end
    elseif (f1*f3>0)
     % one zero, to the right of u3
     if (u3>=u4) % no solution within [u0,u4]
      if (idebug>=5), disp('  case 4.2.a'); end
      ul = u3+(u3-u4)+(u3-u2); f_u = a3 * ul.^3 +   a2 * ul.^2 + a1 * ul + a0 - yev;
     else
      if (idebug>=5), disp('  case 4.2.b'); end
      ul = u4; f_u = f4;
     end
    elseif (i1_in_u04 & i2_in_u04)
     % three zeros, middle one between u0 and u4
     if (idebug>=5), disp('  case 4.3'); end
     ul = (u1+u3) / 2;
     f_u = a3 * ul^3 + a2 * ul^2 + a1 * ul + a0 - yev;
    elseif ( (i1_in_u04 & f0*f1>0) | (i2_in_u04 & f0*f3>0))
     if (idebug>=5), disp('  case 4.4'); end
     % three zeros, right one between u0 and u4
     ul = u4; f_u = f4;
    else
     if (idebug>=5), disp('  case 4.5'); end
     % three zeros, left one between u0 and u4
     ul = u0; f_u = f0;
    end
    df_du = 3 * a3 * ul^2 + 2 * a2 * ul + a1;
   end

   % Newton-Raphson: f(u+du) = f(u) + du * f'(u) = 0
   %                      -->  du = - f(u) / f'(u)

   iter  = 0; tol = 1e-8; maxit = 20;
   if (~isnan(ul) & idebug>=2)
    disp(sprintf(' It %2d: ul = %8.4f in [%8.4f,%8.4f], f(u) = %12.8f, f''(u) = %8.4f', ...
     iter, ul, u0, u4, f_u, df_du))
    if (idebug>=3)
     plot(ul, f_u, '*')
     text(ul, f_u, num2str(iter))
    end
   end

   while(~isnan(ul) & abs(f_u)>=tol & iter<maxit)
    iter  = iter + 1;
    du    = - f_u / df_du;

    % restrict step esp. in first iterations, cf. doubling of step until bracket is found
    du_max = (u4 - u0) * 2^iter;
    if (du>du_max)
     du = du_max;
    elseif (du<-du_max)
     du = -du_max;
    end

    % prevent jumping over u2 in first iteration
    if (iter<=1 & discr>0 & exist('u2'))
     if (ul<u2)
      ul = min(u2, ul+du);
     else
      ul = max(u2, ul+du);
     end
    else
     ul = ul + du;
    end

    % evaluate function and derivative at new position
    f_u   =   a3 * ul^3 +   a2 * ul^2 + a1 * ul + a0 - yev;
    df_du = 3*a3 * ul^2 + 2*a2 * ul   + a1;
    if (idebug>=2)
     disp(sprintf(' It %2d: ul = %8.4f in [%8.4f,%8.4f], f(u) = %12.8f, f''(u) = %8.4f', ...
      iter, ul, u0, u4, f_u, df_du))
    end
    if (idebug>=3)
     plot(ul, f_u, '*')
     text(ul, f_u, num2str(iter))
    end
   end

   if (isnan(ul))
    if (idebug>=1)
     disp(sprintf('solve_newton:  no solution for jout=%3d, jseg=%3d, yev=%8.3f', jout, jseg, yev));
    end
    found = 0;
   elseif (abs(f_u)>=tol)
    disp(sprintf('solve_newton:  no convergence for jout=%3d, jseg=%3d, yev=%8.3f', jout, jseg, yev));
    found = 0;
   elseif (ul<u0 | ul>u4)
    if (idebug>=1)
     disp(sprintf('solve_newton:  ul=%8.4f outside segment [%8.3f,%8.3f] for jout=%3d, yev=%8.3f', ...
      ul, u0, u4, jout, yev));
    end
    found = 0;
   else
    found = 1;
   end

   u_out = useg(1) + ul;

  end % solve_cubic_newton

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [ u_out, found ] = solve_cubic_cardano(ppcoef, useg, yev, idebug, iout, jout, jseg)

   % [ u_out, found ] = #solve_cubic_cardano(ppcoef, useg, yev, idebug, iout, jout, jseg)
   %
   % solve 3rd degree equation f(u) = y, returning one root in interval [useg(1),useg(2)]

   % solve f(u) = a ul^3 + b ul^2 + c ul + d = 0,  d = a3-yev,
   % with ul the local coordinate in the segment, using Cardano's method.
   % #solve_cubic_cardano

   d = ppcoef(1) - yev;
   c = ppcoef(2);
   b = ppcoef(3);
   a = ppcoef(4);

   if (idebug>=4)
    ul = [-25:0.01:15];
    u4 = useg(2) - useg(1);
    f_u = a * ul.^3 + b * ul.^2 + c * ul + d;

    figure(12); clf; hold on;
    plot(ul, f_u);
    v = axis;
    plot([0 0 NaN u4 u4], [v(3:4) NaN v(3:4)], '--');
    plot(v(1:2), 0*[1 1], '--');
    grid on;
    xlabel(sprintf('u_l in seg %d', jseg))
    ylabel('y_l = f(u_l)')
    title(sprintf('Cardano for (i,j)=(%d,%d)',iout,jout));
    l = legend(sprintf('cubic, jseg=%d',jseg));
    % set(l, 'autoupdate','off');
   end

   if (idebug>=4)
    disp(sprintf('  jseg = %d: a = %6.3e, b = %6.3e, c = %6.3e, d = %6.3e', jseg, a, b, c, d));
   end

   tiny_a   = 1e-14;
   tiny_b   = 1e-10;
   tiny_cfc = 1e-20;
   u0     = 0;                  % interval [u0,u4] in local coordinates
   u4     = useg(2) - useg(1);

   if (abs(b)<tiny_b & abs(a)<tiny_a)

    % segment function is linear, within roundoff precision

    if (idebug>=6), disp('  case 1'); end
    ul    = -d / c;

   elseif (abs(a)<tiny_a)

    % segment function is quadratic, within roundoff precision

    if (idebug>=6), disp('  case 2'); end

    a = b; b = c; c = d;
    discr = b^2 - 4*a*c;
    if (discr<0)
     ul = NaN;
    else
     x1 = (-b - sqrt(discr)) / (2*a);
     x2 = (-b + sqrt(discr)) / (2*a);
     if (x1>=u0 & x1<=u4)
      ul = x1;
     else
      ul = x2;
     end
    end

   else

    % Cardano's method

    delta0 = b^2 - 3*a*c;
    delta1 = 2*b^3 - 9*a*b*c + 27*a^2*d;
    discr  = (4*delta0^3 -  delta1^2) / (27*a^2);

    if (abs(discr)>tiny_cfc)
     coefc  = ( ( delta1 + sqrt(delta1^2 - 4*delta0^3) ) / 2 )^(1/3);
     if (abs(coefc)<tiny_cfc)
      disp('Impossible coefc=0?')
      coefc  = ( ( delta1 - sqrt(delta1^2 - 4*delta0^3) ) / 2 )^(1/3);
     end
     ksi  = (-1 + sqrt(-3)) / 2;
    end

    if (abs(discr)<tiny_cfc & abs(delta0)<tiny_cfc)

     % triple root

     x1 = -b / (3*a);
     x2 = x1;
     x3 = x1;
     ul = x1;
     if (idebug>=3)
      disp('cardano: triple root');
     end

    elseif (abs(discr)<tiny_cfc)

     % double root + single root
     disp(discr)

     x1 = (4*a*b*c - 9*a^2*d - b^3) / (a*delta0); % single
     x2 = (9*a*d - b*c) / (2*delta0);             % double
     x3 = x2;
     if (x1>=u0 & x1<=u4)
      ul = x1;
     else
      ul = x2;
     end
     if (idebug>=3)
      disp('cardano: double root');
     end

    elseif (discr<0)

     % one real root + two non-real complex conjugate roots
     % coefc can be complex especially when delta1 < 0.

     x1 = -( b + ksi  *coefc + delta0 / (ksi  *coefc)) / (3*a);
     x2 = -( b + ksi^2*coefc + delta0 / (ksi^2*coefc)) / (3*a);
     x3 = -( b +       coefc + delta0 /        coefc ) / (3*a);

     % select root with smallest imaginary part (theoretically zero)

     x      = [x1, x2, x3];
     [~, i] = min(abs(imag(x)));
     ul = x(i);

     if (idebug>=3)
      disp('cardano: one real root & two non-real complex conjugate roots');
     end

    else

     % three distinct real roots

     x1 = -( b + ksi  *coefc + delta0 / (ksi  *coefc)) / (3*a);
     x2 = -( b + ksi^2*coefc + delta0 / (ksi^2*coefc)) / (3*a);
     x3 = -( b +       coefc + delta0 /        coefc ) / (3*a);

     x = [x1 x2 x3];
     if (max(abs(imag(x)))>1e-12*max(abs(real(x))))
      disp(sprintf('Impossible: nonzero imag.part (%3.1e)???', max(abs(imag(x))) ))
      disp(x)
     end
     x = real(x);

     % select a root in interval [u0,u4] or closest to interval [u0,u4]

     xm = (u0 + u4) / 2;
     [~,i] = min(abs( x-xm ));
     ul = x(i);

     if (idebug>=3)
      disp('cardano: three distinct real roots');
     end

    end

    f1     = a * x1.^3 + b * x1.^2 + c * x1 + d;
    f2     = a * x2.^3 + b * x2.^2 + c * x2 + d;
    f3     = a * x3.^3 + b * x3.^2 + c * x3 + d;

    if (idebug>=4)
     disp('solutions x1, x2, x3:')
     disp([x1, x2, x3; f1, f2, f3]);
    end

   end

   if (abs(imag(ul))>1e-10)
    disp('Impossible: complex ul???')
   end
   u_out = useg(1) + real(ul);

   if (isnan(ul))
    if (idebug>=3)
     disp(sprintf('  solve_cardano: no solution for jout=%3d, jseg=%3d, yev=%8.3f', jout, jseg, yev));
    end
    found = 0;
   elseif (ul<u0 | ul>u4)
    if (idebug>=3)
     disp(sprintf('  solve_cardano: u=%8.4f outside segment [%8.3f,%8.3f] for jout=%3d, yev=%8.3f', ...
      real(ul), u0, u4, jout, yev));
    end
    found = 0;
   else
    found = 1;
   end

  end % solve_cubic_cardano

  %% #write

  function [ ] = write_miniprof(p, fname, is_rail, use_yz, idebug)

   % [ ] = write_miniprof(p, fname, [is_rail], [use_yz], [idebug])
   %
   % write profile to file in Miniprof format
   %  - default (use_yz=0) : write data given in OrigData
   %  - option  (use_yz>0) : overwrite OrigData with ProfileY and ProfileZ
   %
   % use e.g. p = struct('OrigData', [y, z], 'ColumnDef', 'X,Y');

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
   % #write_miniprof

   if (nargin<3 | isempty(is_rail))
    is_rail = -1;
   end
   if (nargin<4 | isempty(use_yz))
    use_yz = 0;
   end
   if (nargin<5 | isempty(idebug))
    idebug = 0;
   end

   if (exist(fname,'file'))
    disp(['ERROR: file ',fname,' already exists, please remove']);
    return
   end

   % use ProfileY, ProfileZ if so requested, else OrigData will be used

   if (use_yz)
    p.OrigData = [ p.ProfileY, p.ProfileZ ];
    p.ColumnDef = 'X,Y';
   end

   f = fopen(fname, 'w');

   known_fields = struct( ...
    'ProfileData',                   'internal', ...
    'ProfileY',                      'internal', ...
    'ProfileZ',                      'internal', ...
    'ProfileS',                      'internal', ...
    'ProfileAngle',                  'internal', ...
    'ProfileCurvature',              'internal', ...
    'ProfileKYield',                 'internal', ...
    'OrigData',                      'internal', ...
    'Mask',                          'internal', ...
    'Fname',                         'internal', ...
    'FileName',                      'str', ...
    'ProgramName',                   'str', ...
    'ProgramDate',                   'str', ...
    'ProgramVer',                    'str', ...
    'Date',                          'str', ...
    'Time',                          'str', ...
    'MPCalDate',                     'str', ...
    'MPCalTime',                     'str', ...
    'UserName',                      'str', ...
    'EmplNum',                       'int', ...
    'Track',                         'str', ...
    'Line',                          'str', ...
    'Direction',                     'str', ...
    'Chainage',                      'str', ...
    'Curve',                         'str', ...
    'Position',                      'str', ...
    'Location_Latitude',             'real', ...
    'Location_Longitude',            'real', ...
    'Location_Altitude',             'real', ...
    'Location_HorizontalAccuracy',   'real', ...
    'Location_VerticalAccuracy',     'real', ...
    'Rail',                          'str', ...
    'TrackSection',                  'unknown', ...
    'KP',                            'unknown', ...
    'Stock',                         'str', ...
    'CarNo',                         'int', ...
    'Bogie',                         'int', ...
    'SurfaceCondition',              'unknown', ...
    'Side',                          'str', ...
    'AxleNo',                        'int', ...
    'Mileage',                       'str', ...
    'WheelID',                       'int', ...
    'RCF',                           'unknown', ...
    'Flat',                          'unknown', ...
    'ReferenceProfile',              'str', ...
    'ProfileAlignment',              'unknown', ...
    'Instrument_Firmware',           'str', ...
    'Instrument_Connection',         'int', ...
    'Windows_Version',               'str', ...
    'Windows_ServicePack',           'str', ...
    'Windows_NetFramework',          'str', ...
    'Windows_ComputerName',          'str', ...
    'Windows_UserName',              'str', ...
    'Battery_Level',                 'real', ...
    'Battery_Capacity',              'real', ...
    'Comment',                       'str', ...
    'Flag',                          'unknown', ...
    'MPTypeNo',                      'int', ...
    'MPSerNo',                       'int', ...
    'Xoffset',                       'real', ...
    'Yoffset',                       'real', ...
    'Elevation',                     'real', ...
    'Tilt',                          'unknown', ...
    'AreaRef',                       'unknown', ...
    'AreaGain',                      'unknown', ...
    'AreaLoss',                      'unknown', ...
    'BaseTilt',                      'unknown', ...
    'BaseOffset',                    'unknown', ...
    'RefPoint1',                     'real', ...
    'RefPoint2',                     'real', ...
    'RefPoint3',                     'real', ...
    'RefPoint4',                     'real', ...
    'RefPoint5',                     'real', ...
    'RefPoint6',                     'real', ...
    'RefPoint7',                     'real', ...
    'RefPoint8',                     'real', ...
    'RefPoint9',                     'real', ...
    'SectionNo',                     'unknown', ...
    'CrownR',                        'real', ...
    'CrownRadius',                   'real', ...
    'TransformState',                'str', ...
    'Transformation',                'str', ...
    'ProfileShaping',                'unknown', ...
    'OriginalHeader',                'unknown', ...
    'ColumnDef',                     'str', ...
    'XYPoints',                      'int', ...
    'RailAlignType',                 'int', ...
    'RailAngle',                     'real', ...
    'RailAngle_AlarmStatus',         'int', ...
    'AlignDX',                       'unknown', ...
    'AlignDY',                       'unknown', ...
    'AlignDA',                       'unknown', ...
    'AlignOX',                       'unknown', ...
    'AlignOY',                       'unknown', ...
    'W1',                            'real', ...
    'W2',                            'real', ...
    'W3',                            'real', ...
    'W1_AlarmLevel',                 'real', ...
    'W1_AlarmStatus',                'int', ...
    'W2_AlarmStatus',                'int', ...
    'W3_AlarmStatus',                'int', ...
    'UOF1',                          'str', ...
    'UOF2',                          'str', ...
    'UOF3',                          'str', ...
    'UOF4',                          'str', ...
    'Chars',                         'str', ...
    'Temperature',                   'real', ...
    'GradeRatio',                    'real', ...
    'Grade',                         'real', ...
    'SuperElevationHeight',          'real', ...
    'SuperElevation',                'real', ...
    'RailGaugePoint',                'str', ...
    'RailRodLength',                 'str', ...
    'RodLength',                     'real', ...
    'GaugeRodLength',                'real', ...
    'GaugeSensorDistance',           'real', ...
    'GaugeSensorStdDev',             'real', ...
    'GaugeSensorSamples',            'real', ...
    'WheelAdjustPoint',              'real', ...
    'WheelDiameterTaperline',        'real', ...
    'AlarmWarningCount',             'int', ...
    'AlarmFailureCount',             'int', ...
    'Gauge',                         'real', ...
    'Gauge_AlarmStatus',             'int', ...
    'DiameterFlange',                'real', ...
    'DiameterTaperline',             'real', ...
    'DiameterFlange_AlarmStatus',    'int', ...
    'DiameterTaperline_AlarmStatus', 'int', ...
    'Sd',                            'real', ...
    'Sh',                            'real', ...
    'qR',                            'real', ...
    'Sd_AlarmStatus',                'int', ...
    'Sh_AlarmStatus',                'int', ...
    'qR_AlarmStatus',                'int'  ...
    );

   % first write all the metadata to the file: keyword=value

   fields_p = fieldnames(p);
   num_fields = length(fields_p);

   for ifld = 1 : num_fields

    my_field_name = cntc.check_field(known_fields, fields_p(ifld));

    if (~isempty(my_field_name))
     % interpret numerical value
     if (strcmp(known_fields.(my_field_name), 'internal'))
      % disp(sprintf('Generated field "%s"', fields_p{ifld}));
     elseif (strcmp(known_fields.(my_field_name), 'int'))
      nval = length(p.(my_field_name));
      fmt = '%s=%d';
      for ival = 2 : nval
       fmt = [fmt, ', %d'];
      end
      fmt = [fmt, '\n'];
      fprintf(f, fmt, my_field_name, p.(my_field_name) );
     elseif (strcmp(known_fields.(my_field_name), 'real') | ...
       strcmp(known_fields.(my_field_name), 'num'))
      % TODO: arrays!
      nval = length(p.(my_field_name));
      fmt = '%s=%6.4f';
      for ival = 2 : nval
       fmt = [fmt, ', %6.4f'];
      end
      fmt = [fmt, '\n'];
      fprintf(f, fmt, my_field_name, p.(my_field_name) );
      % store string value
     elseif (strcmp(known_fields.(my_field_name), 'str') | ...
       strcmp(known_fields.(my_field_name), 'unknown'))
      fprintf(f, '%s=%s\n', my_field_name, p.(my_field_name) );
      % unknown type: error
     else
      disp(sprintf('Field %s has unknown field-type %s.',my_field_name, ...
       known_fields.(my_field_name) ));
     end
     % unknown keyword: report
    else
     disp(sprintf('Unknown field "%s"',fields_p{ifld}));
    end
   end % all fields

   % create format on basis of ColumnDef

   coldef = upper(p.ColumnDef);

   fmt_X = ' %8.4f';
   fmt_Y = ' %8.4f';
   fmt_A = ' %7.4f';
   fmt_C = ' %18.15f';
   fmt_N = ' %3d';
   fmt_K = ' %7.3f';
   str   = ['fmt_', strrep(coldef, ',', ', fmt_') ];
   fmt   = eval([ '[',str,']' ]);
   fmt   = [fmt, ' \n'];

   % Write the profiledata to the file

   fprintf(f, '\n');
   for ipoint = 1 : size( p.OrigData,1 )
    fprintf(f, fmt, p.OrigData(ipoint,:));
   end

   fclose(f);


  end % write_miniprof


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function write_simpack(fname, comment, y, z)

   % #write_simpack
   %
   % Write wheel or rail profile in Simpack prw/prr format

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
   % #write_simpack

   if (isempty(strfind(fname, 'prw')))
    is_wheel = 0; % rail
   else
    is_wheel = 1; % wheel
   end

   % make column vectors

   if (size(y,2)>1), y=y'; end
   if (size(z,2)>1), z=z'; end
   if (max(size(y,2), size(z,2))>1),
    disp('ERROR: profile y, z should be 1D vectors')
    disp(size(y))
    disp(size(z))
    return
   end
   if (length(y)~=length(z))
    disp(sprintf('ERROR: profile y (%d), z (%d) should be the same length', ...
     length(y), length(z)))
    return
   end

   % make y-data decreasing for wheel

   if (is_wheel & y(end)-y(1)>0)
    y = flipud(y);  z = flipud(z);
   end

   f = fopen(fname, 'w');
   fprintf(f,'header.begin\n');
   fprintf(f,'  version  = 1\n');
   fprintf(f,'  type     = %d\n', is_wheel);
   fprintf(f,'header.end\n');
   fprintf(f,'spline.begin\n');
   fprintf(f,'  approx.smooth = 0.000000\n');
   fprintf(f,'  comment       = ''%s''\n', comment);
   fprintf(f,'  units.len.f   = 1000.\n');
   fprintf(f,'  units.ang.f   =    1.\n');
   fprintf(f,'  point.begin\n');
   for iy = 1 : length(y)
    fprintf(f,'   %11.6f %11.6f\n',y(iy), z(iy));
   end
   fprintf(f,'  point.end\n');
   fprintf(f,'spline.end\n');
   fclose(f);

  end % write_simpack

  %------------------------------------------------------------------------------------------------------------

  function [ ierror]=calculate(ire, icp, idebug)
   % [ ierror ] = cntc.calculate(ire, icp, idebug)
   %
   % perform actual CONTACT calculation for a contact problem
   % in:   integer    idebug       - show warnings (2), errors (1) or hide all messages (0)

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 5: m=*, wtd    - default icp=-1
   %% #calculate

   cntc_err_allow  = -12;
   cntc_err_profil = -32;

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = -1; % default: W/R contact, all patches
   end
   if (nargin<3 | isempty(idebug))
    idebug = 1;
   end

   p_ierr = libpointer('int32Ptr',-1);

   cntc.call('cntc.calculate', ire, icp, p_ierr);

   ierror = double(p_ierr.value);
   if (idebug>=1 & ierror==cntc_err_allow)
    disp(sprintf('cntc.calculate: no valid license found for CONTACT library (%d).',ierror));
   elseif (idebug>=1 & ierror==cntc_err_profil)
    disp(sprintf('cntc.calculate: an error is found in the rail or wheel profile specification (%d).',ierror));
   elseif (idebug>=1 & ierror<0)
    disp(sprintf('cntc.calculate: an error occurred in the CONTACT library (%d).',ierror));
   elseif (idebug>=2 & ierror>0)
    disp(sprintf('cntc.calculate: potential contact may be too small: there are %d points adjacent to the boundaries.',ierror));
   end

  end % cntc.calculate

  %% #setcmd : set/get related commands -------------------------------------------

  function result = mygetfield(struct,field);
   %
   % [result] = mygetfield(struct,field);
   %
   % returns the field <field> from the struct <struct>,
   % even if this field contains special characters.
   %
   % INTENDED USE:
   %      a = mygetfield(struct,'substruct.field');
   %
   % Other possible uses (not intended):
   %      a = mygetfield(struct,'substruct.field(3,4)');
   %

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
   % #mygetfield -2

   result = eval(['struct.',field]);
  end

  function contains = myisfield(STRUC,fieldname)
   %
   % [contains] = myisfield(STRUC,fieldname)
   %
   % myisfield extension of isfield, for which the fieldname may contain '.'
   % Example:
   %    STRUC.substruc.field = [0 1];
   %    myisfield(STRUC,'substruc.field')   returns 1
   %    myisfield(STRUC,'substruc')         returns 1,
   %    myisfield(STRUC,'')                 returns 1,
   %    myisfield(STRUC,any_other_string)   returns 0.
   %
   %   1.0     2004-07-28   BvtH (VORtech)   initial version

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
   % #myisfield -2

   % Initialize: assume that the field is in te struct.

   contains   = 1;
   fieldname  = deblank(fieldname);
   structname ='STRUC';

   % Check all the 'words' in the fieldname

   while (contains & ~strcmp(fieldname,''))
    [word,fieldname] = strtok(fieldname,'.');

    % if the field name contains words after this one,
    %  remove the '.' at the start (it is not removed by strtok)
    if (length(fieldname)>0)
     fieldname = fieldname(2:end);
    end

    fieldname = deblank(fieldname);
    word      = deblank(word);

    % Check if the struct contains a field whose name corresponds to
    % the word
    contains = contains & eval(['isfield(',structname,',word)']);

    % Continue checking inside the field.
    structname = [structname,'.',word];
   end
  end

  function s = mysetfield(s,varargin)
   %
   % [s] = mysetfield(s, varargin)
   %
   % Set structure field contents.
   %
   %   S = MYSETFIELD(S,'field',V) sets the contents of the specified
   %   field to the value V.  This is equivalent to the syntax S.field = V.
   %   S must be a 1-by-1 structure.  The changed structure is returned.
   %
   %   S = MYSETFIELD(S,{i,j},'field',{k},V) is equivalent to the syntax
   %       S(i,j).field(k) = V;
   %   In other words, S = SETFIELD(S,sub1,sub2,...,V) sets the
   %   contents of the structure S to V using the subscripts or field
   %   references specified in sub1,sub2,etc.  Each set of subscripts in
   %   parentheses must be enclosed in a cell array and passed to
   %   SETFIELD as a separate input.  Field references are passed as
   %   strings.
   %
   %   See also MYGETFIELD, GETFIELD, RMFIELD, MYISFIELD, ISFIELD, FIELDNAMES.

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
   % #mysetfield -2

   if (isempty(varargin) | length(varargin) < 2)
    error('MYSETFIELD:InsufficientInputs', 'Not enough input arguments.');
   end

   % The most common case
   arglen = length(varargin);
   strField = varargin{1};
   if (arglen==2)
    eval(['s.',deblank(strField),' = varargin{end};']);
    return
   end

   subs = varargin(1:end-1);
   for i = 1:arglen-1
    index = varargin{i};
    if (isa(index, 'cell'))
     types{i} = '()';
    elseif isstr(index)
     types{i} = '.';
     subs{i} = deblank(index); % deblank field name
    else
     error('MYSETFIELD:InvalidType','Inputs must be either cell arrays or strings.');
    end
   end

   % Perform assignment
   try
    s = builtin('subsasgn', s, struct('type',types,'subs',subs), varargin{end});
   catch
    error('MYSETFIELD', lasterr)
   end
  end

  %------------------------------------------------------------------------------------------------------------

  %% #get
  
  function [ tcpu, twall]=getcalculationtime(ire, icp)
   % [ tcpu, twall ] = cntc.getcalculationtime(ire, icp)
   %
   % return accumulated cpu-time and wall-clock-time used since last timer reset for a contact problem
   %  tcpu, twall   - cpu- and wall-clock times used [time]


   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 6: m=*, cp     - require icp>0, default 1

   % #getcalculationtime

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.getcalculationtime: not available for icp=%d',icp));
    return
   end

   p_tcpu  = libpointer('doublePtr',-1);
   p_twall = libpointer('doublePtr',-1);

   cntc.call('cntc.getcalculationtime', ire, icp, p_tcpu, p_twall);

   tcpu  = p_tcpu.value;
   twall = p_twall.value;

  end % cntc.getcalculationtime

  %------------------------------------------------------------------------------------------------------------

  function [fn,tx,ty,mz]=getcontactforces(ire, icp, LI)
   % [fn,tx,ty,mz] = cntc.getcontactforces(ire, icp)
   %
   % return the total forces and torsional moment for a contact problem in contact local coordinates
   %
   %  fn            - total normal force [force]
   %  tx, ty        - total tangential forces [force]
   %  mz            - total torsional moment [force.length]


   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 6: m=*, cp     - require icp>0, default 1
   % # getcontactforces

   l1={1,'fn','total normal force'
    2,'tx','total tangential forces'
    3,'ty','total tangential forces'
    4,'mz','total torsional moment'};

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.getcontactforces: not available for icp=%d',icp));
    return
   end

   p_fn = libpointer('doublePtr',-1);
   p_tx = libpointer('doublePtr',-1);
   p_ty = libpointer('doublePtr',-1);
   p_mz = libpointer('doublePtr',-1);

   cntc.call('cntc.getcontactforces', ire, icp, p_fn, p_tx, p_ty, p_mz);

   fn = p_fn.value;
   tx = p_tx.value;
   ty = p_ty.value;
   mz = p_mz.value;

   if nargin ==3
    values=[fn,tx,ty,mz];
    i2=LI.Cmacro.rowM(l1(:,2));
    LI.Cmacro.Y(i2,ire,LI.cur.j1)=values(vertcat(l1{:,1}));
   end
  end % cntc.getcontactforces

  function [ rvalues]=getcontactlocation(ire, icp,LI)
   % [ rvalues ] = cntc.getcontactlocation(ire, icp)
   %
   % return the contact reference location for a contact problem (module 1 only)
   %
   %  rvalues(lenarr)   - contact location parameters.
   %
   %  Values are listed with cntc.getcontactlocation
   %
   %  The "contact reference point" is the origin of the contact local coordinate system. It is determined by
   %  a heuristic rule and is centered within the contact patch in a weighted sense.
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
   % #getcontactlocation

   l1={
    1,'XCP_TR','x-position of the contact reference point in track coordinates'
    2,'YCP_TR','y-position of the contact reference point in track coordinates'
    3,'ZCP_TR','z-position of the contact reference point in track coordinates'
    4,'DELTCP_TR','contact reference angle: rotation about track x-axis from the track positive z-axis to the contact positive n-axis, with sign according the right-hand rule'
    %
    5,'XCP_R','x-position of the contact reference point in rail profile coordinates'
    6,'YCP_R','y-position of the contact reference point in rail profile coordinates'
    7,'ZCP_R','z-position of the contact reference point in rail profile coordinates'
    8,'SCP_R','s-parameter of the contact reference point, measured along the rail profile'
    9,'DELTCP_R','rotation about rail x-axis from rail positive z-axis to contact positive n-axis'

    10,'XCP_W','x-position of the contact reference point in wheel profile coordinates'
    11,'YCP_W','y-position of the contact reference point in wheel profile coordinates'
    12,'ZCP_W','z-position of the contact reference point in wheel profile coordinates'
    13,'SCP_W','s-parameter of the contact reference point, measured along the wheel profile'
    14,'DELTCP_W','rotation about wheel x-axis from wheel positive z-axis to contact positive n-axis'

    15,'XPN_TR','x-position of the pressure center of gravity in track coordinates'
    16,'YPN_TR','y-position of the pressure center of gravity in track coordinates'

    21,'XW_TR','x-position of wheel profile marker in track coordinates'
    22,'YW_TR','y-position of wheel profile marker in track coordinates'
    23,'ZR_TR','z-position of wheel profile marker in track coordinates'
    24,'ROLLW_TR','roll angle of wheel profile marker in track coordinates'
    25,'YAWW_TR','yaw angle of wheel profile marker in track coordinates'

    26,'YR_TR','y-position of rail profile marker in track coordinates'
    27,'ZR_TR','z-position of rail profile marker in track coordinates'
    28,'ROLLR_TR','roll angle of rail profile marker in track coordinates' % #roll_tr -2

    31,'DY_DEFL','lateral rail shift according to massless rail deflection'
    32,'DZ_DEFL','vertical rail shift according to massless rail deflection'
    };l1(:,2)=lower(l1(:,2));
   if nargin==0
    if nargout==0;disp(l1);else; rvalues=l1; end
    return
   end

   % category 3: m=1, cp     - require icp>0, default 1

   if (nargin<1 | isempty(ire));ire = 1;end
   if (nargin<2 | isempty(icp)); icp = 1;end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.getcontactlocation: not available for icp=%d',icp));
    return
   end

   lenarr = 32;
   p_values = libpointer('doublePtr',zeros(lenarr,1));

   cntc.call('cntc.getcontactlocation', ire, icp, lenarr, p_values);

   rvalues = p_values.value;
   if nargin==3
    i2=LI.Cmacro.rowM(l1(:,2));
    LI.Cmacro.Y(i2,ire,LI.cur.j1)=rvalues(vertcat(l1{:,1}));
    %rvalues=sdtm.toStruct([l1(:,2) num2cell()]);
   end

  end % cntc.getcontactlocation

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [ carea, harea, sarea, parea]=getcontactpatchareas(ire, icp)
   % [ carea, harea, sarea, parea ] = cntc.getcontactpatchareas(ire, icp)
   %
   % return the area of contact for a contact problem
   %  carea   - area of contact patch [area]
   %  harea   - area of adhesion area [area]
   %  sarea   - area of slip area [area]
   %  parea   - area of plasticity area [area]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 6: m=*, cp     - require icp>0, default 1
   % #getcontactpatchareas

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.getcontactpatchareas: not available for icp=%d',icp));
    return
   end

   p_carea  = libpointer('doublePtr',-1);
   p_harea  = libpointer('doublePtr',-1);
   p_sarea  = libpointer('doublePtr',-1);

   cntc.call('cntc.getcontactpatchareas', ire, icp, p_carea, p_harea, p_sarea);

   carea = p_carea.value;
   harea = p_harea.value;
   sarea = p_sarea.value;
   parea = carea - harea - sarea;

  end % cntc.getcontactpatchareas

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [ sol]=getcpresults(ire, icp)
   % [ sol ] = cntc.getcpresults(ire, icp)
   %
   % get results of a single contact patch in the form used by loadcase and plot3d
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 6: m=*, cp     - require icp>0, default 1
   % #getcpresults


   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.getcpresults: not available for icp=%d',icp));
    return
   end

   sol         = struct();

   CNTC = cntc.getmagicnumbers();

   % retrieve control digits needed for plotting

   iparam(1) = CNTC.ic_config;
   iparam(2) = CNTC.ic_tang;
   iparam(3) = CNTC.ic_frclaw;
   iparam(4) = CNTC.ic_discns;
   iparam(5) = CNTC.ic_mater;
   iparam(6) = CNTC.ic_heat;
   values = cntc.getflags(ire, icp, iparam);

   sol.config         = values(1);
   sol.kincns         = struct();
   sol.kincns.t_digit = values(2);
   % sol.fric           = struct();
   % sol.fric.frclaw    = values(3);
   sol.d_digit        = values(4);
   sol.mater          = struct();
   sol.mater.m_digit  = values(5);
   sol.h_digit        = values(6);

   % get material / kinematic parameters needed

   values = cntc.getparameters(ire, icp);
   sol.kincns.veloc  = values.veloc;
   sol.kincns.chi    = values.chi;
   sol.kincns.dq     = values.dq;
   sol.meta          = struct();
   sol.meta.spinxo   = values.spinxo;
   sol.meta.spinyo   = values.spinyo;
   sol.mater.tau_c0  = values.tau_c0;
   use_plast = (sol.mater.m_digit==4 & (sol.mater.tau_c0>1e-10 & sol.mater.tau_c0<1e10));

   % retrieve wheel and rail profile data

   itask   = 1;
   iswheel = 0;
   %isampl  = 0; %Spline representation
   isampl = -1; %Original data
   iparam  = [iswheel, isampl];
   rparam  = [];
   val = cntc.getprofilevalues(ire, itask, iparam, rparam);
   sol.prr = struct('ProfileY',val(:,1), 'ProfileZ',val(:,2));
   itask   = 4;
   sol.prr.ProfileS = cntc.getprofilevalues(ire, itask, iparam, rparam);

   itask   = 1;
   iswheel = 1;
   iparam  = [iswheel, isampl];
   val = cntc.getprofilevalues(ire, itask, iparam, rparam);
   sol.prw = struct('ProfileY',val(:,1), 'ProfileZ',val(:,2));
   itask   = 4;
   sol.prw.ProfileS = cntc.getprofilevalues(ire, itask, iparam, rparam);

   % retrieve wheel, rail and contact position data

   ws_pos = cntc.getwheelsetposition(ire);
   cp_pos = cntc.getcontactlocation(ire, icp);
   sol.meta.s_ws     = ws_pos( 1);
   sol.meta.x_w      = cp_pos(21);
   sol.meta.y_w      = cp_pos(22);
   sol.meta.z_w      = cp_pos(23);
   sol.meta.roll_w   = ws_pos( 4);
   sol.meta.yaw_w    = ws_pos( 5);
   sol.meta.y_r      = cp_pos(26);
   sol.meta.z_r      = cp_pos(27);
   sol.meta.roll_r   = cp_pos(28);
   sol.meta.xcp_r    = cp_pos( 5);
   sol.meta.ycp_r    = cp_pos( 6);
   sol.meta.zcp_r    = cp_pos( 7);
   sol.meta.deltcp_r = cp_pos( 9);
   sol.meta.xcp_w    = cp_pos(10);
   sol.meta.ycp_w    = cp_pos(11);
   sol.meta.zcp_w    = cp_pos(12);
   sol.meta.rnom_whl =  500;
   sol.meta.rnom_rol = 1000;
   sol.meta.npatch   = cntc.getnumcontactpatches(ire);
   sol.meta.ipatch   = icp;

   % get discretization parameters

   [ mx, my, xc1, yc1, dx, dy ]  = cntc.getpotcontact(ire, icp);
   sol.mx = mx; sol.my = my;
   sol.xl = xc1 - 0.5*dx; sol.yl = yc1 - 0.5*dy;
   sol.dx = dx; sol.dy = dy;

   sol.x_offset = []; sol.y_offset = [];
   sol.x = sol.xl + ([1:mx]-0.5) * dx;
   sol.y = sol.yl + ([1:my]-0.5) * dy;

   % get grid-data

   sol.eldiv  = cntc.getelementdivision(ire, icp);
   sol.h      = cntc.getfielddata(ire, icp, CNTC.fld_h);
   sol.mu     = cntc.getfielddata(ire, icp, CNTC.fld_mu);
   sol.pn     = cntc.getfielddata(ire, icp, CNTC.fld_pn);
   sol.px     = cntc.getfielddata(ire, icp, CNTC.fld_px);
   sol.py     = cntc.getfielddata(ire, icp, CNTC.fld_py);
   sol.un     = cntc.getfielddata(ire, icp, CNTC.fld_un);
   sol.ux     = cntc.getfielddata(ire, icp, CNTC.fld_ux);
   sol.uy     = cntc.getfielddata(ire, icp, CNTC.fld_uy);
   sx         = cntc.getfielddata(ire, icp, CNTC.fld_sx);
   sy         = cntc.getfielddata(ire, icp, CNTC.fld_sy);
   sol.srel   = sqrt(sx.^2 + sy.^2);
   if (sol.kincns.t_digit>=2)
    sol.shft   = sol.srel * sol.kincns.dq;
   else
    sol.shft   = sol.srel;
   end
   if (use_plast)
    sol.taucrt     = cntc.getfielddata(ire, icp, CNTC.fld_taucrt);
    sol.uplsx      = cntc.getfielddata(ire, icp, CNTC.fld_uplsx);
    sol.uplsy      = cntc.getfielddata(ire, icp, CNTC.fld_uplsy);
    sol.trcbnd     = min(sol.mu .* sol.pn, sol.taucrt);
   else
    sol.trcbnd     = sol.mu .* sol.pn;
   end
   if (sol.h_digit>0)
    sol.temp1      = cntc.getfielddata(ire, icp, CNTC.fld_temp1);
    sol.temp2      = cntc.getfielddata(ire, icp, CNTC.fld_temp2);
   end

  end % cntc.getcpresults

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [vx, vy, phi]=getcreepages(ire, icp, LI)
   % [vx, vy, phi] = cntc.getcreepages(ire, icp)
   %
   % get the kinematic constants (creepages) for a contact problem
   %
   %  vx, vy, phi    - in rolling, T=2,3: long/lat/spin creepages [-, -, angle/length]
   %                   in shifts,  T=1:   long/lat/spin shift [length, length, angle]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 6: m=*, cp     - require icp>0, default 1
   % #getcreepages

   l1={1,'vx','Creepage in x-direction'
    2,'vy','Creepage in y-direction'
    3,'phi','Spin creepage'};

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.getcreepages: not available for icp=%d',icp));
    return
   end

   p_vx  = libpointer('doublePtr',-1);
   p_vy  = libpointer('doublePtr',-1);
   p_phi = libpointer('doublePtr',-1);

   cntc.call('cntc.getcreepages', ire, icp, p_vx, p_vy, p_phi);

   vx  = p_vx.value;
   vy  = p_vy.value;
   phi = p_phi.value;

   if nargin ==3
    values=[vx,vy,phi];
    i2=LI.Cmacro.rowM(l1(:,2));
    LI.Cmacro.Y(i2,ire,LI.cur.j1)=values(vertcat(l1{:,1}));
   end

  end % cntc.getcreepages

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [ un, ux, uy]=getdisplacements(ire, icp)
   % [ un, ux, uy ] = cntc.getdisplacements(ire, icp)
   %
   % return the displ.differences for all elements in the potential contact area for a contact problem
   %  mx, my                            - number of elements in potential contact area
   %  un(my,mx), ux(my,mx), uy(my,mx)   - displacement difference [length]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 6: m=*, cp     - require icp>0, default 1
   % #getdisplacements

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.getdisplacements: not available for icp=%d',icp));
    return
   end

   [ mx, my ] = cntc.getnumelements(ire, icp);

   if (mx*my <= 0)
    un = [0];
    ux = [0];
    uy = [0];
   else
    lenarr = mx * my;
    p_un = libpointer('doublePtr',zeros(lenarr,1));
    p_ux = libpointer('doublePtr',zeros(lenarr,1));
    p_uy = libpointer('doublePtr',zeros(lenarr,1));

    cntc.call('cntc.getdisplacements', ire, icp, lenarr, p_un, p_ux, p_uy);

    un = reshape(p_un.value, mx, my)';
    ux = reshape(p_ux.value, mx, my)';
    uy = reshape(p_uy.value, mx, my)';
   end

  end % cntc.getdisplacements

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [ eldiv]=getelementdivision(ire, icp)
   % [ eldiv ] = cntc.getelementdivision(ire, icp)
   %
   % return flags for all elements in the potential contact area for a contact problem
   % indicating whether the element is in Exerior (0), Adhesion (1), Slip (2) or Plasticity (3).
   %
   %  mx, my       - number of elements in potential contact area
   %  eldiv(my,mx) - element division of contact area, 0=E, 1=H, 2=S, 3=P.
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 6: m=*, cp     - require icp>0, default 1
   % #getelementdivision

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.getelementdivision: not available for icp=%d',icp));
    return
   end

   [ mx, my ] = cntc.getnumelements(ire, icp);

   if (mx*my <= 0)
    eldiv = [ 0 ];
   else
    lenarr = mx * my;
    p_eldiv = libpointer('int32Ptr',zeros(lenarr,1));

    cntc.call('cntc.getelementdivision', ire, icp, lenarr, p_eldiv);

    % convert to double for consistency with loadcase.m
    % note that for int32-arrays, NaN == 0
    eldiv = double(reshape(p_eldiv.value, mx, my)');
   end

  end % cntc.getelementdivision

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [ fld]=getfielddata(ire, icp, ifld)
   % [ fld ] = cntc.getfielddata(ire, icp, ifld)
   %
   % return the values of field 'ifld' for all elements in the potential contact area for a contact problem
   %
   %  mx, my        - number of elements in potential contact area
   %  fld(my,mx)    - output value of field 'ifld' for all elements of contact area
   %
   %                  ifld  =  cntc.fld_h      =    1 ! for retrieving array h      [length]
   %                           cntc.fld_mu     =    2 ! for retrieving array mu     [-]
   %                           cntc.fld_px     =    3 ! for retrieving array px     [force/area]
   %                           cntc.fld_py     =    4 ! for retrieving array py     [force/area]
   %                           cntc.fld_pn     =    5 ! for retrieving array pn     [force/area]
   %                           cntc.fld_ux     =    7 ! for retrieving array ux     [length]
   %                           cntc.fld_uy     =    8 ! for retrieving array uy     [length]
   %                           cntc.fld_un     =    9 ! for retrieving array un     [length]
   %                           cntc.fld_taucrt =   11 ! for retrieving array taucrt [force/area]
   %                           cntc.fld_uplsx  =   12 ! for retrieving array uplsx  [length]
   %                           cntc.fld_uplsy  =   13 ! for retrieving array uplsy  [length]
   %                           cntc.fld_sx     =   15 ! for retrieving array sx     [-]
   %                           cntc.fld_sy     =   16 ! for retrieving array sy     [-]
   %                           cntc.fld_temp1  =   20 ! for retrieving array temp1  [C]
   %                           cntc.fld_temp2  =   21 ! for retrieving array temp2  [C]
   %                           cntc.fld_wx     =   22 ! for retrieving array wx     [-]
   %                           cntc.fld_wy     =   23 ! for retrieving array wy     [-]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 6: m=*, cp     - require icp>0, default 1
   % #getfielddata

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.gettractions: not available for icp=%d',icp));
    return
   end
   if (nargin<3 | isempty(ifld))
    ifld = CNTC.fld_pn;
   end

   [ mx, my ] = cntc.getnumelements(ire, icp);

   if (mx*my <= 0)
    fld = [0];
   else
    lenarr = mx * my;
    p_val = libpointer('doublePtr',zeros(lenarr,1));

    cntc.call('cntc.getfielddata', ire, icp, ifld, lenarr, p_val);

    fld = reshape(p_val.value, mx, my)';
   end

  end % cntc.getfielddata

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [ values]=getflags(ire, icp, params)
   % [ values ] = cntc.getflags(ire, icp, params)
   %
   % used for retrieving various configuring flags from a contact problem
   %
   %  lenflg         - length of params/values arrays
   %  params(lenflg) - codes of the parameters to be obtained from CONTACT
   %  values(lenflg) - values of the parameters obtained from CONTACT
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 5: m=*, wtd    - default icp=-1
   % #getflags

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = -1; % default: W/R contact, all patches
   end
   if (nargin<3 | isempty(params))
    disp('ERROR in cntc.getflags: params is mandatory.');
    return
   end

   lenarr = length(params);
   p_values = libpointer('int32Ptr',zeros(lenarr,1));

   cntc.call('cntc.getflags', ire, icp, lenarr, params, p_values);

   % convert to double for consistency with loadcase.m
   % note that for int32-arrays, NaN == 0
   values = double(p_values.value);

  end % cntc.getflags

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [values]=getglobalforces(ire, icp, LI)
   % [ values ] = cntc.getglobalforces(ire, icp)
   %
   % return the overall forces for a w/r contact problem (module 1 only)
   %
   % rvalues(lenarr)   - contact forces [force] and moments [force.length]
   %
   % For returned values use cntc.getglobalforces
   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 4: m=1, wtd/cp - default icp=-1
   % #getglobalforces

   l1={
    1,'FX_TR','total force on the output body, component in track longitudinal x-direction'
    2,'FY_TR','total force on the output body, component in track lateral y-direction'
    3,'FZ_TR','total force on the output body, component in track vertical z-direction'
    4,'MX_R_TR','total moment on output body about rail profile marker, component in track x-direction'
    5,'MY_R_TR','total moment on output body about rail profile marker, component in track y-direction'
    6,'MZ_R_TR','total moment on output body about rail profile marker, component in track z-direction'

    7,'FX_WS','total force on the output body, component in wheelset longitudinal x-direction'
    8,'FY_WS','total force on the output body, component in wheelset lateral y-direction'
    9,'FZ_WS','total force on the output body, component in wheelset vertical z-direction'
    10,'MX_W_WS','total moment on output body about wheel profile marker, component in wheelset x-direction'
    11,'MY_W_WS','total moment on output body about wheel profile marker, component in wheelset y-direction'
    12,'MZ_W_WS','total moment on output body about wheel profile marker, component in wheelset z-direction'
    %       Note that the 'output body' is the rail when using CONTACTs unit convention'
    13,'FX_R','total force on the output body, component in rail profile x-direction'
    14,'FY_R','total force on the output body, component in rail profile y-direction'
    15,'FZ_R','total force on the output body, component in rail profile z-direction'
    16,'MX_R_R','total moment on output body about rail profile marker, component in rail x-direction'
    17,'MY_R_R','total moment on output body about rail profile marker, component in rail y-direction'
    18,'MZ_R_R','total moment on output body about rail profile marker, component in rail z-direction'

    19,'FX_W','total force on the output body, component in wheel profile x-direction'
    20,'FY_W','total force on the output body, component in wheel profile y-direction'
    21,'FZ_W','total force on the output body, component in wheel profile z-direction'
    22,'MX_W_W','total moment on output body about wheel profile marker, component in wheel x-direction'
    23,'MY_W_W','total moment on output body about wheel profile marker, component in wheel y-direction'
    24,'MZ_W_W','total moment on output body about wheel profile marker, component in wheel z-direction'
    }; l1(:,2)=lower(l1(:,2));

   LI=cntc.call;

   if nargout==0&&nargin==0
    disp(l1);return
   end
   if (nargin<1 | isempty(ire)); ire = 1; end
   if (nargin<2 | isempty(icp))
    icp = -1;   % default: sum over all contact patches
   end

   lenarr = 30;
   p_values = libpointer('doublePtr',zeros(lenarr,1));
   cntc.call('cntc.getglobalforces', ire, icp, lenarr, p_values);

   values = p_values.value;

   if nargin==3
    i2=LI.Cmacro.rowM(l1(:,2));
    LI.Cmacro.Y(i2,ire,LI.cur.j1)=values(vertcat(l1{:,1}));
   end

   if 1==0
    [i1,i2]=ismember(upper(LI.Out.X{3}(:,1)),l1(:,2));
    'xxxeb'
    % LI.Out.Y()
    dbstack;
   end

  end % cntc.getglobalforces

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [ dx, dy]=getgriddiscretization(ire, icp)
   % [ dx, dy ] = cntc.getgriddiscretization(ire, icp)
   %
   % get the grid discretization step sizes dx,dy for a contact problem
   %
   %  dx, dy        - grid discretization step sizes [length]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 6: m=*, cp     - require icp>0, default 1
   % #getgriddiscretization

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.getgriddiscretization: not available for icp=%d',icp));
    return
   end

   p_dx = libpointer('doublePtr',-1);
   p_dy = libpointer('doublePtr',-1);

   cntc.call('cntc.getgriddiscretization', ire, icp, p_dx, p_dy);

   dx = p_dx.value;
   dy = p_dy.value;

  end % cntc.getgriddiscretization

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [ rvalues]=gethertzcontact(ire, icp)
   % [ rvalues ] = cntc.gethertzcontact(ire, icp)
   %
   % get the parameters from a Hertzian contact problem
   %
   % The following values are returned in rvalues:
   %   1 - A1       - curvature in rolling direction [1/length]
   %   2 - B1       - curvature in lateral direction [1/length]
   %   3 - AA       - semi-axis in rolling direction [length]
   %   4 - BB       - semi-axis in lateral direction [length]
   %   5 - RHO      - effective radius of curvature, 2 / (A1 + B1) [length]
   %   6 - CP       - effective semi-axis, sqrt(AA * BB) [length]
   %   7 - SCALE    - potential contact scale factor [-]
   %   8 - BNEG     - semi-axis of negative half-ellipse in lateral direction [length]
   %   9 - BPOS     - semi-axis of positive half-ellipse in lateral direction [length]
   %  10 - AOB      - ellipticity AA/BB [-]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 1: m=3, cp     - require icp>0, default 1
   % #gethertzcontact

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.gethertzcontact: not available for icp=%d',icp));
    return
   end

   lenarr = 10;
   p_rvalues = libpointer('doublePtr',zeros(lenarr,1));

   cntc.call('cntc.gethertzcontact', ire, icp, lenarr, p_rvalues);

   rvalues = p_rvalues.value;

  end % cntc.gethertzcontact

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [ CNTC ] = getmagicnumbers();
   % [ CNTC ] = cntc.getmagicnumbers();
   %
   % get the struct with 'magic numbers' for configuring CONTACT
   %
   %  CNTC         - struct with 'magic numbers' for configuring CONTACT
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 0: m=*, glob   - no icp needed
   % #getmagicnumbers


   % return a struct with 'magic numbers' for setting flags later on

   CNTC.if_units   = 1933; % code for setting the units convention
   % note that w/r profiles are given in [mm] in all cases
   CNTC.un_cntc    = 1934; % code for selecting CONTACT's unit convention,
   %          length in [mm], veloc in [mm/s], force in [N], acting on body 1
   CNTC.un_spck    = 1935; % code for selecting SIMPACK's unit convention,
   %          length in [m], veloc in [m/s], force in [N], acting on body 2
   CNTC.un_si      = 1936; % code for selecting SI unit convention,
   %          length in [m], veloc in [m/s], force in [N], acting on body 1
   CNTC.un_imper   = 1937; % code for selecting imperial units (n.y.a.),
   %          length in [in], veloc in [???], force in [lbf], on body 1

   CNTC.inp_vamp   = 1984; % code for cntc.readinpfile for input to CONTACT add-on to Vampire (n.y.a.)
   CNTC.inp_spck   = 1990; % code for cntc.readinpfile for input to Simpack user-subroutine
   CNTC.inp_gensys = 1992; % code for cntc.readinpfile for input to CONTACT add-on to GENSYS (n.y.a.)

   CNTC.ic_config  = 1967; % code for setting the control digit C1, CONFIG:
   %   0: wheelset on track, left rail
   %   1: wheelset on track, right rail
   %   4: wheelset on rollers, left side
   %   5: wheelset on rollers, right side
   CNTC.ic_pvtime  = 1970; % code for setting the control digit P, PVTIME:
   %   0: continuation, sequence of cases (n.y.a. in module 1)
   %   2: initiation of contact, first case
   CNTC.ic_bound   = 1971; % code for setting the control digit B, BOUND:
   %   0: compute traction bound from normal problem
   %   2, 3: set parabolic/elliptical traction bound (Fastsim)
   %   4: set elliptical traction bound according to the SDEC approach
   %   5: use non-Hertzian traction bound according to the KPEC approach
   %   6: use non-Hertzian traction bound according to modified ANALYN approach
   CNTC.ic_tang    = 1972; % code for setting the control digit T, TANG:
   %   0: frictionless, no tangential problem
   %   1: shift or rolling with (material-) fixed coordinates
   %   2: transient rolling with (contact-) moving coordinates
   %   3: steady rolling with direct method, moving coordinates
   CNTC.ic_norm    = 1973; % code for setting the control digits N1 and N3, NORM:
   % module 1: N1:  0: vertical position Z_WS prescribed
   %                1: vertical force FZ prescribed
   % module 3: N3:  0: approach PEN prescribed
   %                1: normal force FN prescribed
   CNTC.ic_force   = 1974; % code for setting the control digit F, FORCE:
   %   0: creepages CKSI and CETA prescribed (default)
   %   1: tangential force FX and creepage CETA prescribed
   %   2: tangential forces FX and FY prescribed
   CNTC.ic_frclaw  = 1976; % code for retrieving the control digit L, FRCLAW
   CNTC.ic_discns  = 1977; % code for setting the control digit D, DISCNS (module 1):
   %   2: planar contact, using creep calculation
   %   3: planar contact, extended rigid slip calculation
   %   4: conformal contact on curved surface
   CNTC.ic_inflcf  = 1978; % code for setting the control digit C3, INFLCF:
   %   2: analytical IF for halfspace, piecewise constant
   %   3: analytical IF for halfspace, bilinear elements
   %   4: IF for conformal contact, using angle correction
   CNTC.ic_mater   = 1979; % code for retrieving the control digit M, MATER
   CNTC.ic_xflow   = 1981; % code for setting extended debug print output (X, XFLOW)
   % <=0: no additional debug output
   %  >0: codeword PSFLRIN for profil, smooth, force, locate,
   %      readln, inflcf, nmdbg
   CNTC.ic_heat    = 1982; % code for retrieving the control digit H, HEAT (cntc.getflags)
   CNTC.ic_iestim  = 1983; % code for setting the control digit I, IESTIM:
   %   0: start from zero initial estimate (default)
   %   3: use el. division and tractions of previous case (n.y.a. in module 1)
   CNTC.ic_output  = 1984; % code for setting the O-control digit O, OUTPUT:
   %   0-4: amount of output on the results of the problem
   CNTC.ic_flow    = 1985; % code for setting the W-control digit W, FLOW:
   %   0-9: amount of output on the flow of the computations
   CNTC.ic_return  = 1986; % code for setting the R-control digit R, RETURN:
   %   0-1: perform actual calculation
   %   2-3: perform checks but skip actual calculation
   CNTC.ic_matfil  = 1987; % code for setting the A-control digit A, MATFIL:
   %   0: the mat-file is not created
   %   1: write detailed results for contact area
   %   2: write detailed results for potential contact area
   CNTC.ic_sens    = 1988; % code for setting the control digit S2, SENS:
   %   0: sensitivities are not computed/printed
   %   2: compute sensitivities for normal problem
   %   3: compute sensitivities for normal+tangential problems
   CNTC.ic_ifmeth  = 1989; % code for setting IF_METH for the Blanco-approach (C3=4)
   %   0: fast, constant curvature (default),
   %   1: detailed, varying curvature
   CNTC.ic_ifvari  = 1990; % code for setting VARIANT for the Blanco-approach (C3=4)
   %   variants: 1--4, default 4.
   CNTC.ic_sbsout  = 1991; % code for setting the control digit O_s, OUTPUT_SUBS:
   %   extent of output on subsurf.stresses to out-file:
   %   0: no results are printed
   %   1: print maximum values of primary stress invariants
   %   2: more maximum values: Tresca, principal stresses
   %   3: not used
   %   4: print detailed results: displacements, invariants
   CNTC.ic_sbsfil  = 1992; % code for setting the control digit A_s, MATFIL_SUBS:
   %   0: the subs-file is not created
   %   1: displacements and stress-invariants are written
   %   2: the full stress tensor is written as well
   CNTC.ic_npomax  = 1993; % code for setting the max. #elements in the potential contact area

   CNTC.if_idebug  = 2000; % code for level of print-output of the addon itself
   %     0 = none, 1 = overview (default), 2 = full, 3+ = debug
   CNTC.if_licdbg  = 2001; % code for tracing which license is used
   %     0 = none, 1 = overview (default), 2 = full, 3+ = debug
   CNTC.if_wrtinp  = 2002; % code for activating writing of a CONTACT input-file
   %     0 = no, 1 = write .inp-file (default)
   CNTC.if_openmp  = 2003; % >0: number of threads to use per contact patch in CONTACT
   %     -1 = automatic (#cores), 0 = no multi-threading (default)
   % note: per-patch multi-threading is disabled if multi-
   %       threading is used to compute different contact
   %       patches concurrently
   CNTC.if_timers  = 2004; % code for level of output w.r.t. performance timers
   %     0 = none, 1 = overview, 2 = full

   % codes for cntc.getFieldData:

   CNTC.fld_h      =   1; % for retrieving array h      [length]
   CNTC.fld_mu     =   2; % for retrieving array mu     [-]
   CNTC.fld_px     =   3; % for retrieving array px     [force/area]
   CNTC.fld_py     =   4; % for retrieving array py     [force/area]
   CNTC.fld_pn     =   5; % for retrieving array pn     [force/area]
   CNTC.fld_ux     =   7; % for retrieving array ux     [length]
   CNTC.fld_uy     =   8; % for retrieving array uy     [length]
   CNTC.fld_un     =   9; % for retrieving array un     [length]
   CNTC.fld_taucrt =  11; % for retrieving array taucrt [force/area]
   CNTC.fld_uplsx  =  12; % for retrieving array uplsx  [length]
   CNTC.fld_uplsy  =  13; % for retrieving array uplsy  [length]
   CNTC.fld_sx     =  15; % for retrieving array sx     [-]
   CNTC.fld_sy     =  16; % for retrieving array sy     [-]
   CNTC.fld_temp1  =  20; % for retrieving array temp1  [C]
   CNTC.fld_temp2  =  21; % for retrieving array temp2  [C]
   CNTC.fld_wx     =  22; % for retrieving array wx     [-]
   CNTC.fld_wy     =  23; % for retrieving array wy     [-]

   % in module 3, the calling environment may provide additional information to be stored in the mat-file:

   CNTC.mt_tim    = 2011; % code for setting the simulation time
   CNTC.mt_xr     = 2012; % code for setting xr, long. pos. of contact point on rail
   CNTC.mt_yr     = 2013; % code for setting yr, lat. pos. of contact point on rail
   CNTC.mt_xw     = 2014; % code for setting xw, long. pos. of contact point on wheel
   CNTC.mt_yw     = 2015; % code for setting yw, lat. pos. of contact point on wheel
   CNTC.mt_sw     = 2016; % code for setting sw, long. wheel position in track coordinates
   CNTC.mt_run    = 2021; % code for setting irun, run number
   CNTC.mt_axle   = 2022; % code for setting iax, axle number
   CNTC.mt_side   = 2023; % code for setting iside, side number

   CNTC.err_allow  =  -12; % no license found or invalid license
   CNTC.err_search =  -25; % failure in contact search
   CNTC.err_ftot   =  -26; % no solution in total force iteration
   CNTC.err_norm   =  -27; % no convergence in NORM
   CNTC.err_tang   =  -28; % no convergence in TANG
   CNTC.err_tol    =  -29; % iteration stopped with residual > tolerance
   CNTC.err_icp    =  -31; % invalid ire/icp combination
   CNTC.err_profil =  -32; % error reading w/r profile file(s)
   CNTC.err_frclaw =  -33; % invalid friction parameters
   CNTC.err_discr  =  -34; % invalid discretisation parameters
   CNTC.err_input  =  -39; % error reading inp-file
   CNTC.err_other  =  -99; % other error, unspecified
   CNTC.err_broydn =  CNTC.err_ftot; % obsolete, kept for backward compatibility

  end % cntc.getmagicnumbers

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [ pmax]=getmaximumpressure(ire, icp)
   % [ pmax ] = cntc.getmaximumpressure(ire, icp)
   %
   % return the maximum normal pressure in a contact problem
   %
   %  pnmax         - maximum pressure in contact patch [force/area]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 6: m=*, cp     - require icp>0, default 1
   % #getmaximumpressure

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.getmaximumpressure: not available for icp=%d',icp));
    return
   end

   p_pmax = libpointer('doublePtr',-1);

   cntc.call('cntc.getmaximumpressure', ire, icp, p_pmax);

   pmax = p_pmax.value;

  end % cntc.getmaximumpressure

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [ t1max, t2max]=getmaximumtemperature(ire, icp)
   % [ t1max, t2max ] = cntc.getmaximumtemperature(ire, icp)
   %
   % return the maximum contact temperature in a contact problem
   %
   %  t1max, t2max  - maximum surface temperatures in bodies 1 and 2 in contact patch [C]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 6: m=*, cp     - require icp>0, default 1
   % #getmaximumtemperature

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.getmaximumtemperature: not available for icp=%d',icp));
    return
   end

   p_t1max = libpointer('doublePtr',-1);
   p_t2max = libpointer('doublePtr',-1);

   cntc.call('cntc.getmaximumtemperature', ire, icp, p_t1max, p_t2max);

   t1max = p_t1max.value;
   t2max = p_t2max.value;

  end % cntc.getmaximumtemperature

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [ ptmax]=getmaximumtraction(ire, icp)
   % [ ptmax ] = cntc.getmaximumtraction(ire, icp)
   %
   % return the maximum tangential traction in a contact problem
   %
   %  ptmax         - maximum traction |pt| in contact patch [force/area]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 6: m=*, cp     - require icp>0, default 1
   % #getmaximumtraction

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.getmaximumtraction: not available for icp=%d',icp));
    return
   end

   p_ptmax = libpointer('doublePtr',-1);

   cntc.call('cntc.getmaximumtraction', ire, icp, p_ptmax);

   ptmax = p_ptmax.value;

  end % cntc.getmaximumtraction

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [ sx, sy]=getmicroslip(ire, icp)
   % [ sx, sy ] = cntc.getmicroslip(ire, icp)
   %
   % return the relative micro-slip velocity for all elements in the potential contact area for
   %            a contact problem
   %
   %  mx, my               - number of elements in potential contact area
   %  sx(my,mx), sy(my,mx) - in rolling, T=2,3: relative micro-slip velocity [-]
   %                         in shifts,  T=1:   shift distance [length]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 6: m=*, cp     - require icp>0, default 1
   % #getmicroslip

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.getmicroslip: not available for icp=%d',icp));
    return
   end

   [ mx, my ] = cntc.getnumelements(ire, icp);

   if (mx*my <= 0)
    sx = [ 0 ];
    sy = [ 0 ];
   else
    lenarr = mx * my;
    p_sx = libpointer('doublePtr',zeros(lenarr,1));
    p_sy = libpointer('doublePtr',zeros(lenarr,1));

    cntc.call('cntc.getmicroslip', ire, icp, lenarr, p_sx, p_sy);

    sx = reshape(p_sx.value, mx, my)';
    sy = reshape(p_sy.value, mx, my)';
   end

  end % cntc.getmicroslip

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [ npatch]=getnumcontactpatches(ire)
   % [ npatch ] = cntc.getnumcontactpatches(ire)
   %
   % return the number of contact patches used in a w/r contact problem
   %
   %  npatch        - number of separate contact patches
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 2: m=1, wtd    - no icp needed
   % #getnumcontactpatches

   if (nargin<1 | isempty(ire))
    ire = 1;
   end

   p_npatch = libpointer('int32Ptr',-1);

   cntc.call('cntc.getnumcontactpatches', ire, p_npatch);

   npatch = double(p_npatch.value);

  end % cntc.getnumcontactpatches

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [ mx, my]=getnumelements(ire, icp, LI)
   % [ mx, my ] = cntc.getnumelements(ire, icp)
   %
   % return the number of elements in the potential contact area used for a contact problem,
   %            length of tractions arrays
   %
   %  mx, my        - number of discretization elements in long/lat dirs
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 6: m=*, cp     - require icp>0, default 1
   % #getnumelements

   l1={1,'mx','number of discretization elements in long/lat dirs'
    2,'my','number of discretization elements in long/lat dirs'};
   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.getnumelements: not available for icp=%d',icp));
    return
   end

   p_mx = libpointer('int32Ptr',0);
   p_my = libpointer('int32Ptr',0);

   cntc.call('cntc.getnumelements', ire, icp, p_mx, p_my);

   mx = double(p_mx.value);
   my = double(p_my.value);

   if nargin==3
    values=[mx,my];
    i2=LI.Cmacro.rowM(l1(:,2));
    LI.Cmacro.Y(i2,ire,LI.cur.j1)=values(vertcat(l1{:,1}));
   end

  end % cntc.getnumelements

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [ values]=getparameters(ire, icp, LI)
   % [ values ] = cntc.getparameters(ire, icp)
   %
   % used for retrieving various parameters from a contact problem needed for cntc.getcpresults
   %
   %  lenarr         - length of values array
   %  values(lenarr) - values of the parameters obtained from CONTACT
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 5: m=*, wtd    - default icp=-1
   % #getparameters

   l1={
    1,'veloc'
    2,'chi'
    3,'dq'
    4,'spinxo'
    5,'spinyo'
    6,'tau_c0'
    };

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = -1; % default: W/R contact, all patches
   end

   lenarr = 10;
   p_values = libpointer('doublePtr',zeros(lenarr,1));

   cntc.call('cntc.getparameters', ire, icp, lenarr, p_values);

   tmp    = p_values.value;
   values = struct();
   values.veloc  = tmp(1);
   values.chi    = tmp(2);
   values.dq     = tmp(3);
   values.spinxo = tmp(4);
   values.spinyo = tmp(5);
   values.tau_c0 = tmp(6);

   if nargin==3
    values=[values.veloc, values.chi, values.dq, values.spinxo, ...
     values.spinyo, values.tau_c0];
    i2=LI.Cmacro.rowM(l1(:,2));
    LI.Cmacro.Y(i2,ire,LI.cur.j1)=values(vertcat(l1{:,1}));
   end

  end % cntc.getparameters

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [ pen]=getpenetration(ire, icp, LI)
   % [ pen ] = cntc.getpenetration(ire, icp)
   %
   % return the penetration (approach) for a contact problem
   %
   %  pen           - penetration [length]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 6: m=*, cp     - require icp>0, default 1
   % #getpenetration

   l1={1,'pen','penetration for the contact problem'};
   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.getpenetration: not available for icp=%d',icp));
    return
   end

   p_pen = libpointer('doublePtr',-1);
   cntc.call('cntc.getpenetration', ire, icp, p_pen);
   pen=p_pen.value;

   if nargin==3
    i2=LI.Cmacro.rowM(l1(:,2));
    LI.Cmacro.Y(i2,ire,LI.cur.j1)=pen;
   end

  end % cntc.getpenetration

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [ mx, my, xc1, yc1, dx, dy ]=getpotcontact(ire, icp, LI)
   % [ mx, my, xc1, yc1, dx, dy ] = cntc.getpotcontact(ire, icp)
   %
   % get the parameters of the potential contact area for a contact problem
   %    3: first center + grid sizes,        params = [ mx, my, xc1, yc1, dx, dy ]
   % values= { mx, my, xc1, yc1, dx, dy }
   %  mx, my        - number of elements in x- and y-directions [-]
   %  xc1, yc1      - position of first element center [length]
   %  dx, dy        - grid discretization step sizes [length]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 6: m=*, cp     - require icp>0, default 1
   % #getpotcontact

   l1={
    1, 'mx', 'number of elements in x-directions [-]'
    2, 'my', 'number of elements in y-directions [-]'
    3, 'xc1', 'position of first element center [length]'
    4, 'yc1', ' position of first element center [length]'
    5, 'dx', 'grid discretization step sizes [length]'
    6, 'dy', 'grid discretization step sizes [length]'
    };


   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc_getpotcontact: not available for icp=%d',icp));
    return
   end

   lenarr = 6;
   p_values = libpointer('doublePtr',zeros(lenarr,1));

   cntc.call('cntc.getpotcontact', ire, icp, lenarr, p_values);

   values = p_values.value;

   mx  = round(values(1));
   my  = round(values(2));
   xc1 = values(3);
   yc1 = values(4);
   dx  = values(5);
   dy  = values(6);

   if nargin==3
    i2=LI.Cmacro.rowM(l1(:,2));
    LI.Cmacro.Y(i2,ire,LI.cur.j1)=values(vertcat(l1{:,1}));
   end


  end % cntc_getpotcontact

  %------------------------------------------------------------------------------------------------------------

  %------------------------------------------------------------------------------------------------------------
  function [ values]=getprofilevalues(ire, itask, iparam, rparam)
   % [ values ] = cntc.getprofilevalues(ire, itask, iparam, rparam)
   %
   % Get a wheel or rail profile for a wheel-rail contact problem as a table of values
   %
   %  itask          - select type of outputs:
   %                    -1: ierr   <0: error codes, 0=profile loaded ok, 1=not set
   %                     0: npnt   number of points used in requested sampling method
   %                     1: r/w    get (yr,zr) values for rail or (yw,zw) for wheel profile
   %                     2: trk    get (ytr,ztr) values for rail or wheel profile (principal profile)
   %                     3: gaug   get left-most point within gauge height, offset, point at gauge height
   %                               [ygauge1, zgauge1, yoffs, zoffs, ygauge2, zgauge2]  [length]
   %                     4: arc    get arc-length parameter s along profile
   %                     5: angl   get surface inclination atan2(dz, dy) [angle]
   %  iparam         - integer configuration parameters
   %                     1: itype     0 = rail, 1 = wheel profile
   %                     2: isampl   -1 = sampling cf. original input data;
   %                                  0 = sampling cf. spline representation (default);
   %                                  1 = sampling cf. spline representation at spacing ds_out
   %                            kchk>=2 = sampling cf. spline representation with integer refinement factor
   %  rparam         - real configuration parameters
   %                     1: ds_out  step-size ds used with sampling method isampl=1, default 1mm
   %                     2: x_out   for variable profiles: longitudinal position x_r [mm] (task=1) or
   %                                x_tr [mm] (task=2) or circumferential position th_w [rad] on wheel
   %
   %  tasks 1,2,4: no unit conversion or scaling are applied for profile values
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 2: m=1, wtd    - no icp needed
   % #getprofilevalues

   cntc_err_profil = -32;

   values = [];
   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(itask))
    disp('ERROR in cntc.getprofilevalues: itask is mandatory.');
    return
   end
   if (itask<-1 | itask>5)
    disp(sprintf('ERROR in cntc.getprofilevalues: itask=%d not available.', itask));
    return
   end
   if (nargin<3 | isempty(iparam))
    itype = 0; isampl = 0;
    iparam = [ itype, isampl ];
   end
   if (nargin<4 | isempty(rparam))
    ds_out = 1.0;
    rparam = [ ds_out ];
   end
   nints  = length(iparam);
   nreals = length(rparam);

   if (itask==-1)       % task -1: get error code

    lenarr = 2;
    p_val  = libpointer('doublePtr',zeros(lenarr,1));

    cntc.call('cntc.getprofilevalues_new', ire, itask, nints, iparam, nreals, rparam, 1, p_val);
    values = round(p_val.value);

   elseif (itask>=0 & itask<=5)

    % task 0--5: first get the number of points in the profile at requested s sampling positions

    p_npnt = libpointer('doublePtr',-1);
    cntc.call('cntc.getprofilevalues_new', ire, 0, nints, iparam, nreals, rparam, 1, p_npnt);

    npnt = round(p_npnt.value);

    if (npnt==cntc_err_profil | npnt<=0)
     disp(sprintf(['cntc.getprofilevalues: the rail and/or wheel profile could not be found or ',...
      'processed (%d).'], npnt));
     return;
    end

    % create output-array dependent on task to be performed

    if (itask==0)

     values = [ npnt ];

    elseif (itask>=1 & itask<=5)

     lenarr = [ npnt*2, npnt*2, 6, npnt, npnt ];
     lenarr = lenarr(itask);
     p_val  = libpointer('doublePtr',zeros(lenarr,1));

     cntc.call('cntc.getprofilevalues_new', ire, itask, nints, iparam, nreals, rparam, ...
      lenarr, p_val);
     values = p_val.value;

     if (itask==1 | itask==2)
      values = reshape(values, npnt, 2);
     end

    end % (itask==0 | 1--5)

   else

    %values = [];

   end % if (itask==-1)

  end % cntc.getprofilevalues

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [ veloc]=getreferencevelocity(ire, icp, LI)
   % [ veloc ] = cntc.getreferencevelocity(ire, icp)
   %
   % get the rolling velocity for a contact problem
   %
   %  veloc          - absolute rolling velocity [veloc]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 6: m=*, cp     - require icp>0, default 1
   % #getreferencevelocity

   l1={1,'Vel','Absolute rolling velocity'};

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.getreferencevelocity: not available for icp=%d',icp));
    return
   end

   p_veloc = libpointer('doublePtr',-1);

   cntc.call('cntc.getreferencevelocity', ire, icp, p_veloc);
   veloc=p_veloc.value;

   if nargin==3 % 'LI'
    i2=LI.Cmacro.rowM(l1(:,2));
    LI.Cmacro.Y(i2,ire,LI.cur.j1)=veloc;
   end

  end % cntc.getreferencevelocity

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [ sens]=getsensitivities(ire, icp, lenout, lenin)
   % [ sens ] = cntc.getsensitivities(ire, icp, lenout, lenin)
   %
   % return the sensitivities of the total forces for a contact problem
   %
   %  lenout, lenin      - requested number of outputs (forces) and inputs (creepages or shifts)
   %  sens(lenout,lenin) - matrix of sensitivities -- zeros if not computed
   %
   %  calculation of sensitivities is requested using the flag 'cntc.ic_sens'.
   %  accuracy is configured using cntc.setsolverflags with G=6.
   %
   %  the inputs are ordered   1: pen, 2: cksi, 3: ceta, 4: cphi
   %  in rolling, T=2,3, the units are pen [length], cksi, ceta [-],      cphi [angle/length]
   %  in shifts,  T=1,   the units are pen [length], cksi, ceta [length], cphi [angle]
   %
   %  the outputs are ordered  1: fn,  2: fx,   3: fy,   4: mz
   %  the units are fn [force], fx, fy [-], mz [force.length]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 6: m=*, cp     - require icp>0, default 1
   % #getsensitivities



   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.getsensitivities: not available for icp=%d',icp));
    return
   end
   if (nargin<3 | isempty(lenout))
    lenout = 4;
   end
   if (nargin<4 | isempty(lenin))
    lenin  = 4;
   end

   lenarr = lenout * lenin;
   p_sens = libpointer('doublePtr',zeros(lenarr,1));

   cntc.call('cntc.getsensitivities', ire, icp, lenout, lenin, p_sens);

   sens = reshape(p_sens.value, lenout, lenin);

  end % cntc.getsensitivities

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [ pn, px, py]=gettractions(ire, icp)
   % [ pn, px, py ] = cntc.gettractions(ire, icp)
   %
   % return the tractions for all elements in the potential contact area for a contact problem
   %            note the order of the arguments, with pn (z-direction) occurring before px,py.
   %
   %  mx, my                            - number of elements in potential contact area
   %  pn(my,mx), px(my,mx), py(my,mx)   - surface tractions for all elements of contact area, [force/area]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 6: m=*, cp     - require icp>0, default 1
   % #gettractions



   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.gettractions: not available for icp=%d',icp));
    return
   end

   [ mx, my ] = cntc.getnumelements(ire, icp);

   if (mx*my <= 0)
    pn = [ 0 ];
    px = [ 0 ];
    py = [ 0 ];
   else
    lenarr = mx * my;
    p_pn = libpointer('doublePtr',zeros(lenarr,1));
    p_px = libpointer('doublePtr',zeros(lenarr,1));
    p_py = libpointer('doublePtr',zeros(lenarr,1));

    cntc.call('cntc.gettractions', ire, icp, lenarr, p_pn, p_px, p_py);

    pn = reshape(p_pn.value, mx, my)';
    px = reshape(p_px.value, mx, my)';
    py = reshape(p_py.value, mx, my)';
   end

  end % cntc.gettractions

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [ rvalues]=getwheelsetposition(ire,LI)
   % [ rvalues ] = cntc.getwheelsetposition(ire)
   %
   % return the wheelset position a w/r contact problem (module 1 only)
   %
   %  rvalues(lenarr)   - wheelset position parameters.
   %
   %  The following values are returned, if permitted by the length of rvalues:
   %   1 - S_WS      - s-position of the wheelset center of mass along the track center line
   %   2 - Y_WS      - lateral y-position of the wheelset center of mass in track coordinates
   %   3 - Z_WS      - vertical z-position of the wheelset center of mass in track coordinates
   %   4 - ROLL_WS   - wheelset roll angle with respect to the track plane
   %   5 - YAW_WS    - wheelset yaw angle with respect to the track center line x_tr
   %   6 - PITCH_WS  - wheelset pitch angle, i.e. rotation about the wheelset axle
   %
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 2: m=1, wtd    - no icp needed
   % #getwheelsetposition

   l1={
    1, 'S_WS', 's-position of the wheelset center of mass along the track center line'
    2, 'Y_WS', 'lateral y-position of the wheelset center of mass in track coordinates'
    3, 'Z_WS', ' vertical z-position of the wheelset center of mass in track coordinates'
    4, 'ROLL_WS', 'wheelset roll angle with respect to the track plane'
    5, 'YAW_WS', 'wheelset yaw angle with respect to the track center line x_tr'
    6, 'PITCH_WS', 'wheelset pitch angle, i.e. rotation about the wheelset axle'
    };l1(:,2)=lower(l1(:,2));

   if (nargin<1 | isempty(ire))
    ire = 1;
   end

   lenarr = 6;
   p_values = libpointer('doublePtr',zeros(lenarr,1));

   cntc.call('cntc.getwheelsetposition', ire, lenarr, p_values);

   rvalues = p_values.value;
   if nargin==2
    i2=LI.Cmacro.rowM(l1(:,2));
    LI.Cmacro.Y(i2,ire,LI.cur.j1)=rvalues(vertcat(l1{:,1}));
   end

  end % cntc.getwheelsetposition

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [ rvalues]=getwheelsetvelocity(ire)
   % [ rvalues ] = cntc.getwheelsetvelocity(ire)
   %
   % return the wheelset velocity a w/r contact problem (module 1 only)
   %
   %  rvalues(lenarr)   - wheelset velocity parameters.
   %
   %  The following values are returned, if permitted by the length of rvalues:
   %   1 - VS_WS      - s-velocity of the wheelset center of mass along the track center line
   %   2 - VY_WS      - velocity in lateral y-direction of the wheelset center of mass (track coordinates)
   %   3 - VZ_WS      - velocity in vertical z-direction of the wheelset center of mass (track coordinates)
   %   4 - VROLL_WS   - wheelset roll velocity, angular velocity with respect to the track plane
   %   5 - VYAW_WS    - wheelset yaw velocity, angular velocity with respect to the track center line x_tr
   %   6 - VPITCH_WS  - wheelset pitch velocity, i.e. angular velocity about the wheelset axle
   %
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 2: m=1, wtd    - no icp needed
   % #getwheelsetvelocity



   if (nargin<1 | isempty(ire))
    ire = 1;
   end

   lenarr = 6;
   p_values = libpointer('doublePtr',zeros(lenarr,1));

   cntc.call('cntc.getwheelsetvelocity', ire, lenarr, p_values);

   rvalues = p_values.value;

  end % cntc.getwheelsetvelocity

  %------------------------------------------------------------------------------------------------------------
  function []=set(list)
   %% #SetSDT : SDT distribution to CNTC set commands -1
   if ~iscell(list);list={list};end
   LI=cntc.call;
   for j1=1:size(list,1)
    if ischar(list{j1})
     st1=lower(regexprep(list{j1},'^([^{]*).*','$1'));evt=[];
    else; evt=list{j1}; st1=lower(evt.type);
    end
    switch st1
     case 'solver'
      % #set.solver configuration -2
      cntc.setsolverflags(list{j1});
     case 'calculate'
      %% #set.calculate the contact problem for current wheel -2

      disp(sprintf('Starting case %2d for wheel %d', icase, iwhe));
      ierror = cntc.calculate(iwhe);
      if (ierror~=0), return; end

     case 'stresscalc'
      %% #set.stresscalc compute subsurface stresses -2
      error('Implement')

      ierror = cntc.subs_calculate(iwhe);
      if (ierror~=0), return; end

     case 'initout'
      %% #set.initOut Initialize the output -2
      % Fields stored
      if isempty(evt);iwhe=evalin('caller','iwhe');icase=evalin('caller','j1');
      else;iwhe=evt.iwhe;icase=evt.j1;  end

      l1={'eldiv';'h';'mu';'pn';'px';'py';'un';'ux';'uy';'sx';'sy'};
      % Fields SDT Storage
      Cfield=struct('X',{{[],[],l1,[1,1,1;2,2,1], LI.Traj.X{1}}},'Y',[], ...
       'Xlab',{{'Ngx','Ngy','Comp',{'ire';'iwhe';'icp'},'iTime'}});
      LI.Cfield=Cfield;

      % get detailed results per contact patch, if there is more than one
      % contact between the wheel and the rail
      LI.Cmacro.rowM=sdtu.ivec('ColList', {}); %dDQ
      LI.Cmacro.Xlab={'comp',{'ire';'iwhe';'icp'},'iTime'};
      LI.Cmacro.X={[],[],LI.Traj.X{1}};

     case 'getbasis'
      ire=1;
      values = cntc.getwheelsetposition(ire);

     case 'getout'
      %% #set.GetOut Get results -2
      if isempty(evt);iwhe=evalin('caller','iwhe');icase=evalin('caller','j1');
      else;iwhe=evt.iwhe;icase=evt.j1;  end

      % Initialize the loop
      if icase==1
       cntc.set('initout{}'); 
      end

      if icase==LI.Traj.X{1}(end)
       LI.Cmacro.X{1}=LI.Cmacro.rowM.prop.name;
      end

      if ~isempty(LI.Cmacro.X{2})&&any(LI.Cmacro.X{2}(:,2)==iwhe)
      else      % Need to extend contact dimension
       i2=(1:cntc.getnumcontactpatches(iwhe))'*[0 0 1];
       i2(:,2)=iwhe; i2(:,1)=size(LI.Cmacro.X{2},1)+(1:size(i2,1));
       LI.Cmacro.X{2}=[LI.Cmacro.X{2};i2];
      end
      ind=LI.Cmacro.X{2};ind=ind(ind(:,2)==iwhe,:);

      for jre=1:size(ind,1); % icp = 1 : LI.results{iwhe}.npatch(icase)
       ire=ind(jre,1); iwhe=ind(jre,2);icp=ind(jre,3);

       %% get macro parameters
       cntc.getglobalforces(ire,icp,LI);
       cntc.getcontactlocation(ire,icp, LI);
       cntc.getreferencevelocity(ire,icp, LI);
       cntc.getpenetration(ire, icp, LI);
       cntc.getcreepages(ire, icp, LI);
       cntc.getcontactforces(ire, icp,LI);
       cntc.getwheelsetposition(ire,LI);
       cntc.getparameters(ire, icp, LI);% get material / kinematic parameters needed
       cntc.getpotcontact(ire, icp,LI);
       cntc.getnumelements(ire,icp,LI);


       % itask   = 1; iswheel = 0; isampl=0; iparam=[iswheel, isampl]; rparam=[];
       % RailProfil = cntc.getprofilevalues(ire, itask, iparam, rparam);
       % itask   = 4;
       % RProfileS = cntc.getprofilevalues(ire, itask, iparam, rparam);
       %
       % itask   = 1;iswheel = 1; iparam  = [iswheel, isampl];
       % WheelProfil = cntc.getprofilevalues(ire, itask, iparam, rparam);
       % itask   = 4;
       % WProfileS = cntc.getprofilevalues(ire, itask, iparam, rparam);% a regarder


       %% get micro parameters
       C1=LI.Cfield;
       C1.X{1}=[C1.X{1} LI.Cmacro.Y(72,1,icase)];
       C1.X{2}=[C1.X{2} LI.Cmacro.Y(73,1,icase)];
       % Cfield.Y=zeros(cellfun(@(x)size(x,1),Cfield.X));

       % get grid-data
       r1=cntc.getelementdivision(ire, icp);
       C1.Y(1:size(r1,1),1:size(r1,2),1,ire,icase)=r1;
       r1=cntc.getfielddata(ire, icp, LI.CNTC.fld_h);
       C1.Y(1:size(r1,1),1:size(r1,2),2,ire,icase)=r1;
       r1=cntc.getfielddata(ire, icp, LI.CNTC.fld_mu);
       C1.Y(1:size(r1,1),1:size(r1,2),3,ire,icase)=r1;
       r1=cntc.getfielddata(ire, icp, LI.CNTC.fld_pn);
       C1.Y(1:size(r1,1),1:size(r1,2),4,ire,icase)=r1;
       r1=cntc.getfielddata(ire, icp, LI.CNTC.fld_px);
       C1.Y(1:size(r1,1),1:size(r1,2),5,ire,icase)=r1;
       r1=cntc.getfielddata(ire, icp, LI.CNTC.fld_py);
       C1.Y(1:size(r1,1),1:size(r1,2),6,ire,icase)=r1;
       r1=cntc.getfielddata(ire, icp, LI.CNTC.fld_un);
       C1.Y(1:size(r1,1),1:size(r1,2),7,ire,icase)=r1;
       r1=cntc.getfielddata(ire, icp, LI.CNTC.fld_ux);
       C1.Y(1:size(r1,1),1:size(r1,2),8,ire,icase)=r1;
       r1=cntc.getfielddata(ire, icp, LI.CNTC.fld_uy);
       C1.Y(1:size(r1,1),1:size(r1,2),9,ire,icase)=r1;
       r1=cntc.getfielddata(ire, icp, LI.CNTC.fld_sx);
       C1.Y(1:size(r1,1),1:size(r1,2),10,ire,icase)=r1;
       r1=cntc.getfielddata(ire, icp, LI.CNTC.fld_sy);
       C1.Y(1:size(r1,1),1:size(r1,2),11,ire,icase)=r1;
       LI.Cfield=C1;

       %Need to see the goal of what follow
       if 1==0
        use_plast = (sol.mater.m_digit==4 & (sol.mater.tau_c0>1e-10 & sol.mater.tau_c0<1e10));
        sol.srel   = sqrt(sx.^2 + sy.^2);
        if (sol.kincns.t_digit>=2)
         sol.shft   = sol.srel * sol.kincns.dq;
        else
         sol.shft   = sol.srel;
        end
        if (use_plast)
         sol.taucrt     = cntc.getfielddata(ire, icp, LI.CNTC.fld_taucrt);
         sol.uplsx      = cntc.getfielddata(ire, icp, LI.CNTC.fld_uplsx);
         sol.uplsy      = cntc.getfielddata(ire, icp, LI.CNTC.fld_uplsy);
         sol.trcbnd     = min(sol.mu .* sol.pn, sol.taucrt);
        else
         sol.trcbnd     = sol.mu .* sol.pn;
        end
        if (sol.h_digit>0)
         sol.temp1      = cntc.getfielddata(ire, icp, LI.CNTC.fld_temp1);
         sol.temp2      = cntc.getfielddata(ire, icp, LI.CNTC.fld_temp2);
        end

        % get maximum von mises stress
        iblk = 1;
        save('LI.mat','LI');
        table = subs_getresults(ire, icp, iblk, [1,2,3,8]);
        [vm_max, ii_max] = max(table(:,4));

        results{ire}.cp_force.sigvm(icase,icp) = vm_max;
        results{ire}.cp_force.vm_x(icase,icp)  = table(ii_max,1);
        results{ire}.cp_force.vm_y(icase,icp)  = table(ii_max,2);
        results{ire}.cp_force.vm_z(icase,icp)  = table(ii_max,3);
       end

      end
     case 'getsol'
      %% Get sol to understand plot3d (need to disappear for store in SDT directly)
      if isempty(evt);iwhe=evalin('caller','iwhe');icase=evalin('caller','j1');
      else;iwhe=evt.iwhe;icase=evt.j1;  end

      sol=cntc.getcpresults(iwhe, 1);
      LI.sol{icase}=sol;


     case 'mat'
      %% Set material property
      cntc.setmaterialparameters(list{j1}); %config material
     case 'friction'
      %% friction definition : Coulomb's law
      cntc.setfrictionmethod(list{j1});
     case 'bound'
      %% Boundary condition definition
      cntc.setboundcond(list{j1});
     case 'potcntc'
      %% Potential contact area definition and discretization
      cntc.setpotcontact(list{j1});
     case 'rolling'
      %% Rolling stepsize
      cntc.setrollingstepsize(list{j1});
     case 'track'
      %% Set Track dimensions
      cntc.settrackdimensions(list{j1});
     case 'setprofile'
      %% Set wheel/rail profile
      cntc.setprofileinputfname(list{j1});
     case 'wheelsetdim'
      %% Set wheelset dimension and deviation
      cntc.setwheelsetdimensions(list{j1});
     case 'traj'
      %% #set.traj set wheel position and velocity associated to a trajectory point -2
      % increment icase = j1 in fe_time
      [~,RO]=sdtm.urnPar(['{s_ws,y_ws,z_ws,roll_ws,yaw_ws,pitch_ws,vs,vy,vz,vroll,' ...
       'vyaw,vpitch,dxwhl, dywhl, dzwhl,drollw, dyaww, dpitchw, vxwhl,vywhl, vzwhl,' ...
       ' vrollw, vyaww, vpitchw}'],'{}{}');

      if (~isfield(evt,'j1')&&~isfield(evt,'iwhe'));iwhe=1;icase=1;
      else;iwhe=evt.iwhe;icase=evt.j1;  end

      Ewheel=LI.wheelsetDim.Ewheel;

      % xxxgae unl vnl
      [i1,i2]=ismember(LI.Traj.X{2},RO.Other(1:6));
      pos=zeros(1,6);
      pos(i2(i1))=LI.Traj.Y(icase,i1);
      cntc.setwheelsetposition(iwhe,Ewheel,pos); % set wheel set pos

      [i1,i2]=ismember(LI.Traj.X{2},RO.Other(7:12));
      vel=zeros(1,6);vel(i2(i1))=LI.Traj.Y(icase,i1);
      cntc.setwheelsetvelocity(iwhe,Ewheel,vel);
      LI.cur=evt;

      [i1,i2]=ismember(LI.Traj.X{2},RO.Other(12:end));
      flex=zeros(1,6);flex(i2(i1))=LI.Traj.Y(icase,i1);
      cntc.setwheelsetflexibility(iwhe,Ewheel,flex);

    end
   end
  end

    function []= TempLoop(list,LI)
   %% #TempLoop -2
   iwhe=evalin('caller','iwhe');
   for j1 = 1:size(LI.Traj.Y,1) % Loop on traj, In fe_time step is called j1 here icase
    st1=list;
    if comstr(lower(list{1}),'traj');
     st1{1}=struct('type','Traj','j1',j1,'iwhe',iwhe);
    else
     error('Trajectory must be define first')
    end
    cntc.set(st1);
   end % icase j1

  end
  %% #SDTUI user interface  ----------------------------------------------- -2
  function []=tab()
   % #tab -3
   r1={'Name','level';
       'Library call',1;'CNTC',2 ; 'ProjectWd',2 ; 'CallLog',2 ; 'libname',2;...
    'Experiment Definition',1;'flags',2;'global',2;'wheels',2;'Material parameter', ...
    1;'Bound' ,2; 'Friction',2 ;'Mat',2 ;'Solver Definition',1;'PotCntc',2 ; ...
    'Rolling',2 ; 'Solver',2;'Geometry Definition',1;'Track',2 ; 'prr' ,2; ...
    'prw',2 ; 'wheelsetDim' ,2; 'RailProfile',2 ; 'WheelProfile',2;'Loop parameter'...
    ,1;'cur' ,2; 'Traj',2;'Output',1;'Cfield',2 ; 'Cmacro',2 ; 'sol',2};
   LI=cntc.call;
   r2=vhandle.tab(r1,LI);r2.name='CNTC';
   asTab(r2)
   %ua=struct('ColumnName',{{'Name','value','ToolTip'}},)

  
  end

  %------------------------------------------------------------------------------------------------------------

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% #Set : CONTACT functions Set -1

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %------------------------------------------------------------------------------------------------------------
  function []=setcreepages(ire, icp, vx, vy, phi)
   % [ ] = cntc.setcreepages(ire, icp, vx, vy, phi)
   %
   % set the kinematic constants (creepages) for a contact problem
   % note: vx is ignored when the F-digit is 1 or 2, vy is ignored when F=1.
   %
   %  vx, vy, phi    - in rolling, T=2,3: long/lat/spin creepages [-, -, angle/length]
   %                   in shifts,  T=1:   long/lat/spin shift     [length, length, angle]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 1: m=3, cp     - require icp>0, default 1
   % #setcreepages



   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.setcreepages: not available for icp=%d',icp));
    return
   end
   if (nargin<5 | isempty(vx) | isempty(vy) | isempty(phi))
    disp('ERROR in cntc.setcreepages: creepages vx,vy,phi are mandatory.');
    return
   end

   cntc.call('cntc.setcreepages', ire, icp, vx, vy, phi);

  end % cntc.setcreepages

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function []=setextrarigidslip(ire, icp, wx, wy)
   % [ ] = cntc.setextrarigidslip(ire, icp, wx, wy)
   %
   % set the extra term of the rigid slip for all elements in the potential contact area for a contact problem
   %
   %  mx, my               - number of elements in potential contact area
   %  wx(my,mx), wy(my,mx) - in rolling, T=2,3: extra relative rigid slip  [-]
   %                         in shifts,  T=1:   extra rigid shift distance [length]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 1: m=3, cp     - require icp>0, default 1
   % #setextrarigidslip



   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.setextrarigidslip: not available for icp=%d',icp));
    return
   end
   if (nargin<4 | isempty(wx) | isempty(wy))
    disp('ERROR in cntc.setextrarigidslip: wx, wy are mandatory.');
    return
   end

   [ mx, my ] = cntc.getnumelements(ire, icp);

   if (any(size(wx)~=[my mx]) | any(size(wy)~=[my mx]))
    disp('ERROR in cntc.setextrarigidslip: arrays wx, wy must have size (my,mx).');
    return
   end

   wx = reshape(wx', lenarr, 1);
   wy = reshape(wy', lenarr, 1);

   cntc.call('cntc.setextrarigidslip', ire, icp, lenarr, wx, wy);

  end % cntc.setextrarigidslip

  %------------------------------------------------------------------------------------------------------------

  %------------------------------------------------------------------------------------------------------------
  function []=setflags(ire, icp, params, values)
   % [ ] = cntc.setflags(ire, icp, params, values)
   %
   %
   %  lenflg         - length of params/values arrays
   %  params(lenflg) - codes of the parameters to be communicated to CONTACT
   %  values(lenflg) - values of the parameters to be communicated to CONTACT
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 5: m=*, wtd    - default icp=-1
   % #setflags(pdfsection.2.3)  configure flags for a contact problem

   if nargin==0
    LI=cntc.call; CNTC=LI.CNTC;
    flagM=vhandle.nmap([fieldnames(CNTC) struct2cell(CNTC)]);
    li=cell(LI.flags);
    flags=flagM(li(:,1));flags=vertcat(flags{:});% skip iwhe given first always
    % xxx check if scalar value for each row
    values = vertcat(li{:,2});
    cntc.setflags(LI.global.ire, [], flags, values); %config
    return
   end

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = -1; % default: W/R contact, all patches
   end
   if (nargin<4 | isempty(params) | isempty(values))
    disp('ERROR in cntc.setflags: params and values are mandatory.');
    return
   end
   if (length(params)~=length(values))
    disp('ERROR in cntc.setflags: inconsistent length of params and values.');
    return
   end

   cntc.call('cntc.setflags', ire, icp, length(params), params, values);

  end % cntc.setflags

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function []=setfrictionmethod(ire, icp, imeth, params)
   % [ ] = cntc.setfrictionmethod(ire, icp, imeth, params)
   %
   %  imeth          - type of friction law used: 10 * V-digit + 1 * L-digit
   %  params         - depending on method that is used
   %
   %  L = 0: Coulomb friction,             lparam = [fstat, fkin] (pdf.4.2.1)
   %      1: keep previous parameters,     lparam = [ ]
   %      2: linear falling friction,      lparam = [fkin, flin1, sabsh1, flin2, sabsh2, memdst, mem_s0] (pdf.4.2.2)
   %      3: rational falling friction,    lparam = [fkin, frat1, sabsh1, frat2, sabsh2, memdst, mem_s0] (pdf.4.2.2)
   %      4: exponential falling friction, lparam = [fkin, fexp1, sabsh1, fexp2, sabsh2, memdst, mem_s0] (pdf.4.2.2)
   %      5: exponential falling friction, lparam = [fstat, polach_a, polach_b] (pdf.4.2.2)
   %      6: temperature dep. friction,    lparam = [fref, tref, dfheat, dtheat, memdst, mem_s0](pdf.4.2.3/4.2.4)
   %
   %  When V = 0, params = lparam(1:n), with n=2 for L=0, n=7 for L=2--4, n=3 for L=5, n=6 for L=6
   %  When V = 1, params = [ nvf, ...
   %                         alphvf(1),   lparam( 1 ,1:n) ],  ...
   %                                   ...
   %                         alphvf(nvf), lparam(nvf,1:n) ]
   %
   %  dimensions:  alphvf: [angle],  fstat, fkin, flin1,2, frat1,2, fexp1,2, polach_a, fref, dfheat: [-]
   %               sabsh1,2, mem_s0: [veloc],  memdst: [length],  polach_b [1/veloc],  tref, dtheat [C]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 5: m=*, wtd    - default icp=-1
   % #setfrictionmethod parameters for the friction law for a contact problem -2

   % parameters
   if nargin==0||ischar(ire)
    % package SDT options
    DoOpt=['FrcLaw(0#vd{0,coul,1,Maintain,2,LinFall,3,RatFall,4,ExpFall1,5,ExpFall2,6,TempDep}#"Friction law options")' ...
     'fstat(#%g#"static friction coefficient")' ...
     'fkin(#%g#"kinetic friction coefficient")' ...
     'flin1(#%g#"Coefficient for the linear curve")' ...
     'flin2(#%g#"Coefficient for the linear curve")' ...
     'memdst(#%g#"Characteristic distance for the friction memory effect (pdf4.2.4)")' ...
     'mem_s0(#%g#"Minimum velocity for the friction memory effect(pdf4.2.4)")' ...
     'frat1(#%g#"Rational falling friction coefficient")' ...
     'frat2(#%g#"rational falling friction coefficient")' ...
     'sabsh1(#%g#"Absolute slip velocity for which the size of a term is halved")' ...
     'sabsh2(#%g#"Absolute slip velocity for which the size of a term is halved")' ...
     'fexp1(#%g#"Coefficient for the exponential curve")' ...
     'fexp2(#%g#"Coefficient for the exponential curve")' ...
     'polach_a(#%g#"Coefficient de Polach A")' ...
     'polach_b(#%g#"Coefficient de Polach B")' ...
     'fref(#%g#"Reference value for the coefficient of friction at low surface temperatures.")' ...
     'tref(#%g#"Reference temperature Tre f of equation")' ...
     'dfheat(#%g#"Delta du coefficient de friction dependant de la temperature")' ...
     'dtheat(#%g#"Delta de la temperature")' ...
     'nvf(#%g#"Number of control points for friction variation")' ...
     'alphvf(#%g#"Rail surface inclinations at the control points")'];

    if nargin==0; CAM='';else; CAM=ire;end
    LI=cntc.call;
    LI.Friction=cingui('paramedit -doclean2',DoOpt,{struct,CAM});

    switch LI.Friction.FrcLaw
     case 0; params={'fstat','fkin'}; % 0 Coul
     case 1; params={};% 1 Maintain
     case 2; params={'fkin', 'flin1', 'sabsh1', 'flin2', 'sabsh2', 'memdst', 'mem_s0'};% 2 LinFall
     case 3; params={'nu1','nu2','g1','g2','k0_mf','alfamf','betamf'};% 3 RatFall
     case 4; params={'nu1','nu2','g1','g2','g3','laythk','tau_c0','k_tau'};% 4 ExpFall1
     case 5; params={'nu1','nu2','g1','g2','k0_mf','alfamf','betamf'};% 5 ExpFall2
     case 6 ; params={'nu1','nu2','g1','g2','g3','laythk','tau_c0','k_tau'};% 6 TempDep
      %XXXGae Var
     case 10 ; params={'fstat','fkin'}; % 10 VarCoul
     case 11; params={};% 11 VarMaintain
     case 12; params={'fkin', 'flin1', 'sabsh1', 'flin2', 'sabsh2', 'memdst', 'mem_s0'};% 12 VarLinFall
     case 13; params={'nu1','nu2','g1','g2','k0_mf','alfamf','betamf'};% 13 VarRatFall
     case 14; params={'nu1','nu2','g1','g2','g3','laythk','tau_c0','k_tau'};% 4 VarExpFall1
     case 15; params={'nu1','nu2','g1','g2','k0_mf','alfamf','betamf'};% 15 VarExpFall2
     case 16; params={'nu1','nu2','g1','g2','g3','laythk','tau_c0','k_tau'};% 16 VarTempDep
     otherwise
      error('Wrong parameter')
    end
    params=sdth.sfield('addselected',struct,LI.Friction,params);
    params=struct2cell(params);params=horzcat(params{:});
    cntc.setfrictionmethod(LI.global.ire, [], LI.Friction.FrcLaw, params);
    return
    % end SDT packaging
   end

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = -1; % default: W/R contact, all patches, doesn't support icp >0
   end
   if (nargin<4 | isempty(imeth) | isempty(params))
    disp('ERROR in cntc.setfrictionmethod: imeth, params are mandatory.');
    return
   end

   cntc.call('cntc.setfrictionmethod', ire, icp, imeth, length(params), params);

  end % cntc.setfrictionmethod

  %------------------------------------------------------------------------------------------------------------



  %------------------------------------------------------------------------------------------------------------
  function []=sethertzcontact(ire, icp, ipotcn, params)
   % [ ] = cntc.sethertzcontact(ire, icp, ipotcn, params)
   %
   % set the parameters for a Hertzian contact problem
   %  ipotcn    - type of specification of the Hertzian geometry
   %  params    - depending on method that is used
   %
   %   -6: SDEC approach, union of two half ellipses                  params = [ mx, my, aa , bneg, bpos, scale ]
   %   -5: Hertzian rectangular contact, half-sizes prescribed,       params = [ mx, my, aa , bb , scale ]
   %   -4: Hertzian rectangular contact, curv+half width prescribed,  params = [ mx, my, a1 , bb , scale ]
   %   -3: Hertzian elliptical contact, semi-axes prescribed,         params = [ mx, my, aa , bb , scale ]
   %   -2: Hertzian elliptical contact, ellipticity prescribed,       params = [ mx, my, a1 , aob, scale ]
   %   -1: Hertzian elliptical contact, curvatures prescribed,        params = [ mx, my, a1 , b1 , scale ]
   %
   % dimensions:  mx, my, aob, scale [-],    a1, b1: [1/length],    aa, bneg, bpos, bb: [length]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 1: m=3, cp     - require icp>0, default 1
   % #sethertzcontact



   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.sethertzcontact: not available for icp=%d',icp));
    return
   end
   if (nargin<4 | isempty(ipotcn) | isempty(params))
    disp('ERROR in cntc.sethertzcontact: ipotcn, params are mandatory.');
    return
   end

   cntc.call('cntc.sethertzcontact', ire, icp, ipotcn, length(params), params);

  end % cntc.sethertzcontact

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function []=setinterfaciallayer(ire, icp, imeth, params)
   % [ ] = cntc.setinterfaciallayer(ire, icp, imeth, params)
   %
   % set parameters for the interfacial layer for a contact problem.
   %       *Obsolete, replaced by cntc.setmaterialparameters.
   %  imeth     - type of interfacial layer used
   %  params    - depending on method that is used
   %
   %    0: clean interface, no layer     params = [ ]
   %  2,3: Modified FASTSIM algorithm    params = [ k0_mf, alfamf, betamf ]
   %    4: elasto-plastic layer          params = [ G3, laythk, tau_c0, k_tau ]
   %
   %  dimensions: k0_mf, alfamf, betamf: [-],
   %              G3, tau_c0: [force/area],  laythk: [length],  k_tau: [force/area/length]

   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 5: m=*, wtd    - default icp=-1
   % #setinterfaciallayer



   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = -1; % default: W/R contact, all patches
   end
   if (nargin<3 | isempty(imeth))
    disp('ERROR in cntc.setinterfaciallayer: imeth is mandatory.');
    return
   end
   if (imeth==0)
    params = [ 0 ];
   end
   if (nargin<4 | isempty(params))
    disp('ERROR in cntc.setinterfaciallayer: invalid params provided.');
    return
   end

   cntc.call('cntc.setinterfaciallayer', ire, icp, imeth, length(params), params);

  end % cntc.setinterfaciallayer

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function []=setmaterialparameters(ire, icp, m_digit, rparam)
   % [ ] = cntc.setmaterialparameters(ire, icp, m_digit, rparam)
   %
   % set the M-digit and material parameters for a contact problem
   % values < 10 are used to configure the M-digit
   %    0: purely elastic material,            params = [ nu1, nu2, g1, g2 ]
   %    1: visco-elastic material,             params = [ nu1, nu2, g1, g2, fg1,   fg2,    vt1,    vt2    ] (pdf.4.1.3)
   %    2: modified Fastsim, 1 flexibility,    params = [ nu1, nu2, g1, g2, flx,   k0_mf,  alfamf, betamf ] (pdf.4.1.4)
   %    3: modified Fastsim, 3 flexibilities,  params = [ nu1, nu2, g1, g2, k0_mf, alfamf, betamf         ] (pdf.4.1.4)
   %    4: elastic + elasto-plastic 3rd body,  params = [ nu1, nu2, g1, g2, g3,    laythk, tau_c0, k_tau  ] (pdf.4.1.5)
   % values >= 10 are used to configure the M2-digit, M2 = m_digit - 10
   %   12: force-based proportional damping,   params = [ cdampn, cdampt, dfnmax, dftmax ] (pdf.4.1.7)
   %
   % dimensions: nu1, nu2        [-],        g1, g2    [force/area],
   %             fg1, fg2        [-],        vt1, vt2  [time],
   %             flx         [volume/force], k0_mf, alfamf, betamf [-],
   %             g3, tau_c0   [force/area],  laythk    [length],        k_tau [force/volume]
   %             cdampn, cdampt  [time],     dfnmax, dftmax   [force/time]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 5: m=*, wtd    - default icp=-1
   % #setmaterialparameters(pdfsection.4.1)  configure the material -2

   % parameters
   if nargin==0||ischar(ire)
    % package SDT options
    DoOpt=['Mater1(0#vd{0,Elas,1,VisElas,2,ModFast1Flex,3,ModFast3Flex,4,InterLayer}#"Mater1 options")' ...
     'nu1(#%g#"Body''s 1 poisson''s ratio")' ...
     'nu2(#%g#"Body''s 2 poisson''s ratio")' ...
     'g1(#%g#"Body''s 1 shear modulus")' ...
     'g2(#%g#"Body''s 2 shear modulus")' ...
     'g3(#%g#"Interface layer shear modulus")' ...
     'fg1(#%g#"Ratio of Kelvin-Voigt spring compliance constants")' ...
     'fg2(#%g#"Ratio of Kelvin-Voigt spring compliance constants")' ...
     'vt1(#%g#"Creep relaxation time")' ...
     'vt2(#%g#"Creep relaxation time")' ...
     'flx(#%g#"Kalker flexibility parameter")' ...
     'k0_mf(#%g#"Initial value of Kalker’s reduction factor")' ...
     'alphamf(#%g#"Fraction α_{mf}=k∞/k0")' ...
     'betamf(#%g#"Governs how quickly k changes with increasing slip area")' ...
     'laythk(#%g#"Thickness of the interface layer")' ...
     'tau_c0(#%g#"Initial yield limit at which plasticity starts to occur")' ...
     'k_tau(#%g#"Rate of increase of the yield limit with accumulated plastic deformation.")' ...
     'cdampn(#%g#"Damping coefficient for normal contact forces")' ...
     'cdampt(#%g#"Damping coefficient for tangential contact forces.")' ...
     'dfnmax(#%g#"Maximum absolute value dFn,max used for dFnel/dt.")' ...
     'dftmax(#%g#"Maximum absolute value dFτ,max used for dFτel/dt.")'];

    if nargin==0; CAM='';end
    LI=cntc.call;
    LI.Mat=cingui('paramedit -doclean2',DoOpt,{struct,ire});

    switch LI.Mat.Mater1
     case 0; params={'nu1','nu2','g1','g2'};% 0 Elas
     case 1; params={'nu1','nu2','g1','g2','fg1','fg2','vt1','vt2'};% 1 VisElas
     case 2; params={'nu1','nu2','g1','g2','flx','k0_mf','alfamf','betamf'};% 2 ModFast1Flex
     case 3; params={'nu1','nu2','g1','g2','k0_mf','alfamf','betamf'};% 3 ModFast3Flex
     case 4; params={'nu1','nu2','g1','g2','g3','laythk','tau_c0','k_tau'};% 4 InterLayer
     otherwise; error('Wrong parameter selected')
    end
    params=sdth.sfield('addselected',struct,LI.Mat,params);
    params=struct2cell(params);params=horzcat(params{:});
    cntc.setmaterialparameters(LI.global.ire, [], LI.Mat.Mater1, params);
    return
    % end SDT packaging

   end
   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = -1; % default: W/R contact, all patches
   end
   if (nargin<3 | isempty(m_digit))
    m_digit = 0; % default: linearly elastic material
   end
   if (nargin<4 | isempty(rparam))
    disp('ERROR in cntc.setmaterialparameters: rparam is mandatory.');
    return
   end

   cntc.call('cntc.setmaterialparameters', ire, icp, m_digit, length(rparam), rparam);

  end % cntc.setmaterialparameters

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function []=setmaterialproperties(ire, icp, g1, nu1, g2, nu2)
   % [ ] = cntc.setmaterialproperties(ire, icp, g1, nu1, g2, nu2)
   %
   % set the material properties for a contact problem.    Obsolete, replaced by cntc.setmaterialparameters.
   %
   %  g1, g2         - modulus of rigidity for body 1, 2 [force/area]
   %  nu1, nu2       - Poisson's ratio for body 1, 2 [-]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 5: m=*, wtd    - default icp=-1
   % #setmaterialproperties



   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = -1; % default: W/R contact, all patches
   end
   if (nargin<6 | isempty(g1) | isempty(nu1) | isempty(g2) | isempty(nu2))
    disp('ERROR in cntc.setmaterialproperties: g1, nu1, g2, nu2 are mandatory.');
    return
   end

   cntc.call('cntc.setmaterialproperties', ire, icp, g1, nu1, g2, nu2);

  end % cntc.setmaterialproperties

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function []=setmetadata(ire, icp, params, values)

   % [ ] = cntc.setmetadata(ire, icp, params, values)
   %
   % used for configuring various metadata for a contact problem
   %
   %  params   - codes of the metadata to be communicated to CONTACT
   %  values   - values of the metadata to be communicated to CONTACT
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 1: m=3, cp     - require icp>0, default 1
   % #setmetadata


   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.setmetadata: not available for icp=%d',icp));
    return
   end
   if (nargin<4 | isempty(params) | isempty(values))
    disp('ERROR in cntc.setmetadata: params and values are mandatory.');
    return
   end
   if (length(params) ~= length(values))
    disp('ERROR in cntc.setmetadata: params and values must provide same number of items.');
    return
   end

   cntc.call('cntc.setmetadata', ire, icp, length(params), params, values);

  end % cntc.setmetadata

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function []=setnormalforce(ire, icp, fn)
   % [ ] = cntc.setnormalforce(ire, icp, fn)
   %
   % set the total normal force of the bodies as a whole for a contact problem
   % Note: this function sets control digit N = 1
   %
   %  fn             - total normal force between the two bodies [force]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 1: m=3, cp     - require icp>0, default 1
   % #setnormalforce



   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.setnormalforce: not available for icp=%d',icp));
    return
   end
   if (nargin<3 | isempty(fn))
    disp('ERROR in cntc.setnormalforce: fn is mandatory.');
    return
   end

   cntc.call('cntc.setnormalforce', ire, icp, fn);

  end % cntc.setnormalforce

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function []=setboundcond(ire, icp, val)
   % #setboundcond boundary conditions -2

   if nargin==0||ischar(ire)
    % package SDT options
    DoOpt=['Cond(0#vd{0,Force,1,Penetration}#"Boundary conditions options")' ...
     'fz(#%g#"Vertical effort, train weight")' ...
     'pen(#%g#"Penetration between the 2 bodies")'];

    if nargin==0; CAM='';else; CAM=ire; end
    LI=cntc.call;
    LI.Bound=cingui('paramedit -doclean2',DoOpt,{struct,CAM});

    switch LI.Bound.Cond;
     case 0;  params={'fz'}; % 0 Force
     case 1;  params={'pen'}; % 1 Penetration
     otherwise; error('Wrong parameter selected')
    end
    params=sdth.sfield('addselected',struct,LI.Bound,params);
    Global=LI.global; if iscell(Global);Global=sdtm.toStruct(Global(:,1:2));end
    params=struct2cell(params);params=horzcat(params{:});
    %selection of the boundary condition
    if LI.Bound.Cond==0; cntc.setverticalforce(Global.ire,params);
    else; cntc.setpenetration(Global.ire, Global.icp,params);end
    return
    % end SDT packaging

   end
  end
  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function []=setpenetration(ire, icp, pen)
   % [ ] = cntc.setpenetration(ire, icp, pen)
   %
   % set the approach (penetration) of the bodies as a whole for a contact problem
   % Note: this function sets control digit N = 0
   %
   %  pen            - penetration/approach of the two bodies [length]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 1: m=3, cp     - require icp>0, default 1
   % #setpenetration


   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.setpenetration: not available for icp=%d',icp));
    return
   end
   if (nargin<3 | isempty(pen))
    disp('ERROR in cntc.setpenetration: pen is mandatory.');
    return
   end

   cntc.call('cntc.setpenetration', ire, icp, pen);

  end % cntc.setpenetration

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function []=setverticalforce(ire, fz)
   % [ ] = cntc.setverticalforce(ire, icp, fz)
   %
   % set the total vertical force between the contacting bodies for a w/r contact problem (module 1)
   % Note: this function sets control digit N = 1
   %
   %  fz             - total vertical force between the two bodies [force]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 2: m=1, wtd    - no icp needed
   % #setverticalforce -2

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(fz))
    disp('ERROR in cntc.setverticalforce: fz is mandatory.');
    return
   end

   cntc.call('cntc.setverticalforce', ire, fz);

  end % cntc.setverticalforce

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function []=setpotcontact(ire, icp, ipotcn, params)
   % [ ] = cntc.setpotcontact(ire, icp, ipotcn, params)
   %
   % set the parameters of the potential contact area for a contact problem
   %  ipotcn   - type of specification for the potential contact area
   %  params   - depending on method that is used
   %
   %  for w/r contact, ipotcn = -1:
   %    0: N/A w/r contact, fixed grid sizes,  params = [ dx, ds, n.a. ]
   %   -1: w/r contact, fixed grid sizes,  params = [ dx, ds, a_sep, d_sep, d_comb, [d_turn] ]
   %
   %  for generic contact, ipotcn > 0:
   %    1: lower-left + grid sizes,        params = [ mx, my, xl , yl , dx , dy  ]
   %    2: lower-left + upper right,       params = [ mx, my, xl , yl , xh , yh  ]
   %    3: 1st center + grid sizes,        params = [ mx, my, xc1, yc1, dx , dy  ]
   %    4: 1st center + last center,       params = [ mx, my, xc1, yc1, xcm, ycm ]
   %
   %  dimensions: mx, my [-],   a_sep [angle],
   %              dx, ds, dy, d_sep, d_comb, d_turn, xl, yl, xh, yh, xc1, yc1, xcm, ycm [length]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 5: m=*, wtd    - default icp=-1
   % #setpotcontact(pdf.3.5 & 4.3) -2

   % parameters
   if nargin==0||ischar(ire)
    % package SDT options
    DoOpt=['PosGrid(-1#vd{-1,WR,0,WR_N/A,1,LL_grid,2,LL_UR,3,CentGrid,4,2Cent}#"Potential contact area position and grid options")' ...
     'mx(#%g#"The number of discretisation elements in x-directions")' ...
     'my(#%g#"The number of discretisation elements in y-directions")' ...
     'a_sep(#%g#"Threshold angle[rad] for the angle differenceabove which patches are processed separately (can be given as 90deg)")' ...
     'dx(#%g#"Size of each grid element in longitudinal direction")' ...
     'ds(#%g#"Size of each grid element in lateral direction tangent to the rail surface")' ...
     'dy(#%g#"Size of each grid element in lateral direction")' ...
     'd_sep(#%g#"Threshold distance dsep above which patches are processed separately.")' ...
     'd_comb(#%g#"Threshold distance below which patches are combined fully, if not prohibited by the angle threshold")' ...
     'd_turn(#%g#"Threshold distance below which reference angle turning is used for patches that are processed separately.")' ...
     'xl(#%g#"x coordinate of the lower left vertex of the PCA")' ...
     'yl(#%g#"y coordinate of the lower left vertex of the PCA")' ...
     'xh(#%g#"x coordinate of the high right vertex of the PCA")' ...
     'yh(#%g#"y coordinate of the high right vertex of the PCA")' ...
     'xc1(#%g#"Center x coordinate of the (1,1) mesh of the PCA")' ...
     'yc1(#%g#"Center y coordinate of the (1,1) mesh of the PCA")' ...
     'xcm(#%g#"Center x coordinate of the (mx,my) mesh of the PCA")' ...
     'ycm(#%g#"Center y coordinate of the (mx,my) mesh of the PCA")' ];

    if nargin==0; CAM='';end
    [i1,i2]=regexp(ire,'a_sep[^,]*');
    st2=strrep(strrep(ire(i1+5:i2),'deg','*pi/180'),'rad','');% should be rad
    st2=sprintf('a_sep %.15g',eval(st2));ire=strrep(ire,ire(i1:i2),st2);
    LI=cntc.call;
    LI.PotCntc=cingui('paramedit -doclean2',DoOpt,{struct,ire});

    switch LI.PotCntc.PosGrid
     case -1; params={'dx','ds','a_sep','d_sep','d_comb','[d_turn]'};   % -1 WR
     case 0; params={'dx','ds', 'n.a'};                                 % 0 WR_N/A
     case 1; params={'mx', 'my', 'x1', 'y1', 'dx', 'dy'};               % 1 LL_grid
     case 2; params={'mx', 'my', 'x1', 'y1', 'xh', 'yh'};               % 2 LL_UR
     case 3; params={'mx', 'my', 'xc1', 'yc1', 'dx', 'dy'};             % 3 CentGrid
     case 4; params={'mx', 'my', 'xc1', 'yc1', 'xcm', 'ycm'};           % 4 2Cent
     otherwise; error('Wrong parameter selected');
    end

    params=sdth.sfield('addselected',struct,LI.PotCntc,params);
    params=struct2cell(params);params=horzcat(params{:});
    cntc.setpotcontact(LI.global.ire, -1, LI.PotCntc.PosGrid, params);%XXXGAE icp =-1 bizarre
    return
    % end SDT packaging
   end

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = -1;
   end
   if (nargin<4 | isempty(ipotcn) | isempty(params))
    disp('ERROR in cntc.setpotcontact: ipotcn, params are mandatory.');
    return
   end

   cntc.call('cntc.setpotcontact', ire, icp, ipotcn, length(params), params);

  end % cntc.setpotcontact

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function []=setprofileinputfname(ire, fname, iparam, rparam)

   % [ ] = cntc.setprofileinputfname(ire, fname, iparam, rparam)
   %
   % set a wheel or rail profile filename for a wheel-rail contact problem
   %
   %  fname          - string: name of profile file, absolute path or relative wrt effective working folder
   %  iparam         - integer configuration parameters
   %                     1: iswheel     0 = rail, 1 = wheel profile, -1 = taken from file extension (default)
   %                     2:  -        not used
   %                     3: mirrory   0 or -1 = no mirroring (default), 1 = mirror y coordinate values
   %                     4: mirrorz   0 = autodetect (default), -1 = no mirroring, 1 = mirror z values
   %                     5: errhndl   configuration of error handling.
   %                                   -2 = continue as much as possible, suppress error messages;
   %                                   -1 = suppress warnings; 0: warn and continue (default);
   %                                    1 = signal errors and abort
   %                     6: ismooth   selection of smoothing method. 0 = original smoothing spline (default),
   %                                    1 = weighted PP smoothing spline, 2 = weighted smoothing B-spline (best)
   %  rparam         - real configuration parameters
   %                     1: sclfac    scaling factor for conversion to [mm], e.g. 1e3 for data given in [m]
   %                                  default (sclfac<=0): using the active unit convention
   %                     2: smooth    smoothing parameter lambda for non-weighted spline or l_filt for
   %                                  weighted spline smoothing
   %                     3: maxomit   fraction: signal error if more than maxomit of profile points are
   %                                  discarded after cleanup of profile. Default 0.5, use 1 to disable check.
   %                     4: zigthrs   angle threshold for zig-zag detection. Default 5/6*pi, >=pi to disable.
   %                     5: kinkhigh  angle threshold for kink detection. Default pi/6, >=pi to disable.
   %                     6: kinklow   angle threshold for neighbouring points in kink detection.
   %                                  default kinkhigh/5.
   %                     7: kinkwid   half-width of window used for kink detection, [len], default 2 mm
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 2: m=1, wtd    - no icp needed
   % #setprofileinputfname(pdfsubsection.2.1.2) set a profile FileName -2

   % parameters
   if nargin==0||ischar(ire)
    % package SDT options
    DoOpt=['iswheel(#%g#"1 if is wheel, 0 if is rail")' ...
     'fname(#%s#"File name")' ...
     'mirrory(#%g#"Mirror y coordinate values")' ...
     'mirrorz(#%g#"Mirror z coordinate values")' ...
     'notuse(#%g#"Not use parameter, filled with 0")' ...
     'errhndl(#%g#"Configuration of error handling")' ...
     'ismooth(#%g#"Selection of smoothing method")' ...
     'sclfac(#%g#"Scaling factor for conversion to [mm]")' ...
     'smooth(#%g#"Smoothing parameter lambda")' ...
     'maxomit(#%g#"Error if more than maxomit of profile points are discarded after cleanup of profile")' ...
     'zigthrs(#%g#"Angle threshold for zig-zag detection")' ...
     'kinkhigh(#%g#"Angle threshold for kink detection")' ...
     'kinklow(#%g#"Angle threshold for kink detection")' ...
     'kinkwid(#%g#"Half-width of window used for kink detection")'];

    if nargin==0; CAM='';else; CAM=ire;end
    LI=cntc.call;

    ProfileFname=cingui('paramedit -doclean2',DoOpt,{struct,CAM});
    if ProfileFname.iswheel == 1
     LI.WheelProfile=ProfileFname;
    else
     LI.RailProfile=ProfileFname;
    end
    fname=sdth.sfield('addselected',struct,ProfileFname,{'fname'});
    iparams=sdth.sfield('addselected',struct,ProfileFname,{'iswheel','notuse','mirrory','mirrorz','errhndl','ismooth'});
    rparams=sdth.sfield('addselected',struct,ProfileFname,{'sclfac','smooth','maxomit','zigthrs','kinkhigh','kinklow','kinkwid'});

    iparams.notuse=0;iparams=struct2cell(iparams);iparams=horzcat(iparams{:});
    rparams=struct2cell(rparams);rparams=horzcat(rparams{:});
    fname=struct2cell(fname);fname=fname{1};
    cntc.setprofileinputfname(LI.global.ire, fname, iparams, rparams);
    return

   end
   % end SDT packaging

   if (nargin<1 | isempty(ire));ire = 1;end
   if (nargin<2 | isempty(fname))
    disp('ERROR in cntc.setprofileinputfname: fname is mandatory.');
    return
   end
   if (nargin<3 | isempty(iparam)); iparam = [ -1 ];end
   if (nargin<4 | isempty(rparam));rparam = [ 1., 0. ];end

   LI=cntc.call;
   if ~exist(fname,'file'); fname=fullfile(LI.ProjectWd,fname);end
   fname=sdtu.f.safe(fname);

   len_fname = length(fname);nints= length(iparam);nreals=length(rparam);

   cntc.call('cntc.setprofileinputfname', ire, fname, len_fname, nints, iparam, nreals, rparam);

  end % cntc.setprofileinputfname

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function []=setprofileinputvalues(ire, values, iparam, rparam)
   % [ ] = cntc.setprofileinputvalues(ire, values, iparam, rparam)
   %
   % Set a wheel or rail profile for a wheel-rail contact problem using a table of values
   %
   %  values         - lists of profile points. Can be (npoint,2) or (2,npoint), or
   %                        1d array of length (2*npoint), ordered [y1,z1, ... yi,zi, ... yn,zn]
   %  iparam         - integer configuration parameters
   %                     1: itype     0 = rail, 1 = wheel profile, -1 = taken from file extension (default)
   %                     2:  -        not used
   %                     3: mirrory   0 or -1 = no mirroring (default), 1 = mirror y coordinate values
   %                     4: mirrorz   0 = autodetect (default), -1 = no mirroring, 1 = mirror z values
   %                     5: errhndl   configuration of error handling.
   %                                   -2 = continue as much as possible, suppress error messages;
   %                                   -1 = suppress warnings; 0: warn and continue (default);
   %                                    1 = signal errors and abort
   %                     6: ismooth   selection of smoothing method. 0 = original smoothing spline (default),
   %                                    1 = weighted PP smoothing spline, 2 = weighted smoothing B-spline (best)
   %  rparam         - real configuration parameters
   %                     1: sclfac    scaling factor for conversion to [mm], e.g. 1e3 for data given in [m]
   %                                  default (sclfac<=0): using the active unit convention
   %                     2: smooth    smoothing parameter lambda for non-weighted spline or l_filt for
   %                                  weighted spline smoothing
   %                     3: maxomit   fraction: signal error if more than maxomit of profile points are
   %                                  discarded after cleanup of profile. Default 0.5, use 1 to disable check.
   %                     4: zigthrs   angle threshold for zig-zag detection. Default 5/6*pi, >=pi to disable.
   %                     5: kinkhigh  angle threshold for kink detection. Default pi/6, >=pi to disable.
   %                     6: kinklow   angle threshold for neighbouring points in kink detection.
   %                                  default kinkhigh/5.
   %                     7: kinkwid   half-width of window used for kink detection, [len], default 2 mm
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 2: m=1, wtd    - no icp needed
   % #setprofileinputvalues -2

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(values))
    disp('ERROR in cntc.setprofileinputvalues: values is mandatory.');
    return
   end
   if (nargin<3 | isempty(iparam))
    iparam = [  0 ];
   end
   if (nargin<4 | isempty(rparam))
    rparam = [ -1, 0. ];
   end

   % check sizes, reshape input values to 1-d array of 2*npoint elements

   arrsiz = size(values);
   if (min(arrsiz)==2) % 2-d array
    npoint = max(arrsiz);
    % transpose to (2,npoint) array
    if (arrsiz(1)>2)
     values = values';
    end
    % reshape to 1D array
    values = reshape(values, 2*npoint, 1);
   elseif (min(arrsiz)==1)
    npoint = max(arrsiz) / 2;
   else
    disp('ERROR in cntc.setprofileinputvalues: values must have 2*npoint values.');
    return
   end

   nints     = length(iparam);
   nreals    = length(rparam);
   cntc.call('cntc.setprofileinputvalues', ire, npoint, values, nints, iparam, nreals, rparam);

  end % cntc.setprofileinputvalues

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function []=setreferencevelocity(ire, icp, veloc)
   % [ ] = cntc.setreferencevelocity(ire, icp, veloc)
   %
   % set the rolling velocity for a contact problem
   %
   %  veloc          - absolute rolling velocity [veloc]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 1: m=3, cp     - require icp>0, default 1
   % #setreferencevelocity



   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.setreferencevelocity: not available for icp=%d',icp));
    return
   end
   if (nargin<3 | isempty(veloc))
    disp('ERROR in cntc.setreferencevelocity: veloc is mandatory.');
    return
   end

   cntc.call('cntc.setreferencevelocity', ire, icp, veloc);

  end % cntc.setreferencevelocity

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function []=setrollingstepsize(ire, icp, chi, dq)
   % [ ] = cntc.setrollingstepsize(ire, icp, chi, dq)
   %
   % set the rolling direction and step size for a contact problem
   %
   % in w/r contact,    icp = -1
   %     chi            - ignored
   %     dqrel          - rolling step size relative to grid size dx [-]
   % in generic contact, icp > 0,
   %     chi            - rolling direction [angle]
   %     dq             - rolling step size [length]
   %
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 5: m=*, wtd    - default icp=-1
   % #setrollingstepsize(pdf.4.5) -2

   if nargin==0||ischar(ire)
    % package SDT options
    DoOpt=['StepSize(0#vd{-1,WRCn,1,otherCn}#"Rolling direction and step size")' ...
     'chi(#%g#"Rolling direction angle")' ...
     'dqrel(#%g#"Rolling step size relative to the grid size")' ...
     'dq(#%g#"Rolling step size")'];

    if nargin==0; CAM='';else; CAM=ire ;end
    LI=cntc.call;
    LI.Rolling=cingui('paramedit -doclean2',DoOpt,{struct,CAM});

    switch LI.Rolling.StepSize
     case -1
      % -1 WRCb
      chi=[ ];
      val=sdth.sfield('addselected',struct,LI.Rolling,{'dqrel'});
     case 0
      % 0 OtherCn
      chi=sdth.sfield('addselected',struct,LI.Rolling,{'chi'});
      val=sdth.sfield('addselected',struct,LI.Rolling,{'dq'});
     otherwise
      error('Wrong parameter selected')
    end
    Global=LI.global; if iscell(Global);Global=sdtm.toStruct(Global(:,1:2));end
    if isstruct(chi); chi=struct2cell(chi); end
    val=struct2cell(val);val=val{:};
    cntc.setrollingstepsize(Global.ire, LI.Rolling.StepSize, chi, val);
    return
    % end SDT packaging

   end

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = -1;
    chi =  0;
   end
   if (nargin<4 | isempty(dq))
    disp('ERROR in cntc.setrollingstepsize: dq is mandatory.');
    return
   end

   cntc.call('cntc.setrollingstepsize', ire, icp, chi, dq);

  end % cntc.setrollingstepsize

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function [out]=setsolverflags(ire, icp, gdigit, iparam, rparam)
   % [ ] = cntc.setsolverflags(ire, icp, gdigit, iparam, rparam)
   %
   % set parameters for the iterative solution algorithms
   %  gdigit          - G-digit, which solvers to use
   %  iparams, rparam   - depending on method that is used
   %
   %  0: default           - iparam = [maxgs, maxin, maxnr, maxout],         rparam = [eps]
   %  1: maintain          - iparam = [ ],                                   rparam = [ ]
   %  2: ConvexGS          - iparam = [maxgs, maxin, maxnr, maxout, inislp], rparam = [eps, omegah, omegas, omgslp]
   %  3: SteadyGS          - iparam = [maxgs, maxin, maxnr, maxout, inislp], rparam = [eps, omegah, omegas, omgslp]
   %  4: SlipVel           - iparam = [maxgs, maxin, maxnr, maxout, inislp], rparam = [eps, omgslp]
   %  5: GDsteady          - iparam = [maxgs, maxin, maxnr, maxout]          rparam = [eps, pow_s, omg_s]
   %  6: SensCalc          - iparam = [mxsens],                              rparam = [epsens]

   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 5: m=*, wtd    - default icp=-1
   % #setsolverflags(pdfsection.4.6) Solver definition -2

   if nargin==0||ischar(ire)
    % package SDT options
    DoOpt=['GauSei(0#vd{0,default,1,maintain,2,ConvexGS,3,SteadyGS,4,SlipVel,5,GDSteady, 6,SensCalc}#"Gauss Seidel options")' ...
     'maxgs(10#%g#"Intern loop maximum number of iterations")' ...
     'maxin(10#%g#"NORM/TANG algorithms maximum number of iterations")' ...
     'maxnr(10#%g#"Newton-raphson Maximum number of iterations ")' ...
     'maxout(10#%g#"Outer loop maximum iterations")' ...
     'inislp(#%g#"Initial estimate for slip velocity")' ...
     'eps(#%g#"Stop-criteria")' ...
     'omegah(#%g#"Adhesion area relaxation parameter for ConvexGS/SteadyGS")' ...
     'omegas(#%g#"Slip area relaxation parameter for ConvexGS/SteadyGS")' ...
     'omgslp(#%g#"Relaxation parameter for the slip velocity")' ...
     'mxsens(#%g#"Sensitivity approach maximum iteration")' ...
     'pow_s(#%g#"Gradient Descent parameter")' ...
     'omg_s(#%g#"Gradient Descent parameter")' ...
     'epsens(#%g#"Sensivity approach tolerance")'];

    if nargin==0; CAM='';else; CAM=ire ;end
    LI=cntc.call;
    LI.Solver=cingui('paramedit -doclean2',DoOpt,{struct,CAM});

    switch LI.Solver.GauSei
     case 0; iparam={'maxgs','maxin','maxnr','maxout'};rparam={'eps'}; % 0 Default
     case 1; iparam=[ ];rparam=[ ];                             % 1 Maintain
     case 2; iparam={'maxgs','maxin','maxnr','maxout','inislp'};% 2 ConvexGS
      rparam={'eps','omegah','omegas','omgslp'};
     case 3; iparam={'maxgs','maxin','maxnr','maxout','inislp'};% 3 SteadyGS
      rparam={'eps','omegah','omegas','omgslp'};
     case 4; iparam={'maxgs','maxin','maxnr','maxout','inislp'};% 4 SlipVel
      rparam={'eps','omgslp'};
     case 5; iparam={'maxgs','maxin','maxnr','maxout'};         % 5 GDSteady
      rparam={'eps','pow_s','omg_s'};
     case 6; iparam={'mxsens'};rparam={'epsens'};               % Flags Sensitivities
     otherwise; error('Wrong parameter selected')
    end

    Global=LI.global; if iscell(Global);Global=sdtm.toStruct(Global(:,1:2));end
    iparam=sdth.sfield('addselected',struct,LI.Solver,iparam);
    iparam=struct2cell(iparam);iparam=horzcat(iparam{:});
    rparam=sdth.sfield('addselected',struct,LI.Solver,rparam);
    rparam=struct2cell(rparam);rparam=horzcat(rparam{:});
    cntc.setsolverflags(Global.ire, [], LI.Solver.GauSei, iparam, rparam); % config solver Gauss-Seidel
    return
    % end SDT packaging
   end

   if (nargin<1 | isempty(ire));ire = 1;end
   if (nargin<2 | isempty(icp))
    icp = -1; % default: W/R contact, all patches
   end
   if (nargin<3 | isempty(gdigit))
    disp('ERROR in cntc.setsolverflags: gdigit is mandatory.');
    return
   end
   if (gdigit==1)
    nints = 0; iparam = 0; nreals = 0; rparam = 0;
   end
   if (nargin<5 | isempty(iparam) | isempty(rparam))
    disp('ERROR in cntc.setsolverflags: iparam, rparam are mandatory.');
    return
   end
   if (gdigit~=1)
    nints = length(iparam); nreals = length(rparam);
   end

   cntc.call('cntc.setsolverflags', ire, icp, gdigit, nints, iparam, nreals, rparam);

  end % cntc.setsolverflags

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function []=settangentialforces(ire, icp, fx, fy)
   % [ ] = cntc.settangentialforces(ire, icp, fx, fy)
   %
   % set the total tangential forces a contact problem
   % note: fx is ignored when the F-digit is 0, fy is ignored when F=0 or 1.
   %
   %  fx, fy       - total tangential forces relative to fstat*fn [-]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 1: m=3, cp     - require icp>0, default 1
   % #settangentialforces



   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (nargin<4 | isempty(fx) | isempty(fy))
    disp('ERROR in cntc.settangentialforces: total forces fx, fy are mandatory.');
    return
   end

   cntc.call('cntc.settangentialforces', ire, icp, fx, fy);

  end % cntc.settangentialforces

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function []=settemperaturedata(ire, icp, imeth, params)
   % [ ] = cntc.settemperaturedata(ire, icp, imeth, params)
   %
   % set parameters for the temperature calculation for a contact problem
   %  imeth    - type of temperature model used (H-digit)
   %  params   - depending on method that is used
   %    0: no temperature calculation,     params = []
   %    1: keep old parameters,            params = []
   %    3: calculate temperature based on new parameters and steady rolling,
   %       params = [bktemp1, heatcp1, lambda1, dens1, bktemp2, heatcp2, lambda2, dens2]
   %
   % dimensions:  bktemp: [C],  heatcp: [J/kg-C],  lambda: [W/length-C],  dens: [kg/length^3]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 5: m=*, wtd    - default icp=-1
   % #settemperaturedata



   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = -1; % default: W/R contact, all patches
   end
   if (nargin<3 | isempty(imeth))
    disp('ERROR in cntc.settemperaturedata: imeth is mandatory.');
    return
   end
   if (imeth<=1)
    params = [ 0 ];
   end
   if (nargin<4 | isempty(params))
    disp('ERROR in cntc.settemperaturedata: invalid params provided.');
    return
   end

   cntc.call('cntc.settemperaturedata', ire, icp, imeth, length(params), params);

  end % cntc.settemperaturedata

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function []=settimestep(ire, icp, dt)
   % [ ] = cntc.settimestep(ire, icp, dt)
   %
   % set the time step size dt for a contact problem, in particular for T = 0 or 1 (shifts)
   %
   %  dt          - time step size [time]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 1: m=3, cp     - require icp>0, default 1
   % #settimestep



   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.settimestep: not available for icp=%d',icp));
    return
   end
   if (nargin<3 | isempty(dt))
    disp('ERROR in cntc.settimestep: dt is mandatory.');
    return
   end

   cntc.call('cntc.settimestep', ire, icp, dt);

  end % cntc.settimestep

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function []=settrackdimensions(ire, ztrack, params)
   % [ ] = cntc.settrackdimensions(ire, ztrack, params)
   %
   % set the track or roller-rig description for a wheel-rail contact problem
   %  ztrack    - control digit ZTRACK
   %  params    - depending on method that is used
   %
   %    0: maintain track dimensions     params = [ ]
   %    1: new design track dimensions   params = [gaught, gaugsq, gaugwd, cant, nomrad],   if gaught >  0,
   %                                         or   [gaught, raily0, railz0, cant, nomrad],   if gaught <= 0.
   %    2: new track deviations          params = [dyrail, dzrail, drollr, vyrail, vzrail, vrollr]
   %    3: new dimensions & track deviations for current side of the track
   %                                     params = params(1:5) cf. Z=1 followed by params(6:11) cf. Z=2;
   %                                              additionally, [kyrail, fyrail, kzrail, fzrail]
   %
   % dimensions: gaught, gaugwd, raily0, railz0, nomrad, dyrail, dzrail [length],    cant, drollr [angle]
   %                                                     vyrail, vzrail [veloc],           vrollr [ang.veloc]
   %                                                kyrail, kzrail [force/length], fyrail, fzrail [force]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 2: m=1, wtd    - no icp needed
   % #settrackdimensions(pdf.3.3) -2

   if nargin==0||ischar(ire)
    % package SDT options
    DoOpt=['Design(0#vd{0,Maintain,1,NewDim,2,NewDevia,3,NewBoth}#"Track dimension and deviation options")' ...
     'gaught(#%g#"Height below the track plane at which the gauge width is measured.")' ...
     'gaugwd(#%g#"Distance between the inner faces of the two rails.")' ...
     'gaugsq(#%g#"Reserved for selecting second or further gauge faces instead of the leftmost one")' ...
     'nomrad(#%g#"Nominal radius of the rollers")' ...
     'dyrail(#%g#"Offset Δyrail of the rail profile reference from the design")' ...
     'dzrail(#%g#"Offset Δzrail of the rail profile reference from the design")' ...
     'cant(#%g#"Rail cant")' ...
     'drollr(#%g#"Rotation Δφrail of the rail profile")'...
     'vyrail(#%g#"Velocity vrail  y of the rail origin")' ...
     'vzrail(#%g#"Velocity vrail  z of the rail origin")' ...
     'raily0(#%g#"Position y of the rail origin with respect to track coordinates")' ...
     'railz0(#%g#"Position z of the rail origin with respect to track coordinates")' ...
     'vrollr(#%g#"Angular velocity vrailφ of the rail origin")' ...
     'kyrail(#%g#"xxx")' ...
     'kzrail(#%g#"xxx")' ...
     'fyrail(#%g#"xxx")' ...
     'fzrail(#%g#"xxx")' ];

    if nargin==0; CAM='';else; CAM=ire ;end
    LI=cntc.call;
    LI.Track=cingui('paramedit -doclean2',DoOpt,{struct,CAM});
    % XXXGAE EB
    switch LI.Track.Design
     case 0; params=[];% 0 Maintain
     case 1
      % 1 NewDim
      gaught=sdth.sfield('addselected',struct,LI.Track,{'gaught'});
      gaught=struct2cell(gaught); gaught=gaught{:};
      if gaught>0
       params=sdth.sfield('addselected',struct,LI.Track,{'gaught','gaugsq', 'gaugwd','cant', 'nomrad'});   % values
      else
       params=sdth.sfield('addselected',struct,LI.Track,{'gaught', 'raily0', 'railz0', 'cant', 'nomrad'});   % values
      end
     case 2
      % 2 NewDevia
      params=sdth.sfield('addselected',struct,LI.Track,{'dyrail', 'dzrail', 'drollr', 'vyrail', 'vzrail', 'vrollr'});   % values
     case 3
      % 3 NewBoth
      gaught=sdth.sfield('addselected',struct,LI.Track,{'gaught'});
      gaught=struct2cell(gaught); gaught=gaught{:};
      if gaught>0
       params=sdth.sfield('addselected',struct,LI.Track,{'gaught','gaugsq', ...
        'gaugwd','cant', 'nomrad','dyrail', 'dzrail', 'drollr', 'vyrail', 'vzrail', 'vrollr'});   % values
      else
       params=sdth.sfield('addselected',struct,LI.Track,{'gaught', 'raily0', ...
        'railz0', 'cant', 'nomrad','dyrail', 'dzrail', 'drollr', 'vyrail', 'vzrail', 'vrollr'});   % values
      end
     otherwise
      error('Wrong parameter selected')
    end
    Global=LI.global; if iscell(Global);Global=sdtm.toStruct(Global(:,1:2));end
    params=struct2cell(params);params=horzcat(params{:});
    cntc.settrackdimensions(Global.ire, LI.Track.Design, params); % config solver Gauss-Seidel
    return
    % end SDT packaging

   end

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(ztrack))
    disp('ERROR in cntc.settrackdimensions: ztrack is mandatory.');
    return
   end
   if (nargin<3 | isempty(params))
    disp('ERROR in cntc.settrackdimensions: invalid params provided.');
    return
   end

   cntc.call('cntc.settrackdimensions', ire, ztrack, length(params), params);

  end % cntc.settrackdimensions

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function []=setundeformeddistc(ire, icp, ibase, prmudf)
   % [ ] = cntc.setundeformeddistc(ire, icp, ibase, prmudf)
   %
   % set the undeformed distance function through a formula or by element-wise specification
   %
   %  ibase          - type of undeformed distance specification
   %  nparam         - number of parameters provided
   %  prmudf(nparam) - parameters of undef.dist, depending on method that is used
   %
   %    1: quadratic function            6    params = [b1, b2, b3, b4, b5, b6]
   %    2: circular-x, piecewise-lin-y   5+nn params = [nn, xm, rm, y1, dy1], [b(k), k=1..nn]
   %    3: quadratic plus two sines      8    params = [b1, b2, b3, b4, b5, b6, b7, b8]
   %    9: elementwise specification     npot params = [h(i), i=1..npot] - undeformed distance per elem. [length]
   %                                                   note: positive values == separation between profiles
   %
   % when ibase=9, prmudf may be of size (my,mx) as well, with nparam=mx*my.
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 1: m=3, cp     - require icp>0, default 1
   % #setundeformeddistc



   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in cntc.setundeformeddistc: not available for icp=%d',icp));
    return
   end
   if (nargin<4 | isempty(ibase) | isempty(prmudf))
    disp('ERROR in cntc.setundeformeddistc: ibase, prmudf are mandatory.');
    return
   end

   % check size of prmudf-array

   nparam = prod(size(prmudf));
   if (ibase==9)
    [ mx, my ] = cntc.getnumelements(ire, icp);
    if (nparam~=mx*my)
     disp(sprintf('ERROR in cntc.setundeformeddistc: expecting nparam = mx*my = %d',mx*my));
     return
    end
    if (min(size(prmudf))>1 & any(size(prmudf)~=[my mx]))
     disp('ERROR in cntc.setundeformeddistc:');
     disp(sprintf('      expecting array prmudf with size (%d,%d), got size (%d,%d)',my,mx,size(prmudf)))
     return
    end
   end

   % if prmudf is 2-dimensional, reshape to 1-d array

   if (ibase==9 & min(size(prmudf))>1)
    prmudf = reshape(prmudf', nparam, 1);
   end

   cntc.call('cntc.setundeformeddistc', ire, icp, ibase, nparam, prmudf);

  end % cntc.setundeformeddistc

  %------------------------------------------------------------------------------------------------------------


  %------------------------------------------------------------------------------------------------------------
  function []=setwheelsetdimensions(ire, ewheel, params)
   % [ ] = cntc.setwheelsetdimensions(ire, ewheel, params)
   %
   % set the wheelset description for a wheel-rail contact problem
   %  ewheel         - control digit EWHEEL
   %  nparam         - number of parameters provided
   %  params(nparam) - depending on method that is used
   %
   %  E=1,2,4: keep old wheelset dimensions, ignore params provided
   %  E=3,5 : new wheelset geometry   params = [fbdist, fbpos, nomrad]
   %
   %  dimensions:  fbdist, fbpos, nomrad [length]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 2: m=1, wtd    - no icp needed
   % #setwheelsetdimensions(pdf.3.4) -2

   if nargin==0||ischar(ire)
    % package SDT options
    DoOpt=['Ewheel(#vd{0,Maintain,1,NewPos,2,NewPosVel,3,NewDimProfPosVel,4' ...
     ',NewDimPosVelFlex,5,NewAll}#"Wheelset options")' ...
     'fbdist(#%g#"Lateral distance between the flange backs of the two wheels of the wheelset")' ...
     'fbpos(#%g#"Lateral position of the flange back with respect to the wheel profile origin Ow")' ...
     'nomrad(#%g#"Nominal radius of the wheel")'];

    if nargin==0; CAM='';else; CAM=ire ;end
    LI=cntc.call;
    LI.wheelsetDim=cingui('paramedit -doclean2',DoOpt,{struct,CAM});

    switch LI.wheelsetDim.Ewheel
     case {0,1,2,4}; params=[];    % 0 Maintain
     case {3,5}; params=sdth.sfield('addselected',struct,LI.wheelsetDim, ...
       {'fbdist','fbpos','nomrad'}); % 1 NewDim
      params=struct2cell(params);params=horzcat(params{:});
     otherwise; error('Wrong parameter selected')
    end
    Global=LI.global; if iscell(Global);Global=sdtm.toStruct(Global(:,1:2));end
    cntc.setwheelsetdimensions(Global.ire, LI.wheelsetDim.Ewheel, params);
    return
    % end SDT packaging

   end

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<3 | isempty(ewheel) | isempty(params))
    disp('ERROR in cntc.setwheelsetdimensions: ewheel and params are mandatory.');
    return
   end

   cntc.call('cntc.setwheelsetdimensions', ire, ewheel, length(params), params);

  end % cntc.setwheelsetdimensions

  %------------------------------------------------------------------------------------------------------------

  function []=setwheelsetflexibility(ire, ewheel, params)
   %------------------------------------------------------------------------------------------------------------
   % [ ] = cntc.setwheelsetflexibility(ire, ewheel, params)
   %
   % set the description of wheelset flexibilities for a wheel-rail contact problem
   %    ewheel         - type of wheelset flexibilities specification (E-digit)
   %    nparam         - number of parameters provided
   %    params(nparam) - depending on method that is used
   %
   %    E=0  : keep wheelset flexibility from previous specification, ignore params provided
   %    E=1-3: no wheelset flexibility               params = []
   %    E=4,5: new wheelset flexibility parameters   params = [dxwhl, dywhl, dzwhl, drollw, dyaww, dpitchw,
   %                                                           vxwhl, vywhl, vzwhl, vrollw, vyaww, vpitchw],
   %
   %    dimensions:   dxwhl, dywhl, dzwhl [length],  drollw, dyaww, dpitchw [angle]
   %                  vxwhl, vywhl, vzwhl [veloc],   vrollw, vyaww, vpitchw [ang.veloc]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 2: m=1, wtd    - no icp needed
   % #setwheelsetflexibility  -2

   if nargin==0||ischar(ire)
    % package SDT options
    DoOpt=['Ewheel(#vd{0,Maintain,1,NewPos,2,NewPosVel,3,NewDimProfPosVel,4,' ...
     'NewDimPosVelFlex,5,NewAll}#"Wheelset options")' ...
     'dxwhl(#%g#"xxx")' ...
     'dywhl(#%g#"xxx")' ...
     'dzwhl(#%g#"xxx")' ...
     'drollw(#%g#"xxx")' ...
     'dyaww(#%g#"xxx")' ...
     'dpitchw(#%g#"xxx")'...
     'vxwhl(#%g#"xxx")' ...
     'vywhl(#%g#"xxx")' ...
     'vzwhl(#%g#"xxx")' ...
     'vrollw(#%g#"xxx")' ...
     'vyaww(#%g#"xxx")' ...
     'vpitchw(#%g#"xxx")'];

    if nargin==0; CAM='';else; CAM=ire ;end
    LI=cntc.call;
    LI.wheelsetFlex=cingui('paramedit -doclean2',DoOpt,{struct,CAM});

    switch LI.wheelsetFlex.Ewheel;
     case 0; params=[];% Maintain
     case {1:3}; params=zeros(1,12); % NoFlex
     case {4,5}; % New flexibilities with symmetry
      params={'dxwhl', 'dywhl', 'dzwhl','drollw', 'dyaww', 'dpitchw', 'vxwhl', ...
       'vywhl', 'vzwhl', 'vrollw', 'vyaww', 'vpitchw'};
     otherwise; error('Wrong parameter')
      %% xxxGAE
    end
    Global=LI.global; if iscell(Global);Global=sdtm.toStruct(Global(:,1:2));end
    params=sdth.sfield('addselected',struct,LI.wheelsetFlex,params);
    params=struct2cell(params);params=horzcat(params{:});
    cntc.setwheelsetflexibility(Global.ire, LI.wheelsetFlex.Ewheel, params);
    return
    % end SDT packaging

   end
   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(ewheel))
    ewheel = 0;
   end
   if (nargin<3 | isempty(params))
    disp('ERROR in cntc.setwheelsetflexibility: invalid params provided.');
    return
   end

   cntc.call('cntc.setwheelsetflexibility', ire, ewheel, length(params), params);

  end % cntc.setwheelsetflexibility

  %------------------------------------------------------------------------------------------------------------


  function []=setwheelsetposition(ire, ewheel, params)
   %------------------------------------------------------------------------------------------------------------
   % [ ] = cntc.setwheelsetposition(ire, ewheel, params)
   %
   % set the wheelset position state data for a wheel-rail contact problem
   %  ewheel         - type of position specification (E-digit)
   %  nparam         - number of parameters provided
   %  params(nparam) - depending on method that is used
   %
   %  E=0  : keep wheelset position from previous specification, ignore params provided
   %  E=1-5: new wheelset position   params = [s_ws, y_ws, z_ws, roll_ws, yaw_ws, pitch_ws]
   %
   %  dimensions:   s_ws, y_ws, z_ws [length],       roll, yaw, pitch [angle]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 2: m=1, wtd    - no icp needed
   % #setwheelsetposition(pdf.3.4) -2

   if nargin==0||ischar(ire)
    % package SDT options
    DoOpt=['Ewheel(#vd{0,Maintain,1,NewPos,2,NewPosVel,3,NewDimProfPosVel,4,' ...
     'NewDimPosVelFlex,5,NewAll}#"Wheelset options")' ...
     's_ws(#%g#"Wheelset position along the track center line")' ...
     'y_ws(#%g#"Lateral position of the wheelset center in track coordinates.")' ...
     'z_ws(#%g#"Vertical position zws of the wheelset center")' ...
     'roll_ws(#%g#"Wheelset roll angle with respect to the track plane.")' ...
     'yaw_ws(#%g#"Wheelset yaw angle with respect to the track center line")' ...
     'pitch_ws(#%g#"Wheelset pitch angle, i.e. rotation about the wheelset axle")'];

    if nargin==0; CAM='';else; CAM=ire ;end
    LI=cntc.call;
    LI.wheelsetPos=cingui('paramedit -doclean2',DoOpt,{struct,CAM});

    switch LI.wheelsetPos.Ewheel
     case 0 % Maintain
      params=[];
     case {1:5}; % 1 NewPos
      params=sdth.sfield('addselected',struct,LI.wheelsetPos,{'s_ws', 'y_ws', 'z_ws', ...
       'roll', 'yaw', 'pitch'});
     otherwise
      error('Wrong parameter')
      %% xxxGAE
    end
    Global=LI.global; if iscell(Global);Global=sdtm.toStruct(Global(:,1:2));end
    params=struct2cell(params);params=horzcat(params{:});
    cntc.setwheelsetposition(Global.ire, LI.wheelsetPos.Ewheel, params);
    return
    % end SDT packaging

   end

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(ewheel))
    ewheel = 1;
   end
   if (nargin<3 | isempty(params))
    disp('ERROR in cntc.setwheelsetposition: invalid params provided.');
    return
   end

   cntc.call('cntc.setwheelsetposition', ire, ewheel, length(params), params);

  end % cntc.setwheelsetposition

  %------------------------------------------------------------------------------------------------------------


  function []=setwheelsetvelocity(ire, ewheel, params)
   %------------------------------------------------------------------------------------------------------------
   % [ ] = cntc.setwheelsetvelocity(ire, ewheel, params)
   %
   % continue the wheelset description for a wheel-rail contact problem, set wheelset velocity data
   %  ewheel         - type of velocity specification (E-digit)
   %  nparam         - number of parameters provided
   %  params(nparam) - depending on method that is used
   %
   %  E=0-1: no new wheelset velocity    params = [ ]
   %  E=2-5:    new wheelset velocity    params = [ vs, vy, vz, vroll, vyaw, vpitch ]
   %            for roller-rigs vx is replaced by rpitch (C1=4,5)
   %            position increments v * dt are used in transient shifts (T=1)
   %
   %  dimensions:     vs, vy, vz    : [veloc]   vroll, vyaw, vpitch, rpitch  : [ang.veloc]
   %                  shft_sws--zws : [length]  shft_rol, shft_yaw, shft_pit : [angle]
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 2: m=1, wtd    - no icp needed
   % #setwheelsetvelocity(pdf.3.4) -2

   if nargin==0||ischar(ire)
    % package SDT options
    DoOpt=['Ewheel(#vd{0,Maintain,1,NewPos,2,NewPosVel,3,NewDimProfPosVel,4,' ...
     'NewDimPosVelFlex,5,NewAll}#"Wheelset options")' ...
     'vs(#%g#"Wheelset forward velocity")' ...
     'vy(#%g#"Wheelset lateral velocity")' ...
     'vz(#%g#"Wheelset vertical velocity")' ...
     'vroll(#%g#"Wheelset rate of roll")' ...
     'vyaw(#%g#"Wheelset yaw rate")' ...
     'vpitch(#%g#"Wheelset angular velocity")'...
     'shft_sws(#%g#"Wheelset forward position increment")' ...
     'shft_yws(#%g#"Wheelset lateral position increment")' ...
     'shft_zws(#%g#"Wheelset vertical position increment")' ...
     'shft_rol(#%g#"Wheelset roll angle increment")' ...
     'shft_yaw(#%g#"Wheelset yaw angle increment")' ...
     'shft_pit(#%g#"Wheelset pitch angle icrement")'...
     'rpitch(#%g#"Roller pitch angle increment")'];

    if nargin==0; CAM='';else; CAM=ire ;end
    LI=cntc.call;
    LI.wheelsetVel=cingui('paramedit -doclean2',DoOpt,{struct,CAM});

    switch LI.wheelsetVel.Ewheel
     case {0,1} % Maintain
      params=[];
     case {2:5} % New Velocity
      params=sdth.sfield('addselected',struct,LI.wheelsetVel,{'vs', 'vy', 'vz', ...
       'vrol', 'vyaw', 'vpit'});
     otherwise
      error('Wrong parameter')
      %% xxxGAE
    end
    Global=LI.global; if iscell(Global);Global=sdtm.toStruct(Global(:,1:2));end
    params=struct2cell(params);params=horzcat(params{:});
    cntc.setwheelsetvelocity(Global.ire, LI.wheelsetVel.Ewheel, params);
    return
    % end SDT packaging

   end


   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(ewheel))
    ewheel = 2;
   end
   if (nargin<3 | isempty(params))
    disp('ERROR in cntc.setwheelsetvelocity: invalid params provided.');
    return
   end

   cntc.call('cntc.setwheelsetvelocity', ire, ewheel, length(params), params);

  end % cntc.setwheelsetvelocity

  %------------------------------------------------------------------------------------------------------------
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% #Subs
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [ ] = subs_addblock(ire, icp, iblk, isubs, xparam, yparam, zparam)
   %------------------------------------------------------------------------------------------------------------
   % [ ] = subs_addblock(ire, icp, iblk, isubs, xparam, yparam, zparam)
   %
   % set the parameters for a block of points for the subsurface stress calculation for a contact problem
   %  iblk           - block number; all blocks with this number or higher will be discarded
   %  isubs          - type of block specification
   %  x/y/zparam     - parameters describing x/y/z-coordinates of the block
   %
   %  isubs
   %    1: xparam = [                 ]  yparam = [                 ]  zparam = [ NZ, ZL, DZ   ]
   %    2: xparam = [ IXL, IXINC, IXH ]  yparam = [ IYL, IYINC, IYH ]  zparam = [ NZ, ZL, DZ   ]
   %    3: xparam = [ IX(i), i=1:nx   ]  yparam = [ IY(j), j=1:ny   ]  zparam = [ NZ, ZL, DZ   ]
   %    5: xparam = [                 ]  yparam = [                 ]  zparam = [ Z(k), k=1:nz ]
   %    6: xparam = [ IXL, IXINC, IXH ]  yparam = [ IYL, IYINC, IYH ]  zparam = [ Z(k), k=1:nz ]
   %    7: xparam = [ IX(i), i=1:nx   ]  yparam = [ IY(j), j=1:ny   ]  zparam = [ Z(k), k=1:nz ]
   %    9: xparam = [ X(i),  i=1:nx   ]  yparam = [ Y(j),  j=1:ny   ]  zparam = [ Z(k), k=1:nz ]
   %       ix, iy [-],    x, y, z, zl, dz [length]
   %
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 7: m=*, wtd or cp - default icp=-1
   % #subs_addblock

   % xxxgae WIP
   if nargin==0||ischar(ire)
    % package SDT options
    DoOpt=['BlockSpecif(0#vd{0,Maintain,1,NewVel}#"Wheelset velocity options")' ...
     'iblk(#%g#"Wheelset forward velocity")' ...
     'isubs(#%g#"Wheelset lateral velocity")' ...
     'vz(#%g#"Wheelset vertical velocity")' ...
     'vroll(#%g#"Wheelset rate of roll")' ...
     'vyaw(#%g#"Wheelset yaw rate")' ...
     'vpitch(#%g#"Wheelset angular velocity")'...
     'shft_sws(#%g#"Wheelset forward position increment")' ...
     'shft_yws(#%g#"Wheelset lateral position increment")' ...
     'shft_zws(#%g#"Wheelset vertical position increment")' ...
     'shft_rol(#%g#"Wheelset roll angle increment")' ...
     'shft_yaw(#%g#"Wheelset yaw angle increment")' ...
     'shft_pit(#%g#"Wheelset pitch angle icrement")'...
     'rpitch(#%g#"Roller pitch angle increment")'];

    if nargin==0; CAM='';else; CAM=ire ;end
    LI=cntc.call;
    LI.wheelsetVel=cingui('paramedit -doclean2',DoOpt,{struct,CAM});

    Vel=LI.wheelsetVel.Vel;
    switch Vel
     case 0
      % 0 Maintain
      params=[];
     case 1
      % 1 NewVel xxxgae roller rigs ???
      params=sdth.sfield('addselected',struct,LI.wheelsetVel,{'vx', 'vy', 'vz', ...
       'vrol', 'vyaw', 'vpit'});
     otherwise
      error('Wrong parameter')
      %% xxxGAE
    end
    Global=LI.global; if iscell(Global);Global=sdtm.toStruct(Global(:,1:2));end
    params=struct2cell(params);params=horzcat(params{:});
    cntc.setwheelsetvelocity(Global.ire, Vel, params);
    return
    % end SDT packaging

   end


   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = -1; % default: fall-back for all contact patches
   end
   if (nargin<3 | isempty(iblk))
    iblk = 1;
   end
   if (nargin<7 | isempty(isubs) | isempty(zparam) )
    disp('ERROR in subs_addblock: isubs, zparam are mandatory.');
    return
   end
   if (any(isubs==[1,5]))
    xparam = [ -999 ]; yparam = [ -999 ];
   elseif (isempty(xparam) | isempty(yparam))
    disp(sprintf('ERROR in subs_addblock: xparam, yparam are mandatory for isubs=%d.',isubs));
    return
   end

   npx = length(xparam);
   npy = length(yparam);
   npz = length(zparam);

   cntc.call('subs_addblock', ire, icp, iblk, isubs, npx, npy, npz, xparam, yparam, zparam);

  end % subs_addblock

  %------------------------------------------------------------------------------------------------------------


  function [ ierror ] = subs_calculate(ire, icp, idebug)
   %------------------------------------------------------------------------------------------------------------
   % [ ierror ] = subs_calculate(ire, icp, idebug)
   %
   % perform subsurface stress calculation for a contact problem
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 7: m=*, wtd or cp - default icp=-1
   % #subs_calculate



   cntc_err_allow  = -12;

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = -1; % default: W/R contact, compute for all patches
   end
   if (nargin<3 | isempty(idebug))
    idebug = 1;
   end

   p_ierr = libpointer('int32Ptr',-1);

   cntc.call('subs_calculate', ire, icp, p_ierr);

   ierror = double(p_ierr.value);
   if (idebug>=1 & ierror==cntc_err_allow)
    disp(sprintf('subs_calculate: no valid license found for CONTACT library (%d).',ierror));
   elseif (idebug>=1 & ierror<0)
    disp(sprintf('subs_calculate: an error occurred in the CONTACT library (%d).',ierror));
   end

  end % subs_calculate

  %------------------------------------------------------------------------------------------------------------


  function [ nx, ny, nz ] = subs_getblocksize(ire, icp, iblk)
   %------------------------------------------------------------------------------------------------------------
   % [ nx, ny, nz ] = subs_getblocksize(ire, icp, iblk)
   %
   % get the number of points in a block used for subsurface stress calculation
   %
   %  nx, ny, nz    - number of points used in x-, y- and z-directions
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   %
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 6: m=*, cp     - require icp>0, default 1
   % #subs_getblocksize



   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in subs_getresults: not available for icp=%d',icp));
    return
   end
   if (nargin<3 | isempty(iblk))
    iblk = 1;
   end

   p_nx = libpointer('int32Ptr',-1);
   p_ny = libpointer('int32Ptr',-1);
   p_nz = libpointer('int32Ptr',-1);

   cntc.call('subs_getblocksize', ire, icp, iblk, p_nx, p_ny, p_nz);

   nx = double(p_nx.value);
   ny = double(p_ny.value);
   nz = double(p_nz.value);

  end % subs_getblocksize

  %------------------------------------------------------------------------------------------------------------

  function [ table ] = subs_getresults(ire, icp, iblk, icol)
   %------------------------------------------------------------------------------------------------------------
   % [ table ] = subs_getresults(ire, icp, iblk, icol)
   %   or
   % [ blk ] = subs_getresults(ire, icp, iblk, 'all')
   %
   % return selected columns of the results of the subsurface stress calculation for one block for a
   % contact problem.
   %
   %  iblk              - number of elements in potential contact area
   %  icol(ncol)        - requested columns of result data
   %                       1-- 3: x,y,z        - positions of points for subsurface stress calculation
   %                       4-- 6: ux,uy,uz     - elastic displacements
   %                       7-- 9: sighyd,vm,tr - hydrostatic, von Mises and Tresca stresses
   %                      10--12: sigma1,2,3   - principal stresses
   %                      13--15: sigxx,yx,zx  - components of the full stress tensor
   %                      16--18: sigxy,yy,zy
   %                      19--21: sigxz,yz,zz
   %  table(npnt,ncol)  - result data, with npnt = nx*ny*nz entries filled per column
   %  blk               - structure with subsurface results, as used by loadstrs / plotstrs
   %------------------------------------------------------------------------------------------------------------

   % Copyright 2008-2023 by Vtech CMCC.
   % Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % category 6: m=*, cp     - require icp>0, default 1
   % #subs_getresults



   totcol = 21;

   if (nargin<1 | isempty(ire))
    ire = 1;
   end
   if (nargin<2 | isempty(icp))
    icp = 1;
   end
   if (icp<=0)
    disp(sprintf('ERROR in subs_getresults: not available for icp=%d',icp));
    return
   end

   if (nargin<3 | isempty(iblk))
    iblk = 1;
   end
   [ nx, ny, nz ] = subs_getblocksize(ire, icp, iblk);
   if (any([nx, ny, nz]<=0))
    disp(sprintf('ERROR in subs_getresults: no data for iblk = %d', iblk));
    disp([nx, ny, nz])
    table = [];
    return
   end

   make_struct = 0;
   if (nargin>=4 & ischar(icol))
    make_struct = 1;
    icol = [1:totcol];
   end
   if (nargin<4 | isempty(icol) | any(icol<=0) | any(icol>totcol))
    disp(sprintf('ERROR in subs_getresults: columns must be 1 <= icol <= %d',totcol));
    disp(icol)
    return
   end

   npnt = nx * ny * nz;
   ncol = length(icol);

   p_table = libpointer('doublePtr',zeros(npnt,ncol));

   cntc.call('subs_getresults', ire, icp, iblk, npnt, ncol, icol, p_table);

   table = p_table.value;
   % table = reshape(p_table.value, npnt, ncol);

   if (make_struct)
    blk = struct('nx',nx, 'ny',ny, 'nz',nz);

    % sort data such that x runs fastest, then y, then z
    [~, iperm] = sort(table(:,1)); table = table(iperm,:);
    [~, iperm] = sort(table(:,2)); table = table(iperm,:);
    [~, iperm] = sort(table(:,3)); table = table(iperm,:);

    % the data-lines are given with iz running fastest, ix slowest.
    blk.npoints = nx * ny * nz;
    blk.x       = reshape(table(:, 1), nx, ny, nz); blk.x = squeeze(blk.x(:,1,1));
    blk.y       = reshape(table(:, 2), nx, ny, nz); blk.y = squeeze(blk.y(1,:,1)); blk.y = blk.y';
    blk.z       = reshape(table(:, 3), nx, ny, nz); blk.z = squeeze(blk.z(1,1,:));
    blk.ux      = reshape(table(:, 4), nx, ny, nz);
    blk.uy      = reshape(table(:, 5), nx, ny, nz);
    blk.uz      = reshape(table(:, 6), nx, ny, nz);
    blk.sighyd  = reshape(table(:, 7), nx, ny, nz);
    blk.sigvm   = reshape(table(:, 8), nx, ny, nz);
    blk.sigtr   = reshape(table(:, 9), nx, ny, nz);
    blk.sigma1  = reshape(table(:,10), nx, ny, nz);
    blk.sigma2  = reshape(table(:,11), nx, ny, nz);
    blk.sigma3  = reshape(table(:,12), nx, ny, nz);
    blk.sigxx   = reshape(table(:,13), nx, ny, nz);
    blk.sigxy   = reshape(table(:,16), nx, ny, nz);
    blk.sigyy   = reshape(table(:,17), nx, ny, nz);
    blk.sigxz   = reshape(table(:,19), nx, ny, nz);
    blk.sigyz   = reshape(table(:,20), nx, ny, nz);
    blk.sigzz   = reshape(table(:,21), nx, ny, nz);
    table = blk;

   end % make_struct

  end % subs_getresults

  %------------------------------------------------------------------------------------------------------------


  function out=autoGen(varargin)
   %% #autoGen automated generation from initial files


   if 1==1
    %% attempts at extracting pop from data
    f1='D:\APP\win64\contact_v24.1\doc\user-guide.txt';
    IN=sdtu.f.asChar8(f1);
    IN=regexprep(IN,'Vtech CMCC[^\n]*\n[^\n]*\n[^\d]*[^\n]*\n','');
    sdtu.f.asChar8(strrep(f1,'user-guide','xxx'),IN);

    [st,i1]=sdtm.findTok(IN,'([a-zA-Z]\d*\s+-[^:]{4,20}:[^\x{A}\x{D}]*)');
    ua=struct('ColumnName',{{'Tag','choice','value'}},'table',{{}}, ...
     'ColumnWidth',[80 30 -1],'name','param','setSort',4);
    i1(end+1,1)=i1(end)+1000;
    for j1=1:size(st,1)
     r2='\s*(?<short>[^\s]*)\s+-\s+(?<tag>[^\s:]*)\s+:\s+(?<DisplayName>.*)';
     r2=regexp(st{j1,1},r2,'names');
     r2.DisplayName=regexprep(r2.DisplayName,':\s*$','');
     r3=struct2cell(r2);
     ua.table(end+1,[2 1 3])=r3;ua.level(size(ua.table,1),1:2)=1;
     st2=IN(i1(j1,2):i1(j1+1,1));
     r3='\n(?<choice>\d[^\s]*)\s+-\s+(?<value>[^\x{A}\x{D}]*)';r3=regexp(st2,r3,'names')';
     if isempty(r3);continue;end
     r3=struct2cell(r3)';
     ua.table(end+(1:size(r3,1)),[2 3])=r3;
     ua.level(end+1:size(ua.table,1),1)=2;
    end
    %[~,i1]=sort(ua.table(:,1));ua.table=ua.table(i1,:);
    ua=vhandle.tab(ua);
    f1=fullfile(sdtweb('_onedrive','SNCF*\exchange\'),'CntcParam.mat');
    sdtm.save(f1,'ua');

    asTab(ua)

    assignin('base','ua',ua)


   elseif 1==3
    %% Initial inclusion of functions
    %[~,r1]=sdtu.grep('-s -R','^\s*function','D:\APP\win64\contact_v24.1\matlab');
    [~,r1]=sdtu.grep('-s -R','\nfunction','D:\APP\win64\contact_v24.1\matlab_intfc');
    r1.match
    fname='d:/balmes/xxx2.m';fid=fopen(fname,'wb');
    for j1=1:length(r1.files)
     IN=sdtu.f.asChar8(r1.files{j1});
     fwrite(fid,uint8(IN),'uint8');continue;
     r2.block=sdtu.f.findBlocks(IN,'mathelp');
     r2.fun=regexp(r2.block{1},'(function[^=\(]*)(\([^)]*\)|=[^\(]*)*','tokens');

     r2.fun=[regexprep(r2.fun{1}{end},'=\s*','%% #') newline];
     r2.block=[r2.block((1:2)');{r2.fun};r2.block(3)];
     if 1==2
      st1={'function a',' function out=b','function [a,b]=c(d)','function a(b,c)'
       'function a(b ...\n, c)'}
      regexp(st1,'(function[^=\(]*)(\([^)]*\)|=[^\(]*)*','tokens');sdtm.toString(ans)
     end
     fwrite(fid,uint8(horzcat(r2.block{:})),'uint8');
     %regexp('function a(b,c)','(function[^=\(]*)(\([^)]*\))','tokens');sdtm.toString(ans)
    end
    fclose(fid);edit(fname)

   end

   %% decompress zip



  end % autoGen

 end % Static
end

