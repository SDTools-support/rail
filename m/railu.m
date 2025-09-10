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


function   m_rail=rail_section(RO);
%% #dv_rail_section : rail section mesh -3

RA=struct;
if ischar(RO); 
  nmap=d_rail('nmap');
  if isKey(nmap,RO);RA=nmap(RO);
  else; nmap=nmap('Map:Sections');
   if isKey(nmap,RO);RA=nmap(RO);end
  end
end
if isfield(RO,'Arm')
     RA=RO.Arm;
end
if isfield(RA,'lar1')
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

    % Node with mass element on the middle-top of the rail section  
    [NodeID,nodes]=feutil(sprintf('findnode y==%.15g & z==%.15g',...
        0.5*(min(m_rail.Node(:,6))+max(m_rail.Node(:,6))),...
        max(m_rail.Node(:,7))),m_rail);
    
    if isempty(nodes)
      sdtw('_err','there must be a node on the middle-top of the rail section');
    end


end

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


end % Static
end
