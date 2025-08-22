function [out,out1]=t_rail(varargin)

% t_rail: non regression testing for rail applications
%
% <a href="matlab:d_rail('nmap')">List configurations</a>
% <a href="matlab:d_rail('MeshDb')">Mesh database</a>

%  Etienne Balmes, Guillaume Vermot des Roches
% Copyright SDTools & SNCF I&R, 2019-2025

%#ok<*ASGLU,*NASGU,*NOSEM,*NBRAK>
obj=[];evt=[];
if nargin==0; help t_rail
    return; 
else; [CAM,Cam]=comstr(varargin{1},1);carg=2;
end

if comstr(Cam,'ref');[CAM,Cam]=comstr(CAM,4);
%% 

if comstr(Cam,'mesh')
%% #RefMesh check examples

r2=sdtu.log('CmdDisp{Listen{diag,warn},stack}'); % Turn logViewer on
RT=d_rail('nmap.Trk21ref');nmap=RT.nmap; 
%sdtweb d_rail 'MeshDOE' / sdtweb d_rail nameToMeshRO
% nmap('Ref21')
RT.nmap('CurExp')={'MeshCfg{d_rail(Ref21):21ref}','RunCfg{feplot}'};
sdtm.range(RT)
fecom showFiPro

 railu.rail_section('60-E1')

matgui('jpl.asMd')

elseif comstr(Cam,'period')
%% #RefPeriod : reactivate periodic computations on one sleeper

RT=rail19('nmap.Trk23Period');nmap=RT.nmap; 
sdtm.range(RT); mo2=RT.nmap('CurModel');SE=mo2;
fe_case(mo2,'info')

RD=fe_range('BuildGrid',struct('ncx',linspace(100,2.01,30)));
RD.freq=logspace(3,log10(3e5),100);mo2=stack_set(mo2,'info','EigOpt',[5 40 0]);
[def,hist]=fe_homo('dftDisp',mo2,RD);
hist.Ylab={'Frequency','Hz'};
cdm.urnVec(hist,'{x,1}{y,:}{gf30}');set(gca,'xscale','log')
def.Range.param.ak.LabFcn='';def.Range.param.kx.LabFcn='';
def.Range.param.ncx.LabFcn='sprintf('' %i@lx=%.2f m'',nnz(def.data(1:ch,2)==def.data(ch,2)),val*.6)';
[~,i1]=sortrows(def.data,[2 1]); d2=fe_def('subdef',fe_homo('dftDefLong',def),i1);

RD=struct('mno',(0:4)','cf',102,'ci',112,'ViewHist',hist);
%fe_homo('dfpViewHist',RD,RD.ViewHist);

fe_homo('dfpInitSelDef',SE,d2,RD);fecom colordataevala-edgealpha.1;
fecom('animcolorinstant');
a=comgui('iminfo',102);%

else;
        sdtweb('_link','sdtweb(''t_exp19'',''ref'')');
end


elseif comstr(Cam,'cvs'); out=sdtcheck('revision'); 
elseif comstr(Cam,'@');out=eval(CAM);
elseif ~isempty(CAM); error('%s',CAM);
end
end



