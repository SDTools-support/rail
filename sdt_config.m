
if sd('_user','balmes') % BALMES
 pw1=pwd;
 if ~exist('dyn_utils','file');sd _pdynavoie;end
 if ~exist('rail19','file');sd _psncf_ir;end
 addpath(fullfile(pw1,'m'));clear pw1;
 pw0=comgui('cd /o/sdtdata/rail19');
 if ~exist(pw0,'dir');pw0='d:/balmes';end
end
