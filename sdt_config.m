
if sd('_user','balmes') % BALMES
 pw1=pwd;
 sd _pdynavoie
 addpath(fullfile(pw1,'m'));clear pw1;
 pw0=comgui('cd /o/sdtdata/rail19');
 if ~exist(pw0,'dir');pw0='d:/balmes';end
end
