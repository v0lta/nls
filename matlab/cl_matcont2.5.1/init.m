
[list,val] = spparms;
spparms('default');
for i=1:length(list(:,1))
    if strcmp(list(i,:),'umfpack')
        spparms('umfpack',0); %switch umfpack off, use v4solver
    end
end

addpath([cd '/Continuer/']);
addpath([cd '/Equilibrium/']);
addpath([cd '/LimitCycle/']);
addpath([cd '/PeriodDoubling/']);
addpath([cd '/Systems/']);
addpath([cd '/LimitPoint/']);
addpath([cd '/Hopf/']);
addpath([cd '/Help/']);
addpath([cd '/LimitPointCycle/']);
addpath([cd '/NeimarkSacker/']);
addpath([cd '/GUI/']);
addpath([cd '/Testruns/']);
addpath([cd '/BranchPoint/']);
addpath([cd '/BranchPointCycle/']);
addpath([cd '/Homoclinic/']);
addpath([cd '/HomoclinicSaddleNode/']);
addpath([cd '/MultilinearForms/']);

% Find and go to correct directory (class directory)
p = mfilename('fullpath');
p = p(1:length(p)-length(mfilename));
p = strcat(p,'/LimitCycle');
curdir = cd;
cd(p);
 
% Compile the c-files (optimized)
plf=mexext;
if (strcmp(plf,'mexw64') || strcmp(plf,'mexa64') || strcmp(plf,'mexs64'))
  if ~(exist('BVP_LC_jac.mexa64','file') ...
      || exist('BVP_LC_jac.mexw64','file') ...
      || exist('BVP_LC_jac.mexs64','file'))
    mex -largeArrayDims -O BVP_LC_jac64.c    -output BVP_LC_jac; 
    mex -largeArrayDims -O BVP_PD_jac64.c    -output BVP_PD_jac;
    mex -largeArrayDims -O BVP_BPC_jacC64.c  -output BVP_BPC_jacC;
    mex -largeArrayDims -O BVP_BPC_jacCC64.c -output BVP_BPC_jacCC;
    mex -largeArrayDims -O BVP_LPC_jac64.c   -output BVP_LPC_jac;
    mex -largeArrayDims -O BVP_NS_jac64.c    -output BVP_NS_jac;
    mex -largeArrayDims -O BVP_LCX_jac64.c   -output BVP_LCX_jac;
  end
elseif ~(exist('BVP_BPC_jacC.mexsol','file') || exist('BVP_BPC_jacC.mexglx','file')...
        || exist('BVP_BPC_jacC.mexmac','file') || exist('BVP_BPC_jacC.mexw32','file')...
        || exist('BVP_BPC_jacC.mexmaci','file'))
  mex -O BVP_LC_jac.c;
  mex -O BVP_PD_jac.c;
  mex -O BVP_BPC_jacC.c;
  mex -O BVP_BPC_jacCC.c;
  mex -O BVP_LPC_jac.c;
  mex -O BVP_NS_jac.c;
  mex -O BVP_LCX_jac.c;  
end
 
% Return to directory we started in
cd (curdir);
