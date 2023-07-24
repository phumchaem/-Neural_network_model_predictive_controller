function [sys,States0,str,tsampling,simStateCompliance] = snnmpc(t,States,Inputs,flag,net,ui,yi,u0,ts,w,nc,np,ub,lb)
global nID nI
nID = net.numInputDelays;
nI = net.numInputs;
switch flag  
  case 0
    [sys,States0,str,tsampling,simStateCompliance] = mdlInitializeSizes(net,ui,yi,u0,ts);
  case 1
    sys=mdlDerivatives;
  case 2
    sys=mdlUpdate(States,Inputs,net,w,nc,np,ub,lb);
  case 3
    sys=mdlOutputs(t,States,u0,ts);
  case 4
    sys=mdlGetTimeOfNextVarHit;
  case 9
    sys=mdlTerminate;
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end

function [sys,States0,str,tsampling,simStateCompliance]=mdlInitializeSizes(net,ui,yi,u0,ts)
global nID nI
sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = nID*(nI+1)+2;
sizes.NumOutputs     = 1;
sizes.NumInputs      = nI+2;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;  
sys = simsizes(sizes);

% initialize the initial conditions

Xi1 = reshape(repmat(ui,1,nID)',nID*nI,1);
Xi2 = repmat(yi,1,nID)';
States0 = [Xi1; Xi2; 0; u0];

% str is always an empty matrix
str = [];

tsampling  = [ts 0];
simStateCompliance = 'UnknownSimState';

function sys=mdlDerivatives
sys = [];

function sys=mdlUpdate(States,Inputs,net,w,nc,np,ub,lb)
global nID nI Tin
u = Inputs(1:nI);
y = Inputs(nI+1);
Tsp = Inputs(nI+2);
Xi1 = reshape(States(1:nI*nID),nID,nI);
Xi2 = States((1:nID)+nI*nID);
ierr = States(end-1);

Xs = cell(2,1);
Xs{1} = u; 
Xs{2} = y;
Xi = [num2cell(Xi1',1); num2cell(Xi2')];
Ai = cell(2,0);
[ym,newXi,newAi] = net(Xs,Xi,Ai);
err = Tsp-y;
ierr = err+ierr;
q1 = (1-u(1)/100)*10;
q2 = u(1)/100*10;
Tmix = (ym{1}.*q1+q2*303.15)/(q1+q2);
TspStar = Tmix(1)+err+ierr;
[netc,Xic,Aic] = closeloop(net,newXi,newAi);
lastu = u(1);
u0 = ones(1,nc)*u(1);
opt = optimset('MaxFunEvals',20,'Display','Off');
% uOpt = fminsearch(@objfun,u0,opt,netc,Xic,Aic,Tsp+(ym-y),lastu,w,np,nc);
%{
if abs(err)>0.1
    we = w;
else
    we = abs(err)/0.1*w;
end
%}
uOpt = fmincon(@objfun,u0,[],[],[],[],u0*0+lb,u0*0+ub,[],opt,netc,Xic,Aic,TspStar,lastu,w,np,nc,u(2:end));

newStates1 = reshape([newXi{1,:}]',nI*nID,1);
newStates2 = [newXi{2,:}]';
States = [newStates1; newStates2; err; uOpt(1)];
sys = States;

function sys=mdlOutputs(t,States,u0,ts)
global nID

if t<nID*ts
   
    sys = u0;
elseif States(end)>0
    sys = States(end);
          
else
     sys=0;  

   
    
end

function sys=mdlGetTimeOfNextVarHit
sampleTime = 1;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;

function sys=mdlTerminate
sys = [];

function [J,Tf] = objfun(u,netc,Xic,Aic,Tsp,lastu,w,np,nc,uOther)
Xf = num2cell([u,ones(1,np-nc)*u(end); repmat(uOther,1,np)],1);
Tf = netc(Xf,Xic,Aic);
q1 = (1-u(1)/100)*10;
q2 = u(1)/100*10;
Tmix = ([Tf{:}].*q1+q2*303.15)/(q1+q2);
er = Tmix - Tsp;

du = u-[lastu,u(1:end-1)];
J =er*er'+du*du'*w;