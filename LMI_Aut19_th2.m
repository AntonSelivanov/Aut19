function feas=LMI_Aut19_th2(D,cf,F,K,l,cb,alpha)
% This MATLAB program checks the feasibility of LMIs from Theorem 2 of the paper 
% A. Selivanov and E. Fridman, "Delayed H-infinity control of 2D diffusion systems under delayed pointlike measurements," Automatica, 2019. 

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)
% and SeDuMi solver (http://sedumi.ie.lehigh.edu/)

% Input: 
% D         - diffusion matrix from (8)
% cf, F     - parameters of nonlinearity from (9)
% K         - controller gain from (15) 
% l         - defined in (22), characterizes sensors
% cb        - defined in (12), characterizes actuators 
% alpha     - decay rate 

% Output: 
% feas =1 if feasible, =0 otherwise
%% Decision variables and notations 
sdpvar p1 p2 p3 mu0 mu1 mu2 mu3 mu4 mu5 mu6 mu7 mu8
P=[p1 p2; p2 p3]; 
pbar=[p1; 2*p2; p3]; 
dbar=[D(1,1); D(1,2)+D(2,1); D(2,2)]; 
%% LMIs 
Phi=blkvar; 
Phi(1,1)=-2*(K-alpha)-(mu7+mu8)*pi^2+mu5*cf+mu4*cb*K^2; 
Phi(1,3)=1; 
Phi(1,4)=1-mu4*cb*K; 
Phi(1,5)=1; 
Phi(2,2)=-pbar*dbar'-dbar*pbar'+[0 0 mu6; 0 mu3*(2*l/pi)^4-2*mu6 0; mu6 0 0]; 
Phi(2,3)=-pbar; 
Phi(2,4)=-pbar; 
Phi(2,5)=-pbar; 
Phi(3,3)=-mu5; 
Phi(4,4)=-mu0/K^2+mu4*cb; 
Phi(5,5)=-mu4; 
Phi=sdpvar(Phi); 

PhiNabla=-2*D-2*(K-alpha)*P+mu5*F+(2*l/pi)^2*diag([mu1,mu2])+diag([mu7,mu8]); 
%% Solution of LMIs 
LMIs=[P>=0, mu7>=0, mu8>=0, diag([mu1,mu2,mu3])>=mu0*ones(3,3), Phi<=0, PhiNabla<=0]; 
options=sdpsettings('solver','sedumi','verbose',0);
sol=optimize(LMIs,[],options); 

feas=0; 
if sol.problem==0
    [primal,~]=check(LMIs); 
    feas=min(primal)>=0; 
else
    yalmiperror(sol.problem) 
end