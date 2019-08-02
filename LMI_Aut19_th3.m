function feas=LMI_Aut19_th3(D,cf,F,K,l,cb,alpha,OmegaM,cinf,tauM,du,gamma)
% This MATLAB program checks the feasibility of LMIs from Theorem 3 of the paper 
% A. Selivanov and E. Fridman, "Delayed H-infinity control of 2D diffusion systems under delayed pointlike measurements," Automatica, 2019. 

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)
% and SeDuMi solver (http://sedumi.ie.lehigh.edu/)

% Input: 
% D         - diffusion matrix from (8)
% cf, F     - parameters of nonlinearity from (9)
% K         - controller gain from (34) 
% l         - defined in (22), characterizes sensors
% cb        - defined in (12), characterizes actuators 
% alpha     - decay rate 
% OmegaM    - maximum subdomain square 
% cinf      - maximum infinity-norm of c_i
% tauM      - delay bound from (37) 
% du, gamma - H-inf parameters from (44) 

% Output: 
% feas =1 if feasible, =0 otherwise
%% Decision variables and notations 
sdpvar p1 p2 p3 mu0 mu1 mu2 mu3 mu4 mu5 mu6 mu7 mu8 gamma1 r
P=[p1 p2; p2 p3]; 
pbar=[p1; 2*p2; p3]; 
dbar=[D(1,1); D(1,2)+D(2,1); D(2,2)]; 
gamma2=gamma1*gamma^2; 
gamma3=gamma1*du^2; 
%% LMIs 
PhiTilde=blkvar; 
PhiTilde(1,1)=-2*(K-alpha)-(mu7+mu8)*pi^2+mu5*cf+mu4*cb*K^2+gamma1+gamma3*K^2; 
PhiTilde(1,3)=1; 
PhiTilde(1,4)=1-mu4*cb*K-gamma3*K; 
PhiTilde(1,5)=1-gamma3*K; 
PhiTilde(1,6)=1-mu4*K*cb-gamma3*K; 
PhiTilde(1,7)=1-gamma3*K; 
PhiTilde(1,8)=1; 
PhiTilde(1,9)=-tauM*r*K; 
PhiTilde(2,2)=-pbar*dbar'-dbar*pbar'+[0 0 mu6; 0 mu3*(2*l/pi)^4-2*mu6 0; mu6 0 0]; 
PhiTilde(2,3)=-pbar; 
PhiTilde(2,4)=-pbar; 
PhiTilde(2,5)=-pbar; 
PhiTilde(2,6)=-pbar; 
PhiTilde(2,7)=-pbar; 
PhiTilde(2,8)=-pbar; 
PhiTilde(2,9)=tauM*r*dbar;
PhiTilde(3,3)=-mu5; 
PhiTilde(3,9)=tauM*r; 
PhiTilde(4,4)=-mu0/K^2+mu4*cb+gamma3; 
PhiTilde(4,5)=gamma3; 
PhiTilde(4,6)=mu4*cb+gamma3; 
PhiTilde(4,7)=gamma3; 
PhiTilde(4,9)=tauM*r; 
PhiTilde(5,5)=-mu4+gamma3; 
PhiTilde(5,6)=gamma3; 
PhiTilde(5,7)=gamma3; 
PhiTilde(5,9)=tauM*r; 
PhiTilde(6,6)=-r+mu4*cb+gamma3; 
PhiTilde(6,7)=gamma3; 
PhiTilde(6,9)=tauM*r; 
PhiTilde(7,7)=-gamma2/K^2+gamma3; 
PhiTilde(7,9)=tauM*r; 
PhiTilde(8,8)=-gamma2; 
PhiTilde(8,9)=tauM*r; 
PhiTilde(9,9)=-r*exp(-2*alpha*tauM)/(K^2*OmegaM*cinf); 
PhiTilde=sdpvar(PhiTilde); 

PhiNabla=-2*D-2*(K-alpha)*P+mu5*F+(2*l/pi)^2*diag([mu1,mu2])+diag([mu7,mu8]); 
%% Solution of LMIs 
LMIs=[P>=0, mu7>=0, mu8>=0, diag([mu1,mu2,mu3])>=mu0*ones(3,3), PhiTilde<=0, PhiNabla<=0]; 
options=sdpsettings('solver','sedumi','verbose',0);
sol=optimize(LMIs,[],options); 

feas=0; 
if sol.problem == 0
    [primal,~]=check(LMIs); 
    feas=min(primal)>=0; 
else
    yalmiperror(sol.problem) 
end