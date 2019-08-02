% This MATLAB program checks the feasibility of LMIs from Theorems 2 and 3 of the paper 
% A. Selivanov and E. Fridman, "Delayed H-infinity control of 2D diffusion systems under delayed pointlike measurements," Automatica, 2019.
%% System parameters 
D=eye(2)/(2*pi^2);                                  % diffusion matrix from (8)
betaU=2; betaT=50; gammaa=4;                        % parameters of the nonlinearity 
cf=(4*betaT/gammaa*exp(-2)-betaU)^2; F=zeros(2);    % nonlinearity bounds from (9)
%% Control without delayes and disturbances 
K=10;       % controller gain from (15) 
alpha=.01;  % decay rate 
l=.0785;    % defined in (22), characterizes sensors
cb=.01;     % defined in (12), characterizes actuators 

if LMI_Aut19_th2(D,cf,F,K,l,cb,alpha) 
    disp('Theorem 2: feasible') 
else
    disp('Theorem 2: not feasible') 
end
%% H-inf control with delays 
K=10;                       % controller gain from (34) 
alpha=.01;                  % decay rate 
N=8^2;                      % number of subdomains 
epsilon=.0125;              % parameter of measurements from (14) 
l=1/(2*sqrt(N))+epsilon/2;  % defined in (22), characterizes sensors
cb=.01;                     % defined in (12), characterizes actuators 
OmegaM=1/N;                 % maximum subdomain square 
cinf=1/epsilon^2;           % maximum infinity-norm of c_i
tauM=.001;                  % delay bound from (37) 
du=.1;                      % | H-inf parameters from (44) 
gamma=100;                  % |

if LMI_Aut19_th3(D,cf,F,K,l,cb,alpha,OmegaM,cinf,tauM,du,gamma)
    disp('Theorem 3: feasible') 
else
    disp('Theorem 3: not feasible') 
end