%% Overview and Authorship
%
% Some bullshit I'll put later. Please see blah blah blah
%
% Author: Garrett Higginbotham
% Institute: University of Alabama at Birmingham
% Date: June 14, 2018

%% System Information

% Simulation Parameters
L = 10;                              % Volume of system
T = 4000;                            % Run-time of system
NL = 20;                             % Spatial slicing
NT = 40000;                          % Temporal slicing

% Particle A: Unbleached (Au)
NAu = 1000;                          % Initial population of Au in volume
SAu = 0;                             % Initial population of Au on surface
DAu = 1.0;                           % Diffusion constant of Au
konAu = 0.1;                         % Association rate of Au
koffAu = 0.05;                       % Disassociation rate of Au
mAu = 1.0;                           % Mass of Au

% Particle A: Bleached (Ab)
NAb = 0;                             % Initial population of Ab in volume
SAb = 0;                             % Initial population of Ab on surface
DAb = DAu;                           % Diffusion constant of Ab
konAb = konAu;                       % Association rate of Ab
koffAb = koffAu;                     % Disassociation rate of Ab
mAb = mAu;                           % Mass of Ab

% Particle B: Unbleached (Bu)
NBu = 1000;                          % Initial population of Bu in volume
SBu = 0;                             % Initial population of Bu on surface
DBu = 1.0;                           % Diffusion constant of Bu
konBu = 0.1;                         % Association rate of Bu
koffBu = 0.05;                       % Disassociation rate of Bu
mBu = 1.0;                           % Mass of Bu

% Particle B: Bleached (Bb)
NBb = 0;                             % Initial population of Bb in volume
SBb = 0;                             % Initial population of Bb on surface
DBb = DBu;                           % Diffusion constant of Bb
konBb = konBu;                       % Association rate of Bb
koffBb = koffBu;                     % Disassociation rate of Bb
mBb = mBu;                           % Mass of Bb

% Particle Cuu: Au + Bu -> Cuu
NCuu = 0;                            % Initial Population of Cuu in volume
SCuu = 0;                            % Initial Population of Cuu on surface
DCuu = 1.0;                          % Diffusion Constant of Cuu
konCuu = 0;                          % Association rate of Cuu
koffCuu = 0;                         % Disassociation rate of Cuu
mCuu = mAu + mBu;                    % Mass of Cuu
kbinduu = 0.01;                      % Binding rate of Cuu
kdisuu = 0.005;                      % Dissipation rate of Cuu

% Particle Cub: Au + Bb -> Cub
NCub = 0;                            % Initial Population of Cub in volume
SCub = 0;                            % Initial Population of Cub on surface
DCub = DCuu;                         % Diffusion Constant of Cub
konCub = 0;                          % Association rate of Cub
koffCub = 0;                         % Disassociation rate of Cub
mCub = mAu + mBb;                    % Mass of Cub
kbindub = kbinduu;                   % Binding rate of Cub
kdisub = kdisuu;                     % Dissipation rate of Cub

% Particle Cbu: Ab + Bu -> Cbu
NCbu = 0;                            % Initial population of Cbu in volume
SCbu = 0;                            % Initial population of Cbu on surface
DCbu = DCuu;                         % Diffusion constant of Cbu
konCbu = 0;                          % Association rate of Cub
koffCbu = 0;                         % Disassociation rate of Cub
mCbu = mAb + mBu;                    % Mass of Cub
kbindbu = kbinduu;                   % Binding rate of Cub
kdisbu = kdisuu;                     % Dissipation rate of Cub

% Particle Cbb: Ab + Bb -> Cbb
NCbb = 0;                            % Initial population of Cbb in volume
SCbb = 0;                            % Initial population of Cbb on surface
DCbb = DCuu;                         % Diffusion constant of Cbb
konCbb = 0;                          % Association rate of Cbb
koffCbb = 0;                         % Disassociation rate of Cbb
mCbb = mAb + mBb;                    % Mass of Cbb
kbindbb = kbinduu;                   % Binding rate of Cbb
kdisbb = kdisuu;                     % Dissipation rate of Cbb

% FRAP Information
frapT = NT / 2;                      % Time at which FRAP occurs
killA = 0.9;                         % Percentage of A bleached
killB = 0.9;                         % Percentage of B bleached
killC = 0.9;                         % Percentage of C bleached

%% storage Setup

dx = L / NL;                         % Spatial step
dt = T / NT;                         % Temporal step
x = 0 : dx : L;                      % Spatial vector
t = 0 : dt : T;                      % Temporal vector
nAu = zeros(NT + 1, NL + 1);         % Volume population vector for Au
nAb = zeros(NT + 1, NL + 1);         % Volume population vector for Ab
nBu = zeros(NT + 1, NL + 1);         % Volume population vector for Bu
nBb = zeros(NT + 1, NL + 1);         % Volume population vector for Bb
nCuu = zeros(NT + 1, NL + 1);        % Volume population vector for Cuu
nCub = zeros(NT + 1, NL + 1);        % Volume population vector for Cub
nCbu = zeros(NT + 1, NL + 1);        % Volume population vector for Cbu
nCbb = zeros(NT + 1, NL + 1);        % Volume population vector for Cbb
sAu = zeros(NT + 1, 1);              % Surface population vector for Au
sAb = zeros(NT + 1, 1);              % Surface population vector for Ab
sBu = zeros(NT + 1, 1);              % Surface population vector for Bu
sBb = zeros(NT + 1, 1);              % Surface population vector for Bb
sCuu = zeros(NT + 1, 1);             % Surface population vector for Cuu
sCub = zeros(NT + 1, 1);             % Surface population vector for Cub
sCbu = zeros(NT + 1, 1);             % Surface population vector for Cbu
sCbb = zeros(NT + 1, 1);             % Surface population vector for Cbb
mass = zeros(NT + 1, 1);             % Total system mass

%% Set the Initial Distributions

% Set the initial volume distributions
muAu = L / 2;
stdAu = L / 2;
muAb = L / 2;
stdAb = L / 2;
muBu = L / 2;
stdBu = L / 2;
muBb = L / 2;
stdBb = L / 2;
muCuu = L / 2;
stdCuu = L / 2;
muCub = L / 2;
stdCub = L / 2;
muCbu = L / 2;
stdCbu = L / 2;
muCbb = L / 2;
stdCbb = L / 2;

for i = 2 : NL
    nAu(1, i) = NAu * exp(-(x(i)-muAu)^2/(2*stdAu^2)) / sqrt(2*pi*stdAu^2);
    nAb(1, i) = NAb * exp(-(x(i)-muAb)^2/(2*stdAb^2)) / sqrt(2*pi*stdAb^2);
    nBu(1, i) = NBu * exp(-(x(i)-muBu)^2/(2*stdBu^2)) / sqrt(2*pi*stdBu^2);
    nBb(1, i) = NBb * exp(-(x(i)-muBb)^2/(2*stdBb^2)) / sqrt(2*pi*stdBb^2);
    nCuu(1, i) = NCuu * exp(-(x(i)-muCuu)^2/(2*stdCuu^2)) / sqrt(2*pi*stdCuu^2);
    nCub(1, i) = NCub * exp(-(x(i)-muCub)^2/(2*stdCub^2)) / sqrt(2*pi*stdCub^2);
    nCbu(1, i) = NCbu * exp(-(x(i)-muCbu)^2/(2*stdCbu^2)) / sqrt(2*pi*stdCbu^2);
    nCbb(1, i) = NCbb * exp(-(x(i)-muCbb)^2/(2*stdCbb^2)) / sqrt(2*pi*stdCbb^2);
end

% for i = 2 : NL
%     nAu(1, i) = NAu / NL;
%     nAb(1, i) = NAb / NL;
%     nBu(1, i) = NBu / NL;
%     nBb(1, i) = NBb / NL;
%     nCuu(1, i) = NCuu / NL;
%     nCub(1, i) = NCub / NL;
%     nCbu(1, i) = NCbu / NL;
%     nCbb(1, i) = NCbb / NL;
% end

% Preserve the boundary at x = 0
nAu(1,1) = nAu(1,2);
nAb(1,1) = nAb(1,2);
nBu(1,1) = nBu(1,2);
nBb(1,1) = nBb(1,2);
nCuu(1,1) = nCuu(1,2);
nCub(1,1) = nCub(1,2);
nCbu(1,1) = nCbu(1,2);
nCbb(1,1) = nCbb(1,2);

% Set the initial surface distribution
sAu(1) = SAu;
sAb(1) = SAb;
sBu(1) = SBu;
sBb(1) = SBb;
sCuu(1) = SCuu;
sCub(1) = SCub;
sCbu(1) = SCbu;
sCbb(1) = SCbb;

% We must preserve the boundary at x = L
nAu(1,NL+1) = (DAu/(DAu+konAu*dx))*nAu(1,NL) + ((koffAu*dx)/(DAu+konAu*dx)) * sAu(1);
nAb(1,NL+1) = (DAb/(DAb+konAb*dx))*nAb(1,NL) + ((koffAb*dx)/(DAb+konAb*dx)) * sAb(1);
nBu(1,NL+1) = (DBu/(DBu+konBu*dx))*nBu(1,NL) + ((koffBu*dx)/(DBu+konBu*dx)) * sBu(1);
nBb(1,NL+1) = (DBb/(DBb+konBb*dx))*nBb(1,NL) + ((koffBb*dx)/(DBb+konBb*dx)) * sBb(1);
nCuu(1,NL+1) = (DCuu/(DCuu+konCuu*dx))*nCuu(1,NL) + ((koffCuu*dx)/(DCuu+konCuu*dx)) * sCuu(1);
nCub(1,NL+1) = (DCub/(DCub+konCub*dx))*nCub(1,NL) + ((koffCub*dx)/(DCub+konCub*dx)) * sCub(1);
nCbu(1,NL+1) = (DCbu/(DCbu+konCbu*dx))*nCbu(1,NL) + ((koffCbu*dx)/(DCbu+konCbu*dx)) * sCbu(1);
nCbb(1,NL+1) = (DCbb/(DCbb+konCbb*dx))*nCbb(1,NL) + ((koffCbb*dx)/(DCbb+konCbb*dx)) * sCbb(1);

% Compute the mass of the system
mr=mAu*(sum(nAu(1,:))+sAu(1))+mAb*(sum(nAb(1,:))+sAb(1))+mBu*(sum(nBu(1,:))+sBu(1))+mBb*(sum(nBb(1,:))+sBb(1));
mp=mCuu*(sum(nCuu(1,:))+sCuu(1))+mCub*(sum(nCub(1,:))+sCub(1))+mCbu*(sum(nCbu(1,:))+sCbu(1))+mCbb*(sum(nCbb(1,:))+sCbb(1));
mass(1) = mr + mp;

%% Calculate the Time Evolution Until FRAP

for j = 1 : frapT - 1
    % Compute the evolution of the bulk of the volume
    for i = 2 : NL
        nAu(j+1,i)=nAu(j,i)+((DAu*dt)/(dx^2))*(nAu(j,i+1)-2*nAu(j,i)+nAu(j,i-1));
        nAb(j+1,i)=nAb(j,i)+((DAb*dt)/(dx^2))*(nAb(j,i+1)-2*nAb(j,i)+nAb(j,i-1));
        nBu(j+1,i)=nBu(j,i)+((DBu*dt)/(dx^2))*(nBu(j,i+1)-2*nBu(j,i)+nBu(j,i-1));
        nBb(j+1,i)=nBb(j,i)+((DBb*dt)/(dx^2))*(nBb(j,i+1)-2*nBb(j,i)+nBb(j,i-1));
        nCuu(j+1,i)=nCuu(j,i)+((DCuu*dt)/(dx^2))*(nCuu(j,i+1)-2*nCuu(j,i)+nCuu(j,i-1));
        nCub(j+1,i)=nCub(j,i)+((DCub*dt)/(dx^2))*(nCub(j,i+1)-2*nCub(j,i)+nCub(j,i-1));
        nCbu(j+1,i)=nCbu(j,i)+((DCbu*dt)/(dx^2))*(nCbu(j,i+1)-2*nCbu(j,i)+nCbu(j,i-1));
        nCbb(j+1,i)=nCbb(j,i)+((DCbb*dt)/(dx^2))*(nCbb(j,i+1)-2*nCbb(j,i)+nCbb(j,i-1));
    end
    
    % Preserve the boundary at x = 0
    nAu(j+1,1) = nAu(j+1,2);
    nAb(j+1,1) = nAb(j+1,2);
    nBu(j+1,1) = nBu(j+1,2);
    nBb(j+1,1) = nBb(j+1,2);
    nCuu(j+1,1) = nCuu(j+1,2);
    nCub(j+1,1) = nCub(j+1,2);
    nCbu(j+1,1) = nCbu(j+1,2);
    nCbb(j+1,1) = nCbb(j+1,2);
    
    % Compute the evolution of the surface densities
    sAu(j+1)=sAu(j)+dt*(konAu*nAu(j,NL+1)-(koffAu+kbinduu*sBu(j)+kbindub*sBb(j))*sAu(j)+kdisuu*sCuu(j)+kdisub*sCub(j));
    sAb(j+1)=sAb(j)+dt*(konAb*nAb(j,NL+1)-(koffAb+kbindbu*sBu(j)+kbindbb*sBb(j))*sAb(j)+kdisbu*sCbu(j)+kdisbb*sCbb(j));
    sBu(j+1)=sBu(j)+dt*(konBu*nBu(j,NL+1)-(koffBu+kbinduu*sAu(j)+kbindbu*sAb(j))*sBu(j)+kdisuu*sCuu(j)+kdisbu*sCbu(j));
    sBb(j+1)=sBb(j)+dt*(konBb*nBb(j,NL+1)-(koffBb+kbindub*sAu(j)+kbindbb*sAb(j))*sBb(j)+kdisub*sCub(j)+kdisbb*sCbb(j));
    sCuu(j+1)=sCuu(j)+dt*(konCuu*nCuu(j,NL+1)-(koffCuu+kdisuu)*sCuu(j)+kbinduu*sAu(j)*sBu(j));
    sCub(j+1)=sCub(j)+dt*(konCub*nCub(j,NL+1)-(koffCub+kdisub)*sCub(j)+kbindub*sAu(j)*sBb(j));
    sCbu(j+1)=sCbu(j)+dt*(konCbu*nCbu(j,NL+1)-(koffCbu+kdisbu)*sCbu(j)+kbindbu*sAb(j)*sBu(j));
    sCbb(j+1)=sCbb(j)+dt*(konCbb*nCbb(j,NL+1)-(koffCbb+kdisbb)*sCbb(j)+kbindbb*sAb(j)*sBb(j));
 
    % Preserve the boundary at x = L
    nAu(j+1,NL+1) = (DAu*nAu(j+1,NL) + koffAu*dx*sAu(j+1))/(DAu+konAu*dx);
    nAb(j+1,NL+1) = (DAb*nAb(j+1,NL) + koffAb*dx*sAb(j+1))/(DAb+konAb*dx);
    nBu(j+1,NL+1) = (DBu*nBu(j+1,NL) + koffBu*dx*sBu(j+1))/(DBu+konBu*dx);
    nBb(j+1,NL+1) = (DBb*nBb(j+1,NL) + koffBb*dx*sBb(j+1))/(DBb+konBb*dx);
    nCuu(j+1,NL+1) = (DCuu*nCuu(j+1,NL) + koffCuu*dx*sCuu(j+1))/(DCuu+konCuu*dx);
    nCub(j+1,NL+1) = (DCub*nCub(j+1,NL) + koffCub*dx*sCub(j+1))/(DCub+konCub*dx);
    nCbu(j+1,NL+1) = (DCbu*nCbu(j+1,NL) + koffCbu*dx*sCbu(j+1))/(DCbu+konCbu*dx);
    nCbb(j+1,NL+1) = (DCbb*nCbb(j+1,NL) + koffCbb*dx*sCbb(j+1))/(DCbb+konCbb*dx);
    
    % Compute the mass of the system
    mr=mAu*(sum(nAu(j,:))+sAu(j))+mAb*(sum(nAb(j,:))+sAb(j))+mBu*(sum(nBu(j,:))+sBu(j))+mBb*(sum(nBb(j,:))+sBb(j));
    mp=mCuu*(sum(nCuu(j,:))+sCuu(j))+mCub*(sum(nCub(j,:))+sCub(j))+mCbu*(sum(nCbu(j,:))+sCbu(j))+mCbb*(sum(nCbb(j,:))+sCbb(j));
    mass(j) = mr + mp;
end

%% Execute FRAP

killedA = killA * sAu(frapT);
killedB = killB * sBu(frapT);
killedCuu = killC * sCuu(frapT);
killedCub = killC * sCub(frapT);
killedCbu = killC * sCbu(frapT);

sAu(frapT) = sAu(frapT) - killedA;
sAb(frapT) = sAb(frapT) + killedA;
sBu(frapT) = sBu(frapT) - killedB;
sBb(frapT) = sBb(frapT) + killedB;
sCuu(frapT) = sCuu(frapT) - killedCuu;
sCub(frapT) = sCub(frapT) - killedCub;
sCbu(frapT) = sCbu(frapT) - killedCbu;
sCbb(frapT) = sCbb(frapT) + killedCuu + killedCub + killedCbu;

% Compute the mass of the system
mr=mAu*(sum(nAu(frapT,:))+sAu(frapT))+mAb*(sum(nAb(frapT,:))+sAb(frapT))+mBu*(sum(nBu(frapT,:))+sBu(frapT))+mBb*(sum(nBb(frapT,:))+sBb(frapT));
mp=mCuu*(sum(nCuu(frapT,:))+sCuu(frapT))+mCub*(sum(nCub(frapT,:))+sCub(frapT))+mCbu*(sum(nCbu(frapT,:))+sCbu(frapT))+mCbb*(sum(nCbb(frapT,:))+sCbb(frapT));
mass(frapT) = mr + mp;

%% Calculate the Remaining Time Evolution

for j = frapT : NT
    % Compute the evolution of the bulk of the volume
    for i = 2 : NL
        nAu(j+1,i)=nAu(j,i)+((DAu*dt)/(dx^2))*(nAu(j,i+1)-2*nAu(j,i)+nAu(j,i-1));
        nAb(j+1,i)=nAb(j,i)+((DAb*dt)/(dx^2))*(nAb(j,i+1)-2*nAb(j,i)+nAb(j,i-1));
        nBu(j+1,i)=nBu(j,i)+((DBu*dt)/(dx^2))*(nBu(j,i+1)-2*nBu(j,i)+nBu(j,i-1));
        nBb(j+1,i)=nBb(j,i)+((DBb*dt)/(dx^2))*(nBb(j,i+1)-2*nBb(j,i)+nBb(j,i-1));
        nCuu(j+1,i)=nCuu(j,i)+((DCuu*dt)/(dx^2))*(nCuu(j,i+1)-2*nCuu(j,i)+nCuu(j,i-1));
        nCub(j+1,i)=nCub(j,i)+((DCub*dt)/(dx^2))*(nCub(j,i+1)-2*nCub(j,i)+nCub(j,i-1));
        nCbu(j+1,i)=nCbu(j,i)+((DCbu*dt)/(dx^2))*(nCbu(j,i+1)-2*nCbu(j,i)+nCbu(j,i-1));
        nCbb(j+1,i)=nCbb(j,i)+((DCbb*dt)/(dx^2))*(nCbb(j,i+1)-2*nCbb(j,i)+nCbb(j,i-1));
    end
    
    % Preserve the boundary at x = 0
    nAu(j+1,1) = nAu(j+1,2);
    nAb(j+1,1) = nAb(j+1,2);
    nBu(j+1,1) = nBu(j+1,2);
    nBb(j+1,1) = nBb(j+1,2);
    nCuu(j+1,1) = nCuu(j+1,2);
    nCub(j+1,1) = nCub(j+1,2);
    nCbu(j+1,1) = nCbu(j+1,2);
    nCbb(j+1,1) = nCbb(j+1,2);
    
    % Compute the evolution of the surface densities
    sAu(j+1)=sAu(j)+dt*(konAu*nAu(j,NL+1)-(koffAu+kbinduu*sBu(j)+kbindub*sBb(j))*sAu(j)+kdisuu*sCuu(j)+kdisub*sCub(j));
    sAb(j+1)=sAb(j)+dt*(konAb*nAb(j,NL+1)-(koffAb+kbindbu*sBu(j)+kbindbb*sBb(j))*sAb(j)+kdisbu*sCbu(j)+kdisbb*sCbb(j));
    sBu(j+1)=sBu(j)+dt*(konBu*nBu(j,NL+1)-(koffBu+kbinduu*sAu(j)+kbindbu*sAb(j))*sBu(j)+kdisuu*sCuu(j)+kdisbu*sCbu(j));
    sBb(j+1)=sBb(j)+dt*(konBb*nBb(j,NL+1)-(koffBb+kbindub*sAu(j)+kbindbb*sAb(j))*sBb(j)+kdisub*sCub(j)+kdisbb*sCbb(j));
    sCuu(j+1)=sCuu(j)+dt*(konCuu*nCuu(j,NL+1)-(koffCuu+kdisuu)*sCuu(j)+kbinduu*sAu(j)*sBu(j));
    sCub(j+1)=sCub(j)+dt*(konCub*nCub(j,NL+1)-(koffCub+kdisub)*sCub(j)+kbindub*sAu(j)*sBb(j));
    sCbu(j+1)=sCbu(j)+dt*(konCbu*nCbu(j,NL+1)-(koffCbu+kdisbu)*sCbu(j)+kbindbu*sAb(j)*sBu(j));
    sCbb(j+1)=sCbb(j)+dt*(konCbb*nCbb(j,NL+1)-(koffCbb+kdisbb)*sCbb(j)+kbindbb*sAb(j)*sBb(j));

    % Preserve the boundary at x = L
    nAu(j+1,NL+1) = (DAu*nAu(j+1,NL) + koffAu*dx*sAu(j+1))/(DAu+konAu*dx);
    nAb(j+1,NL+1) = (DAb*nAb(j+1,NL) + koffAb*dx*sAb(j+1))/(DAb+konAb*dx);
    nBu(j+1,NL+1) = (DBu*nBu(j+1,NL) + koffBu*dx*sBu(j+1))/(DBu+konBu*dx);
    nBb(j+1,NL+1) = (DBb*nBb(j+1,NL) + koffBb*dx*sBb(j+1))/(DBb+konBb*dx);
    nCuu(j+1,NL+1) = (DCuu*nCuu(j+1,NL) + koffCuu*dx*sCuu(j+1))/(DCuu+konCuu*dx);
    nCub(j+1,NL+1) = (DCub*nCub(j+1,NL) + koffCub*dx*sCub(j+1))/(DCub+konCub*dx);
    nCbu(j+1,NL+1) = (DCbu*nCbu(j+1,NL) + koffCbu*dx*sCbu(j+1))/(DCbu+konCbu*dx);
    nCbb(j+1,NL+1) = (DCbb*nCbb(j+1,NL) + koffCbb*dx*sCbb(j+1))/(DCbb+konCbb*dx);

    % Compute the mass of the system
    mr=mAu*(sum(nAu(j,:))+sAu(j))+mAb*(sum(nAb(j,:))+sAb(j))+mBu*(sum(nBu(j,:))+sBu(j))+mBb*(sum(nBb(j,:))+sBb(j));
    mp=mCuu*(sum(nCuu(j,:))+sCuu(j))+mCub*(sum(nCub(j,:))+sCub(j))+mCbu*(sum(nCbu(j,:))+sCbu(j))+mCbb*(sum(nCbb(j,:))+sCbb(j));
    mass(j) = mr + mp;
end