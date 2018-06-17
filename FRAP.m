%% Overview and Authorship
%
% This code will implement a reactive-diffusive process in which two
% particles bind together to form a third particle, which may then unbind
% to give the old particle. Furthermore, a FRAP procedure is implemented at
% a pre-defined time. For more information on the model, please see
% FRAP_Note.pdf.
%
% Author: Garrett Higginbotham
% Institute: University of Alabama at Birmingham
% Date: June 7, 2018

%% System Information

% Simulation Parameters
L = 10;                              % Volume of system
T = 4000;                            % Run-time of system
NL = 20;                             % Spatial slicing
NT = 40000;                          % Temporal slicing

% Particle A: Unbleached
NAU = 1000;                          % Initial population of AU in volume
SAU = 0;                             % Initial population of AU on surface
DAU = 1.0;                           % Diffusion constant of AU
konAU = 0.1;                         % Association rate of AU
koffAU = 0.05;                       % Disassociation rate of AU
mAU = 1.0;                           % Mass of AU

% Particle A: Bleached
NAB = 1000;                          % Initial population of AB in volume
SAB = 0;                             % Initial population of AB on surface
DAB = DAU;                           % Diffusion constant of AB
konAB = konAU;                       % Association rate of AB
koffAB = koffAU;                     % Disassociation rate of AB
mAB = mAU;                           % Mass of AB

% Particle B: Unbleached
NBU = 1000;                          % Initial population of BU in volume
SBU = 0;                             % Initial population of BU on surface
DBU = 1.0;                           % Diffusion constant of BU
konBU = 0.1;                         % Association rate of BU
koffBU = 0.05;                       % Disassociation rate of BU
mBU = 1.0;                           % Mass of BU

% Particle B: Bleached
NBB = 1000;                          % Initial population of BB in volume
SBB = 0;                             % Initial population of BB on surface
DBB = DBU;                           % Diffusion constant of BB
konBB = konBU;                       % Association rate of BB
koffBB = koffBU;                     % Disassociation rate of BB
mBB = mBU;                           % Mass of BB

% Particle C: Unbleached
SCU = 0;                             % Initial population of CU
kbindU = 0.0005;                       % Binding rate of CU
kdisU = 0.01;                       % Disassociation rate of CU

% Particle C: Bleached
SCB = 0;                             % Initial population of CB
kbindB = 0.0005;                       % Binding rate of CB
kdisB = 0.01;                       % Disassociation rate of CB

% FRAP Information
frapT = NT / 2;                      % Time at which FRAP occurs
killA = 0.9;                         % Percentage of A bleached
killB = 0.9;                         % Percentage of B bleached
killC = 0.9;                         % Percentage of C bleached

%% Storage Setup

dx = L / NL;                         % Spatial step
dt = T / NT;                         % Time step
x = 0 : dx : L;                      % Spatial vector
t = 0 : dt : T;                      % Temporal vector
nAU = zeros(NT + 1, NL + 1);         % Volume population vector for AU
nAB = zeros(NT + 1, NL + 1);         % Volume population vector for AB
nBU = zeros(NT + 1, NL + 1);         % Volume population vector for BU
nBB = zeros(NT + 1, NL + 1);         % Volume population vector for BB
sigAU = zeros(NT + 1, 1);            % Surface population vector for AU
sigAB = zeros(NT + 1, 1);            % Surface population vector for AB
sigBU = zeros(NT + 1, 1);            % Surface population vector for BU
sigBB = zeros(NT + 1, 1);            % Surface population vector for BB
sigCU = zeros(NT + 1, 1);            % Surface population vector for CU
sigCB = zeros(NT + 1, 1);            % Surface population vector for CB

%% Set the Initial Distribution

% Set the initial volume distribution
muAU = L / 2;
stdAU = L / 2;
muAB = L / 2;
stdAB = L / 2;
muBU = L / 2;
stdBU = L / 2;
muBB = L / 2;
stdBB = L / 2;

for i = 2 : NL
    nAU(1, i) = NAU * exp(-(x(i)-muAU)^2/(2*stdAU^2)) / sqrt(2*pi*stdAU^2);
    nAB(1, i) = NAB * exp(-(x(i)-muAB)^2/(2*stdAB^2)) / sqrt(2*pi*stdAB^2);
    nBU(1, i) = NBU * exp(-(x(i)-muBU)^2/(2*stdBU^2)) / sqrt(2*pi*stdBU^2);
    nBB(1, i) = NBB * exp(-(x(i)-muBB)^2/(2*stdBB^2)) / sqrt(2*pi*stdBB^2);
end

% Set the initial surface distribution
sigAU(1) = SAU;
sigAB(1) = SAB;
sigBU(1) = SBU;
sigBB(1) = SBB;
sigCU(1) = SCU;
sigCB(1) = SCB;

% We must preserve the boundary at x = 0
nAU(1,1) = nAU(1,2);
nAB(1,1) = nAB(1,2);
nBU(1,1) = nBU(1,2);
nBB(1,1) = nBB(1,2);

% We must preserve the boundary at x = L
nAU(1,NL+1) = (DAU/(DAU+konAU*dx))*nAU(1,NL) + ((koffAU*dx)/(DAU+konAU*dx)) * sigAU(1);
nAB(1,NL+1) = (DAB/(DAB+konAB*dx))*nAB(1,NL) + ((koffAB*dx)/(DAB+konAB*dx)) * sigAB(1);
nBU(1,NL+1) = (DBU/(DBU+konBU*dx))*nBU(1,NL) + ((koffBU*dx)/(DBU+konBU*dx)) * sigBU(1);
nBB(1,NL+1) = (DBB/(DBB+konBB*dx))*nBB(1,NL) + ((koffBB*dx)/(DBB+konBB*dx)) * sigBB(1);

%% Calculate the Time Evolution

for j = 1 : frapT - 1
    
    % Compute the evolution of the bulk of the volume
    for i = 2 : NL
        nAU(j+1,i)=nAU(j,i)+((DAU*dt)/(dx^2))*(nAU(j,i+1)-2*nAU(j,i)+nAU(j,i-1));
        nAB(j+1,i)=nAB(j,i)+((DAB*dt)/(dx^2))*(nAB(j,i+1)-2*nAB(j,i)+nAB(j,i-1));
        nBU(j+1,i)=nBU(j,i)+((DBU*dt)/(dx^2))*(nBU(j,i+1)-2*nBU(j,i)+nBU(j,i-1));
        nBB(j+1,i)=nBB(j,i)+((DBB*dt)/(dx^2))*(nBB(j,i+1)-2*nBB(j,i)+nBB(j,i-1));
    end
    
    % Preserve the boundary at x = 0
    nAU(j+1,1) = nAU(j+1,2);
    nAB(j+1,1) = nAB(j+1,2);
    nBU(j+1,1) = nBU(j+1,2);
    nBB(j+1,1) = nBB(j+1,2);
    
    % Compute the evolution of the surface densities
    chunk1 = kbindU*(sigAB(j)*sigBU(j)+sigAU(j)*sigBB(j)+sigAU(j)*sigBU(j));
    chunk2 = kbindB*sigAB(j)*sigBB(j);
    chunk3 = kdisB*sigCB(j)+kdisU*sigCU(j);
    sigCU(j+1) = (1-dt*kdisU)*sigCU(j)+dt*chunk1;
    sigCB(j+1) = (1-dt*kdisB)*sigCB(j)+dt*chunk2;
    sigAU(j+1) = (1-dt*koffAU)*sigAU(j)+konAU*dt*nAU(j,NL+1)-dt*(chunk1+chunk2-chunk3);
    sigAB(j+1) = (1-dt*koffAB)*sigAB(j)+konAB*dt*nAB(j,NL+1)-dt*(chunk1+chunk2-chunk3);
    sigBU(j+1) = (1-dt*koffBU)*sigBU(j)+konBU*dt*nBU(j,NL+1)-dt*(chunk1+chunk2-chunk3);
    sigBB(j+1) = (1-dt*koffBB)*sigBB(j)+konBB*dt*nBB(j,NL+1)-dt*(chunk1+chunk2-chunk3);
    
    % Preserve the boundary at x = L
    nAU(j+1,NL+1) = (DAU/(DAU+konAU*dx))*nAU(j+1,NL) + ((koffAU*dx)/(DAU+konAU*dx)) * sigAU(j+1);
    nAB(j+1,NL+1) = (DAB/(DAB+konAB*dx))*nAB(j+1,NL) + ((koffAB*dx)/(DAB+konAB*dx)) * sigAB(j+1);
    nBU(j+1,NL+1) = (DBU/(DBU+konBU*dx))*nBU(j+1,NL) + ((koffBU*dx)/(DBU+konBU*dx)) * sigBU(j+1);
    nBB(j+1,NL+1) = (DBB/(DBB+konBB*dx))*nBB(j+1,NL) + ((koffBB*dx)/(DBB+konBB*dx)) * sigBB(j+1);
end

%% Execute FRAP

sigAU(frapT) = (1 - killA)*sigAU(frapT);
sigAB(frapT) = (1 + killA)*sigAB(frapT);
sigBU(frapT) = (1 - killB)*sigBU(frapT);
sigBB(frapT) = (1 + killB)*sigBB(frapT);
sigCU(frapT) = (1 - killC)*sigCU(frapT);
sigCB(frapT) = (1 + killC)*sigCB(frapT);

% Do we need to readjust the boundary condition? I don't think so. It
% seems reasonable to expect that FRAP would produce a discontinuity in the
% results. But it's also the case that the lasering process isn't
% instantaneous. I suppose it's possible that the system might have time to
% readjust before the next time step is calculated.

%% Calculate the Remaining Time Evolution

for j = frapT : NT
        
    % Compute the evolution of the bulk of the volume
    for i = 2 : NL
        nAU(j+1,i)=nAU(j,i)+((DAU*dt)/(dx^2))*(nAU(j,i+1)-2*nAU(j,i)+nAU(j,i-1));
        nAB(j+1,i)=nAB(j,i)+((DAB*dt)/(dx^2))*(nAB(j,i+1)-2*nAB(j,i)+nAB(j,i-1));
        nBU(j+1,i)=nBU(j,i)+((DBU*dt)/(dx^2))*(nBU(j,i+1)-2*nBU(j,i)+nBU(j,i-1));
        nBB(j+1,i)=nBB(j,i)+((DBB*dt)/(dx^2))*(nBB(j,i+1)-2*nBB(j,i)+nBB(j,i-1));
    end
    
    % Preserve the boundary at x = 0
    nAU(j+1,1) = nAU(j+1,2);
    nAB(j+1,1) = nAB(j+1,2);
    nBU(j+1,1) = nBU(j+1,2);
    nBB(j+1,1) = nBB(j+1,2);
    
    % Compute the evolution of the surface densities
    chunk1 = kbindU*(sigAB(j)*sigBU(j)+sigAU(j)*sigBB(j)+sigAU(j)*sigBU(j));
    chunk2 = kbindB*sigAB(j)*sigBB(j);
    chunk3 = kdisB*sigCB(j)+kdisU*sigCU(j);
    sigCU(j+1) = (1-dt*kdisU)*sigCU(j)+dt*chunk1;
    sigCB(j+1) = (1-dt*kdisB)*sigCB(j)+dt*chunk2;
    sigAU(j+1) = (1-dt*koffAU)*sigAU(j)+konAU*dt*nAU(j,NL+1)-dt*(chunk1+chunk2-chunk3);
    sigAB(j+1) = (1-dt*koffAB)*sigAB(j)+konAB*dt*nAB(j,NL+1)-dt*(chunk1+chunk2-chunk3);
    sigBU(j+1) = (1-dt*koffBU)*sigBU(j)+konBU*dt*nBU(j,NL+1)-dt*(chunk1+chunk2-chunk3);
    sigBB(j+1) = (1-dt*koffBB)*sigBB(j)+konBB*dt*nBB(j,NL+1)-dt*(chunk1+chunk2-chunk3);
    
    % Preserve the boundary at x = L
    nAU(j+1,NL+1) = (DAU/(DAU+konAU*dx))*nAU(j+1,NL) + ((koffAU*dx)/(DAU+konAU*dx)) * sigAU(j+1);
    nAB(j+1,NL+1) = (DAB/(DAB+konAB*dx))*nAB(j+1,NL) + ((koffAB*dx)/(DAB+konAB*dx)) * sigAB(j+1);
    nBU(j+1,NL+1) = (DBU/(DBU+konBU*dx))*nBU(j+1,NL) + ((koffBU*dx)/(DBU+konBU*dx)) * sigBU(j+1);
    nBB(j+1,NL+1) = (DBB/(DBB+konBB*dx))*nBB(j+1,NL) + ((koffBB*dx)/(DBB+konBB*dx)) * sigBB(j+1);
end

% figure(1);
% hold on;
% plot(t, sigAU);
% 
% figure(2);
% hold on;
% plot(t, sigAB);
% 
% figure(3);
% hold on;
% plot(t, sigBU);
% 
% figure(4);
% hold on;
% plot(t, sigBB);
% 
% figure(5);
% hold on;
% plot(t, sigCU);
% 
% figure(6);
% hold on;
% plot(t, sigCB);
% 
% figure(7);
% hold on;
% plot(t, nAU(:,NL+1));
% 
% figure(8);
% hold on;
% plot(t, nAB(:,NL+1));
% 
% figure(9);
% hold on;
% plot(t, nBU(:,NL+1));
% 
% figure(10);
% hold on;
% plot(t, nBB(:,NL+1));