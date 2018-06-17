%% Overview and Authorship
%
% This code will implement a reactive-diffusive process in which two
% particles form together to create a third particle. The details of the
% model and the numerical technique implemented here is outlined in the
% associated file ReactionDiffusion.pdf.
%
% Author: Garrett Higginbotham
% Date: June 3, 2018
% Institute: University of Alabama at Birmingham

%% System Parameters
L = 10; % Volume of system
T = 1000; % Run-time of simulation
NL = 200; % Spatial slicing
NT = 1000000; % Temporal slicing
% Particle A
NA = 1000; % Initial total population of A in volume
SA = 0; % Initial population of A on surface
DA = 1.0; % Diffusion constant A
konA = 0.1; % Association rate A
koffA = 0.05; % Disassociation rate A
% Particle B
NB = 1000; % Initial total population of B in volume
SB= 0; % Initial population of B on surface
DB = 1.0; % Diffusion constant B
konB = 0.1; % Association rate B
koffB = 0.05; % Disassociation rate B
% Particle C
SC = 0; % Initial population of C surface
kbind = 0.01; % Binding rate of particles A and B
kdis = 0.001; % Disassociation rate of particle C

%% Setup Calculations and Convergence Check
dx = L / NL; % Spatial spacing
dt = T / NT; % Temporal spacing
% My convergence testing isn't working
% % Check for convergence. We will use the convergence criterion
% %    dt << min{dx^2/(2*Dl), 1}
% % Other criteria may be needed to compensate for the reaction processes
% ratA = dx^2 / (2 * DA);
% ratB = dx^2 / (2 * DB);
% if dt > (ratA | ratB | 1)
%     error('Convergence not achieved.');
% end
% Define arrays
x = 0 : dx : L; % Position vector
t = 0 : dt : T; % Time vector
nA = zeros(NT + 1, NL + 1); % Volume population vector for particle A
nB = zeros(NT + 1, NL + 1); % Volume population vector for particle B
sigA = zeros(1, NT + 1); % Boundary population vector for pasrticle A
sigB = zeros(1, NT + 1); % Boundary population vector for pasrticle B
sigC = zeros(1, NT + 1); % Boundary population vector for pasrticle C

%% Set the Initial Distributions
% Few details are given here because the details would be pertinent only to
% the specific initial distribution. The only details are given for the
% surface distribution.
muA = L / 2;
muB = L / 2;
stdA = L / 2;
stdB = L / 2;
for i = 2 : NL
    nA(1, i) = NA * exp(-(x(i)-muA)^2/(2*stdA^2)) / sqrt(2*pi*stdA^2);
    nB(1, i) = NB * exp(-(x(i)-muB)^2/(2*stdB^2)) / sqrt(2*pi*stdB^2);
end
sigA(1) = SA; % Store the initial number of A particles
sigB(1) = SB; % Store the initial number of B particles
sigC(1) = SC; % Store the initial number of C particles
% We must preserve equation 6
nA(1,1) = nA(1,2);
nB(1,1) = nB(1,2);
% We must preserve equation 7
nA(1,NL+1) = (DA/(DA+konA*dx))*nA(1,NL) + ((koffA*dx)/(DA+konA*dx)) * sigA(1);
nB(1,NL+1) = (DB/(DB+konB*dx))*nB(1,NL) + ((koffB*dx)/(DB+konB*dx)) * sigB(1);

%% Calculate the Time Evolution
for j = 1 : NT
    % Use equation 5 to compute the evolution of the bulk of the volume
    for i = 2 : NL
        nA(j+1,i)=nA(j,i)+((DA*dt)/(dx^2))*(nA(j,i+1)-2*nA(j,i)+nA(j,i-1));
        nB(j+1,i)=nB(j,i)+((DB*dt)/(dx^2))*(nB(j,i+1)-2*nB(j,i)+nB(j,i-1));
    end
    % Preserve left boundary using equation 6
    nA(j+1,1) = nA(j+1,2);
    nB(j+1,1) = nB(j+1,2);
    % Use equations 8 and 9 to compute the evolution of the surface
    % densities
    sigA(j+1) = sigA(j) + dt*(konA*nA(j,NL+1)-koffA*sigA(j)-kbind*sigA(j)*sigB(j)+kdis*sigC(j));
    sigB(j+1) = sigB(j) + dt*(konB*nB(j,NL+1)-koffB*sigB(j)-kbind*sigA(j)*sigB(j)+kdis*sigC(j));
    sigC(j+1) = sigC(j) + dt*(kbind*sigA(j)*sigB(j) - kdis*sigC(j));
    % We must preserve equation 7
    nA(j+1,NL+1) = (DA/(DA+konA*dx))*nA(j+1,NL) + ((koffA*dx)/(DA+konA*dx)) * sigA(j+1);
    nB(j+1,NL+1) = (DB/(DB+konB*dx))*nB(j+1,NL) + ((koffB*dx)/(DB+konB*dx)) * sigB(j+1);
end
