set(0, 'DefaultFigureWindowStyle', 'Docked')

clc
clear
close all

% ELEC 4707 - Assignment 1
% Brandon Schelhaas
% 101036851

global const

% Add constants to the constants struct
const.m_0 = 9.10938356e-31; % kg
const.m_n = 0.26*const.m_0;
const.k = 1.38064852e-23; 

%------------------------------------------
%                  PART 1
%------------------------------------------

% Thermal velocity at T = 300
T = 300;
v_th = sqrt((2 * const.k * T)/const.m_n);

% Mean free path for tau_mn = 0.2ps
tau_mn = 0.2e-12;
lambda = tau_mn * v_th;

% Define active region's perimeter 
regionLength = 200e-9; % nm
regionWidth = 100e-9; % nm

% Simulation setup
numElectrons = 5000; % Set the amount of electrons to simulate
numToPlot = 10; % Set the amount of electrons to plot 
dt = ((0.01)*regionLength)/v_th; % set a time step to move electrons 1/100 of region length
numSteps = 200;
t = zeros(1, numSteps);

% Setup initial electron positions
pos.x = zeros(numElectrons, 2); % two columns for old and new positions
pos.y = zeros(numElectrons, 2); % two columns for old and new positions
pos.x(:,1) = regionLength .* rand(numElectrons, 1); % fill initial position with random val
pos.y(:,1) = regionWidth .* rand(numElectrons, 1); % fill initial position with random val
eColours = hsv(numToPlot); % get colours for each electron

% Initialize electron velocity
vel.x = zeros(numElectrons, 1);
vel.y = zeros(numElectrons, 1);
angles = 2*pi*rand(numElectrons, 1);
vel.x(:,1) = v_th * cos(angles(:, 1));
vel.y(:,1) = v_th * sin(angles(:, 1));

% Run modelling script
electron_modelling

%------------------------------------------
%                PART 2
%------------------------------------------

% Reset positions and velocities
pos.x = zeros(numElectrons, 2);
pos.y = zeros(numElectrons, 2);
vel.x = zeros(numElectrons, 1);
vel.y = zeros(numElectrons, 1);

% Place some electrons
pos.x(:,1) = regionLength * rand(numElectrons, 1); % fill initial position with random val
pos.y(:,1) = regionWidth * rand(numElectrons, 1); % fill initial position with random val

% Assign each electron a random velocity by Maxwell-Boltzmann distribution
vel.x(:,1) = v_th/sqrt(2) .* randn(numElectrons, 1);
vel.y(:,1) = v_th/sqrt(2) .* randn(numElectrons, 1);
%            v_th/sqrt(2) = sqrt(2kT/m)/sqrt(2) = sqrt(kT/m) = sigma

% Calculate the probability for scattering
p_scat = 1 - exp(-dt/tau_mn);

% Run collisions with MFP script
collisions_with_mfp

%------------------------------------------
%                PART 3
%------------------------------------------

% Reset positions and velocities
pos.x = zeros(numElectrons, 2);
pos.y = zeros(numElectrons, 2);
vel.x = zeros(numElectrons, 1);
vel.y = zeros(numElectrons, 1);

% Create boxes
boxes.box1.coords = [80e-9 60e-9 40e-9 40-9];
boxes.box1.type = 'Specular';
boxes.box2.coords = [80e-9 0e-9 40e-9 40e-9];
boxes.box2.type = 'Diffusive';

% Place electrons
for i = 1:numElectrons
    randx = regionLength * rand();
    randy = regionWidth * rand();
    flag1 = (randx >= 80e-9) && (randx <= 120e-9);
    flag2 = (randy >= 60e-9) || (randy <= 40e-9);
    
    while (flag1 && flag2)
        randx = regionLength * rand();
        randy = regionWidth * rand();
        flag1 = (randx >= 80e-9) && (randx <= 120e-9); % check box1 and box2 x coordinate
        flag2 = (randy >= 60e-9) || (randy <= 40e-9); % check box1 and box2 y coordinate
    end
    
    pos.x(i,1) = randx;
    pos.y(i,1) = randy;
end

% Assign each electron a random velocity by Maxwell-Boltzmann distribution
vel.x(:,1) = v_th/sqrt(2) .* randn(numElectrons, 1);
vel.y(:,1) = v_th/sqrt(2) .* randn(numElectrons, 1);

% Run enhancements script
enhancements



