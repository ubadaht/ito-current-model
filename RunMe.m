clear all; close all; clc;

GtoSlow = 0.0376; % Endo [mS/µF]
FoRT = 3597.59; % 1/mV
% Define constants
Ki = 120; % mM
Ko = 5.4; % mM
R = 8314; % J/(kmol · K)
T = 310; % K
F = 96485; % C/mol
%EK = (1/FoRT)* log(Ko / Ki); % Nernst potential
EK = (R*T)/F*log(Ko / Ki);

% range of Vm values
Vm_values = -40:10:80;


%% Ito_slow from HH based simulink model
% Initialize arrays to store peak Itoslow values
peak_Itos_values = zeros(size(Vm_values));
% Initialize array to store steady-state currents
steady_state_currents = zeros(size(Vm_values));

% Loop through each Vm value
for i = 1:length(Vm_values)
    Vm = Vm_values(i);
    
    % Simulate the model
    sim('Ito.slx');
    
    % Extract time and Itoslow data
    Itos = ans.ItoSlow;
    actv = ans.xtoss;
    inact = ans.ytoss;

    % Find peak Itoslow
    peak_Itos = max(Itos);
    peak_Itos_values(i) = peak_Itos;

    peak_actv(i)=max(actv);
    peak_inact(i)=max(inact);
end

% Plot peak Itoslow against Vm
figure(1);
subplot(2,3,3);
plot(Vm_values, peak_Itos_values, 'bo-');
xlabel('Vm (mV)');
ylabel('(A/F)');
title('IV relationship (Grandi, Pasqualini, and Bers 2010)');
grid on;
ylim([0 18]);

subplot(2,3,1);
plot(Vm_values, peak_actv, 'bo-');
xlim([-40, 80]);  % Limit x-axis
ylim([0, 1]);  % Limit y-axis
title('Ito SS activation');
xlabel('mV');
ylabel('I/Imax');
grid On

subplot(2,3,2);
plot(Vm_values, peak_inact, 'bo-');
xlim([-40, 80]);  % Limit x-axis
ylim([0, 1]);  % Limit y-axis
title('Ito SS Inactivation');
xlabel('mV');
ylabel('I/Imax');
grid On



% Define time span
tspan = [0 100];

% Initial conditions
y0 = [1 0 0 0]; % Initial values of C, O, IC, IO
G = 0.0376;

% Preallocate arrays to store peaks
peaks = zeros(13, 1); % Since we are going from -40 to 80 with an increment of 10, we will have 13 values

%% Ito Markov model simulink
peaks_val = zeros(13, 1); % Since we are going from -40 to 80 with an increment of 10, we will have 13 values
Vm_values = -40:10:80;

% Loop through Vm values
for i = 1:numel(Vm_values)
    % Calculate Nernst potential for current Vm
    Vm = Vm_values(i);
    EK = (R*T)/F*log(Ko / Ki); % Nernst potential
    
    % Simulate the model
    sim('MMSimulink.slx');
    
    % Get the time and current data
    t = ans.tout;
    Itos = ans.Itomm;
    
    % Find and store the peak current
    peaks_val(i) = max(Itos);
end

% Plot the peak current against Vm
subplot(2,3,4);
plot(Vm_values, peaks_val, 'o-');
xlabel('Vm (mV)');
ylabel('Current (A/F)');
title('IV relationship - Markov Model (Simulink)');
grid On;
ylim([0 18]);

%% Markov Model with ODE45
% Loop through Vm values from -40 to 80 with an increment of 10
for i = 1:13
    % Define Vm value
    Vm_mm = -40 + (i - 1) * 10;
    
    % Nernst potential
    EK = (R*T)/F*log(Ko / Ki);
    
    % Solve the system of differential equations
    [t, y] = ode45(@(t, y) myODEs(t, y, Vm_mm), tspan, y0);
    
    % Calculate current
    I = G * y(:,2) .* (Vm_mm - EK);
    
    % Find peak current and store it
    peaks(i) = max(I);
end

% Plot the results
subplot(2,3,5);
plot(Vm_values, peaks, 'o-');
xlabel('Vm (mV)');
ylabel('Current (A/F)');
title('IV relationship (ODE45)');
grid On
ylim([0 18]);

%% MARKOV MODEL FUNCTION
function dydt = myODEs(t, y, Vm)
    xtoss = 1 / (1 + exp(-(Vm - 19.0) / 13));
    ytoss = 1 / (1 + exp((Vm + 19.5) / 5));
    tau_xtos = (9 / (1 + exp((Vm + 3) / 15))) + 0.5;
    tau_ytos = (800 / (1 + exp((Vm + 60) / 10))) + 30;
    
    alpha = xtoss / tau_xtos;
    beta = (1 - xtoss) / tau_xtos; 
    gamma = ytoss / tau_ytos; 
    delta = (1 - ytoss) / tau_ytos; 
    
    C = y(1);
    O = y(2);
    IC = y(3);
    IO = y(4);
    
    dydt = zeros(4, 1);
    dydt(1) = beta * O + delta * IC - (alpha + gamma) * C;
    dydt(2) = alpha * C + delta * IO - (beta + gamma) * O;
    dydt(3) = beta * IO + gamma * C - (alpha + delta) * IC;
    dydt(4) = alpha * IC + gamma * O - (beta + delta) * IO;
end
