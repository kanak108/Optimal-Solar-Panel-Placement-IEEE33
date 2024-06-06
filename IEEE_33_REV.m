clc;
clear all;

% Load the Z_bus and Loads data from .m files
LD = Z_bus();
BD = load_data();

Vbase = 11;
Sbase = 100;
Zbase = (Vbase^2)/Sbase ;
LD(:,4:5) = LD(:,4:5  ) / Zbase;
BD(:,2:3) = BD(:,2:3)/(1000*Sbase);
% Display the data
disp('Z_bus data:');
disp('From Bus No.  To Bus No.  Resistance  Reactance');
disp(LD);

disp('Loads data:');
disp('Bus No.  Real Load (kW)  Reactive Load (kVAR)');
disp(BD);

N=max(max(LD(:,2:3)));
Sload = complex(BD(:,2),BD(:,3));

% Initialize variables
num_buses = 33; % Total number of buses

V = ones(size(BD,1), 1); % Assume initial voltage at each bus is 1 p.u.
Z= complex(LD(:,4),LD(:,5))

Iline= zeros(size(LD,1),1);
I=2000;

for i = 1:I
    Iload = conj(Sload./V);
    for j = size(LD,1):-1:1
        c=[];
        e=[];
        [c e]=find(LD(:,2:3)==LD(j,3));
        if size(c,1)==1
            Iline(LD(j,1))==Iload(LD(j,3));
        else

        Iline(LD(j,1))=Iload(LD(j,3))+sum(Iline(LD(c,1)))-Iline(LD(j,1));
        end
    end
        % Forward Sweep
    for j = 1:size(LD, 1)
        V(LD(j, 3)) = V(LD(j, 2)) - Iline(LD(j, 1)) * Z(j);
    end
end

disp('Bus voltages (in per unit):');
disp(V);

% Plot the voltage profile
figure;
subplot(2,1,1);
plot(1:num_buses, abs(V), '-o', 'LineWidth', 2);
xlabel('Bus Number');
ylabel('Voltage Magnitude (p.u.)');
title('Voltage Profile of the 33-Bus System');
grid on;
xlim([1 num_buses]);
xticks(1:num_buses);

% Calculate losses in each line
Ploss = zeros(size(LD, 1), 1);
for j = 1:size(LD, 1)
    Ploss(j) = abs(Iline(j))^2 * real(Z(j));
end

disp('Line losses (in p.u.):');
disp(Ploss);

% Plot the line losses
subplot(2,1,2);
plot(1:size(LD, 1), Ploss, '-o', 'LineWidth', 2);
xlabel('Line Number');
ylabel('Losses (p.u.)');
title('Line Losses of the 33-Bus System');
grid on;
xlim([1 size(LD, 1)]);
xticks(1:size(LD, 1));
%% 
%----solar panel-----
% Initialize variables to store results
power_losses = zeros(num_buses, 1);
voltage_improvement = zeros(num_buses, 1);

% Loop through each bus and simulate the addition of the solar panel
for bus = 1:num_buses
    % Simulate the addition of the solar panel at the current bus
    Sload_with_solar = Sload;
    Sload_with_solar(bus) = Sload_with_solar(bus) + 1; % Add 1MW solar panel
    V_with_solar = V; % Initialize voltage vector
    
    % Perform power flow calculation with the modified load
    for iter = 1:I
        Iload = conj(Sload_with_solar./V_with_solar);
        for j = size(LD, 1):-1:1
            c = [];
            e = [];
            [c, e] = find(LD(:, 2:3) == LD(j, 3));
            if size(c, 1) == 1
                Iline(LD(j, 1)) == Iload(LD(j, 3));
            else
                Iline(LD(j, 1)) = Iload(LD(j, 3)) + sum(Iline(LD(c, 1))) - Iline(LD(j, 1));
            end
        end
        
        % Forward Sweep
        for j = 1:size(LD, 1)
            V_with_solar(LD(j, 3)) = V_with_solar(LD(j, 2)) - Iline(LD(j, 1)) * Z(j);
        end
    end
    
    % Calculate power losses and voltage improvement with the solar panel
    power_losses(bus) = sum(abs(Iline) .^ 2); % Total power losses
    voltage_improvement(bus) = sum(abs(V_with_solar)) - sum(abs(V)); % Total voltage improvement
end

% Find the bus with the maximum voltage improvement or minimum power losses
[max_voltage_improvement, optimal_bus_voltage] = max(voltage_improvement);
[min_power_losses, optimal_bus_power] = min(power_losses);

disp('Optimal Location for Solar Panel Placement:');
disp(['Based on voltage improvement: Bus Number - ', num2str(optimal_bus_voltage)]);
disp(['Based on power losses: Bus Number - ', num2str(optimal_bus_power)]);


%% 
%----adding the solar panel------
% Optimal Location for Solar Panel Placement
optimal_bus = optimal_bus_voltage; % Choose the bus with maximum voltage improvement

% Add the solar panel power to the load data at the optimal bus
BD(optimal_bus, 2) = BD(optimal_bus, 2) + 1000; % Add 1MW in kW

% Re-run the power flow calculation with the updated load data
Sload = complex(BD(:, 2), BD(:, 3)); % Update load data
V = ones(size(BD, 1), 1); % Re-initialize voltage vector

% Perform power flow calculation
for i = 1:I
    Iload = conj(Sload ./ V);
    for j = size(LD, 1):-1:1
        c = [];
        e = [];
        [c, e] = find(LD(:, 2:3) == LD(j, 3));
        if size(c, 1) == 1
            Iline(LD(j, 1)) == Iload(LD(j, 3));
        else
            Iline(LD(j, 1)) = Iload(LD(j, 3)) + sum(Iline(LD(c, 1))) - Iline(LD(j, 1));
        end
    end
    
    % Forward Sweep
    for j = 1:size(LD, 1)
        V(LD(j, 3)) = V(LD(j, 2)) - Iline(LD(j, 1)) * Z(j);
    end
end

% Display the voltage profile of the system with the integrated solar panel
disp('Bus voltages (in per unit) with integrated solar panel:');
disp(V);

% Plot the updated voltage profile
figure;
plot(1:num_buses, abs(V), '-o', 'LineWidth', 2);
xlabel('Bus Number');
ylabel('Voltage Magnitude (p.u.)');
title('Voltage Profile of the 33-Bus System with Integrated Solar Panel');
grid on;
xlim([1 num_buses]);
xticks(1:num_buses);
