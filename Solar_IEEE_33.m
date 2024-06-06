clc;
clear all;

% Load the Z_bus and Loads data from .m files
LD = Z_bus();
BD = load_data();

Vbase = 11; % in kV
Sbase = 100; % in MVA
Zbase = (Vbase^2)/Sbase; % in Ohms
LD(:,4:5) = LD(:,4:5) / Zbase; % Convert line impedances to per unit
BD(:,2:3) = BD(:,2:3)/(1000*Sbase); % Convert loads to per unit

% Display the data
disp('Z_bus data:');
disp('From Bus No.  To Bus No.  Resistance  Reactance');
disp(LD);

disp('Loads data:');
disp('Bus No.  Real Load (p.u.)  Reactive Load (p.u.)');
disp(BD);

N = max(max(LD(:,2:3)));
Sload = complex(BD(:,2), BD(:,3));

% Initialize variables
num_buses = 33; % Total number of buses

V = ones(size(BD,1), 1); % Assume initial voltage at each bus is 1 p.u.
Z = complex(LD(:,4), LD(:,5));

Iline = zeros(size(LD,1), 1);
I = 2000;

% Calculate initial losses (without solar panel)
V_initial = V;
Iline_initial = Iline;

for iter = 1:I
    Iload = conj(Sload ./ V_initial);
    for j = size(LD,1):-1:1
        c = [];
        e = [];
        [c, e] = find(LD(:,2:3) == LD(j,3));
        if size(c,1) == 1
            Iline_initial(LD(j,1)) = Iload(LD(j,3));
        else
            Iline_initial(LD(j,1)) = Iload(LD(j,3)) + sum(Iline_initial(LD(c,1))) - Iline_initial(LD(j,1));
        end
    end
    
    % Forward Sweep
    for j = 1:size(LD, 1)
        V_initial(LD(j, 3)) = V_initial(LD(j, 2)) - Iline_initial(LD(j, 1)) * Z(j);
    end
end

% Calculate initial losses
Ploss_initial = zeros(size(LD, 1), 1);
for j = 1:size(LD, 1)
    Ploss_initial(j) = abs(Iline_initial(j))^2 * real(Z(j));
end
total_Ploss_initial = sum(Ploss_initial);

% Solar panel parameters
solar_panel_power = 2.5; % Real power generation of the solar panel in MW
% Convert solar panel power to per unit
solar_panel_power_pu = solar_panel_power / Sbase;

% Variable to store the best placement result
best_bus = 0;
min_loss = total_Ploss_initial;

% Iterate over all buses to find the best placement
for bus = 1:num_buses
    % Reset load data to initial state
    Sload_temp = Sload;

    % Adjust load data to include solar panel injection (only real power)
    Sload_temp(bus) = Sload_temp(bus) - complex(solar_panel_power_pu, 0);
    
    % Calculate losses after adding the solar panel
    V_temp = V;
    Iline_temp = Iline;

    for iter = 1:I
        Iload = conj(Sload_temp ./ V_temp);
        for j = size(LD,1):-1:1
            c = [];
            e = [];
            [c, e] = find(LD(:,2:3) == LD(j,3));
            if size(c,1) == 1
                Iline_temp(LD(j,1)) = Iload(LD(j,3));
            else
                Iline_temp(LD(j,1)) = Iload(LD(j,3)) + sum(Iline_temp(LD(c,1))) - Iline_temp(LD(j,1));
            end
        end
        
        % Forward Sweep
        for j = 1:size(LD, 1)
            V_temp(LD(j, 3)) = V_temp(LD(j, 2)) - Iline_temp(LD(j, 1)) * Z(j);
        end
    end

    % Calculate losses after adding the solar panel
    Ploss_temp = zeros(size(LD, 1), 1);
    for j = 1:size(LD, 1)
        Ploss_temp(j) = abs(Iline_temp(j))^2 * real(Z(j));
    end
    total_Ploss_temp = sum(Ploss_temp);
    
    % Check if the current bus results in a lower loss
    if total_Ploss_temp < min_loss
        min_loss = total_Ploss_temp;
        best_bus = bus;
    end
end

% Calculate percentage reduction in losses at the best bus
loss_reduction_percentage = ((total_Ploss_initial - min_loss) / total_Ploss_initial) * 100;

% Display results
disp(['Best bus for solar panel placement: ', num2str(best_bus)]);
disp(['Total real power loss before solar panel addition (in p.u.): ', num2str(total_Ploss_initial)]);
disp(['Total real power loss after solar panel addition at best bus (in p.u.): ', num2str(min_loss)]);
disp(['Percentage reduction in losses: ', num2str(loss_reduction_percentage), '%']);

% Plot the voltage profile before and after solar panel addition at the best bus
figure;
subplot(2,1,1);
plot(1:num_buses, abs(V_initial), '-o', 'LineWidth', 2);
xlabel('Bus Number');
ylabel('Voltage Magnitude (p.u.)');
title('Voltage Profile of the 33-Bus System (Before Solar Panel)');
grid on;
xlim([1 num_buses]);
xticks(1:num_buses);

% Recalculate voltage profile after solar panel addition at the best bus
Sload(best_bus) = Sload(best_bus) - complex(solar_panel_power_pu, 0);

V_after_solar_best = V;
Iline_after_solar_best = Iline;

for iter = 1:I
    Iload = conj(Sload ./ V_after_solar_best);
    for j = size(LD,1):-1:1
        c = [];
        e = [];
        [c, e] = find(LD(:,2:3) == LD(j,3));
        if size(c,1) == 1
            Iline_after_solar_best(LD(j,1)) = Iload(LD(j,3));
        else
            Iline_after_solar_best(LD(j,1)) = Iload(LD(j,3)) + sum(Iline_after_solar_best(LD(c,1))) - Iline_after_solar_best(LD(j,1));
        end
    end
    
    % Forward Sweep
    for j = 1:size(LD, 1)
        V_after_solar_best(LD(j, 3)) = V_after_solar_best(LD(j, 2)) - Iline_after_solar_best(LD(j, 1)) * Z(j);
    end
end

subplot(2,1,2);
plot(1:num_buses, abs(V_after_solar_best), '-o', 'LineWidth', 2);
xlabel('Bus Number');
ylabel('Voltage Magnitude (p.u.)');
title(['Voltage Profile of the 33-Bus System (After Solar Panel at Bus ', num2str(best_bus), ')']);
grid on;
xlim([1 num_buses]);
xticks(1:num_buses);
