% Title: Ballisitic Trajectory Calculator
% Author: James Emerson Parkus
% Date: June 23rd, 2017
% Purpose: Take initial conditions of a projectile motion problem and
% analyze the situation extensively and provide a plot of motion, fit a
% quadratic curve, and present the data.
%
% Initial Conditions
% 1. velocity vector
%   - Vector magnitude and launch angle
% 2. object mass
%
% Assumptions
% 1. negligible air resistance
% 2. ground is flat
% 3. sunny weather
% 4. no wind

clc
clear

%% Constants
gravity = -1*9.8066; % m/s^2 - Standard gravity at sea level
drag_coefficient = 2.0;
rho_air = 1.225; % density of air

%% Units
distance_unit = '[m]';
time_unit = '[s]';
velocity_unit = '[m/s]';
acceleration_unit = '[m/s^2]';
force_unit = '[N]';
energy_unit = '[J]';
mass_unit = '[kg]';
angle_unit = '[deg]';
area_unit = '[m^2]';
density_unit = '[kg/m^3]';

%% Initial Conditions
initial_velocity_vector_magnitude = input('Enter the initial velocity vector magnitude in m/s\n');
initial_angle = input('Enter the initial angle of launch in degrees\n');
initial_height(1) = input('Enter the initial height in meters\n'); % m - Height from ground
mass = input('Enter the objects mass in kg\n');
drag_truth = input('Do you want to account for atmospheric drag or not? 1 = Yes, 2 = No \n');
if drag_truth == 1
    cross_area = input('Enter the cross-sectional area in m^2 \n');
    wind_velocity = input('Enter the velocity of the wind in m/s \n');
end
initial_velocity_x(1) = initial_velocity_vector_magnitude*cos(deg2rad(initial_angle));
initial_velocity_y(1) = initial_velocity_vector_magnitude*sin(deg2rad(initial_angle));
initial_distance_x(1) = 0;
initial_kinetic_energy(1) = 1/2*mass*initial_velocity_y(1)^2; % Joules
initial_gravitational_potential_energy(1) = mass*abs(gravity)*initial_height(1); % Joules
%% Iterative Sequence
t(1) = 0;
dt = 1*10^-3; % s - Time Step 0.001s
i = 1;
if drag_truth == 1
    while initial_height >= 0
        % Iterative Calculations
        final_velocity_y(i) = initial_velocity_y(i) + gravity*dt;
        final_velocity_x(i) = initial_velocity_x(i) - 1/(2*mass)*rho_air*drag_coefficient*cross_area*(initial_velocity_x(i) + wind_velocity)^2*dt;
        final_height(i) = initial_height(i) + initial_velocity_y(i)*dt + 1/2*gravity*dt^2;
        final_kinetic_energy(i) = 1/2*mass*final_velocity_y(i)^2;
        final_gravitational_potential_energy(i) = mass*abs(gravity)*final_height(i);
        final_distance_x(i) = initial_distance_x(i) + final_velocity_x(i)*dt;
        
        % Changes between final and initial conditions
        delta_velocity_y(i) = final_velocity_y(i) - initial_velocity_y(i);
        delta_height(i) = final_height(i) - initial_height(i);
        delta_kinetic_energy(i) = final_kinetic_energy(i) - initial_kinetic_energy(i);
        delta_gravitational_potential_energy(i) = final_gravitational_potential_energy(i) - initial_gravitational_potential_energy(i);
        
        % Updating the variables
        initial_velocity_y(i+1) = final_velocity_y(i);
        initial_velocity_x(i+1) = final_velocity_x(i);
        initial_height(i+1) = final_height(i);
        initial_distance_x(i+1) = final_distance_x(i);
        initial_kinetic_energy(i+1) = final_kinetic_energy(i);
        initial_gravitational_potential_energy(i+1) = final_gravitational_potential_energy(i);
        
        t(i+1) = t(i) + dt;
        i = i + 1;
    end
else
    while initial_height >= 0
        % Iterative Calculations
        final_velocity_y(i) = initial_velocity_y(i) + gravity*dt;
        final_velocity_x(i) = initial_velocity_x(i);
        final_height(i) = initial_height(i) + initial_velocity_y(i)*dt + 1/2*gravity*dt^2;
        final_kinetic_energy(i) = 1/2*mass*final_velocity_y(i)^2;
        final_gravitational_potential_energy(i) = mass*abs(gravity)*final_height(i);
        final_distance_x(i) = initial_distance_x(i) + final_velocity_x(i)*dt;
        
        % Changes between final and initial conditions
        delta_velocity_y(i) = final_velocity_y(i) - initial_velocity_y(i);
        delta_height(i) = final_height(i) - initial_height(i);
        delta_kinetic_energy(i) = final_kinetic_energy(i) - initial_kinetic_energy(i);
        delta_gravitational_potential_energy(i) = final_gravitational_potential_energy(i) - initial_gravitational_potential_energy(i);
        
        % Updating the variables
        initial_velocity_y(i+1) = final_velocity_y(i);
        initial_velocity_x(i+1) = final_velocity_x(i);
        initial_height(i+1) = final_height(i);
        initial_distance_x(i+1) = final_distance_x(i);
        initial_kinetic_energy(i+1) = final_kinetic_energy(i);
        initial_gravitational_potential_energy(i+1) = final_gravitational_potential_energy(i);
        
        t(i+1) = t(i) + dt;
        i = i + 1;
    end
end
%% Total, Maximum, Final Conditions
% Total Conditions
total_velocity_y = sum(delta_velocity_y);
total_height = sum(delta_height);
total_distance_x = sum(final_distance_x);
total_kinetic_energy = sum(delta_kinetic_energy);
total_gravitational_potential_energy = sum(delta_gravitational_potential_energy);

% Maximum Conditions
maximum_height = max(final_height);
maximum_velocity_y = max(final_velocity_y);
maximum_kinetic_energy = max(final_kinetic_energy);
maximum_gravitational_potential_energy = max(final_gravitational_potential_energy);

% Final Conditions
% final_velocity_x = initial_velocity_x;
final_angle = rad2deg(atan(final_velocity_y(i-1)/final_velocity_x(i-1)));
final_velocity_vector_magnitude = sqrt(final_velocity_x(i-1)^2 + final_velocity_y(i-1)^2);

%% Data Visualization
p1 = polyfit(t(1:i-1),final_height(1:i-1),2)
y1 = polyval(p1,t(1:i-1));
figure()
hold on
plot(t(1:i-1),final_height(1:i-1));
plot(t(1:i-1),y1,'-');
xlabel('Time [s]');
ylabel('Vertical Distance [m]');
legend('Calculated Data','Quadratic Approximation');
grid on
hold off

% figure()
% plot(final_distance_x(1:i-1),final_velocity_y(1:i-1));
% xlabel('Horizontal Distance [m]');
% ylabel('Vertical Velocity [m/s]');
% grid on
%
% figure()
% plot(t(1:i-1),final_velocity_y(1:i-1));
% xlabel('Time [s]');
% ylabel('Vertical Velocity [m/s]');
% grid on

p2 = polyfit(t(1:i-1),final_kinetic_energy(1:i-1),2)
y2 = polyval(p2,t(1:i-1));
figure()
hold on
plot(t(1:i-1),final_kinetic_energy(1:i-1));
plot(t(1:i-1),y2,'-');
xlabel('Time [s]');
ylabel('Kinetic Energy [J]');
legend('Calculated Data','Quadratic Approximation');

p3 = polyfit(t(1:i-1),final_gravitational_potential_energy(1:i-1),2)
y3 = polyval(p3,t(1:i-1));
plot(t(1:i-1),final_gravitational_potential_energy(1:i-1));
plot(t(1:i-1),y3,'-');
xlabel('Time [s]');
ylabel('Gravitational Potential Energy [J]');
legend('Calculated Data','Quadratic Approximation');
grid on
hold off

figure()
plot(final_distance_x,final_height)
xlabel('Horizontal Distance [m]');
ylabel('Vertical Distance [m]');
grid on

%% Data Concatenation
linedivider = '-----------';
result = {
    'Variable Title','Numeric Value','Unit';
    linedivider,linedivider,linedivider;
    'Initial Conditions','','';
    linedivider,linedivider,linedivider;
    'Initial Velocity Vector Magnitude',initial_velocity_vector_magnitude,velocity_unit;
    'Initial Launch Angle',initial_angle,angle_unit;
    'Initial Vertical Velocity',initial_velocity_y(1),velocity_unit;
    'Initial Horizontal Velocity',initial_velocity_x(1),velocity_unit;
    'Initial Height',initial_height(1),distance_unit;
    'Initial Horizontal Distance',initial_distance_x(1),distance_unit;
    'Initial Kinetic Energy',initial_kinetic_energy(1),energy_unit;
    'Initial Gravitational Potential Energy',initial_gravitational_potential_energy(1),energy_unit;
    'Object Mass',mass,mass_unit;
    linedivider,linedivider,linedivider;
    'Final Conditions','','';
    linedivider,linedivider,linedivider;
    'Final Velocity Vector Magnitude',final_velocity_vector_magnitude,velocity_unit;
    'Final Landing Angle',final_angle,angle_unit;
    'Final Vertical Velocity',final_velocity_y(i-1),velocity_unit;
    'Final Horizontal Velocity',final_velocity_x(i-1),velocity_unit;
    'Final Height',final_height(i-1),distance_unit;
    'Final Horizontal Distance',final_distance_x(i-1),distance_unit;
    'Final Kinetic Energy',final_kinetic_energy(i-1),energy_unit;
    'Final Gravitational Potential Energy',final_gravitational_potential_energy(i-1),energy_unit;
    linedivider,linedivider,linedivider;
    'Total & Maximum Conditions','','';
    linedivider,linedivider,linedivider;
    'Total Vertical Velocity',total_velocity_y,velocity_unit;
    'Total Height',total_height,distance_unit;
    'Total Kinetic Energy',total_kinetic_energy,energy_unit;
    'Total Gravitational Potential Energy',total_gravitational_potential_energy,energy_unit;
    'Maximum Vertical Velocity',maximum_velocity_y,velocity_unit;
    'Maximum Height',maximum_height,distance_unit;
    'Maximum Kinetic Energy',maximum_kinetic_energy,energy_unit;
    'Maximum Gravitational Potential Energy',maximum_gravitational_potential_energy,energy_unit;
    };

display(result);
xlswrite('Projectile Motion Analysis Program Data',result,1);













































































