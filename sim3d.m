%% This program simulates and calculates the centre of mass of a bomb following a projectile motion
% In this program it is assumed that the bomb starts at on piece and then
% explodes into multiple fragments
clear ; clc; %clear previous program
%% Setting up paremters
M = 100; % the mass in kg
phi = pi/3; % the angle that the bomb starts at in deg
theta = 0; % angle from x axis, in x-y plane. shot along x-axis
speed = 200; % the initial speed of the bomb in m/s
x_speed = speed*cos(theta)*cos(phi) ; % x component of the speed
y_speed = speed*sin(theta)*cos(phi);
z_speed = speed*sin(phi) ; % z component of the speed
v = [x_speed, y_speed, z_speed];
g = 9.81; % the acceleraton of gravity in m+/s^2
z_at_peak = z_speed^2/(2*g);
time_at_peak = z_speed/g;
x_at_peak = x_speed*time_at_peak;
y_at_peak = y_speed*time_at_peak; %0
%range of explosion is between quarter and three-quarters through
%doing this here because we need speed and position right before explosion
% explosion_time = time_at_peak + rand*time_at_peak*0.5;
explosion_time = time_at_peak*0.5 +rand*time_at_peak;
velocity_pre_explosion = [x_speed, y_speed, 0]; %z set later, x-y constant
explosion_position = [0, 0, 0];
%% the trajectory of the bomb
tau = 0.01; %timestep
Nsteps = ceil((2*time_at_peak)/tau); % number of steps
explosion_step = ceil(explosion_time/tau);
x_position = zeros(1, Nsteps+1);
y_position = zeros(1, Nsteps+1);
z_position = zeros(1, Nsteps+1);
for i=2:Nsteps+1
% now we will update the velocity and the position since the acceleration
% is constant
    t = i*tau ; % current time
% from the kinematic equations
    x_position(i) = x_speed*t; % x = xo +v(ox)t 
    y_position(i) = y_speed*t; %0
    z_position(i) = z_speed*t - 0.5*g*t^2 ; % y = yo + v(oy)t-gt^2/2
    if i == explosion_step
        velocity_pre_explosion(3) = z_speed - g*t;
        explosion_position = [x_position(i), y_position(i), z_position(i)];
    end
    if y_position(i) < 0 %to counter error due to rounding
        y_position(i) = 0;
        x_position(i) = x_speed*2*time_at_peak;
    end
end
%% plotting 
clf ; % clear and move forward
plot3(x_position(1, 1:explosion_step), y_position(1, 1:explosion_step), ...
    z_position(1, 1:explosion_step),'r-','LineWidth',2)
grid on
hold on
plot3(x_position, y_position, z_position,'g--','Linewidth',1.2)
title('position y against position x')
xlabel('x position (m)')
ylabel('position (m)')
zlabel('height (m)')
grid on
hold on

%part 2
%randomise some aspects of explosion
explosion_impulse_time = 0.1;
n_fragments = 4 + ceil(rand*3); %5-7
fragment_masses = zeros(1, n_fragments);
remaining_mass = M;
for i=1:n_fragments-1
    %largest fragment is less than 60% of remainder, greater than 10%
    %this way we don't have (or 'ignore') minisicule pieces
    fragment_masses(i) = (0.1*remaining_mass) + (0.5*rand*remaining_mass);
    remaining_mass = remaining_mass - fragment_masses(i);
end
fragment_masses(n_fragments) = remaining_mass;

%if each fragment is a circle sector, with sides R = 0.25,
%the angle for each sector is given by
fragment_angles_theta = (2*pi)*(fragment_masses/M);
fragment_angles_phi = (pi)*(fragment_masses/M);
%randomising a starting point, relative to the horizontal
explosion_force = 10000 + (rand*5000); %newtons
impulse = explosion_force*explosion_impulse_time;
x_speed_fragments = zeros(1, n_fragments);
y_speed_fragments = zeros(1, n_fragments);
z_speed_fragments = zeros(1, n_fragments);
starting_angle_theta = rand*2*pi;
starting_angle_phi = rand*pi;%randomised initially, cimulative
%getting initial velocities, conservation of momentum
for i=1:n_fragments
    
    %the centre of mass of each circle sector with angle gamma, is along
    %the radius starting from gamma/2 (due to symmetry).
    %so, the angle from from the horizontal to the direction each fragment
    %will fly off in is given by
    
    alpha = starting_angle_theta + fragment_angles_theta(i)/2;
    beta = starting_angle_phi + fragment_angles_phi(i)/2;
    
    %mv(aft) - mv(bef) = ft
    %v(aft) = (ft/m)+v(bef) or v(aft) = (ft + v(bef))/m
    
    x_speed_fragments(i) = (impulse*cos(alpha)*cos(beta) + ...
        (fragment_masses(i)*velocity_pre_explosion(1)))/fragment_masses(i);
    
    y_speed_fragments(i) = (impulse*sin(alpha)*cos(beta) + ...
        (fragment_masses(i)*velocity_pre_explosion(2)))/fragment_masses(i);
    
    z_speed_fragments(i) = (impulse*sin(beta) + ...
        (fragment_masses(i)*velocity_pre_explosion(3)))/fragment_masses(i);
    
    starting_angle_theta = starting_angle_theta + fragment_angles_theta(i);
    starting_angle_phi = starting_angle_phi + fragment_angles_phi(i);
end
%these get extended, in the following loop
fragment_x_positions = zeros(n_fragments, 1);
fragment_y_positions = zeros(n_fragments, 1);
fragment_z_positions = zeros(n_fragments, 1);

all_midair = ones(1, n_fragments);
reached_ground = 0; %n of fragments that reached ground
midair = true;
k = 0; %while loop counter
while (midair)
    t = tau*k;
    k = k+1;
    for i=1:n_fragments
        if all_midair(i)
            fragment_x_positions(i, k) = (x_speed_fragments(i)*t)+explosion_position(1);
            fragment_y_positions(i, k) = (y_speed_fragments(i)*t)+explosion_position(2);
            fragment_z_positions(i, k) = ((z_speed_fragments(i)*t)-(0.5*g*t^2))+explosion_position(3);
            if (fragment_z_positions(i, k) < 0)
                fragment_z_positions(i, k) = 0;
                t_to_ground_1 = (z_speed_fragments(i)/g) + (((z_speed_fragments(i)^2) + (2*g*explosion_position(3)))^0.5)/g;
%                 t_to_ground_2 = (z_speed_fragments(i)/g) - (((z_speed_fragments(i)^2) + (2*g*explosion_position(3)))^0.5)/g;
                fragment_x_positions(i, k) = explosion_position(1) + (x_speed_fragments(i)*t_to_ground_1);
                fragment_y_positions(i, k) = explosion_position(2) + (y_speed_fragments(i)*t_to_ground_1);
                
                all_midair(i) = 0;
                reached_ground = reached_ground + 1;
                if (reached_ground == n_fragments)
                    midair = false;
                    %DO NOT ADD 'BREAK', still need x-y positions to be set
                end
            end
        else
            % already on the floor, x-y positions constant. z = 0
            fragment_x_positions(i, k) = fragment_x_positions(i, k-1);
            fragment_y_positions(i, k) = fragment_y_positions(i, k-1);
        end
    end
end

%plotting the trajectory of the peices
for i=1:n_fragments
    plot3(fragment_x_positions(i, :), fragment_y_positions(i, :), fragment_z_positions(i, :),'LineWidth',0.8);
    grid on
    hold on
end


centre_of_mass = zeros(3, k);
for j=1:k
    for i=1:n_fragments
        centre_of_mass(1, j) = centre_of_mass(1, j) + (fragment_x_positions(i, j) * fragment_masses(i));
        centre_of_mass(2, j) = centre_of_mass(2, j) + (fragment_y_positions(i, j) * fragment_masses(i));
        centre_of_mass(3, j) = centre_of_mass(3, j) + (fragment_z_positions(i, j) * fragment_masses(i));
    end
end
centre_of_mass = centre_of_mass/M;

plot3(centre_of_mass(1, :),centre_of_mass(2, :),centre_of_mass(3, :), ...
    'b-', 'lineWidth', 2) % plotting the centre of mass trajectory
f_names = strings(1, 3 + n_fragments);
f_names(1) = 'bomb';
f_names(2) = 'expected trajectory';
for i=3:2+n_fragments
    f_names(i) = 'fragment ' + string(i-2);
end
f_names(n_fragments+3) = 'centre of mass';
legend(f_names) %remove this to view graph clearly
