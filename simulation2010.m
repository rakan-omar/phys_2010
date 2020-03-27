%% This program simulates and calculates the centre of mass of a bomb following a projectile motion
% In this program it is assumed that the bomb starts at on piece and then
% explodes into multiple fragments
clear ; clc; %clear previous program
%% Setting up paremters
M = 100; % the mass in kg
theta = pi/3; % the angle that the bomb starts at in deg
speed = 200; % the initial speed of the bomb in m/s
x_speed = speed*cos(theta) ; % x component of the speed
y_speed = speed*sin(theta) ; % y component of the speed
v = [x_speed, y_speed];
g = 9.81; % the acceleraton of gravity in m+/s^2
y_at_peak = y_speed^2/(2*g);
time_at_peak = y_speed/g;
x_at_peak = x_speed*time_at_peak;

%range of explosion is between quarter and three-quarters through
%doing this here because we need speed and position right before explosion
% explosion_time = time_at_peak + rand*time_at_peak*0.5;
explosion_time = time_at_peak*0.5 +rand*time_at_peak;
velocity_pre_explosion = [x_speed,0];
explosion_position = [0, 0];
%% the trajectory of the bomb
tau = 0.01; %timestep
Nsteps = ceil((2*time_at_peak)/tau); % number of steps
explosion_step = ceil(explosion_time/tau);
x_position = zeros(1, Nsteps+1);
y_position = zeros(1, Nsteps+1);
for i=2:Nsteps+1
% now we will update the velocity and the position since the acceleration
% is constant
    t = i*tau ; % current time
% from the kinematic equations
    x_position(i) = x_speed*t ; % x = xo +v(ox)t 
    y_position(i) = y_speed*t - 0.5*g*t^2 ; % y = yo + v(oy)t-gt^2/2
    if i == explosion_step
        velocity_pre_explosion(2) = y_speed - g*t;
        explosion_position = [x_position(i), y_position(i)];
    end
    if y_position(i) < 0 %to counter error due to rounding
        y_position(i) = 0;
        x_position(i) = x_speed*2*time_at_peak;
    end
end
%% plotting 
clf ; % clear and move forward
plot(x_position(1, 1:explosion_step),y_position(1, 1:explosion_step),'r-','LineWidth',2)
grid on
hold on
plot(x_position,y_position,'g--','Linewidth',1.2)
title('position y against position x')
xlabel('range position in m')
ylabel('height position in m')
grid on
hold on

%part 2
%randomise some aspects of explosion
explosion_impulse_time = 0.1;
n_fragments = 5;
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
fragment_angles = 2*pi*fragment_masses/M;
%randomising a starting point, relative to the horizontal
explosion_force = 3000 + (rand*1500); %newtons
impulse = explosion_force*explosion_impulse_time;
x_speed_fragments = zeros(1, n_fragments);
y_speed_fragments = zeros(1, n_fragments);
starting_angle = 2*pi*rand; %randomised initially, cimulative
%getting initial velocities, conservation of momentum
for i=1:n_fragments
    
    %the centre of mass of each circle sector with angle gamma, is along
    %the radius starting from gamma/2 (due to symmetry).
    %so, the angle from from the horizontal to the direction each fragment
    %will fly off in is given by
    
    tilda = starting_angle + fragment_angles(i)/2;
    % mv(aft) - mv(bef) = ft
    
    x_speed_fragments(i) = ((impulse * cos(tilda))+ ...
        (fragment_masses(i)*velocity_pre_explosion(1)))/fragment_masses(i);   
    y_speed_fragments(i) = ((impulse * sin(tilda))+ ...
        (fragment_masses(i)*velocity_pre_explosion(2)))/fragment_masses(i);  
    
    starting_angle = starting_angle + fragment_angles(i);
end
%these get extended, in the following loop
fragment_x_positions = zeros(n_fragments, 1);
fragment_y_positions = zeros(n_fragments, 1);

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
            fragment_y_positions(i, k) = ((y_speed_fragments(i)*t)-(0.5*g*t^2))+explosion_position(2);
            if (fragment_y_positions(i, k) < 0)
                fragment_y_positions(i, k) = 0;
                t_to_ground_1 = (y_speed_fragments(i)/g) + (((y_speed_fragments(i)^2) + (2*g*explosion_position(2)))^0.5)/g;
%                 t_to_ground_2 = (y_speed_fragments(i)/g) - (((y_speed_fragments(i)^2) + (2*g*explosion_position(2)))^0.5)/g;
                fragment_x_positions(i, k) = explosion_position(1) + (x_speed_fragments(i)*t_to_ground_1);
                all_midair(i) = 0;
                reached_ground = reached_ground + 1;
                if (reached_ground == n_fragments)
                    midair = false;
                    %DO NOT ADD 'BREAK', still need x positions to be set
                end
            end
        else
            fragment_x_positions(i, k) = fragment_x_positions(i, k-1);
        end
    end
end

%plotting the trajectory of the peices
for i=1:n_fragments
    plot(fragment_x_positions(i, :),fragment_y_positions(i, :),'LineWidth',0.8);
    grid on
    hold on
end


centre_of_mass = zeros(2, k);
for j=1:k
    for i=1:n_fragments
        centre_of_mass(1, j) = centre_of_mass(1, j) + (fragment_x_positions(i, j) * fragment_masses(i));
        centre_of_mass(2, j) = centre_of_mass(2, j) + (fragment_y_positions(i, j) * fragment_masses(i));
    end
end
centre_of_mass = centre_of_mass/M;

plot(centre_of_mass(1, :),centre_of_mass(2, :),'b-', 'lineWidth', 2) % plotting the centre of mass trajectory
f_names = strings(1, 3 + n_fragments);
f_names(1) = 'bomb';
f_names(2) = 'expected trajectory';
for i=3:2+n_fragments
    f_names(i) = 'fragment ' + string(i-2);
end
f_names(n_fragments+3) = 'centre of mass';
legend(f_names) %remove this to view graph clearly
