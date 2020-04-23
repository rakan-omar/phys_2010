clear ; clc;
M = 20; % the mass in kg
BOMB_RADIUS = 0.25; %in meters
Bomb_Volume = (3/4)*pi*BOMB_RADIUS^3; %2 dimensinal
density = M/Bomb_Volume; % mass per unit volume
explosion_force = 2000 + (rand*1000); %newtons
explosion_impulse_time = 0.1;
impulse = explosion_force*explosion_impulse_time;
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

%if each fragment is a circle sector (with sides R = 0.25)
%the angle for each sector is given by
x_speed_fragments = zeros(1, n_fragments);
y_speed_fragments = zeros(1, n_fragments);
z_speed_fragments = zeros(1, n_fragments);

fragment_angles_theta = (2*pi)*(fragment_masses/M);
fragment_angles_phi = (pi)*(fragment_masses/M);
starting_angle_theta = rand*2*pi;
starting_angle_phi = rand*pi;
for i=1:n_fragments
    alpha = starting_angle_theta + fragment_angles_theta(i)/2;
    beta = starting_angle_phi + fragment_angles_phi(i)/2;
    x_speed_fragments(i) = impulse*cos(alpha)*cos(beta)/fragment_masses(i);
    y_speed_fragments(i) = impulse*sin(alpha)*cos(beta)/fragment_masses(i);
    z_speed_fragments(i) = impulse*sin(beta)/fragment_masses(i);
    starting_angle_theta = starting_angle_theta + fragment_angles_theta(i);
    starting_angle_phi = starting_angle_phi + fragment_angles_phi(i);
end
numberofsteps = 100;
fragment_x_positions = zeros(n_fragments, numberofsteps);
fragment_y_positions = zeros(n_fragments, numberofsteps);
fragment_z_positions = zeros(n_fragments, numberofsteps);

for k=1:numberofsteps
    t = 0.01*k;
    for i=1:n_fragments
        fragment_x_positions(i, k) = x_speed_fragments(i)*t;
        fragment_y_positions(i, k) = y_speed_fragments(i)*t;
        fragment_z_positions(i, k) = (z_speed_fragments(i)*t) - (5*t^2);
    end
end
for i=1:n_fragments
    plot3(fragment_x_positions(i, :), fragment_y_positions(i, :), fragment_z_positions(i, :))
    hold on
end
hold off %erase this to see a bunch of different explosions together
