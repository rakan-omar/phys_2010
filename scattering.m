clear ; clc;
M = 20; % the mass in kg
BOMB_RADIUS = 0.25; %in meters
Bomb_area = pi*BOMB_RADIUS^2; %2 dimensinal
density = M/Bomb_area; % mass per unit area
explosion_force = 400 + (rand*200); %newtons
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

fragment_angles = (2*pi)*(fragment_masses/M);
starting_angle = rand*2*pi;
for i=1:n_fragments
    tilda = starting_angle + fragment_angles(i)/2;
    x_speed_fragments(i) = impulse*cos(tilda)/fragment_masses(i);
    y_speed_fragments(i) = impulse*sin(tilda)/fragment_masses(i);
    starting_angle = starting_angle + fragment_angles(i);
end
numberofsteps = 100;
fragment_x_positions = zeros(n_fragments, numberofsteps);
fragment_y_positions = zeros(n_fragments, numberofsteps);
for k=1:numberofsteps
    t = 0.01*k;
    for i=1:n_fragments
        fragment_x_positions(i, k) = x_speed_fragments(i)*t;
        fragment_y_positions(i, k) = (y_speed_fragments(i)*t) - (5*t^2);
    end
end
for i=1:n_fragments
    plot(fragment_x_positions(i, :), fragment_y_positions(i, :))
    hold on
end
hold off %erase this to see a bunch of different explosions together