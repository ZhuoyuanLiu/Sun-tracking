% Define the vertices of a rectangular plate with thickness
thickness = 0.2;
vertices = [-1 -1 -thickness/2; 1 -1 -thickness/2; 1 1 -thickness/2; -1 1 -thickness/2; 
            -1 -1 thickness/2; 1 -1 thickness/2; 1 1 thickness/2; -1 1 thickness/2];

% Define the faces of the plate
faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];

% Rotation angle (in degrees)
theta = 180; % Sun's motion across the sky 
yaw_angle = 30; % Yaw rotation angle
% Initial yaw angle
initial_yaw = -15;

% Define the declination of the sun, and the right ascension
Dec_sun = 13; % As of 8/17

% Number of frames for animation
numFrames = 100;

% Define the point and direction of the rotation axis
axisPoint = [0, -1.5, 0];
axisDirection = [0, 1, 0.9004]; % Rotation axis, Polar axis elevation 42 degs

% Normalize the direction vector
axisDirection = axisDirection / norm(axisDirection);

% Define a fixed roll angle (in degrees)
roll_angle = -45;   % 45 -> sun noon 51;  80 -> sun noon 86
% Rotation matrix around an arbitrary axis
R_axis = @(theta) [...
    cosd(theta) + axisDirection(1)^2 * (1 - cosd(theta)), ...
    axisDirection(1) * axisDirection(2) * (1 - cosd(theta)) - axisDirection(3) * sind(theta), ...
    axisDirection(1) * axisDirection(3) * (1 - cosd(theta)) + axisDirection(2) * sind(theta); ...
    axisDirection(2) * axisDirection(1) * (1 - cosd(theta)) + axisDirection(3) * sind(theta), ...
    cosd(theta) + axisDirection(2)^2 * (1 - cosd(theta)), ...
    axisDirection(2) * axisDirection(3) * (1 - cosd(theta)) - axisDirection(1) * sind(theta); ...
    axisDirection(3) * axisDirection(1) * (1 - cosd(theta)) - axisDirection(2) * sind(theta), ...
    axisDirection(3) * axisDirection(2) * (1 - cosd(theta)) + axisDirection(1) * sind(theta), ...
    cosd(theta) + axisDirection(3)^2 * (1 - cosd(theta))];

% Yaw rotation matrix around z-axis
R_yaw = @(yaw) [cosd(yaw) -sind(yaw) 0; sind(yaw) cosd(yaw) 0; 0 0 1];

% Roll rotation matrix around the plate's longitudinal axis
R_roll = @(roll) [1 0 0; 0 cosd(roll) -sind(roll); 0 sind(roll) cosd(roll)];

% Create a figure
figure;
hold on;
axis equal;
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
view(3);



% Enable interactive rotation
rotate3d on;

% Plot the ground plane
patch([-10 10 10 -10], [-10 -10 10 10], [0 0 0 0], [0.8 0.8 0.8], 'FaceAlpha', 0.5);

% Plot the horizontal coordinate hemisphere
[az, el] = meshgrid(linspace(-pi, pi, 50), linspace(0, pi/2, 50));
[x, y, z] = sph2cart(az, el, 10);
surf(x, y, z, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
% Add directional labels to the axes
text(-10, 0, 0, 'East', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'blue'); % East to X axis negative
text(10, 0, 0, 'West', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'blue'); % West to X axis positive
text(0, -10, 0, 'North', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'blue'); % North to Y axis negative
text(0, 10, 0, 'South', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'blue'); % South to Y axis positive
% Apply the initial yaw rotation to the vertices
initial_rotatedVertices = (R_yaw(initial_yaw) * vertices')';

% Plot the initial plate
patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'cyan');

% Plot the rotation axis
line([axisPoint(1) axisPoint(1) + axisDirection(1) * 2], ...
     [axisPoint(2) axisPoint(2) + axisDirection(2) * 2], ...
     [axisPoint(3) axisPoint(3) + axisDirection(3) * 2], 'Color', 'red', 'LineWidth', 2);

% Pause to visualize the initial state
pause(1);

% Initialize an array to store the normal vectors and angles
normal_vectors = zeros(numFrames, 3);
azimuth_angles = zeros(numFrames, 1);
elevation_angles = zeros(numFrames, 1);
sun_positions = zeros(numFrames, 3);

% Animate the rotation
for i = 1:numFrames
    % Compute the rotation for the current frame
    currentTheta = 90 + theta * (i / numFrames);
    currentYaw = yaw_angle * (i / numFrames);
    currentR_axis = R_axis(currentTheta);
    currentR_yaw = R_yaw(currentYaw + initial_yaw);

    % Apply the rotation to the vertices with respect to the axis
    rotatedVertices = (currentR_axis * (vertices' - axisPoint') + axisPoint')';
    rotatedVertices = (currentR_yaw * rotatedVertices')';
    % Apply the fixed roll rotation
    rotatedVertices = (R_roll(roll_angle) * rotatedVertices')';
    % Clear the previous plot
    cla;
    
    % Add directional labels to the axes
    text(-10, 0, 0, 'East', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'blue'); % East to X axis negative
    text(10, 0, 0, 'West', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'blue'); % West to X axis positive
    text(0, -10, 0, 'North', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'blue'); % North to Y axis negative
    text(0, 10, 0, 'South', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'blue'); % South to Y axis positive
    % Replot the ground plane
    patch([-10 10 10 -10], [-10 -10 10 10], [0 0 0 0], [0.8 0.8 0.8], 'FaceAlpha', 0.5);
    
    % Replot the horizontal coordinate hemisphere
    surf(x, y, z, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    
    % Plot the rotated plate
    patch('Vertices', rotatedVertices, 'Faces', faces, 'FaceColor', 'cyan');
    
    % Plot the rotation axis
    line([axisPoint(1) axisPoint(1) + axisDirection(1) * 2], ...
         [axisPoint(2) axisPoint(2) + axisDirection(2) * 2], ...
         [axisPoint(3) axisPoint(3) + axisDirection(3) * 2], 'Color', 'red', 'LineWidth', 2);
    
    hold on;
    
    % Compute the normal vector of the plate using cross product of two edges
    edge1 = rotatedVertices(2, :) - rotatedVertices(1, :);
    edge2 = rotatedVertices(4, :) - rotatedVertices(1, :);
    normal = -cross(edge1, edge2);
    normal = normal / norm(normal); % Normalize the normal vector
    % Store the normal vector
    normal_vectors(i, :) = normal;
    
    % Calculate the azimuth and elevation angles
    azimuth = atan2d(normal(2), normal(1));
    elevation = asind(normal(3));
    
    % Store the angles
    azimuth_angles(i) = azimuth;
    elevation_angles(i) = elevation;

    % Simulate the sun's position in the X-Z plane using azimuth and
    % elevation
    sun_azimuth = azimuth_angles(i); % Azimuth angle of the sun
    sun_elevation = elevation_angles(i); % Elevation angle 
    
    % Calculate the sun's position vector
    unit_length = 10;
    sun_x = unit_length*cosd(sun_elevation) * cosd(sun_azimuth);
    sun_y = unit_length*cosd(sun_elevation) * sind(sun_azimuth);
    sun_z = unit_length*sind(sun_elevation);
    sun_position = [sun_x, sun_y, sun_z];
    
    % Move the sun's direction vector to originate from the center of the plate
    plate_center = mean(rotatedVertices); % Calculate the center of the plate
    sun_position = sun_position + plate_center; % Offset sun position by the plate center

    % Store the sun's position
    sun_positions(i, :) = sun_position;
    
    % Plot the sun's position
    plot3(sun_position(1),  sun_position(2), sun_position(3), 'o-', 'Color', 'yellow', 'LineWidth', 2);
    
    % Update the figure and allow interactive rotation
    drawnow;
    
    % Pause to create animation effect
    pause(0.05);
end

% Display the normal vectors over time
figure;
subplot(2, 1, 1);
plot(1:numFrames, azimuth_angles, '-r', 'DisplayName', 'Azimuth');
xlabel('Frame');
ylabel('Azimuth Angle (degrees)');
title('Azimuth Angle of Normal Vector Over Time');
legend show;
grid on;
subplot(2, 1, 2);
plot(1:numFrames, elevation_angles, '-b', 'DisplayName', 'Elevation');
xlabel('Frame');
ylabel('Elevation Angle (degrees)');
title('Elevation Angle of Normal Vector Over Time');
legend show;
grid on;

figure;
plot(1:numFrames, normal_vectors(:, 1), '-r', 'DisplayName', 'X component');
hold on;
plot(1:numFrames, normal_vectors(:, 2), '-g', 'DisplayName', 'Y component');
plot(1:numFrames, normal_vectors(:, 3), '-b', 'DisplayName', 'Z component');
xlabel('Frame');
ylabel('Normal Vector Component');
title('Normal Vector Components Over Time');
legend show;
grid on;
