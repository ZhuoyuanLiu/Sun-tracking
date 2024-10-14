% Load the saved data file
load('sun_data_2021.mat', 'sun_data');

% Convert Right Ascension (RA) from hh:mm:ss to decimal hours
RA_values = sun_data.RA;
RA_num = zeros(size(RA_values)); % Initialize array to hold numeric RA values

for i = 1:length(RA_values)
    parts = str2double(strsplit(RA_values{i}, ':')); % Split into hh, mm, ss
    RA_num(i) = parts(1) + parts(2)/60 + parts(3)/3600; % Convert to decimal hours
end

% Convert Declination (DEC) from dd:mm:ss to decimal degrees
DEC_values = sun_data.DEC;
DEC_num = zeros(size(DEC_values)); % Initialize array to hold numeric DEC values

for i = 1:length(DEC_values)
    parts = str2double(strsplit(DEC_values{i}, ':')); % Split into dd, mm, ss
    DEC_num(i) = abs(parts(1)) + parts(2)/60 + parts(3)/3600; % Convert to decimal degrees
    if parts(1) < 0
        DEC_num(i) = -DEC_num(i); % Adjust sign for negative declination
    end
end

% Number of frames for animation
numFrames = 100; % 0.1 hour per frame

% Define the vertices of a rectangular plate with thickness
thickness = 0.2;
vertices = [-1 -1 -thickness/2; 1 -1 -thickness/2; 1 1 -thickness/2; -1 1 -thickness/2; 
            -1 -1 thickness/2; 1 -1 thickness/2; 1 1 thickness/2; -1 1 thickness/2];

% Define the faces of the plate
faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];

% Loop through each day of the year (365 days)
for day = 1:365
    display(day);
    % Set the Right Ascension and Declination for the current day
    RA_sample = RA_num(day);
    Dec_sample = DEC_num(day);
    
    % Define the point and direction of the rotation axis
    axisPoint = [0, 1.5, 0];
    axisDirection = [0, 1, 0.9004]; % Rotation axis, Polar axis elevation 42 degs

    % Normalize the direction vector
    axisDirection = axisDirection / norm(axisDirection);

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
    R_roll = @(roll) [cosd(roll) -sind(roll) 0; sind(roll) cosd(roll) 0; 0 0 1];

    % Pitch rotation matrix around the x-axis
    R_pitch = @(pitch) [1 0 0; 0 cosd(pitch) -sind(pitch); 0 sind(pitch) cosd(pitch)];

    % Create a figure for each day
    figure;
    hold on;
    axis equal;
    grid on;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    view(3);

    % Plot the ground plane
    patch([-10 10 10 -10], [-10 -10 10 10], [0 0 0 0], [0.8 0.8 0.8], 'FaceAlpha', 0.5);

    % Plot the horizontal coordinate hemisphere
    [az, el] = meshgrid(linspace(-pi, pi, 50), linspace(0, pi/2, 50));
    [x, y, z] = sph2cart(az, el, 10);
    surf(x, y, z, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    
    % Add directional labels to the axes
    text(-10, 0, 0, 'East', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'blue');
    text(10, 0, 0, 'West', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'blue');
    text(0, -10, 0, 'North', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'blue');
    text(0, 10, 0, 'South', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'blue');

    % Define the pitch and roll angles
    pitch_angle = 90 - (180 - (90 - Dec_sample + 42));
    roll_angle = 0;  % Fixed roll angle

    % Rotation angle (in degrees)
    theta_change = 15 * 14; % Sun's motion across the sky

    % Initial yaw angle
    initial_yaw = 0;

    % Animate the rotation for the current day
    for i = 1:numFrames
        % Compute the rotation for the current frame
        Initial_theta = (90 + RA_sample);
        currentTheta = Initial_theta - theta_change * (i / numFrames);
        currentYaw = initial_yaw * (i / numFrames);
        currentR_axis = R_axis(currentTheta);
        currentR_yaw = R_yaw(-currentYaw + initial_yaw);

        % Apply the fixed roll rotation
        rotatedVertices = (R_roll(roll_angle) * vertices')';
    
        % Apply the fixed pitch rotation
        rotatedVertices = (R_pitch(pitch_angle) * rotatedVertices')';

        % Apply the yaw and axis rotation
        rotatedVertices = (currentR_yaw * rotatedVertices')';
        rotatedVertices = (currentR_axis * (rotatedVertices' - axisPoint') + axisPoint')';

        % Clear the previous plot
        cla;
        
        % Replot the ground plane and hemisphere
        patch([-10 10 10 -10], [-10 -10 10 10], [0 0 0 0], [0.8 0.8 0.8], 'FaceAlpha', 0.5);
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

        % Calculate the azimuth and elevation angles
        azimuth = atan2d(normal(2), normal(1));
        elevation = asind(normal(3));

        % Simulate the sun's position in the X-Z plane using azimuth and elevation
        unit_length = 10;
        sun_x = -unit_length*cosd(elevation) * cosd(azimuth);
        sun_y = -unit_length*cosd(elevation) * sind(azimuth);
        sun_z = -unit_length*sind(elevation);
        sun_position = [sun_x, sun_y, sun_z];

        % Move the sun's direction vector to originate from the center of the plate
        plate_center = mean(rotatedVertices); % Calculate the center of the plate
        sun_position = sun_position + plate_center; % Offset sun position by the plate center

        % Plot the sun's position
        plot3(sun_position(1),  sun_position(2), sun_position(3), 'o-', 'Color', 'yellow', 'LineWidth', 2);

        % Update the figure
        drawnow;

        % Pause to create animation effect
        pause(0.05);
    end
    
    % Close the figure after the animation of the current day is done
    close(gcf);
    
end
