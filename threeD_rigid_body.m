% ----------------------------
% Let the XY plane be the horizontal plane
% ----------------------------


% Define the vertices of a rectangular plate with thickness
thickness = 0.2;
vertices = [-1 -1 -thickness/2; 1 -1 -thickness/2; 1 1 -thickness/2; -1 1 -thickness/2; 
            -1 -1 thickness/2; 1 -1 thickness/2; 1 1 thickness/2; -1 1 thickness/2];

% Define the faces of the plate
faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];

% Rotation angle (in degrees)
theta = 180; % Complete rotation

% Number of frames for animation
numFrames = 100;

% Define the point and direction of the rotation axis
axisPoint = [0, -1.5, 0];
axisDirection = [0, 1, 0.9004]; % Arbitrary direction vector (Polar axis elevation 42 degs )

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

% Create a figure
figure;
axis equal;
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
view(3);

% Enable interactive rotation
rotate3d on;

% Plot the initial plate
patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'cyan');

% Plot the rotation axis
line([axisPoint(1) axisPoint(1) + axisDirection(1) * 2], ...
     [axisPoint(2) axisPoint(2) + axisDirection(2) * 2], ...
     [axisPoint(3) axisPoint(3) + axisDirection(3) * 2], 'Color', 'red', 'LineWidth', 2);

% Pause to visualize the initial state
pause(1);

% Animate the rotation
for i = 1:numFrames
    % Compute the rotation for the current frame
    currentTheta = 90 + theta * (i / numFrames);
    currentR_axis = R_axis(currentTheta);
    
    % Apply the rotation to the vertices with respect to the axis
    rotatedVertices = (currentR_axis * (vertices' - axisPoint') + axisPoint')';
    
    % Clear the previous plot
    cla;
    
    % Plot the rotated plate
    patch('Vertices', rotatedVertices, 'Faces', faces, 'FaceColor', 'cyan');
    
    % Plot the rotation axis
    line([axisPoint(1) axisPoint(1) + axisDirection(1) * 2], ...
         [axisPoint(2) axisPoint(2) + axisDirection(2) * 2], ...
         [axisPoint(3) axisPoint(3) + axisDirection(3) * 2], 'Color', 'red', 'LineWidth', 2);
    
    % Update the figure and allow interactive rotation
    drawnow;
    
    % Pause to create animation effect
    pause(0.05);
end
