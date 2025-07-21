function vis_realign_image()
    % Ask user to select the motion parameter file
    path = "/home/xym/Desktop/MRI_ROMT/rOMT_spdup-main/data/ours/test1/DCE_nii_data";
    file = "rp_T1_FLASH_3D_baseline_0000.txt";
    if isequal(file, 0)
        disp('User selected Cancel');
        return;
    end
    
    param_file = fullfile(path, file);
    
    % Load the motion parameters
    motion_params = load(param_file);
    
    % Separate translation and rotation
    translation = motion_params(:, 1:3); % Already in mm
    rotation_rad = motion_params(:, 4:6); % In radians
    
    % Convert rotation from radians to mm displacement
    head_radius_mm = 50;
    rotation_mm = rotation_rad * head_radius_mm;
    
    time_points = 1:size(motion_params, 1);
    
    % Plotting
    figure('Name', "Motion Parameters: " + file, 'NumberTitle', 'off', 'Position', [100, 100, 800, 600]);
    
    % Plot Translation
    subplot(2, 1, 1);
    plot(time_points, translation);
    title('Translational Motion');
    xlabel('Volume Number');
    ylabel('Translation (mm)');
    legend('X', 'Y', 'Z', 'Location', 'best');
    grid on;
    
    % Plot Rotation
    subplot(2, 1, 2);
    plot(time_points, rotation_mm);
    title('Rotational Motion');
    xlabel('Volume Number');
    ylabel('Displacement (mm)');
    legend('Pitch (X)', 'Roll (Y)', 'Yaw (Z)', 'Location', 'best');
    grid on;

    % Save the figure to a file
    [~, name, ~] = fileparts(file);
    output_filename = fullfile(path, name + ".png");
    print(gcf, output_filename, '-dpng', '-r300'); % Save as PNG with 300 DPI
    fprintf('Figure saved to %s\n', output_filename);
    close(gcf); % Close the figure window to prevent it from trying to display