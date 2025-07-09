function visualize_velocity_field(cfg, global_max_speed)
%VISUALIZE_VELOCITY_FIELD Saves and visualizes velocity fields from rOMT results.
%
% Syntax:  visualize_velocity_field(cfg, global_max_speed)
%
% Inputs:
%    cfg - Configuration struct containing all necessary parameters 
%          (u, out_dir, true_size, strid, view_azi_elevation, etc.).
%    global_max_speed - The maximum speed calculated across all time frames,
%                       used for consistent color scaling in visualizations.

fprintf('Saving and visualizing velocity fields...\n');
vis_dir = fullfile(cfg.out_dir, 'velocity_field_visualizations');
if ~exist(vis_dir, 'dir')
    mkdir(vis_dir);
end

n = cfg.true_size;
dim = 3;

% Create a downsampled grid for quiver plot to avoid visual clutter
strid = 4; % Stride: plot one vector every N voxels
[x_grid, y_grid, z_grid] = meshgrid(1:strid:n(2), 1:strid:n(1), 1:strid:n(3));

% If global_max_speed is zero (no motion), magnify_vel would be Inf.
% Handle this by setting a small default, although the plotting loop will likely be skipped.
if global_max_speed > 0
    magnify_vel = 10 / global_max_speed; % Magnification factor for arrow length
else
    magnify_vel = 1; % Default value, no magnification
end


parfor tind = 1:length(cfg.u)
    % Reshape the giant column vector into a (space*dim) x nt matrix
    u_interval_matrix = reshape(cfg.u{tind}, [prod(n) * dim, cfg.nt]);
    
    ti = cfg.first_time + (tind-1) * cfg.time_jump;
    tf = ti + cfg.time_jump;
    
    for k = 1:cfg.nt
        % --- 1. Process Data for the Current Frame ---
        vel_frame = u_interval_matrix(:, k);
        
        % Reshape to (voxels x dim) to separate x,y,z components
        u_xyz = reshape(vel_frame, [prod(n), dim]);
        
        % Reshape individual components into 3D matrices
        u_x = reshape(u_xyz(:,1), n);
        u_y = reshape(u_xyz(:,2), n);
        u_z = reshape(u_xyz(:,3), n);
        
        % Downsample the vector fields for plotting
        u_x_ds = u_x(1:strid:end, 1:strid:end, 1:strid:end);
        u_y_ds = u_y(1:strid:end, 1:strid:end, 1:strid:end);
        u_z_ds = u_z(1:strid:end, 1:strid:end, 1:strid:end);
        
        % Calculate speed (magnitude) for the downsampled vectors
        mags = sqrt(u_x_ds.^2 + u_y_ds.^2 + u_z_ds.^2);
        
        % --- 2. Create Visualization for the Current Frame ---
        fig = figure('Visible', 'off', 'Colormap', jet(256));
        q = quiver3(x_grid, y_grid, z_grid, u_y_ds*magnify_vel, u_x_ds*magnify_vel, u_z_ds*magnify_vel, 'AutoScale', 'off');
        hold on;
        % --- 3. Color arrows by magnitude (Robust Method) ---
        % Flatten the magnitudes of the vectors for coloring
        C = mags(:);
        
        % Get the current colormap
        cmap = colormap;
        % Normalize magnitudes to the range [0, 1] based on global max speed
        C_normalized = C / global_max_speed;
        C_normalized(C_normalized > 1) = 1; % Cap at max
        % Map normalized values to colormap indices
        color_indices = round(C_normalized * (size(cmap, 1) - 1)) + 1;
        % Get the RGB color for each arrow
        arrow_colors = cmap(color_indices, :);
        % Replicate colors for each vertex of the quiver head and tail
        % (This is the part that requires care)
        num_arrows = numel(q.UData);
        % Check if any arrows were drawn to avoid errors
        if num_arrows > 0
            % Create the correctly sized ColorData array.
            % Head has 3 vertices per arrow, Tail has 2.
            head_colors = repelem(arrow_colors, 3, 1);
            head_colors_rgba = uint8([head_colors, ones(size(head_colors,1), 1)]' * 255);
            tail_colors = repelem(arrow_colors, 2, 1);
            tail_colors_rgba = uint8([tail_colors, ones(size(tail_colors,1), 1)]' * 255);
            set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', head_colors_rgba);
            set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', tail_colors_rgba);
        end
        % Add semi-transparent brain mask
        msk_patch = patch(isosurface(1:n(2), 1:n(1), 1:n(3), cfg.msk, 0.5));
        msk_patch.FaceColor = [0.8, 0.8, 0.8];
        msk_patch.EdgeColor = 'none';
        msk_patch.FaceAlpha = 0.1;
        % Set plot properties
        axis equal; % Use 'axis equal' to preserve the aspect ratio
        axis tight;
        view(cfg.view_azi_elevation);
        if cfg.flip_z
            set(gca, 'ZDir', 'reverse');
        end
        title(sprintf('Velocity Field: Interval %d-%d, Frame %d', ti, tf, k));
        xlabel('x-axis'); ylabel('y-axis'); zlabel('z-axis');
        % Add colorbar and set fixed limits
        hcb = colorbar;
        hcb.Label.String = 'Speed (voxels/frame)';
        clim([0, global_max_speed]);
        set(gcf,'unit','normalized','position',[0.1, 0.1, 0.6, 0.7],'Color',[0.85,0.85,0.93], 'InvertHardcopy', 'off');
        % Save the visualization
        png_filename = sprintf('%s/velocity_field_interval_%d_to_%d_frame_%02d.png', vis_dir, ti, tf, k);
        saveas(fig, png_filename);
        close(fig); % Close the figure to free up memory
    end
end
fprintf('Velocity fields and speed maps saved in %s\n', vis_dir);

end
