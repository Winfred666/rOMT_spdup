function lymphatics_percent_image(lymph)



%Create mean image based on the baselines 
for FN=1:size(lymph.bas_loc,1)
    clear v img;

    v=spm_vol(char(lymph.bas_loc(FN)));
    img=spm_read_vols(v);

    if FN==1
        img_bas=img;
    else
        img_bas=img_bas+img;
    end

end


img_bas=img_bas/size(lymph.bas_loc,1);

% --- ask for a mask from lymph.msk ---

clear v_msk img_msk
v_msk=spm_vol(lymph.msk);
mask=spm_read_vols(v_msk);
mask=mask>0.01;

img_sum=zeros(size(img_bas));

% --- Initialize vectors to store total signal change over time ---
num_time_points = size(lymph.src_loc,1);
total_absolute_change = zeros(num_time_points, 1);
total_relative_change = zeros(num_time_points, 1);

% --- PASS 1: Analyze all images to find the target mass ---
disp('Pass 1: Analyzing signal change to find target mass...');
for FN=1:size(lymph.src_loc,1)
    v = spm_vol(char(lymph.src_loc(FN)));
    img = spm_read_vols(v);

    % Calculate absolute change for plotting later
    img_abs_change = img - img_bas;
    total_absolute_change(FN) = sum(img_abs_change(mask));

    % Calculate relative change and its total sum (mass)
    img_rel_change = (img - img_bas + eps) ./ (img_bas + eps) * 100;
    total_relative_change(FN) = sum(img_rel_change(mask));
end

% --- Determine the target mass from the analysis pass ---
target_mass = max(total_relative_change);
fprintf('Target mass for normalization set to: %f\n', target_mass);

% --- PASS 2: Correct images to match target mass and save ---
disp('Pass 2: Applying mass conservation and saving corrected images...');
corrected_relative_change = zeros(num_time_points, 1); % For verification plot

for FN=1:size(lymph.src_loc,1)
    clear v img fdir fname;
    v = spm_vol(char(lymph.src_loc(FN)));
    img = spm_read_vols(v);

    % Recalculate the relative change image
    img_rel_change = (img - img_bas + eps) ./ (img_bas + eps) * 100;
    
    % Get the current mass for this image (already computed)
    current_mass = total_relative_change(FN);
    
    % Calculate the correction factor to match the target mass
    % Add a small epsilon to avoid division by zero if current_mass is 0
    correction_factor = target_mass / (current_mass + eps);
    
    % Apply the correction
    corrected_img = img_rel_change * correction_factor;
    
    % Store the new total mass for the verification plot
    corrected_relative_change(FN) = sum(corrected_img(mask));

    % --- Original file writing with the corrected image ---
    [fdir, fname] = fileparts(char(lymph.src_loc(FN)));
    v.fname = [fdir,'/p',fname,'.nii'];
    spm_write_vol(v, corrected_img);
    fprintf('Corrected and saved %s\n', v.fname);

    % Sum all corrected images
    img_sum = img_sum + corrected_img;
end


v.fname=[fdir,'/pbase_',fname,'.nii'];
spm_write_vol(v,img_bas);
v.fname=[fdir,'/psum_',fname,'.nii'];
spm_write_vol(v,img_sum);

% --- Create and Save Mass Conservation Plot ---
disp('Generating and saving mass conservation plot...');
time_points = 1:num_time_points;

figure('Name', 'Mass Conservation Analysis', 'Visible', 'off');

% Use yyaxis to plot both on different scales
yyaxis left;
plot(time_points, total_absolute_change, '-o', 'LineWidth', 1.5, 'DisplayName', 'Original Absolute Change');
hold on;
ylabel('Total Signal Change');
title(['Mass Conservation Analysis for ' fname]);

yyaxis right;
plot(time_points, total_relative_change, '--s', 'LineWidth', 1.5, 'DisplayName', 'Original Relative Change');
hold on;
plot(time_points, corrected_relative_change, '-d', 'LineWidth', 2.0, 'DisplayName', 'Corrected Relative Change');
ylabel('Total Relative Signal Change (%)');


xlabel('Time Point (Image Index)');
legend('show', 'Location', 'best');
grid on;

% Save the figure to the output directory
plot_fname = fullfile(fdir, 'mass_conservation_analysis.png');
saveas(gcf, plot_fname);
fprintf('Saved mass conservation plot to %s\n', plot_fname);
close(gcf); % Close the invisible figure


end