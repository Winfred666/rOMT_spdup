% Suppose D_tensor is XxYxZx3x3
D_tensor = load(fullfile(base_path, 'DTI_tensor_3_3.mat')).D_tensor;
sz = size(D_tensor);
N = prod(sz(1:3));

MD = zeros(sz(1:3));
FA = zeros(sz(1:3));
NormD = zeros(sz(1:3));

for x = 1:sz(1)
    for y = 1:sz(2)
        for z = 1:sz(3)
            D = squeeze(D_tensor(x,y,z,:,:));
            % Check symmetry
            if norm(D - D', 'fro') > 1e-6
                warning('Non-symmetric tensor at (%d,%d,%d)',x,y,z);
            end
            % Eigen decomposition
            lam = eig((D+D')/2);  % enforce symmetry
            if any(lam < 0)
                % Negative eigenvalues = bad tensor
                lam(lam < 0) = 0;
            end
            % Scalars
            MD(x,y,z)   = mean(lam);
            NormD(x,y,z)= norm(D,'fro');
            
            if any(isnan(lam)) || all(lam==0)
                FA(x,y,z) = 1e-5;
                continue
            end
            lam_mean = mean(lam);
            num = sqrt( sum((lam - lam_mean).^2) );
            den = sqrt( sum(lam.^2) );
            FA(x,y,z) = sqrt(3/2) * num / den;

        end
    end
end



% Quick visualization

VIS = MD;

VIS(VIS < 0) = 0;

volshow(VIS);   % average diffusivity field

figure; histogram(VIS(:), 200);
xlabel('MD (mm^2/s)'); ylabel('Count'); title('MD distribution');




% Assume MD is X×Y×Z (mean diffusivity volume)
[x,y,z] = ndgrid(1:size(VIS,1), 1:size(VIS,2), 1:size(VIS,3));

% pick slice locations
x_slices = round(linspace(1, size(VIS,1), 25));
y_slices = round(linspace(1, size(VIS,2), 25));
z_slices = round(linspace(1, size(VIS,3), 25));

figure

hs = slice(y, x, z, VIS, x_slices, y_slices, z_slices);

% add brain mask

brain_mask=load('template.mat').tmp;
%magnify = 1;%1;
%mskfv = isosurface(x,y,z,brain_mask,0.2);
%mskp = patch(mskfv);
%mskp.FaceColor = [.3,.3,.3];
%mskp.FaceAlpha= 0.2;
%mskp.EdgeColor = 'none';
%mskp.EdgeAlpha= 0;
%mskp.DisplayName = 'mask';





% smooth interpolation and transparency
set(hs, 'EdgeColor','none', 'FaceColor','interp', 'FaceAlpha',0.05);
alpha('color');
alphamap(linspace(0,1,100)); 
% visualization aesthetics
axis image; axis tight;
xlabel('x'); ylabel('y'); zlabel('z');
title('Fractional Anisotropy');
colormap(jet);

% robust range
clim([0, prctile(VIS(:),100)]);
colorbar;

% background color and 3D view
set(gca, 'Color', [0.85,0.85,0.93]);
view(45,30);  % azimuth, elevation
grid on; box off;



diffVIS = MD;
y_slices = round(linspace(1, size(diffVIS,2), 16)); % 16 slices

figure;
for i = 1:numel(y_slices)
    subplot(4,4,i); % 4x4 grid of subplots
    imagesc(squeeze(diffVIS(:,y_slices(i),:)));
    axis image off;
    
    % symmetric colormap for positive/negative differences
    clim([0, max(abs(diffVIS(:)))]);
    colormap(gca, jet);
    
    title(['y = ' num2str(y_slices(i))]);
end

colorbar('Position',[0.93 0.1 0.02 0.8]); % add one global colorbar
sgtitle('slices of FA');

