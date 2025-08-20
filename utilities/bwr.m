function c = bwr(m)
%BWR Blue-White-Red colormap.
%   bwr(M) returns an M-by-3 matrix containing a blue-white-red colormap,
%   which is useful for displaying divergent data.
%   bwr (with no arguments) is the same length as the current colormap.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(bwr)
%
%   See also HSV, HOT, COOL, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

if nargin < 1
   m = size(get(gcf,'colormap'),1);
end

% Define the anchor points in the colormap
% Blue -> White -> Red
map_colors = [0 0 1;  % Blue
              1 1 1;  % White
              1 0 0]; % Red

% The positions of the anchor points
% Blue at the bottom (0), White in the middle (0.5), Red at the top (1)
map_positions = [0; 0.5; 1];

% Linearly interpolate to create a colormap of size 'm'
new_positions = linspace(0, 1, m)';
c = interp1(map_positions, map_colors, new_positions);

end