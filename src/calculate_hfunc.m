function hfunc = calculate_hfunc(input_type, varargin)
% CALCULATE_HFUNC Calculate the helicity function of a given object based
% on the input type provided.
%
% The input type can be one of "image", "volume" or "mesh". The function
% takes a variable number of input arguments based on the input type
% provided.
%
% Inputs:
%   input_type  - string that specifies the type of input object. Must be
%                 one of "image", "volume" or "mesh".
%   varargin    - variable number of input arguments based on the input
%                 type provided.
%
% Outputs:
%   hfunc       - A scalar value representing the helicity function of the
%                 given object.
%
% Usage:
%   hfunc = calculate_hfunc("image", img, mask, ps, orientation, centroid, delta_alpha, delta_rho)
%       Calculates the helicity function for an image IMG with a given
%       pixel size PS, taking only the pixels marked by the MASK into
%       account. ORIENTATION is the angle between the helical axis and the
%       horizontal axis in radians, CENTROID is a 2-element vector (x0,y0)
%       with the center of the object. DELTA_ALPHA and DELTA_RHO specify
%       the bin size for the helicity function.
%   hfunc = calculate_hfunc("volume", vol, mask, vs, orientation, centroid, delta_alpha, delta_rho)
%       calculates the helicity function for volumetric data VOL with voxel
%       size VS, taking only the voxels marked by the MASK into account.
%       ORIENTATION is a vector (u,v,w) parallel to the helical axis,
%       CENTROID is a 3-element vector (x0,y0,z0) with the center of the
%       object. DELTA_ALPHA and DELTA_RHO specify the bin size for the
%       helicity function.
%   hfunc = calculate_hfunc("surface", verts, faces, orientation, centroid, delta_alpha, delta_rho)
%       calculates the helicity function for a surface mesh with vertices
%       VERTS and connectivity list FACES. ORIENTATION is a vector (x,y,z)
%       parallel to the helical axis, CENTROID is a 3-element vector
%       (x0,y0,z0) with the center of the object. DELTA_ALPHA and DELTA_RHO
%       specify the bin size for the helicity function.
%
% See also:
%   HelicityFunction
%
switch input_type
    case "image"
        img         = varargin{1};
        mask        = varargin{2};
        ps          = varargin{3};
        orientation = varargin{4};
        centroid    = varargin{5};
        delta_alpha = varargin{6};
        delta_rho   = varargin{7};
        hfunc = calculate_hfunc_image(img, ps, mask, ...
            orientation, centroid, delta_alpha, delta_rho);
    case "volume"
        vol         = varargin{1};
        mask        = varargin{2};
        vs          = varargin{3};
        orientation = varargin{4};
        centroid    = varargin{5};
        delta_alpha = varargin{6};
        delta_rho   = varargin{7};
        hfunc = calculate_hfunc_volume(vol, vs, mask, ...
            orientation, centroid, delta_alpha, delta_rho);
    case "surface"
        verts       = varargin{1};
        faces       = varargin{2};
        orientation = varargin{3};
        centroid    = varargin{4};
        delta_alpha = varargin{5};
        delta_rho   = varargin{6};
        hfunc = calculate_hfunc_mesh(verts, faces, ...
            orientation, centroid, delta_alpha, delta_rho);
    otherwise
        error("Invalid input_type. Must be 'image', 'volume' or 'surface'");
end

% -------------------------------------------------------------------------
function hfunc = calculate_hfunc_image(img, ps, mask, orientation, centroid, delta_alpha, delta_rho)
arguments
    img (:,:)
    ps (1,1)
    mask (:,:)
    orientation (1,1)
    centroid (2,1)
    delta_alpha (1,1)
    delta_rho (1,1)
end
assert(all(size(img) == size(mask)), "img and mask must have the same dimensions.")

idx_list = find(mask);

% Define the coordinate grid
[x,y] = meshgrid(1:size(img,2), 1:size(img,1));
x = ps * x;
y = ps * flipud(y);

% Use the gradient for edge detection and orientation mapping
[gmag, gdir] = imgradient(img);

% First calculate rho and alpha in each pixel
rho = zeros(size(img));
alpha = zeros(size(img));
for i = 1:length(idx_list)
    idx = idx_list(i);
    rho(idx) = abs(cos(orientation) * (centroid(2) - y(idx)) - sin(orientation) * (centroid(1) - x(idx)));
    alpha(idx) = deg2rad(gdir(idx)) - orientation;  % theta is wrt x-axis, so don't subtract pi/2 to account for
                                                    % the difference between normal and inclination angle
    while alpha(idx) >= pi/2, alpha(idx) = alpha(idx) - pi; end
    while alpha(idx) < -pi/2, alpha(idx) = alpha(idx) + pi; end
end

% Then create a helicity function and fill it with the previously
% calculated values
hfunc = HelicityFunction(delta_alpha, delta_rho, max(rho,[],"all"));
for i = 1:length(idx_list)
    idx = idx_list(i);
    hfunc = add_value(hfunc, alpha(idx), rho(idx), gmag(idx));
end

% -------------------------------------------------------------------------
function hfunc = calculate_hfunc_volume(vol, vs, mask, orientation, centroid, delta_alpha, delta_rho)
arguments
    vol (:,:,:)
    vs (1,1)
    mask (:,:,:)
    orientation (3,1)
    centroid (3,1)
    delta_alpha (1,1)
    delta_rho (1,1)
end
error("Not implemented.")

% -------------------------------------------------------------------------
function hfunc = calculate_hfunc_mesh(verts, faces, orientation, centroid, delta_alpha, delta_rho)
arguments
    verts (:,3)
    faces (:,3)
    orientation (3,1)
    centroid (3,1)
    delta_alpha (1,1)
    delta_rho (1,1)
end

% Transform verts such that the z-axis is the helical axis
R = rotmatrix_vec2vec(orientation, [0 0 1]);
verts = (verts - centroid') * R';

% Get surface description
A = verts(faces(:,1),:);
B = verts(faces(:,2),:);
C = verts(faces(:,3),:);
surface_vectors = 1/2 * cross(B - A, C - A, 2);

positions = (A + B + C) / 3;
areas = vecnorm(surface_vectors, 2, 2);
normals = surface_vectors ./ areas;

% Get rho and alpha for each surface element
% The inclination angle can be calculate using the arctan of the ratio
% of the z-component of the normals and the theta-component
[theta, rho, ~] = cart2pol(positions(:,1), positions(:,2), positions(:,3));
normals_theta = -normals(:,1) .* sin(theta) + normals(:,2) .* cos(theta);
alpha = atan2(normals(:,3), normals_theta) + pi/2;  % difference of pi/2 between surface and normal
alpha(alpha > pi/2) = alpha(alpha > pi/2) - pi;     % remap ]pi/2,3pi/2] to ]-pi/2,pi/2]

% Initialize and fill the helicity function
hfunc = HelicityFunction(delta_alpha, delta_rho, max(rho));
for idx = 1:size(faces, 1)
    hfunc = add_value(hfunc, alpha(idx), rho(idx), areas(idx));
end

% -------------------------------------------------------------------------
function R = rotmatrix_vec2vec(moving, fixed)
arguments
    moving (3,1) double
    fixed  (3,1) double
end

a = moving / vecnorm(moving);
b = fixed / vecnorm(fixed);

if a == b
    R = eye(3);
elseif a == -b
    R = -eye(3);
    R(2,2) = 1;
else
    % following https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    v = cross(a, b);
    c = dot(a, b);
    vx = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
    R = eye(3) + vx + vx^2 / (1 + c);
end
