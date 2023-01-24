% HelicityFunction is a class for handling helicity functions, which
% describe the helicity of a given structure based on its helical
% parameters rho and alpha.
%
% Properties:
%   alpha_center   - The center of the bins for the alpha axis
%   alpha_upper    - The upper bound of the bins for the alpha axis
%   rho_center     - The center of the bins for the rho axis
%   rho_upper      - The upper bound of the bins for the rho axis
%   bin_area       - The area of each bin in the helicity function
%   normalization_factor - The factor used to normalize the helicity function
%   hfunc          - The non-normalized helicity function
%
% Methods:
%   HelicityFunction - Constructs a new helicity function object
%   plus             - Adds two helicity functions together
%   total_helicity   - Calculates the total helicity from the helicity function
%   add_value        - Adds a value to the helicity function at the correct bin
%   plot             - Plots the helicity function
%
% See also:
%   calculate_hfunc
%
classdef HelicityFunction

    properties (SetAccess = protected, GetAccess = public)
        alpha_center
        alpha_upper
        rho_center
        rho_upper
        bin_area
        normalization_factor
        hfunc
    end

    methods
        function obj = HelicityFunction(delta_alpha, delta_rho, max_rho)
            % Initialize a new helicity function
            %
            % Inputs:
            %   delta_alpha - The width of the bins for the alpha axis
            %   delta_rho   - The width of the bins for the rho axis
            %   max_rho     - The maximum value for the rho axis
            n_alpha = round(pi/2 / delta_alpha);
            alpha_edges = linspace(0, pi/2, n_alpha + 1);
            obj.alpha_center = (alpha_edges(1:end-1) + alpha_edges(2:end)) / 2;
            obj.alpha_upper = alpha_edges(2:end);

            n_rho = ceil(max_rho / delta_rho);
            max_rho = n_rho * delta_rho;
            rho_edges = linspace(0, max_rho, n_rho + 1);
            obj.rho_center = (rho_edges(1:end-1) + rho_edges(2:end)) / 2;
            obj.rho_upper = rho_edges(2:end);

            obj.hfunc = zeros(n_rho, n_alpha);

            obj.bin_area = (obj.alpha_center(2) - obj.alpha_center(1)) * (obj.rho_center(2) - obj.rho_center(1));
            obj.normalization_factor = 0;
        end

        function obj = plus(obj1, obj2)
            % Add two helicity functions together
            %
            % Inputs:
            %   obj1 - First helic
            %   obj2 - Second helicity function object
            %
            % Output:
            %   obj - A new helicity function object that is the sum of obj1 and obj2

            % Determine which object has the larger rho axis and assign to obj
            if length(obj1.rho_center) >= length(obj2.rho_center)
                obj = obj1;
                obj0 = obj2;
            else
                obj = obj2;
                obj0 = obj1;
            end

            nrho0 = length(obj0.rho_center);
            assert(all(obj.alpha_center == obj0.alpha_center), "Alpha axes must match exactly")
            assert(all(obj.rho_center(1:nrho0) == obj0.rho_center), "Discretization of rho axes must match")

            obj.normalization_factor = obj.normalization_factor + obj0.normalization_factor;
            obj.hfunc(1:nrho0,:) = obj.hfunc(1:nrho0,:) + obj0.hfunc;
        end

        function htot = total_helicity(obj)
            % TOTAL_HELICITY Calculate the total helicity from the helicity function
            %
            % Output:
            %   htot - The total helicity calculated from the helicity function
            htot = sum(obj.hfunc, "all") / obj.normalization_factor;
        end

        function obj = add_value(obj, alpha, rho, value)
            % ADD_VALUE Add a value to the helicity function at the right location
            %
            % Inputs:
            %   alpha - The value to add to the alpha axis
            %   rho   - The value to add to the rho axis
            %   value - The value to add to the helicity function
            i = find(rho <= obj.rho_upper, 1);
            j = find(abs(alpha) <= obj.alpha_upper, 1);
            if isempty(i)
                warning("Rho is out of bounds")
                return;
            end

            obj.normalization_factor = obj.normalization_factor + value;
            if j > 1 && j < length(obj.alpha_center)
                obj.hfunc(i,j) = obj.hfunc(i,j) + sign(alpha) * value;
            end
        end

        function plot(obj)
            % PLOT Plot the helicity function
            %
            % This function will create a plot of the helicity function with the
            % x-axis as the alpha axis, y-axis as the rho axis and the color map
            % representing the helicity function.
            hfunc_norm = obj.hfunc / obj.bin_area / obj.normalization_factor;
            imagesc(obj.alpha_center, obj.rho_center, hfunc_norm);
            set(gca, "YDir", "normal")

            cmax = max(abs(hfunc_norm), [], "all");
            clim([-cmax, cmax])
            colormap(gca, flipud(brewermap([], "RdBu")))
            set(gca, "Color", brewermap(1, "RdBu"))
            colorbar
        end
    end
end
