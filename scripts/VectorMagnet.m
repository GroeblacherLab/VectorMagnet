%==========================================================================
% VectorMagnet.m
% Author: GroeblacherLab / Gaia Da Prato
% Date: 03/09/2024
%
% Description:
%   Defines a 3D vector magnet composed of three orthogonal superconducting
%   bobbins. Computes the total magnetic field at an arbitrary point in space
%   by summing contributions from each coil using their geometry, position,
%   and polarity.
%
% Related publication:
%   G. Da Prato et al., Rev. Sci. Instrum. 96, 065208 (2025)
%   https://doi.org/10.1063/5.0270187
%==========================================================================

classdef VectorMagnet
    properties
        bobx  % Bobbin aligned along x axis
        boby  % Bobbin aligned along y axis
        bobz  % Bobbin aligned along z axis

        posx  % 3D position of bobx center [x y z]
        posy  % 3D position of boby center
        posz  % 3D position of bobz center

        polx  % Polarity sign of bobx (±1)
        poly  % Polarity sign of boby (±1)
        polz  % Polarity sign of bobz (±1)
    end

    methods
        function obj = VectorMagnet(bobx, boby, bobz, posx, posy, posz, polx, poly, polz)
            % Constructor for VectorMagnet
            obj.bobx = bobx;
            obj.boby = boby;
            obj.bobz = bobz;

            obj.posx = posx;
            obj.posy = posy;
            obj.posz = posz;

            obj.polx = polx;
            obj.poly = poly;
            obj.polz = polz;
        end

        function [Bx, By, Bz] = tot_magnetic_field(obj, pos, I)
            % Computes the total magnetic field at a 3D position.
            % Inputs:
            %   pos - 3D position [x y z] in mm
            %   I   - Currents [Ix Iy Iz] in A through each bobbin
            % Outputs:
            %   Bx, By, Bz - Magnetic field components at pos [T]

            % Relative coordinates to each bobbin center
            x_x = pos(1) - obj.posx(1); x_y = pos(2) - obj.posx(2); x_z = pos(3) - obj.posx(3);
            y_x = pos(1) - obj.posy(1); y_y = pos(2) - obj.posy(2); y_z = pos(3) - obj.posy(3);
            z_x = pos(1) - obj.posz(1); z_y = pos(2) - obj.posz(2); z_z = pos(3) - obj.posz(3);

            x_side = sign(x_x); y_side = sign(y_y); z_side = sign(z_z);

            % Calculate magnetic field from each bobbin
            [bobx_ax, bobx_r, ~] = obj.bobx.magnetic_field_off_axis_bobbin(x_x, sqrt(x_y^2+x_z^2), I(1));
            [boby_ax, boby_r, ~] = obj.boby.magnetic_field_off_axis_bobbin(y_y, sqrt(y_z^2+y_x^2), I(2));
            [bobz_ax, bobz_r, ~] = obj.bobz.magnetic_field_off_axis_bobbin(z_z, sqrt(z_x^2+z_y^2), I(3));
            bobx_r = abs(bobx_r); boby_r = abs(boby_r); bobz_r = abs(bobz_r);

            % Bobbin X
            bobx_x = bobx_ax*obj.polx;
            if abs(x_z) > 1e-8 && abs(x_y) > 1e-8
                bobx_y = bobx_r / sqrt(1 + (x_z / x_y)^2) * obj.polx * x_side * sign(x_y);
                bobx_z = abs(bobx_y * x_z / x_y) * obj.polx * x_side * sign(x_z);
            elseif abs(x_z) <= 1e-8
                bobx_z = 0;
                bobx_y = bobx_r * obj.polx * x_side * sign(x_y);
            elseif abs(x_y) <= 1e-8
                bobx_y = 0;
                bobx_z = bobx_r * obj.polx * x_side * sign(x_z);
            end

            % Bobbin Y
            boby_y = boby_ax * obj.poly;
            if abs(y_x) > 1e-8 && abs(y_z) > 1e-8
                boby_z = boby_r / sqrt(1 + (y_x / y_z)^2) * obj.poly * y_side * sign(y_z);
                boby_x = abs(boby_z * y_x / y_z) * obj.poly * y_side * sign(y_x);
            elseif abs(y_x) <= 1e-8
                boby_x = 0;
                boby_z = boby_r * obj.poly * y_side * sign(y_z);
            elseif abs(y_z) <= 1e-8
                boby_z = 0;
                boby_x = boby_r * obj.poly * y_side * sign(y_x);
            end

            % Bobbin Z
            bobz_z = bobz_ax * obj.polz;
            if abs(z_y) > 1e-8 && abs(z_x) > 1e-8
                bobz_x = bobz_r / sqrt(1 + (z_y / z_x)^2) * obj.polz * z_side * sign(z_x);
                bobz_y = abs(bobz_x * z_y / z_x) * obj.polz * z_side * sign(z_y);
            elseif abs(z_y) <= 1e-8
                bobz_y = 0;
                bobz_x = bobz_r * obj.polz * z_side * sign(z_x);
            elseif abs(z_x) <= 1e-8
                bobz_x = 0;
                bobz_y = bobz_r * obj.polz * z_side * sign(z_y);
            end

            % Sum field components
            Bx = bobx_x + boby_x + bobz_x;
            By = bobx_y + boby_y + bobz_y;
            Bz = bobx_z + boby_z + bobz_z;
        end

        function [points_matrix, field_matrix] = field_map_in_plane(obj, norm, P0, grid_size, grid_range, I)
            % Generates a magnetic field map in a plane
            % Inputs:
            %   norm       - Normal vector to the plane
            %   P0         - Center point of the plane [x y z]
            %   grid_size  - [Nx, Ny] number of points in u and v directions
            %   grid_range - [Rx, Ry] span in each direction (mm)
            %   I          - Currents [Ix Iy Iz] in A
            %
            % Outputs:
            %   points_matrix - Nx*Ny x 3 grid point coordinates
            %   field_matrix  - Nx*Ny x 3 magnetic field vectors at each point

            points_matrix = obj.generatePlanePoints(norm, P0, grid_size, grid_range);
            Bx = zeros(size(points_matrix,1),1); By = Bx; Bz = Bx;

            for i = 1:size(points_matrix,1)
                [Bx_tmp, By_tmp, Bz_tmp] = obj.tot_magnetic_field(points_matrix(i,:), I);
                Bx(i) = Bx_tmp;
                By(i) = By_tmp;
                Bz(i) = Bz_tmp;
            end
            field_matrix = [Bx(:), By(:), Bz(:)];
        end

        function points_matrix = generatePlanePoints(~, normal_vec, P0, grid_size, grid_range)
            % Generates a set of points on a plane given a normal vector.
            %
            % Inputs:
            % normal_vec - A 1x3 vector representing the normal to the plane [nx, ny, nz].
            % P0 - A 1x3 vector representing the central point on the plane [x0, y0, z0].
            % grid_size - A 1x2 vector [rows, cols] specifying the number of grid points in each direction.
            % grid_range - A scalar specifying the range in each direction.
            %
            % Outputs:
            % points_matrix - An (rows*cols)x3 matrix where each row represents a point [x, y, z] on the plane.

            % Find two vectors that lie on the plane
            if normal_vec(1) ~= 0 || normal_vec(2) ~= 0
                v1 = [-normal_vec(2), normal_vec(1), 0];
            else
                v1 = [1, 0, 0];  % Choose a different vector if normal_vec is parallel to z-axis
            end
            v2 = cross(normal_vec, v1);  % Second vector perpendicular to v1 and normal_vec

            % Normalize the vectors (optional, if you want unit vectors)
            v1 = v1 / norm(v1);
            v2 = v2 / norm(v2);

            % Generate a 2D grid of points (u, v) using grid size and range
            u = linspace(-grid_range(1)/2, grid_range(1)/2, grid_size(1));  % u-coordinates
            v = linspace(-grid_range(2)/2, grid_range(2)/2, grid_size(2));  % v-coordinates
            [U, V] = meshgrid(u, v);

            % Calculate the coordinates of the points on the plane
            X = P0(1) + U * v1(1) + V * v2(1);
            Y = P0(2) + U * v1(2) + V * v2(2);
            Z = P0(3) + U * v1(3) + V * v2(3);
            % Store the points in a matrix
            points_matrix = [X(:), Y(:), Z(:)];
        end
    end
end
