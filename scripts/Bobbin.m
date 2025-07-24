%==========================================================================
% Bobbin.m
% Author: GroeblacherLab / Gaia Da Prato
% Date: 03/09/2024
%
% Description:
%   Class for modeling a superconducting bobbin (coil) used in a
%   cryogenic three-axis vector magnet. Includes methods to estimate
%   physical parameters and compute the magnetic field both on-axis
%   and off-axis using analytical formulas.
%
% Related publication:
%   G. Da Prato et al., Rev. Sci. Instrum. 96, 065208 (2025)
%   https://doi.org/10.1063/5.0270187
%==========================================================================

classdef Bobbin
    properties
        r_in               % Internal radius of the bobbin (mm)
        r_ex               % External radius of the bobbin (mm)
        h                  % Height of the bobbin (mm)
        t                  % Wire thickness (mm)
        p                  % Winding pitch (mm)

        layers             % Number of radial layers
        winds_per_layer    % Number of windings per layer
        winds              % Total number of windings
        wire_length        % Total wire length (m)
        volume             % Volume occupied by the windings (mm^3)
        weight             % Estimated bobbin mass (g)
    end

    methods
        function obj = Bobbin(r_in, r_ex, h, t, p)
            % Constructor: initializes bobbin geometry and derived quantities.
            %
            % Inputs:
            %   r_in - internal radius [mm]
            %   r_ex - external radius [mm]
            %   h    - height [mm]
            %   t    - wire thickness [mm]
            %   p    - winding pitch [mm]

            rho = 8.96; % Copper density in g/cm^3
            obj.r_in = r_in;
            obj.r_ex = r_ex;
            obj.h = h;
            obj.t = t;
            obj.p = p;
            obj.layers = floor((r_ex - r_in) / t); % number of layers
            obj.winds_per_layer = floor(h / p); % number of windings in each layer
            obj.winds = obj.layers * obj.winds_per_layer; % number of total windings
            obj.wire_length = obj.winds*pi*(r_in*2+t*(obj.layers+1))*1e-3; % wire length in meters
            obj.volume = pi * (r_ex^2 - r_in^2) * h; % Volume occupied by windings in mm^3
            obj.weight = obj.volume * rho * 1e-3; % Windings mass in grams
        end

        function [B_z, B_rho, B] = magnetic_field_off_axis_loop(~, z, rho, r, I)
            % Computes the magnetic field components of a single loop at
            % an off-axis point (z, rho).
            %
            % Inputs:
            %   z   - axial distance from loop center [mm]
            %   rho - radial distance from axis [mm]
            %   r   - loop radius [mm]
            %   I   - current in loop [A]
            %
            % Outputs:
            %   B_z    - axial component [T]
            %   B_rho  - radial component [T]
            %   B      - total magnitude [T]

            mu0 = 4 * pi * 1e-4; % [T·m/A]

            alpha = rho / r;
            beta = z / r;
            gamma = z / rho;
            Q = (1 + alpha)^2 + beta^2;
            k = sqrt(4 * alpha / Q);

            if I < 1e-8
                B_z = 0; B_rho = 0; B = 0;
                return
            end

            B0 = mu0 * I / (2 * r);
            [K, E] = ellipke(k^2); % Elliptic integrals

            if rho == 0
                % On-axis field
                if z == 0
                    B = B0;
                else
                    B = mu0 * I * r^2 / (2 * (r^2 + z^2)^(1.5));
                end
                B_z = B; B_rho = 0;
            else
                B_z = B0 / (pi * sqrt(Q)) * ...
                    (E * (1 - alpha^2 - beta^2) / (Q - 4 * alpha) + K);
                B_rho = B0 * gamma / (pi * sqrt(Q)) * ...
                    (E * (1 + alpha^2 + beta^2) / (Q - 4 * alpha) - K);
                B = sqrt(B_z^2 + B_rho^2);
            end
        end

        function [B_z, B_rho, B] = magnetic_field_off_axis_bobbin(obj, z, rho, I)
            % Computes total magnetic field at a point (z, rho) from
            % all the windings in the bobbin.
            %
            % Inputs:
            %   z, rho - coordinates relative to bobbin center [mm]
            %   I      - current [A]
            %
            % Outputs:
            %   B_z, B_rho - axial and radial components [T]
            %   B          - total magnitude [T]

            Rin = obj.r_in;
            Rex = obj.r_ex;
            th = obj.t;
            pitch = obj.p;

            B_z = 0;
            B_rho = 0;
            B = 0;

            N = floor((Rex - Rin) / th);    % radial turns
            M = floor(obj.h / pitch);       % axial turns

            % Displace reference frame from center to bobbin top
            z = z - obj.h / 2;

            % Radial and axial loop positions
            R = obj.GenerateFrameList(Rin + th / 2, th, N);
            z_v = obj.GenerateFrameList(z + th / 2, pitch, M);

            for i = 1:N
                for j = 1:M
                    [Bz_tmp, Brho_tmp, B_tmp] = ...
                        obj.magnetic_field_off_axis_loop(z_v(j), rho, R(i), I);
                    B_z = B_z + Bz_tmp;
                    B_rho = B_rho + Brho_tmp;
                    B = B + B_tmp;
                end
            end
        end

        function V = GenerateFrameList(~, A, S, N)
            % Generates a 1D array V starting from A, with step S, and number of
            % elements N
            V = A + S * (0:N-1);
        end

    end
end
