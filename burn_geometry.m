function [Ab,Vb,Vc] = burn_geometry(ri,ro,h,bd)
%  This burnback model only assumes burning with the cylindrical 
%  perforation of the grain.  Modify as needed for your rocket motor. 

    % ri: inner radius
    % ro: outer radius (constant)
    % h: grain height
    % bd burned distance

    if bd >= ro-ri % motor is burnt out
        Ab = 0; % [m^2] 
        
    else % there is grain remaining
        %% BURN AREA
        Ab = 2*pi*(ri + bd)*(h-2*bd) + 2*pi*(ro.^2-(ri+bd).^2); % [m^2] total burn area for cylindrically perforated grain
    
    end
    
    %% VOLUME of PROPELLANT CONSUMED
    Vtot = pi*ro^2*h;
    
    Vb = 3.1416*((ri + bd)^2 - ri^2)*(h-2*bd); % [m^3]
    
    %% CHAMBER VOLUME
    Vc = Vtot-Vb; % [m^3]  ignores the volume between the end of 