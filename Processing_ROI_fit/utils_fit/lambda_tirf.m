function [l11, l22, l33, l12, l23] = lambda_tirf(delta_n)

    l11     = ( 0.336535 + cos( delta_n/2 ).^3 .* ( -0.588936 + 0.252401*cos( delta_n ) ) + 4.2397*cos( delta_n/2 ) .* sin( delta_n/4 ).^6 + ( 4.47524 + 0.706617*cos( delta_n ) ) .* sin( delta_n/4 ).^6 ) ...
                ./ ( 1.68267 - 1.68267*cos( delta_n/2 ).^3 + 7.30441*sin( delta_n/4 ).^4 + 3.65221*cos( delta_n/2 ) .* sin( delta_n/4 ).^4 );

    l22     = ( 0.336535 + cos( delta_n/2 ).^3 .* ( -0.588936 + 0.252401*cos( delta_n ) ) + 8.90824*cos( delta_n/2 ) .* sin( delta_n/4 ).^6 + ( 9.40314 + 1.48471*cos( delta_n ) ) .* sin( delta_n/4 ).^6 ) ...
                ./ ( 1.68267 - 1.68267*cos( delta_n/2 ).^3 + 7.30441*sin( delta_n/4 ).^4 + 3.65221*cos( delta_n/2 ) .* sin( delta_n/4 ).^4 );

    l33     = ( 1.37482 - 1.0096*cos( delta_n/2 ).^5 + cos( delta_n/2 ).^3 .* ( -0.639136 + 0.273915*cos( delta_n ) ) ) ...
                ./ ( 1.68267 - 1.68267*cos( delta_n/2 ).^3 + 7.30441*sin( delta_n/4 ).^4 + 3.65221*cos( delta_n/2 ) .* sin( delta_n/4 ).^4 );
    
    l12     = ( (-2.425 - 2.29736.*cos(delta_n./2) - 0.382894.*cos(delta_n) ) .* sin(delta_n./4).^6 ) ...
                ./ ( 1.68267 - 1.68267.*cos(delta_n./2).^3 + 7.30441.*sin(delta_n./4).^4 + 3.65221.*cos(delta_n./2) .* sin(delta_n./4).^4 );
    
    l23     = ( -0.0780416 + cos(delta_n./2).^3 .* (0.136573 - 0.0585312.*cos(delta_n)) ) ...
                ./ ( 1.68267 - 1.68267.*cos(delta_n./2).^3 + 7.30441.*sin(delta_n./4).^4 + 3.65221.*cos(delta_n./2) .* sin(delta_n./4).^4 );
% %     varargout = {l12,l23};
end