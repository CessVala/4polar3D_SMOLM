function [l_iso] = lambda_iso(delta)
    l_iso = (1-cos(delta./2)).*(2+cos(delta./2))/6;
end