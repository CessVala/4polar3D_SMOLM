% for compat octave/matlab

function fid = display_out()
%   fid = 0 ; % pas de sortie (NULL)
%   fid = 1 ; % = stdout
   fid = 2 ; % = stderr pour Octave / en rouge pour Matlab
end