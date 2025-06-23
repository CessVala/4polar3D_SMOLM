function savetiff(path, im)
% =========================================================================
% Function: savetiff
%
% Description:
%   Saves a 2D or 3D numeric image array to an uncompressed multi-page TIFF.
%   The input image is assumed to be arranged as (x, y, t).
%
%   The image is converted to 16-bit unsigned integers (`uint16`) before
%   saving. If normalization or contrast scaling is required, this must be 
%   done before calling this function.
%
% Inputs:
%   path    - (char) Full path to save the TIFF file, including filename and .tif extension
%   im      - (numeric array) 2D image data of any numeric type 
%             Dimensions should be (height, width, depth), where depth is 
%             usually time stacks
%
% Output:
%   Writes a TIFF file to disk at the specified location
%
% Example usage:
%   savetiff('C:\output\my_stack.tif', my_image_stack);
%
% Authors:
%   Charitra S. Senthil Kumar - Institut Fresnel  
%   Miguel Sison              - Institut Fresnel  
%   Cesar Valades-Cruz        - Institute of Hydrobiology (IHB), CAS
%
% Date: June 2025
% =========================================================================

    % Convert input to uint16 (no normalization here – assume preprocessed)
    im = uint16(im);

    % Display message
    disp(['Writing stack to ' path ' ...'])

    % Write first image slice (creates file)
    imwrite(squeeze(im(:,:,1)), path, 'tiff', 'Compression', 'none'); 

    % Append remaining slices
    for indz = 2:size(im, 3)
        imwrite(squeeze(im(:,:,indz)), path, ...
                'tiff', 'WriteMode', 'append', 'Compression', 'none'); 
    end

    % Completion message
    disp('... done!')
end
