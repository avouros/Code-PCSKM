% DESCRIPTION: Exports the current figure into different formats.
% ARGUMENTS:
% -get_current_figure:    figure for exporting, for current "active" figure 
%                         use gcf.
%
% -get_path:              path in which the figure should be saved.
%
% -get_name:              name of the generated file (without extension).
%
% -get_export_format:     see List of accepted formats.      
%
% -get_export_properties: see List of accepted qualities: 
%                         low: 150dpi (save as)
%                         medium: 200dpi
%                         high: 300dpi
%                         super: 600dpi
%                         screen: screen resultion

% List of accepted formats:
% - .pdf
% - .svg
% - .eps
% - .jpg
% - .png
% - .tif
% - .bmp

% List of accepted qualities:
% - Low Quality
% - Medium Quality
% - High Quality
% - Super Quality
% - Screen Quality

% Example:
% Save the current figure to the path 'C:\Documents\' as 'myFig1.jpg':
% export_figure(gcf, 'C:\Documents', 'myFig1', '.jpg', 'High Quality');

% Author:
% Avgoustinos Vouros
% avouros1@sheffield.ac.uk

function export_figure(get_current_figure, get_path, get_name, get_export_format, get_export_properties)
    
    % File name and extension
    file = strcat(get_name,get_export_format);
    fpath = fullfile(get_path,file);
    
    % If .fig format is requested just save and terminate
    if isequal(get_export_format,'.fig')
        saveas(get_current_figure, fpath);
        return;
    end   
    
    % Arrange the export format
    switch get_export_properties
        case 'Low Quality'
            % Not many changes here
            switch get_export_format
                case '.pdf'
                    get_export_format = '.pdf';
                case '.svg'
                    get_export_format = '.svg';
                case '.eps'  
                    get_export_format = '.epsc2'; %eps level 2 colored
                case '.jpg'
                    get_export_format = '.jpeg';
                case '.png'
                    get_export_format = '.png';
                case '.tif'
                    get_export_format = '.tiff'; %compressed                      
                case '.bmp'      
                    get_export_format = '.bmp';
                otherwise
                    error(strcat(get_export_format,': option unavailable'));
            end
        otherwise
            % Prepare for 'print'
            switch get_export_format
                case '.pdf'
                    get_export_format = '-dpdf';
                case '.svg'
                    get_export_format = '-dsvg';
                case '.eps'  
                    get_export_format = '-depsc'; %eps level 3 colored
                case '.jpg'
                    get_export_format = '-djpeg';
                case '.png'
                    get_export_format = '-dpng';
                case '.tif'
                    get_export_format = '-dtiffn'; %compressed                      
                case '.bmp'      
                    get_export_format = '-dbmp';
                otherwise
                    error(strcat(get_export_format,': option unavailable'));
            end            
    end
    
    % Output file with extension
    file = strcat(get_name,get_export_format);
    fpath = fullfile(get_path,file);  
    
    % Arrange the output
    switch get_export_properties
        case 'Low Quality'
            saveas(get_current_figure, fpath); 
        case 'Medium Quality'
            % use 200dpi
            fpath = fullfile(get_path,get_name);
            print(get_current_figure, fpath, get_export_format,'-r200');            
        case 'High Quality'
            % use 300dpi
            fpath = fullfile(get_path,get_name);
            print(get_current_figure, fpath, get_export_format,'-r300');
        case 'Super Quality'
            % use 600dpi
            fpath = fullfile(get_path,get_name);
            print(get_current_figure, fpath, get_export_format,'-r600');
        case 'Screen Quality'
            % use screen resolution
            fpath = fullfile(get_path,get_name);
            print(get_current_figure, fpath, get_export_format,'-r0');            
        otherwise
            error(strcat(get_export_properties,': option unavailable'));
    end   
end
