function ResetPlotOptions
% Reset Matlab plot options to default values
close all
set( 0, 'DefaultAxesFontSize', 12 )
set( 0, 'DefaultAxesFontName',  'Helvetica' )
set( 0, 'DefaultAxesFontAngle', 'normal' )
gridState = 'off' ;
set( 0, 'DefaultAxesXGrid', gridState )
set( 0, 'DefaultAxesYGrid', gridState )
set( 0, 'DefaultAxesZGrid', gridState )
set( 0, 'DefaultTextFontSize', 12 )
set( 0, 'DefaultTextFontName', 'Helvetica' )
set(0,  'DefaultAxesColor',   'w', ...
        'DefaultAxesXColor',  'k', ...
        'DefaultAxesYColor',  'k', ...
        'DefaultAxesZColor',  'k', ...
        'DefaultTextColor',   'black', ...
        'DefaultLineColor',   'black' )
set( 0, 'DefaultLineLineWidth', 1 )
set( 0, 'DefaultUicontrolBackgroundColor', 'w' )
set( 0, 'DefaultAxesBox', 'off') % doesn't work