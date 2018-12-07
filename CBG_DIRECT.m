function P = CBG_DIRECT(Xlim, Ylim, Nsamples)
% CBG_DIRECT - returns the points equidistant within the area defined by
% Xlim and Ylim.
%
% Syntax:  P = CBG_DIRECT(Xlim, Ylim, Nsamples)
%
% Inputs:
%    Xlim - Limits on the x-axis
%    Ylim - Limits on the y-axis
%
% Outputs:
%    P - Matrix [3 x 2 x NReality]. Contains the information of square
%    center, inferior and superior limits of the new squares. The 3rd
%    dimension shows the number of newly generated points in the area.
%
% Example: 
%    To-do
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: To-do

%------------- BEGIN CODE --------------
n_xdim = ceil(sqrt(Nsamples));
n_ydim = ceil(Nsamples/n_xdim);
NReality = n_xdim * n_ydim;

DeltaX = ( Xlim(2) - Xlim(1) ) / n_xdim;
DeltaY = ( Ylim(2) - Ylim(1) ) / n_ydim;

limits_x = Xlim(1):DeltaX:Xlim(2);
limits_y = Ylim(1):DeltaY:Ylim(2);

centers_x = limits_x(2:end) - diff(limits_x)/2;
centers_y = limits_y(2:end) - diff(limits_y)/2;

[C_x,C_y] = meshgrid(centers_x,centers_y);  % Centers
[Li_x,Li_y] = meshgrid(limits_x(1:end-1),limits_y(1:end-1));  % Inferior limits
[Ls_x,Ls_y] = meshgrid(limits_x(2:end),limits_y(2:end));  % Superior limits

% Store results in matrix for each new possibility. Each point is
% characterized by (1) centers, (2) inferior limits and (3) superior
% limits.
P = zeros(3,2,NReality); 
for i = 1:NReality
    P(1,:,i) = [C_x(i) C_y(i)];  % Center square
    P(2,:,i) = [Li_x(i) Ls_x(i)];  % Limits X-dimension
    P(3,:,i) = [Li_y(i) Ls_y(i)];  % Limits Y-dimension
end


% EOF