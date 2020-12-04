function [w, SNR, repeat] = CBG_sdp_solver(channel, M, N, Pt_max, sinrmin, sigma2, antenna_allocation)
% CBG_sdp_solver - The beamforming algorithm that provides the optimum
%                  power allocation as to minimize the transmit power while 
%                  meeting the SNR demands
%
% Syntax:  [w, SNR] = CBG_sdp_solver(H_abs, M, N, Pt_max, ...
%                                    sinrmin, sigma2, antenna_allocation)
%
% Inputs:
%    H_abs - [N x M] channel matrix
%    M - Number of transmitters
%    N - Number of receivers
%    Pt_max - Maximum transmit power in linear scale
%    sinrmin - [N x 1] vector with requested SNR
%    sigma2 - [N x 1] noise vector
%    antenna_allocation - [M x N] binary matrix revealing the antenna
%                         assignation across users
%
% Outputs:
%    w - [1 x N] weights
%    SNR - [1 x N] achieved SNR in linear scale
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

% Compute basics channel
H_abs = abs(channel);  % Absolute gain
Hrici_angle = angle(channel);  % Angle

% Compute Autocorrelation Matrices for each user
R_user = zeros(M,M,N);
for rxID = 1:N
    R_user(:,:,rxID) = H_abs(rxID,:)'*H_abs(rxID,:);
end

X_user1 = zeros(M,M);  % Inizialize variable
X_user2 = zeros(M,M);  % Inizialize variable
X_user3 = zeros(M,M);  % Inizialize variable

% Creating association matrices
association_A = zeros(M,M,N);
association_B = ones(M,M,N);
for i=1:M
    for j=1:N
        if(antenna_allocation(i,j)==1)
          association_A(i,i,j)=1;
        elseif (antenna_allocation(i,j)==0)
          association_B(:,i,j)=0;   
          association_B(i,:,j)=0; 
        end
    end
end

cvx_begin quiet

    variables X_user1(M,M) X_user2(M,M) X_user3(M,M)

    % Objective function - Power minimization
    minimize(  trace(association_A(:,:,1).*X_user1) + ...
               trace(association_A(:,:,2).*X_user2) + ...
               trace(association_A(:,:,3).*X_user3) )

    % Constraints of the optimization problem
    subject to
    ( ...
    diag(association_A(:,:,1).*X_user1) + ...
    diag(association_A(:,:,2).*X_user2) + ... 
    diag(association_A(:,:,3).*X_user3)   ...
    ) ...
    <= Pt_max*ones(M,1);  %#ok. Constraint on total transmitted power per antenna 
    trace(R_user(:,:,1)*(association_B(:,:,1).*X_user1)) - ...
          (  ...
          (   sinrmin(1) * trace( R_user(:,:,1) * (association_B(:,:,2).*X_user2) )   ) + ...
          (   sinrmin(1) * trace( R_user(:,:,1) * (association_B(:,:,3).*X_user3) )   )   ...
          ) ... 
          >= sigma2(1)*sinrmin(1);  %#ok. Contraint on minimum SINR for user 1
    trace(R_user(:,:,2)*(association_B(:,:,2).*X_user2)) - ...
          (  ...
          (   sinrmin(2) * trace( R_user(:,:,2) * (association_B(:,:,1).*X_user1) )   ) + ...
          (   sinrmin(2) * trace( R_user(:,:,2) * (association_B(:,:,3).*X_user3) )   )   ...
          ) ...
          >= sigma2(2)*sinrmin(2);  %#ok. Contraint on minimum SINR for user 2
    trace(R_user(:,:,3)*(association_B(:,:,3).*X_user3)) - ...
          (  ...
          (   sinrmin(3) * trace( R_user(:,:,3) * (association_B(:,:,1).*X_user1) )   ) + ...
          (   sinrmin(3) * trace( R_user(:,:,3) * (association_B(:,:,2).*X_user2) )   )   ...
          ) ...
          >= sigma2(3)*sinrmin(3);  %#ok. Contraint on minimum SINR for user 3
    X_user1 == hermitian_semidefinite(M);  %#ok
    X_user2 == hermitian_semidefinite(M);  %#ok
    X_user3 == hermitian_semidefinite(M);  %#ok

cvx_end

X_user1 = X_user1.*association_B(:,:,1);
X_user2 = X_user2.*association_B(:,:,2);
X_user3 = X_user3.*association_B(:,:,3);

if any(isnan(X_user1(:))) || any(isnan(X_user2(:))) || any(isnan(X_user3(:)))
    repeat = true;
    w = [];
    SNR = [];
else
    repeat = false;
    SNR =  zeros(1,N);
    w_sdp = [];  % general output
    nPos = [];  % to check if rank 1 (everything alright)
    for rxID = 1:N
        if rxID == 1
            eval=diag(svd(X_user1));
            sol = sqrt(diag(X_user1));  % Weighs user 1
        elseif rxID == 2
            eval=diag(svd(X_user2));
            sol = sqrt(diag(X_user2));  % Weighs user 2
        elseif rxID == 3
            eval=diag(svd(X_user3));
            sol = sqrt(diag(X_user3));  % Weighs user 3
        end
        % Final weights
        w_sdp =[w_sdp , sol];  %#ok<AGROW>
        % To check unitary rank
        nPos = [nPos sum(eval(2:end)<10e-6)];  %#ok<AGROW>
        % Achieved SNR
        SNR(rxID) = ((H_abs(rxID,:)*sol)^2)/(((H_abs(rxID,:)*sol)^2)+sigma2(rxID));
    end

    % Checking whether the solution is rank 1
    if(length(unique(nPos))~=1)
       warning('Randomization required')  % For not robust formulation rank = 1 [1].
    end

    ang_w = -1.*Hrici_angle;  % Equalization
    w = w_sdp.*cos(ang_w')+1i.* w_sdp.*sin(ang_w'); 
end



% EOF