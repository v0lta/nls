%function  [Q0,evl,NSub] =  computeBase(hom, A0, unstable_flag, NSub)
%
% Compute an orthonormal basis Q0 for the invariant subspace
% corresponding to the NSub 
%       most unstable eigenvalues of A0 if unstable_flag = 1
%   and most stable eigenvalues of A0 if unstable_flag = 0.
%
% If resizeable_flag is true, then the size of the space NSub may
% be adjusted. 

function  [Q0,evl,NSub] =  computeBase(A0, unstable_flag, NSub)

if unstable_flag
    
    % Make all eigenvalues positive
    evl = eig(A0);
    if min(evl) < 0
        A0 = A0 + eye(size(A0)) * (abs(min(evl))+1);
    end
    
    % Compute eigenvalues and -vectors, ordered
    % (the ones with largest norm are first)
    [VU, DU] = eigs(A0,size(A0,1),'LR');
    % Select first NSub eigenvectors: unstable eigenspace
    VU = VU(:,1:NSub);
    
    % Compute orthonormal basis for the eigenspace
    [Q0,RU] = qr(VU);
    
else
    % Make all eigenvalues negative
    evl = eig(A0);     
    if max(evl) > 0
        A0 = A0 - eye(size(A0)) * (max(evl)+1);
    end
    
    % Compute eigenvalues and -vectors, ordered
    % (the ones with largest norm are first)
    [VS, DS] = eigs(A0,size(A0,1),'LR');
    % Select first NSub eigenvectors: stable eigenspace
    VS = VS(:,1:NSub);
        
    % Compute orthonormal basis for the eigenspace
    [Q0,RS] = qr(VS);
end
