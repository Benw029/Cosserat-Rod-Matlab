function VisualizeRod(R)
% Parameters
[nTime, nPoints] = size(R);
seg_end = nPoints;    % single segment
r_disk   = 0.001;
r_height = 0.005;
opt_pairs = {'tipframe',1,'segframe',0,'baseframe',1,'projections',0,'baseplate',0};

% Extract last time step
R_last = R(nTime,:);   % 1 x nPoints cell

% Build nPoints x 16 matrix
T_last = zeros(nPoints,16);
for j = 1:nPoints
    Rj = R_last{j};
    T = eye(4);        % 4x4 homogeneous
    T(1:3,1:3) = Rj;
    % Flatten columnwise
    T_last(j,:) = reshape(T,1,16); 
end

% Now T_last is nPoints x 16, ready for draw_tdcr
draw_tdcr(T_last, seg_end, r_disk, r_height, opt_pairs{:});


end