% generate uniformly distributed random vectors,e satifying
% sum(e) = Etotal
% 0 <= e <= d
% Generate only feasible points using null space bounding box calculation
% Jon Simon 2020-03-04
numNeeded = 10000; % number of random vectors that are needed
Etotal = 7;
d = [1; 0.5; 1; 0.9; 0.25; 1; 0.5; 2]; % upper bound on components
N = length(d); % number of components
% make sure that it is possible to find solutions
if sum(d) < Etotal
    error('No solution possible, try increasing components of d')
end
% find particular solution in the interior of the bounding box
% define vector that extends from origin to maximal corner of bounding box
v = d;
% find intersection of vector with hyperplane where sum(e) = Etotal
e0 = Etotal/sum(d)*v;
% find null space of solutions
A = null(ones(1,N));
% translate coordinates of bounding box
% to a frame where origin is at e0
dU = d - e0;
dL = -e0; % original lower bound is zero
Bbnd = [dL dU]; % box bounds in translated coordinates
% generate random feasible point that are on the hyperplane and in the
% interior of the bounding box
E = zeros(N,numNeeded); % preallocate array to hold vectors
for k = 1:numNeeded

    % find random vector in null space
    w = 2*rand(N-1,1)-1; % ranges from -1 to 1
    e = A*w;

    % random point will be found in direction of e, but need to determine
    % allowable scaling so that it will be inside of the bounding box

    % for each component of e find matrix of scaling factors Alpha, where
    % Alpha(i,j)*e = Bbnd(i,j); where Alpha(i,1)  is the scale factor
    % that would bring the ith component of Alpha(i,1)*e to the ith lower bound
    % and Alpha(i,2) is the scale factor that would bring the ith component
    % of Alpha(i,2)*e to the ith upper bound
    Alpha = Bbnd./e;

    % find the allowable range of alpha, note that since 0*e is a point that is
    % on the hyperplane and in the interior of the bounding box, alpha = 0 is always in
    % allowable range so the lower bound on alpha is a negative value and
    % the upper bound is a positive value
    alphaL = max(Alpha(Alpha(:)<=0)); % greatest lower bound on Alpha, use (:) to look at elements as column
    alphaU = min(Alpha(Alpha(:)>=0)); % least upper bound on Alpha

    % randomly choose a value of the scale factor that is within the
    % allowable range
    alpha = alphaL + rand*(alphaU-alphaL);

    % generate point on hyperplane that is inside of bounding box
    E(:,k) = e0 + alpha*e ; % particular solution plus null space vector;

end
% check
colSums = sum(E) % gives column sums, should all equal Etotal
minComp = min(E,[],2) % gives minimum of each component over all of the points in E, should all be >= zero
maxComp = max(E,[],2) % gives maximum component over all of the points in E, should all be <= d
% display sum(d)-Etotal, which gives indication of difficulty of proble
disp({'sum(d)-Etotal = ',sum(d)-Etotal})
% if successful all of these should be true
check1 = all(abs(colSums -Etotal) < 1e-10) % allow some tolerance
check2 = all(minComp>=0)
check3 = all(maxComp<=d)
