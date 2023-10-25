function ShapeDriver()
% do all the required tests in the pdf file.
%
%%----------------for open shapes---------------------------
% *) Use "RandOpenShapes" to generate 40 noisy shapes for the above four
%    functions (10 for each)
% *) Compute the pairwise distances between the 40 shapes
% *) plot the image by 'imagesc' and 'colorbar'
%
% *) Use "classifyOpenCurves" to generate another 100 noisy shapes for
%    above four functions (25 for each) and classify the 100 shapes,
%    compute the percent of correctness.
% *) output the percent
% 
%%----------------for closed shapes--------------------------
% *) load the data
% *) use "classifyClosedCurves" to compute the percent of correctness
% *) output the percent
end

%%----------------for both open and closed shapes-----------
function outCs = centerandrescaling(Cs)
% remove translation and rescaling so that the curves have
% center at the origin and F-norms are one.
%
% input: 
%        Cs: a 2 by n by m array, where n is the number of points on a
%        curve and m is the number of curves
%
% output: 
%        outCs is also a 2 by n by m array.
%
end

%%----------------for open shapes---------------------------
function Cs = RandOpenShapes(myfun, num)
% Randomly generate a few functions with noise
% input: 
%       myfun: a function handle
%       num: the number of noisy curves
% output:
%        Cs: a 2 by n by num array, where n is the number of points in the
%            array 0:0.02:pi, and it is also the number of points on the
%            curves; num is the number of noisy curves. In this project,
%            num = 10
end

function distO = distOpenShapes(C1, C2)
% Compute distance between open shapes
% input:
%      C1: a 2 by n matrix, which represents an open curve
%      C2: a 2 by n matrix, which represents an open curve
% output:
%      distO: the distance between the open shapes of C1 and C2
end

function percent = classifyOpenCurves(trainingdata, fhandles, num)
% generate another 100 noisy shapes for above four functions (25 for each) 
% and classify the 100 shapes, output the percent of correctness.
%
% input:
%      trainingdata: a 2 by n by m array, where n is the number of points
%                    in the curve, m is the number of training curves. In
%                    this project, m = 40.
%      fhandles: a cell of the four function handles, cos, sin, x^2, and x
%      num: the number of curves in each trainingdata. In this project, 
%           num = 10
% output:
%      the percent of correctness
%
end

%%----------------for closed shapes--------------------------
function distC = distClosedShapes(C1, C2)
% Compute distance between closed shapes, the function distOpenShapes can
% be used here.
% input:
%      C1: a 2 by n matrix, which represents a closed curve
%      C2: a 2 by n matrix, which represents a closed curve
% output:
%      distO: the distance between the closed shapes of C1 and C2
end

function percent = classifyClosedCurves(trainingdata, testdata, clusternumtrain, clusternumtest)
% classify the given closed shapes
% input: 
%      trainingdata: a 2 by n by numtrain array, where n is 100 and
%                    numtrain is 300 in this project.
%      testdata: a 2 by n by numtest array, where n is 100 and numtest is
%                100 in this project.
%      clusternumtrain: the number of shapes in each cluster of the training 
%                       data, in this project it is 15
%      clusternumtest: the number of shapes in each clusterof the test 
%                       data, in this project it is 5
% output:
%      the percent of correctness
end

%%----------------for the bonus question--------------------------
function clusterclosedcurves(M)
% you need to figure this out by yourself.
end

% you can add more functions if you want.