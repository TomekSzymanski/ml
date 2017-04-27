function TH = linRegGradVecFull(X, Y, alpha = 0.01)
	CONVERGENCE_TRESHOLD=alpha/100;

	m=size(X)(1);
	n=size(X)(2);
	TH=zeros(n,1);
	thPrev=TH;
	
	diff=realmax;
	steps=0;
	costFn=zeros(n,1);
	while diff > CONVERGENCE_TRESHOLD
		costFn = sumErrorsTimesVar(m, n, thPrev,X,Y);
		TH = thPrev - (alpha/m)*costFn;
		diff=max(abs(TH-thPrev));
		thPrev=TH;
		steps++;
	end
if (steps < 2)
    warning(sprintf('gradient descent ended after %d iterations', steps))
  endif
end

% calculates sum of squared differences times x value at given Xidx
function nextApprox = sumErrorsTimesVar(m, n, thetas, X, Y)
	%disp(sprintf('thetas: %d, %d', thetas(1), thetas(2)))
	sum = X'* (linRegFn(thetas, X) - Y);
	nextApprox=sum;
end


% returns m*1
function Y = linRegFn(thetas, X)
	Y=X*thetas;
end


%!test
%! X=[1,1; 1, 2];
%! Y=[2;4];
%! assert (round(linRegGradVecFull(X, Y)) == [0; 2])

%!test
%! X=[1,1; 1,2; 1,3; 1,4];
%! Y=[3.9;7.1; 9.9; 12.1];
%! assert (round(linRegGradVecFull(X, Y)) == [1; 3])

%!test
%! inputSize = 1e3;
%! X=[ones(inputSize, 1) [1:inputSize]'];
%! Y=[X(:,2)] .*2 .- 0.5 .+ rand(inputSize, 1);
%! assert (round(linRegGradVecFull(X, Y, alpha = 0.000001)) == [0; 2])

% 2-dim linear regression
%!test
%! X=[1,1,1; 1,2,2; 1,3,3; 1,4,4];
%! Y=[3.9;7.1; 9.9; 12.1];
%! assert (round(100 * linRegGrad(X, Y)) == [134; 138; 138])