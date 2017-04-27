function TH = linReg(X, Y)
TH=pinv(X'*X)*X'*Y;