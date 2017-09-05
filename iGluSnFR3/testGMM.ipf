#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.





Function fitGauss(wavename,muHat,aName)
// This is the code to fit guassian CURVES (not distributions) to the 
// temporal average of iGluSnFR linescans
// This code uses Igor's built-in curve-fitting software. If it doesn't work, don't blame me.
string waveName,muHat, aName
string coefName = (waveName[0,strlen(waveName)-8])+"Coef"
string fitName =(waveName[0,strlen(waveName)-8])+"fit"
string spatFilt = (waveName[0,strlen(waveName)-8])+"SF"
string sigmaHatName = (waveName[0,strlen(waveName)-8])+"sigmaHat"
duplicate/o $muHat, tmpMu
duplicate/o $aName, tmpA
duplicate/o $waveName, dat
variable k = dimsize(tmpMu,0)

make/o/n=(3*k) W_coef
// Attempt to make the convergence paramters better below...doesn't really work
Variable/g  V_FitMaxIters =1000
Variable/g V_FitTol = 0.000011
Variable/g V_FitOPtions = 1
W_Coef = 10
variable muBand = 2

//For each possible number of gaussians to fit, initialize values (means by cursors)
// and constraints
if (k==1)
	W_Coef[1]= tmpMu[0]
	w_coef[0] = tmpA[0]
	Make/o/T/N=3 T_constraints
	print tmpa[0]
	string c1 = "k1<"+num2str((tmpMu[0]+muBand))
	string c2 = "k1>"+num2str((tmpMu[0]-muBand))
	T_Constraints[0] = {"K0>0",c1,c2}
	FUncFit/H="000" gauss1 W_Coef DAT /D /C=T_Constraints
	makeData(W_Coef)
elseif (k==2)
	W_Coef[0] = tmpA[0]
	W_coef[1] = tmpA[1]
	W_Coef[2] = tmpMu[0]
	W_coef[3] = tmpMu[1]
	Make/O/T/N=6 T_Constraints
	 c1 = "k2<"+num2str((tmpMu[0]+muBand))
	 c2 = "k2>"+num2str((tmpMu[0]-muBand))
	string c3 = "k3<"+num2str((tmpMu[1]+muBand))
	string c4 = "k3>"+num2str((tmpMu[1]-muBand))
	T_Constraints[0] = {"K0 > 0","K1 > 0",c1,c2,c3,c4}
	FuncFit/H="000000" gauss2 W_coef DAT /D /C=T_Constraints 
	makeData(W_Coef)
elseif (k==3)
	W_coef[0] = tmpA[0]
	W_Coef[1] = tmpA[1]
	W_coef[2] = tmpA[2]
	W_coef[3] = tmpMu[0]
	W_coef[4] = tmpMu[1]
	W_coef[5] = tmpMu[2]
	c1 = "k3<"+num2str((tmpMu[0]+muBand))
	c2 = "k3>"+num2str((tmpMu[0]-muBand))
	c3 = "k4<"+num2str((tmpMu[1]+muBand))
	c4 = "k4>"+num2str((tmpMu[1]-muBand))
	string c5 = "k5<"+num2str((tmpMu[2]+muBand))
	string c6 = "k5>"+num2str((tmpMu[2]-muBand))
	Make/O/T/N=9 T_Constraints
	T_Constraints[0] = {"K0 > 0","K1 > 0","K2 > 0",c1,c2,c3,c4,c5,c6}
	FuncFit/H="000000000" gauss3 W_coef dat /D /C=T_Constraints 
	makedata(w_Coef)
elseif (k==4)
	W_coef[0] = tmpA[0]
	W_Coef[1] = tmpA[1]
	W_coef[2] = tmpA[2]
	W_Coef[3] = tmpA[3]
	W_coef[4] = tmpMu[0]
	W_coef[5] = tmpMu[1]
	W_coef[6] = tmpMu[2]
	W_Coef[7] = tmpMu[3]
	Make/O/T/N=12 T_Constraints
	c1 = "k4<"+num2str((tmpMu[0]+muBand))
	c2 = "k4>"+num2str((tmpMu[0]-muBand))
	c3 = "k5<"+num2str((tmpMu[1]+muBand))
	c4 = "k5>"+num2str((tmpMu[1]-muBand))
	c5 = "k6<"+num2str((tmpMu[2]+muBand))
	c6 = "k6>"+num2str((tmpMu[2]-muBand))
	string c7 = "k7>"+num2str((tmpMu[3]-muBand))
	string c8 = "k7<"+num2str((tmpMu[3]+muBand))
	T_Constraints[0] = {"K0 > 0","K1 > 0","K2 > 0","K3>0",c1,c2,c3,c4,c5,c6,c7,c8}
	FuncFit/H="000000000000" gauss4 W_coef dat /D /C=T_Constraints 
	makedata(W_Coef)
elseif (k==5)
	W_coef[0] = tmpA[0]
	W_Coef[1] = tmpA[1]
	W_coef[2] = tmpA[2]
	W_Coef[3] = tmpA[3]
	W_coef[4] = tmpA[4]
	W_coef[5] = tmpMu[0]
	W_coef[6] = tmpMu[1]
	W_coef[7] = tmpMu[2]
	w_Coef[8] = tmpMu[3]
	w_coef[9] = tmpMu[4]
	Make/O/T/N=15 T_Constraints
	c1 = "k5<"+num2str((tmpMu[0]+muBand))
	c2 = "k5>"+num2str((tmpMu[0]-muBand))
	c3 = "k6<"+num2str((tmpMu[1]+muBand))
	c4 = "k6>"+num2str((tmpMu[1]-muBand))
	c5 = "k7<"+num2str((tmpMu[2]+muBand))
	c6 = "k7>"+num2str((tmpMu[2]-muBand))
	c7 = "k8>"+num2str((tmpMu[3]-muBand))
	c8 = "k8<"+num2str((tmpMu[3]+muBand))
	string c9 = "k9>" +num2str(tmpMu[4]-muBand)
	string c10 = "k9<" +num2str(tmpMu[4]+muBand)
	T_Constraints[0] = {"K0 > 0","K1 > 0","K2 > 0","K3>0","K4>0",c1,c2,c3,c4,c5,c6,c7,c8,c9,c10}
	FuncFit/H="000000000000000" gauss5 W_coef dat /D /C=T_Constraints 
	makeData(W_coef)
elseif (k==6)
	w_coef[6] = tmpMu[0]
	W_coef[7] = tmpMu[1]
	W_coef[8] = tmpMu[2]
	W_coef[9] = tmpMu[3]
	w_Coef[10] = tmpMu[4]
	w_coef[11] = tmpMu[5]
	W_coef[0] = tmpA[0]
	W_Coef[1] = tmpA[1]
	W_coef[2] = tmpA[2]
	W_Coef[3] = tmpA[3]
	W_coef[4] = tmpA[4]
	W_coef[5] = tmpA[5]
	c1 = "k6<"+num2str((tmpMu[0]+muBand))
	c2 = "k6>"+num2str((tmpMu[0]-muBand))
	c3 = "k7<"+num2str((tmpMu[1]+muBand))
	c4 = "k7>"+num2str((tmpMu[1]-muBand))
	c5 = "k8<"+num2str((tmpMu[2]+muBand))
	c6 = "k8>"+num2str((tmpMu[2]-muBand))
	c7 = "k9>"+num2str((tmpMu[3]-muBand))
	c8 = "k9<"+num2str((tmpMu[3]+muBand))
	c9 = "k10>" +num2str(tmpMu[4]-muBand)
	c10 = "k10<" +num2str(tmpMu[4]+muBand)
	string c11= "k11>"+num2str(tmpMu[5]-muBand)
	string c12= "k11<"+num2str(tmpMu[5]+muBand)
	Make/O/T/N=18 T_Constraints
	T_Constraints[0] = {"K0 > 0","K1 > 0","K2 > 0","K3>0","K4>0","K5>0",c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12}
	FuncFit/H="000000000000000000" gauss6 W_coef dat /D /C=T_Constraints 
	
endif
wave fit_dat, sumDat, xS, data, sigma2, mu
duplicate/o W_Coef, $CoefName
duplicate/o fit_dat, $fitName
duplicate/o data, $spatFilt
duplicate/o sigma2, $sigmaHatName
killwaves W_Coef, dat, fit_dat, sumDat, tmpMu, tmpA, xS, sigma2, mu
killwaves/z T_Constraints, data, A

end
	
	




Function makeData(w_Coef)
//Initialize some shit
wave W_Coef
variable k = dimsize(W_Coef,0)/3
variable xLen = 128
make/o/n=(xLen) xS
variable i
for (i=0;i<xLen;i+=1)
	xS[i] = i
endfor
make/o/n=(k) mu, sigma2,A
if (k==1)
	A = {W_Coef[0]}
	mu = {W_Coef[1]}
	sigma2 = {w_Coef[2]}
elseif (k==3)
	mu = {w_Coef[3],W_Coef[4],W_Coef[5]}
	sigma2 = {w_Coef[6],W_Coef[7],W_Coef[8]}
	A = {w_Coef[0],W_Coef[1],W_Coef[2]}
elseif (k==2)
	A = {W_Coef[0],W_Coef[1]}
	mu = {W_Coef[2],W_Coef[3]}
	sigma2 = {w_Coef[4],W_Coef[5]}
elseif (k==4)
	mu = {w_Coef[4],W_Coef[5],W_Coef[6],W_Coef[7]}
	sigma2 = {w_Coef[8],W_Coef[9],W_Coef[10],W_Coef[11]}
	A = {w_Coef[0],W_Coef[1],W_Coef[2],W_Coef[3]}
elseif (k==5)
	mu = {W_Coef[5],W_Coef[6],W_Coef[7],w_Coef[8],w_Coef[9]}
	sigma2 = {w_Coef[10],W_Coef[11],W_Coef[12],W_Coef[13],W_Coef[14]}
	A = {w_Coef[0],W_Coef[1],W_Coef[2],W_Coef[3],W_Coef[4]}
elseif (k==6)
	mu = {W_Coef[6],W_Coef[7],w_Coef[8],w_Coef[9],W_Coef[11],W_Coef[12]}
	sigma2 = {w_Coef[13],W_Coef[14],W_Coef[15],W_Coef[16],W_Coef[17],W_Coef[18]}
	A = {w_Coef[0],W_Coef[1],W_Coef[2],W_Coef[3],W_Coef[4],W_Coef[5]}

endif



normpdf(xS,mu,sigma2,A)
end



Function normpdf(xS,mu,sigma2,A)
// Computes the pdf of the gaussian mixture 
wave xS, mu, sigma2, A
if (dimsize(mu,0) >0)
	make/o/n=(dimsize(xS,0),dimsize(mu,0)) data
else 
	make/o/n=(dimsize(xS,0)) data
endif
duplicate/o xS, sumDat
sumDat = 0
variable xLen = dimsize(xS,0)
variable k = dimsize(mu,0)
variable i,j
for (i=0;i<xLen;i+=1)
	for (j=0;j<k;j+=1)
		data[i][j] = A[j]*(2 * pi * sigma2[j])^-.5 * exp(-(xS[i]-mu[j])^2/(2*sigma2[j]))
		sumDat[i] = sumDat[i] + data[i][j]
	endfor
endfor

end

/////////////////////////////// Fit Functions //////////////////////
//All these functions are just gaussian mixtures with one or more components
// For each gaussian, A, Mu, Sigma2. 
Function gauss1(w,x) : FitFunc
	Wave w
	Variable x

	Variable/g V_FitNumIts

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = (A1*(2*pi*sigmaS1)^-.5 * exp(-(x-mu1)^2/(2*sigmaS1)))
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = A1
	//CurveFitDialog/ w[1] = mu1
	//CurveFitDialog/ w[2] = sigmaS1


	return (w[0]*(2*pi*w[2])^-.5 * exp(-(x-w[1])^2/(2*w[2])))

End





Function gauss2(w,x) : FitFunc
	Wave w
	Variable x

	Variable/g V_FitNumIts

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = (A1*(2*pi*sigmaS1)^-.5 * exp(-(x-mu1)^2/(2*sigmaS1))) + (A2 * (2*pi*sigmaS2)^-.5 * exp(-(x-mu2)^2/(2*sigmaS2)))
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ w[0] = A1
	//CurveFitDialog/ w[1] = A2
	//CurveFitDialog/ w[2] = mu1
	//CurveFitDialog/ w[3] = mu2
	//CurveFitDialog/ w[4] = sigmaS1
	//CurveFitDialog/ w[5] = sigmaS2

	return (w[0]*(2*pi*w[4])^-.5 * exp(-(x-w[2])^2/(2*w[4]))) + (w[1] * (2*pi*w[5])^-.5 * exp(-(x-w[3])^2/(2*w[5])))

End


Function gauss3(w,x) : FitFunc
	Wave w
	Variable x
	variable/g V_FitMaxIters
	variable/g v_FitTol
	Variable/g V_FitNumIts

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = (A1*(2*pi*sigma21)^-.5 * exp(-(x-mu1)^2/(2*sigma21))) + (A2 * (2*pi*sigma22)^-.5 * exp(-(x-mu2)^2/(2*sigma22)))  + (A3 * (2*pi*sigma23)^-.5 * exp(-(x-mu3)^2/(2*sigma23)))
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 9
	//CurveFitDialog/ w[0] = A1
	//CurveFitDialog/ w[1] = A2
	//CurveFitDialog/ w[2] = A3
	//CurveFitDialog/ w[3] = mu1
	//CurveFitDialog/ w[4] = mu2
	//CurveFitDialog/ w[5] = mu3
	//CurveFitDialog/ w[6] = sigma21
	//CurveFitDialog/ w[7] = sigma22
	//CurveFitDialog/ w[8] = sigma23

	return (w[0]*(2*pi*w[6])^-.5 * exp(-(x-w[3])^2/(2*w[6]))) + (w[1] * (2*pi*w[7])^-.5 * exp(-(x-w[4])^2/(2*w[7])))  + (w[2] * (2*pi*w[8])^-.5 * exp(-(x-w[5])^2/(2*w[8])))
	
	variable V_FitNumIters
//print V_FitNumIters
End

Function gauss4(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) =  (A1 * (2 * pi * sigmaS1)^-.5 * exp(-(x-mu1)^2/(2*sigmaS1))) +  A2 * (2 * pi * sigmaS2)^-.5 * exp(-(x-mu2)^2/(2*sigmaS2)) +  A3 * (2 * pi * sigmaS3)^-.5 * exp(-(x-mu3)^2/(2*sigmaS3))  +  A4 * (2 * pi * sigmaS4)^-.5 * exp(-(x-mu4)^2/(2*sigmaS4)) 
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 12
	//CurveFitDialog/ w[0] = A1
	//CurveFitDialog/ w[1] = A2
	//CurveFitDialog/ w[2] = A3
	//CurveFitDialog/ w[3] = A4
	//CurveFitDialog/ w[4] = mu1
	//CurveFitDialog/ w[5] = mu2
	//CurveFitDialog/ w[6] = mu3
	//CurveFitDialog/ w[7] = mu4
	//CurveFitDialog/ w[8] = sigmaS1
	//CurveFitDialog/ w[9] = sigmaS2
	//CurveFitDialog/ w[10] = sigmaS3
	//CurveFitDialog/ w[11] = sigmaS4

	return (w[0] * (2 * pi * w[8])^-.5 * exp(-(x-w[4])^2/(2*w[8]))) +  w[1] * (2 * pi * w[9])^-.5 * exp(-(x-w[5])^2/(2*w[9])) +  w[2] * (2 * pi * w[10])^-.5 * exp(-(x-w[6])^2/(2*w[10]))  +  w[3] * (2 * pi * w[11])^-.5 * exp(-(x-w[7])^2/(2*w[11])) 
End

Function gauss5(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = (A1 * (2 * pi * sigmaS1)^-.5 * exp(-(x-mu1)^2/(2*sigmaS1))) +  A2 * (2 * pi * sigmaS2)^-.5 * exp(-(x-mu2)^2/(2*sigmaS2)) +  A3 * (2 * pi * sigmaS3)^-.5 * exp(-(x-mu3)^2/(2*sigmaS3))  +  A4 * (2 * pi * sigmaS4)^-.5 * exp(-(x-mu4)^2/(2*sigmaS4)) +  A5 * (2 * pi * sigmaS5)^-.5 * exp(-(x-mu5)^2/(2*sigmaS5)) 
	//CurveFitDialog/ 
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 15
	//CurveFitDialog/ w[0] = A1
	//CurveFitDialog/ w[1] = A2
	//CurveFitDialog/ w[2] = A3
	//CurveFitDialog/ w[3] = A4
	//CurveFitDialog/ w[4] = A5
	//CurveFitDialog/ w[5] = mu1
	//CurveFitDialog/ w[6] = mu2
	//CurveFitDialog/ w[7] = mu3
	//CurveFitDialog/ w[8] = mu4
	//CurveFitDialog/ w[9] = mu5
	//CurveFitDialog/ w[10] = sigmaS1
	//CurveFitDialog/ w[11] = sigmaS2
	//CurveFitDialog/ w[12] = sigmaS3
	//CurveFitDialog/ w[13] = sigmaS4
	//CurveFitDialog/ w[14] = sigmaS5

	return (w[0] * (2 * pi * w[10])^-.5 * exp(-(x-w[5])^2/(2*w[10]))) +  w[1] * (2 * pi * w[11])^-.5 * exp(-(x-w[6])^2/(2*w[11])) +  w[2] * (2 * pi * w[12])^-.5 * exp(-(x-w[7])^2/(2*w[12]))  +  w[3] * (2 * pi * w[13])^-.5 * exp(-(x-w[8])^2/(2*w[13])) +  w[4] * (2 * pi * w[14])^-.5 * exp(-(x-w[9])^2/(2*w[14])) 
	
End


