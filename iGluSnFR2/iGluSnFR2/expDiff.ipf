#pragma TextEncoding = "MacRoman"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function expDiff(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = aHat *(1-exp(-x/tauRiseHat)) *exp(-x/tauFallHat) 
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = aHat
	//CurveFitDialog/ w[1] = tauRiseHat
	//CurveFitDialog/ w[2] = tauFallHat

	return w[0] *(1-exp(-x/w[1])) *exp(-x/w[2]) 
End
