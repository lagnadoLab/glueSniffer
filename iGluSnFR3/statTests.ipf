#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.



////////////// Various statistical tests - made as needed




///////////////// Poisson Test - KS ////////////////
function poissTest(data,conditions)
/// KS Test for poisson distribution. Computes Empirical CDF and estimated CDF, and
// runs a KS test for statistical difference. Plots both CDFs.

// INPUT : a sequence of integer-valued samples
// OUTPUTS : cdf - empirical CDF
//				 estCFD - estimated CDF (using MLE estimate of lambda)
wave data
variable conditions
variable maxDat = waveMax(data)
make/o/n=(maxDat+1) pdf, cdf, xS, estCDF, estPDF
pdf = 0
cdf = 0
xS = p

// calculate empirical pdfs and cdfs
variable nDats = dimsize(data,0)
variable i,j
for (i=0;i<nDats;i+=1)
	for (j=0;j<maxDat+1;j+=1)
		if (data[i]==j)
			pdf[j]+=1
		endif
	endfor
endfor
pdf/= nDats

for (i=0;i<maxDat+1;i+=1)
	cdf[i] = sum(pdf,0,i)
endfor

// Calculate estimated CDF
variable lambdaHat = mean(data)

for (i=0;i<maxDat+1;i+=1)
	estPDF[i] = exp(-lambdaHat) * lambdaHat^i / factorial(i)
	estCDF[i] = sum(estPDF,0,i)
endfor







//Display Stuff
if (conditions ==1)
StatsKSTest/q estCDF,cdf
display/k=1 cdf, estcdf
ModifyGraph rgb(estCDF)=(0,0,0)
Legend/C/N=text0/J/F=0/A=MC "\\s(cdf) Empirical CDF\r\\s(estCDF) Estimated CDF"
Label left "Cumulative Probability";DelayUpdate
Label bottom "Value"
else
statsksTest/q estCDF,cdf
endif


killwaves/z pdf, estPDF, estCDF, cdf
end











///////////// CHi2 for nAFC stimulus discrimination //////////////
function chi2nAFC(rightVec,nVec)
wave rightVec,nVec

variable nA = dimsize(rightVec,0)
variable nTot = sum(nVec)
variable expN = nTot/nA
duplicate/o nVec, expVec
expVec = expN

StatsChiTest  rightVec,expVec

end




