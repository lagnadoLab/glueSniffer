#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.






/////////////////////////////////////////////////
///////////////////////Expected vesicles///////
////////////////////////////////////////
function expVes(p,binSize)
wave p
variable binSize
variable i
variable ex = 0
for (i=1;i<dimsize(p,0);i+=1)
	ex+= i * p[i]
endfor

ex /= binSize
print "Bin Width is: " + num2str(binSize)
print "Expected vesicles per second is: " + num2str(ex)
ex*=binSize
if (ex <1)
variable hBino = - ( ex *log(ex)/log(2)) - ((1-ex) * log(1-ex)/log(2))
print "Bino Entropy is: " + num2str(hBino)
endif
end

/////////////////////////////////////////////////////////////
//////////////// Entropy //////////////////////////
////////////////////////////////////////////////////
function entropy(p)
wave p

if (sum(p) < .999999)
	print "ENTROPY error - Probability Vector doesn't sum to 1"
endif

variable nX = dimsize(p,0)
variable i
variable h = 0
for (i=0;i<nX;i+=1)
	h-=  p[i] * log(p[i])/log(2)
endfor

print "Entropy is: " + num2str(h)
end


//////////////////////////////////////////////////////////////////
/////////////////////////// Mutual Information (2 Variables)//////////
////////////////////////////////////////////////////////////////////

function mutualInfo(pMat2Name)
string pMat2Name

string outName = pMat2Name + "_MI"
string outPR = pmat2Name + "_pR"
string outPs = pmat2Name + "_pS"

duplicate/o $pMat2Name, pMat2
/// I(X;Y) = Sum xInX sum yIny p(x,y) log( p(x,y)/p(x)p(y))

///Compute Marginal Densities
matrixop/o pS = sumCols(pMat2)
matrixop/o pR = sumRows(pmat2)

//find n(s) and n(r)
variable nS = dimsize(pS,1)
variable nR = dimsize(pR,0)

// Check it is a true probability distribution
if (abs(1-sum(pS)) >= .0001)
	print "Marginal p(S) doesn't sum to 1!"
elseif (abs(1-sum(pR))>=.001)
	print "Marginal P(R) doesn't sum to 1!"
elseif (abs(1-sum(pMat2))>=.001)
	print "Joint doesn't sum to 1!"
endif


variable r,s
variable mutInf = 0

for (r=0;r<nR;r+=1)
	for (s=0;s<nS;s+=1)
		if (pMat2[r][s] != 0 && pR[r] !=0)
		
		mutInf -= pMat2[r][s] * log(pS[s] * pR[r] / pMat2[r][s])/log(2)
		endif
	endfor
endfor

make/o/n=1 mutInfV
mutInfV = mutInf
duplicate/o mutInfV, $outName

duplicate/o ps, $outPs
duplicate/o pr, $outPr

print "Mutual information for " + pMat2Name + " is: " + num2str(mutInf)
killwaves/z pS,pR, mutInfV
end




////////////////////////////////////////////////////////////
///////////// Specific Information (2 Variables) /////////////
/////////////////////////////////////////////////////////
	
function specInf(pMat2Name)
string pMat2Name


duplicate/o $pmat2name, pMat2

matrixop/o sumSs = sumCols(pMat2)
matrixop/o sumRs = sumrows(pMat2)
if (abs(1-sum(sumSs)) >= .0001)
	print "Marginal p(S) doesn't sum to 1!"
elseif (abs(1-sum(sumRs))>=.001)
	print "Marginal P(R) doesn't sum to 1!"
elseif (abs(1-sum(pMat2))>=.001)
	print "Joint doesn't sum to 1!"
endif





// Number of stimuli and quantal letters (including no event)
variable nQ = dimsize(pmat2,0)
variable nS = dimsize(pmat2,1)


variable quanta, stim
make/o/n=(nQ) i2

i2 = 0
matrixop/o pR = sumRows(pMat2)
variable sumR = sum(pr)
pR/=sumR

duplicate/o pMat2, conditionOnR
variable ps = 1/nS
variable i,j

for (i=0;i<nQ;i+=1)
	for (j=0;j<nS;j+=1)
	
		if (pR[i] >0)
			conditionOnR[i][j] = pMat2[i][j]/pR[i]
		
		endIf
	endfor
endfor
	
for (i=0;i<nQ;i+=1)
	for (j=0;j<nS;j+=1)	
	if (conditionOnR[i][j]>0)
		i2[i] += conditionOnR[i][j] *((log(conditionOnR[i][j]))/log(2))
	endif
	endfor
endfor
variable minPart = 0

for (j=0;j<nS;j+=1)
	minPart+= ps * log(ps)/log(2)
endfor


i2-=minPart
string i2Out = pmat2Name + "_i2"
string i2OutPR = pmat2name + "_i2_pR"


duplicate/o i2, $i2Out
duplicate/o pr, $i2Outpr

killwAVES/Z i2,prps,minpart,conditionOnR,pr,ps
end
