#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


/////////////////////// GMM using the EM algorithm
Function emGMM(waveName)
// This function runs an EM algorithm for gaussian mixture models. This sort of algorithm has
// been done to death, but in breath, it switches between estimating which gaussian cluster
// each data point originated from and estimating the means and variances of the components
// INPUT - input is the name of the wave corresponding to amplitude values extracted
// 	using Wiener deconvolution
// OUTPUTS - Depending on which partition you select in the GUI, the outputs are named by
// apending to the wavename _K, followed by whichever number of guassian clusters you ended up selecting
// as well as _muHat, _aHat, _sigmaHat, _clustMAP, and _simData.
// 	_muHat, _aHat, and _sigmaHat are k-length waves (where k is the number of components)
//    where each indice represent the parameter for the corresponding group
// 	_clustMAP is a vector assigning each data point to a cluster, and _simData is the wave
// 	plotting the mixtures of gaussians computed


string waveName
string/g dataName = (waveName[0,strlen(waveName)-1])

duplicate/o $wavename, data
variable nObj = dimsize(data,0)
variable maxDat = wavemax(data)
variable scale = maxDat/100
duplicate/o data, scaleData
scaleData /= scale // Set range to 0-100 to avoid Igor's shitty underflow problem


// Initialize cluster name strings
string k1Clust = dataName+"_k1"
string k2Clust = dataName+"_k2"
string k3Clust = dataName+"_k3"
string k4Clust = dataName+"_k4"
string k5Clust = dataName+"_k5"
string k6Clust = dataName+"_k6"
string k7Clust = dataName+"_k7"
string k8Clust = dataName+"_k8"
string k9Clust = dataName+"_k9"
string k10Clust = dataName+"_k10"
string k11Clust = dataName+"_k11"
string k12Clust = dataName+"_k12"
string k13Clust = dataName+"_k13"
string k14Clust = dataName+"_k14"


// Run algorithm on each possible number of clusters
variable logP2 = kEMgmm(scaleData,data,2,k2Clust,scale)
variable logP3 = kEMgmm(scaleData,data,3,k3Clust,scale)
variable logP4 = kEMgmm(scaleData,data,4,k4Clust,scale)
variable logP5 = kEMgmm(scaleData,data,5,k5Clust,scale)
variable logP6 = kEMgmm(scaleData,data,6,k6Clust,scale)
variable logP7 = kEMgmm(scaleData,data,7,k7Clust,scale)
variable logP8 = kEMgmm(scaleData,data,8,k8Clust,scale)
variable logP9 = kEMgmm(scaleData,data,9,k9Clust,scale)
variable logP10 = kEMgmm(scaleData,data,10,k10Clust,scale)
variable logP11 = kEMgmm(scaleData,data,11,k11Clust,scale)
variable logP12 = kEMgmm(scaleData,data,12,k12Clust,scale)
variable logP13 = kEMgmm(scaleData,data,13,k13Clust,scale)
variable logP14 = kEMgmm(scaleData,data,14,k14Clust,scale)


//Record log probabilities and some metrics for model selection
make/o/N=(7) kP
variable i
kP = {logP2,logP3,logP4,logP5,logP6,logP7,logP8,logP9,logP10,logP11,logP12,logP13,logP14}
duplicate/o kP, AIC, BIC
AIC = 0
BIC = 0
for (i=0;i<7;i+=1)	
	AIC[i] = -2 * kP[i] + 2 * (i+1) * 3
	BIC[i] = -2 * kP[i] + 3 * (i+1)* log(nObj)
endfor
string logPName = dataName +"_logP"
duplicate/o kp, $logPName


// Clean some stuff up, organize
killwaves kP, AIC, BIC

string s1 = k1Clust + "_simData"
string s2 = k2Clust + "_simData"
string s3 = k3Clust + "_simData"
string s4 = k4Clust + "_simData"
string s5 = k5Clust + "_simData"
string s6 = k6Clust + "_simData"
string s7 = k7Clust + "_simData"
string s8 = k8Clust + "_simData"
string s9 = k9Clust + "_simData"
string s10 = k10Clust + "_simData"
string s11 = k11Clust + "_simData"
string s12 = k12Clust + "_simData"
string s13 = k13Clust + "_simData"
string s14 = k14Clust + "_simData"

killwaves/z simMat
concatenate/NP=1 {$s2,$s3,$s4,$s5,$s6,$s7,$s8,$s9,$s10,$s11,$s12,$s13,$s14}, simMat
string/g simMatName = dataName +"simMat"
duplicate/o simMat, $simMatName
make/o/n=(100) ampHist
Histogram/P/B={0,maxDat/100,100} data, ampHist
string/g histName = dataName + "_AH"
duplicate/o ampHist, $histName
variable/g onClust = 0
gmmWindow()
end




//////////////////////////////
// This builds the GUI for selecting the results
Function gmmWindow()
string/g xSpace
string/g simMatName
variable/g onClust
string/g histName
duplicate/o $histName, aHist
duplicate/o $simMatName, cdata
duplicate/o $xSpace, xS


////Panel
duplicate/o/RMD=[][onClust] cdata, clustDat
Display/k=1/W=(100,0,800,400)/N=clusterThingy aHist as "Cluster GUI"
appendtograph clustDat vs xS
modifygraph margin(bottom)=100
ModifyGraph rgb(clustDat)=(8738,8738,8738)
Label bottom "\\Z14 Deconvolved Amplitude"
Label left "\\Z14 Count/Scaled Amplitude"
variable/g nQuantaShow
nQuantaShow = onClust + 2
///Val Display
ValDisplay cluster value = #"nQuantaShow"

//Button
Button nextClust size = {100,25}, proc = nextClustProc, title="Next"
Button chooseClust size={100,25}, proc = chooseClustProc, title = "Select"
end


//////////// Next button
Function nextClustProc(ba) : ButtonControl
struct WMButtonAction &ba
	switch(ba.eventcode)
		case 2:
			variable/g onClust
			dowindow/k clusterThingy
			if (onClust<13)
				onClust +=1
				gmmWindow()
			elseif (onClust==13)
				onClust=0
				gmmWindow()
			endif
		endswitch
	end




//////////// Save button
//Poorly written code that saves the selected partition selected, and deletes all others
Function chooseClustProc(ba) : ButtonControl
struct WMButtonAction &ba
	switch(ba.eventcode)
		case 2:
			variable/g onClust
			dowindow/k clusterThingy
			string/g dataname
			string k1Clust = dataName+"_k1"
			string k2Clust = dataName+"_k2"
			string k3Clust = dataName+"_k3"
			string k4Clust = dataName+"_k4"
			string k5Clust = dataName+"_k5"
			string k6Clust = dataName+"_k6"
			string k7Clust = dataName+"_k7"
			string k8Clust = dataName+"_k8"
			string k9Clust = dataName+"_k9"
			string k10Clust = dataName+"_k10"
			string k11Clust = dataName+"_k11"
			string k12Clust = dataName+"_k12"
			string k13Clust = dataName+"_k13"
			string k14Clust = dataName+"_k14"
			
			
			string sig1 = k1Clust+"_sigmaHat"
			string mu1 = k1Clust + "_muHat"
			string a1 = k1Clust + "_aHat"
			string map1 = k1Clust + "_clustMAP"
			string sig2 = k2Clust+"_sigmaHat"
			string mu2 = k2Clust + "_muHat"
			string a2 = k2Clust + "_aHat"
			string map2 = k2Clust + "_clustMAP"
			string sig3 = k3Clust+"_sigmaHat"
			string mu3 = k3Clust + "_muHat"
			string a3 = k3Clust + "_aHat"
			string map3 = k3Clust + "_clustMAP"
			string sig4 = k4Clust+"_sigmaHat"
			string mu4 = k4Clust + "_muHat"
			string a4 = k4Clust + "_aHat"
			string map4 = k4Clust + "_clustMAP"
			string sig5 = k5Clust+"_sigmaHat"
			string mu5 = k5Clust + "_muHat"
			string a5 = k5Clust + "_aHat"
			string map5 = k5Clust + "_clustMAP"
			string sig6 = k6Clust+"_sigmaHat"
			string mu6 = k6Clust + "_muHat"
			string a6 = k6Clust + "_aHat"
			string map6 = k6Clust + "_clustMAP"
			string sig7 = k7Clust+"_sigmaHat"
			string mu7 = k7Clust + "_muHat"
			string a7 = k7Clust + "_aHat"
			string sig8 = k8Clust+"_sigmaHat"
			string mu8 = k8Clust + "_muHat"
			string a8 = k8Clust + "_aHat"
			string sig9 = k9Clust+"_sigmaHat"
			string mu9 = k9Clust + "_muHat"
			string a9 = k9Clust + "_aHat"
			string sig10 = k10Clust+"_sigmaHat"
			string mu10 = k10Clust + "_muHat"
			string a10 = k10Clust + "_aHat"
			
			string sig11 = k11Clust+"_sigmaHat"
			string mu11 = k11Clust + "_muHat"
			string a11 = k1Clust + "_aHat"
			
			string sig12 = k12Clust+"_sigmaHat"
			string mu12 = k12Clust + "_muHat"
			string a12 = k12Clust + "_aHat"
			
			string sig13 = k13Clust+"_sigmaHat"
			string mu13 = k13Clust + "_muHat"
			string a13 = k13Clust + "_aHat"
			
			string sig14 = k14Clust+"_sigmaHat"
			string mu14 = k14Clust + "_muHat"
			string a14 = k14Clust + "_aHat"
			
			string map7 = k7Clust + "_clustMAP"
			string map8 = k8Clust + "_clustMAP"
			string map9 = k9Clust + "_clustMAP"
			string map10 = k10Clust + "_clustMAP"
			string map11 = k11Clust + "_clustMAP"
			string map12 = k12Clust + "_clustMAP"
			string map13 = k13Clust + "_clustMAP"
			string map14 = k14Clust + "_clustMAP"
			
			string sim1 = k1Clust + "_simData"
			string sim2 = k2Clust + "_simData"
			string sim3 = k3Clust + "_simData"
			string sim4 = k4Clust + "_simData"
			string sim5 = k5Clust + "_simData"
			string sim6 = k6Clust + "_simData"
			string sim7 = k7Clust + "_simData"
			string sim8 = k8Clust + "_simData"
			string sim9 = k9Clust + "_simData"
			string sim10 = k10Clust + "_simData"
			string sim11 = k11Clust + "_simData"
			string sim12 = k12Clust + "_simData"
			string sim13 = k13Clust + "_simData"
			string sim14 = k14Clust + "_simData"
			string simDataNa = dataName + "_sim"
			string finalSigma = dataName + "_S2"
			string finalMu = dataName + "_Mu"
			string finalA = dataName + "_AVals"

			if (onClust==0)
				duplicate/o $sim2, $simDataNa
				duplicate/o $sig2, $finalSigma
				duplicate/o $mu2, $finalMu
				duplicate/o $a2, $finalA
				duplicate/o $map2, MAP
			elseif (onClust==1)
				duplicate/o $sim3, $simDataNa
				duplicate/o $sig3, $finalSigma
				duplicate/o $mu3, $finalMu
				duplicate/o $a3, $finalA
				duplicate/o $map3, MAP
			elseif (onClust==2)
				duplicate/o $sim4, $simDataNa
				duplicate/o $sig4, $finalSigma
				duplicate/o $mu4, $finalMu
				duplicate/o $a4, $finalA
				duplicate/o $map4, MAP
			elseif (onClust==3)
				duplicate/o $sim5, $simDataNa
				duplicate/o $sig5, $finalSigma
				duplicate/o $mu5, $finalMu
				duplicate/o $a5, $finalA
				duplicate/o $map5, MAP
			elseif (onClust==4)
				duplicate/o $sim6, $simDataNa
				duplicate/o $sig6, $finalSigma
				duplicate/o $mu6, $finalMu
				duplicate/o $a6, $finalA
				duplicate/o $map6, MAP
			elseif (onClust==5)
				duplicate/o $sim7, $simDataNa
				duplicate/o $sig7, $finalSigma
				duplicate/o $mu7, $finalMu
				duplicate/o $a7, $finalA
				duplicate/o $map7, MAP
			elseif (onClust==6)
				duplicate/o $sim8, $simDataNa
				duplicate/o $sig8, $finalSigma
				duplicate/o $mu8, $finalMu
				duplicate/o $a8, $finalA
				duplicate/o $map8, MAP
			elseif (onClust==7)
				duplicate/o $sim9, $simDataNa
				duplicate/o $sig9, $finalSigma
				duplicate/o $mu9, $finalMu
				duplicate/o $a9, $finalA
				duplicate/o $map9, MAP
			elseif (onClust==8)
				duplicate/o $sim10, $simDataNa
				duplicate/o $sig10, $finalSigma
				duplicate/o $mu10, $finalMu
				duplicate/o $a10, $finalA
				duplicate/o $map10, MAP
			elseif (onClust==9)
				duplicate/o $sim11, $simDataNa
				duplicate/o $sig11, $finalSigma
				duplicate/o $mu11, $finalMu
				duplicate/o $a11, $finalA
				duplicate/o $map11, MAP
			elseif (onClust==10)
				duplicate/o $sim12, $simDataNa
				duplicate/o $sig12, $finalSigma
				duplicate/o $mu12, $finalMu
				duplicate/o $a12, $finalA
				duplicate/o $map12, MAP
			elseif (onClust==11)
				duplicate/o $sim13, $simDataNa
				duplicate/o $sig13, $finalSigma
				duplicate/o $mu13, $finalMu
				duplicate/o $a13, $finalA
				duplicate/o $map13, MAP
			elseif (onClust==12)
				duplicate/o $sim14, $simDataNa
				duplicate/o $sig14, $finalSigma
				duplicate/o $mu14, $finalMu
				duplicate/o $a14, $finalA
				duplicate/o $map14, MAP
			
			endif
			
				
				
			killwaves/z $k1Clust,$k2Clust,$k3Clust,$k4Clust,$k5Clust,$k6Clust,$k7Clust,$k8Clust,$k9Clust,$k10Clust,$k11Clust,$k12Clust,$k13Clust,$k14Clust
			killwaves/z $sig1,$sig2,$sig3,$sig4,$sig5,$sig6,$sig7,$sig8,$sig9,$sig10,$sig11,$sig12,$sig13,$sig14
			killwaves/z $mu1,$mu2,$mu3,$mu4,$mu5,$mu6,$mu7,$mu8,$mu9,$mu10,$mu11,$mu12,$mu13,$mu14
			killwaves/z $a1,$a2,$a3,$a4,$a5,$a6,$a7,$a8,$a9,$a10,$a11,$a12,$a13,$a14
		
			killwaves/z $map1,$map2,$map3,$map4,$map5,$map6,$map7,$map8,$map9,$map10,$map11,$map12,$map13,$map14
			killwaves/z $sim1,$sim2,$sim3,$sim4,$sim5,$sim6,$sim7,$sim8,$sim9,$sim10,$sim11,$sim12,$sim13,$sim14
			killwaves/z simDat,simMat,sumR,timeSeriesP,sigma2Hat,scaleData,sd,muHat,mapGrab,levelsWave
			killwaves/z cdata,ampHist,aHat,aHist,data,fit, clustDat
			duplicate/o $dataname, tempQ
			
			string quantaName = dataName + "Q"
			MAP+=1
			duplicate/o MAP, $quantaName
		
			string eventName = (dataName[0,strlen(dataName)-2]) + "E"
			display/k=1 $quantaName vs $eventName
			Label bottom "\Z14 Time (s)"
			Label Left "\Z14 Estimted Quanta"
			ModifyGraph mode=1
			SetAxis left 0,*
			killwaves/z tempQ
			variable nF = dimsize($eventName,0)
			variable i
			variable j
			make/o/n=(onClust+3) pWave, qWave
			string pWaveName = quantaName + "P"
			string qWaveName = quantaName + "Sum"
			pWave = 0
		
			
			for (i=0;i<nF;i+=1)
				for (j=0;j<onClust+3;j+=1)
					if (MAP[i]==j)
						pWave[j]+=1
					endif
				endfor
			endfor
			for (i=0;i<onClust+3;i+=1)
				qWave[i]=pWave[i] * i
			endfor
			variable qSum  = sum(qWave)
			qWave/=qSum
			pWave/=nF
			duplicate/o pWave, $pWaveName
			killwaves/z pWave
			display $pWaveName
			label left "Proportion of Events"
			label bottom "nQuanta"
			duplicate/o qWave, $qWaveName
			killwaves/z qWave
			display $qWaveName
			label left "Proportion of Quanta Distributed"
			label bottom "nQuanta"
		endswitch
	
	end



//////////////////////////EM
// Base of the EM algorithm
Function kEMgmm(scaleData,data,k,kDataName,scale)
wave scaleData,data
variable k
string kDataName
variable scale
variable tol = 0
string whichRun ="Running EM GMM on " + num2str(k) + " Clusters"
print whichRun

// Initialize some theta values - means are evenly spread across range of data and sigmas are 
// just globally estimated. Do some scaling (same as underflow issue above)
variable maxDat = wavemax(scaleData)
variable nObj = dimsize(scaleData,0)
make/o/n=(k) muHat,sigma2Hat,aHat
muHat = (p+1) * maxDat/(k+1)
sigma2Hat = sqrt(variance(scaleData))
aHat = 1/nObj
variable group,obj
make/o/d/n=(nobj,k) w
make/o/n=(nobj) MAP
variable objMAP
////////////////// EM
variable it
variable its = 25
string sigmaName = kDataname+"_sigmaHat"
			string muName = kDataName + "_muHat"
			string aName = kDataname + "_aHat"
			string clustMAPName = kDataName + "_clustMAP"
make/o/n=(its) timeSeriesP
for (it=0;it<its;it+=1)
// E Step
	// Given the parameters, estimate which group each data point belongs to
	
	for (obj=0;obj<nObj;obj+=1)
		for (group=0;group<k;group+=1)
			w[obj][group] = -.5*ln((2*pi*sigma2Hat[group])) -(scaleData[obj]-muHat[group])^2-(2*sigma2Hat[group])
			//if (w[obj][group]==-Inf)
			//	w[obj][group]=0
			//endif
		endfor
		duplicate/o/rmd=[obj][] w, tempRow
		variable logMax = wavemax(tempRow)
		w[obj][]-=logMax
	endfor
	w[][] = exp(w)
	matrixop/o sumR = sumRows(w)
	w[][] /= sumR[p]
/// M Step
	// Using clustering result above, reestimate paramters
	// AHat
	matrixop/o aHat = sumCols(w)
	aHat/=nObj
	//muHat
	muHat = 0
	for (group=0;group<k;group+=1)
		for (obj=0;obj<nObj;obj+=1)
			muHat[group] +=w[obj][group] * scaleData[obj]
		endfor
		muHat[group]/=aHat[group]*nObj
	endfor
	//sigma2Hat
	sigma2Hat = 0
	for (group=0;group<k;group+=1)
		for (obj=0;obj<nObj;obj+=1)
			sigma2Hat[group]+=w[obj][group] * ((scaleData[obj]-muHat[group])^2)
		endfor
		sigma2Hat[group] /=aHat[group]*nObj
		sigma2Hat[group] = sqrt(sigma2Hat[group])
		if (sigma2Hat[group]==0)
			sigma2Hat[group] = sqrt(variance(scaleData))
		endif
	endfor
	///// Likelihood/Convergence Shit
	variable logP = 0
	for (obj=0;obj<nObj;obj+=1)
		duplicate/o/RMD=[obj][] w, mapGrab
		objMAP = wavemax(mapGrab)
		findvalue/t=.001/V=(objMAP) mapGrab
		MAP[obj] =V_Value
		if (MAP[obj] == -1)
			MAP[obj] = 0
		endif
		logP+= -.5 * ln((2*pi*sigma2Hat[MAP[obj]]))   -(scaleData[obj]-muHat[MAP[obj]])^2-(2*sigma2Hat[MAP[obj]]) 
	endfor
	timeSeriesP[it] = logP
	print it
	if (it>20)
		if (timeSeriesP[it] >= sum(timeSeriesP,it-9,it)/10 + tol)
			
			duplicate/o MAP, $clustMAPName
			sigma2Hat *= (scale)^2
			muhat *= scale
			duplicate/o sigma2Hat, $sigmaName
			duplicate/o muhat, $muName
			duplicate/o aHat, $aName
			analEM(scaleData,data,kDataName,muname,sigmaName,aName,scale)
			deletepoints it, its, timeSeriesP
			string itRan = "Converged in " + num2str(it) + " iterations"
				print itRan
			return logP
		elseif (it==its-1)
			duplicate/o MAP, $clustMAPName
			sigma2Hat *= (scale)^2
			muhat *= scale
			duplicate/o sigma2Hat, $sigmaName
			duplicate/o muhat, $muName
			duplicate/o aHat, $aName
			killwaves sigma2Hat,aHat,muHat
			analEM(scaleData,data,kDataName,muname,sigmaName,aName,scale)
			deletepoints it, its, timeSeriesP
			itRan = "Did not converge in " + num2str(its) + " iterations"
			print itRan
			return logP
		endif
		killwaves sigma2Hat,aHat,muHat
		killwaves/z tempRow
	endif

			
	
	
	
endfor


end




///////////////////Analyze and Display Results
Function analEM(scaleData,data,kDataName,muname,sigmaName,aName,scale)

wave scaleData,data
variable scale
string kDataName, muname, sigmaName, aName
duplicate/o $muName, muHat
duplicate/o $sigmaName, sigma2Hat
duplicate/o $aName, aHat
variable k = dimsize(muHat,0)
variable maxDat = wavemax(data)
print maxDat
Make/N=100/O sd;DelayUpdate
Histogram/B={0,maxDat/100,100}/P data,sd;DelayUpdate
variable nS = 100
make/o/n=(ns) xS, simDat
xS = p
XS/=ns
xS*=maxDat
simDat = 0
string/g xSpace = "xSpace"
duplicate/o xS, $xSpace
variable i,j
for (i=0;i<nS;i+=1)
	for (j=0;j<k;j+=1)
		simDat[i] += aHat[j]*(2*pi*sigma2Hat[j])^-.5 * exp(-(xS[i]-muHat[j])^2/(2*sigma2Hat[j]))
	endfor
endfor

string simulatedData = kDataName + "_simData"
duplicate/o simDat, $simulatedData
end
		











/////////////////////// Make some fucking data to test this shit //////////////////////////
Function emMakeData()
variable nSamps = 200
variable k = 3
make/o/n=(nSamps) data
make/o/n=(k) mu,sigma2,A
mu[0] = {20,40,60}
sigma2 = 25
A[0] = {.5,.75,1}
A*= nSamps

variable i,j
for (i=0;i<nSamps;i+=1)
	if (i<=A[0])
		data[i] = gnoise(sqrt(sigma2[0]))+mu[0]
	elseif (i>A[k-2])
		data[i] = gnoise(sqrt(sigma2[k-1]))+mu[k-1]
	else
		for (j=1;j<k;j+=1)
			if (i<=A[j]&&i>A[j-1])
				data[i] = gnoise(sqrt(sigma2[j])) + mu[j]
			endif
		endfor
	endif
endfor

end
