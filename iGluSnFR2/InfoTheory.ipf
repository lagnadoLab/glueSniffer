#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
Function infoTheory(baseName,sumWind)
string baseName
variable sumWind

string eventName = baseName + "_E"
string ampName = baseName + "_A"

//events(t) Create a wave of zeros and event amplitudes as function of time
duplicate/o $baseName, sums
sums = 0
duplicate/o $ampName, amps
duplicate/o $eventName, events
variable nEvents = numpnts(events)
variable i
for (i=0;i<nEvents;i+=1)
	sums[x2pnt(sums,events[i])] =1
endfor

// r(t) : Create a wave of rolling sums per unit sumWindow
duplicate/o sums, rates
variable boxL = sumWind / deltax(rates)
boxsmooth boxL, sums, rates
rates = rates * boxL

//Analyze frequencies
variable biggestW = wavemin(rates)+1
duplicate/o rates, sumE
sumE = sumE/ biggestW

variable nPoints = numpnts(sumE)

for (i=0;i<nPoints;i+=1)
	sumE[i] = round(sumE[i])
endfor

variable biggestWord = waveMax(sumE)+1


Make/N=(biggestWord)/O sumE_Hist;DelayUpdate
Histogram/P/B={0,1,biggestWord} sumE,sumE_Hist;DelayUpdate
display/k=1 sumE_Hist
ModifyGraph mode=5




variable meanR = mean(rates)/sumWind
variable varR = variance(rates)/sumWind
string annString2 = "Mean Rate ="+num2str(meanR)
string annString3 = "Rate Variance =" +num2str(varR)
string annString = "WindowSize ="+ num2str(sumWind)
display/k=1 rates
TextBox/C/N=text0 annString
AppendText/N=text0 annString2
AppendText/N=text0 annString3

end
