#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


Function infoTheory2(baseName,sumWind,trig1,trig2,trig3,tstop)
string basename
variable sumWind, trig1,trig2,trig3,tstop

string eventName = baseName + "_E"
string ampName = baseName + "_A"

string totalSums = baseName + "E_EAll"
string firstSums = baseName + "E_E1"
string secondSums = baseName + "E_E2"
string thirdSums = baseName + "E_E3"
string firstHist = firstSums + "H"
string secondHist = secondSums + "H"
string totalHist = totalSums + "H"
string thirdHist = thirdSums + "H"

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

//
duplicate/o/r=(trig1,trig2) sums, firstEs
duplicate/o/r=(trig2,trig3) sums, secondEs
duplicate/o/r=(trig3,tstop) sums, thirdEs

//r(t) rolling sum
duplicate/o firstEs, firstRs
duplicate/o secondEs, secondRs
duplicate/o thirdEs, thirdRs
duplicate/o sums, rates



variable boxL = sumWind/deltax(rates)
boxsmooth boxL, firstEs, firstRs
boxsmooth boxL, secondEs, secondRs
boxsmooth boxL, thirdEs, thirdRs
boxsmooth boxL, sums, rates
firstRs *= boxL
secondRs *= boxL
thirdRs *=boxL
rates *=boxL

//analyse stuff
variable biggestW = wavemin(rates)+1
rates /=biggestW
firstRs/= biggestW
secondRs/=biggestW
thirdRs/=biggestW

variable nPoints1 = numpnts(firstRs)
variable nPoints2 = numpnts(secondRs)
variable nPoints3 = numpnts(thirdRs)
variable nPoints = numpnts(rates)


for (i=0;i<nPoints1;i+=1)
	firstRs[i] = round(firstRs[i])
endfor

for (i=0;i<nPoints2;i+=1)
	secondRs[i] = round(secondRs[i])
endfor

for (i=0;i<nPoints3;i+=1)
	thirdRs[i] = round(thirdRs[i])
endfor

for (i=0;i<nPoints;i+=1)
	rates[i] = round(rates[i])
endfor

variable biggestWord = waveMax(rates)



Make/N=(biggestWord)/O rates_Hist;DelayUpdate
Histogram/P/B={0,1,biggestWord} rates,rates_Hist;DelayUpdate
display/k=1 rates_Hist
ModifyGraph mode=5

Make/N=(biggestWord)/O firstRs_Hist;DelayUpdate
Histogram/P/B={0,1,biggestWord} firstRs,firstRs_Hist;DelayUpdate
display/k=1 firstRs_Hist
ModifyGraph mode=5

Make/N=(biggestW)/O secondRS_Hist;DelayUpdate
Histogram/P/B={0,1,biggestWord} secondRs,secondRs_Hist;DelayUpdate
display/k=1 secondRs_Hist
ModifyGraph mode=5


Make/N=(biggestW)/O thirdRs_Hist;DelayUpdate
Histogram/P/B={0,1,biggestWord} thirdRs,thirdRs_Hist;DelayUpdate
display/k=1 thirdRs_Hist
ModifyGraph mode=5

duplicate/o firstRs, $firstSums
duplicate/o secondRs, $secondSums
duplicate/o thirdRs, $thirdSums
duplicate/o rates, $totalSums

duplicate/o rates_hist, $totalHist
duplicate/o firstRs_Hist, $firstHist
duplicate/o secondRs_Hist, $secondHist
duplicate/o thirdRs_Hist, $thirdHist


variable mean1 = mean(firstRs)
variable mean2 = mean(secondRs)
variable mean3 = mean(thirdRs)
variable grandMean = mean(rates)

variable var1 = variance(firstRs)
variable var2 = variance(secondRs)
variable var3 = variance(thirdRs)
variable grandVar = variance(rates)


string s1a = "Mean =" + num2str(mean1)
string s2a = "Mean =" + num2str(mean2)
string s3a = "Mean =" + num2str(mean3)
string sta = "Mean =" + num2str(grandMean)
string s1b = "Variance =" + num2str(var1)
string s2b = "Variance =" + num2str(var2)
string s3b = "Variance =" + num2str(var3)
string stb = "Variance =" + num2str(grandVar)


display/k=1 $totalSums
TextBox/C/N=text0 sta
AppendText/N=text0 stb

display/k=1 $firstSums
TextBox/C/N=text0 s1a
AppendText/N=text0 s1b


display/k=1 $secondSums
TextBox/C/N=text0 s2a
AppendText/N=text0 s2b

display/k=1 $thirdSums
TextBox/C/N=text0 s3a
AppendText/N=text0 s3b


killwaves firstRs, secondRs, thirdRs, rates







end
