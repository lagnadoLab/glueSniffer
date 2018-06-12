#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.



function lifModel(eventName,ampName)

string eventName, ampName


duplicate/o $eventName, events
duplicate/o $ampName, amps
variable vThresh = -50
variable vRest = -60
variable nVari = .5
variable tauLeak = 20


variable gMin = 1
variable gMax = 10
variable gStep = .5
variable nG = 1+(gMax-gMin)/gStep

make/o/n=(nG) gVals
gVals = (p+2)/(1+nG)*10


//////////// Create QTSs
variable dt = .001
variable nEv = dimsize(events,0)
variable maxEv = wavemax(events)+1
variable nP = maxEv/dt
make/o/n=(nP) QTS, qts2
make/o/n=(nP,nG) vv, v2, spikes, spikes2
setscale/p x,0,dt,qts,qts2,vv,v2,spikes,spikes2
variable aveQ = sum(amps)/nEv
variable i,j
for (i=0;i<nEv;i+=1)
	qts[x2pnt(qts,events[i])] = amps[i]
	qts2[x2pnt(qts2,events[i])] = aveQ
endfor


vv[0][] = vRest
v2[0][] = vRest
spikes = 0
spikes2 = 0

variable maxSpikes = 1000


variable onWhichSpike = 0
variable onWhichSpike2 = 0
////////////////main sim body
killwaves/z gVS
make/o/n=(2,nG) nSpikes
for (i=0;i<nG;i+=1)
onWhichSpike2 = 0
onWhichSpike = 0
make/o/n=(maxSpikes) spikeTime, spikeTime2
spikeTime = 0
spikeTime2 = 0
	for (j=1;j<nP;j+=1)
	
		vv[j][i] = vv[j-1][i] + gnoise(sqrt(nVari))+ (qts[j] * gVals[i])
		vv[j][i] = vv[j][i] - ((vv[j][i]-vRest) /tauLeak)
		if (vv[j][i] > vThresh)
			spikes[j][i] = 1
			vv[j][i] = vRest
			spikeTime[onWhichSpike] = j * dt
			onWhichSpike += 1
		endif



		v2[j][i] = v2[j-1][i] + gnoise(sqrt(nVari))+ (qts2[j] * gVals[i])
		v2[j][i] = v2[j][i] - ((v2[j][i]-vRest) /tauLeak)
		if (v2[j][i] > vThresh)
			spikes2[j][i] = 1
			v2[j][i] = vRest
			spikeTime2[onWhichSpike2] = j * dt
			onWhichSpike2 += 1
		endif
		
	endfor
	deletepoints onWhichSpike, maxSpikes-onWhichSpike, spikeTime
	deletepoints onWhichSpike2, maxSpikes-onWhichSpike2, spikeTime2
	vecStren(spikeTime,spikeTime2)
	wave vs12
		nSpikes[0][i] = onWhichSpike
		nSpikes[1][i] = onWhichSpike2
	concatenate/np=1 {vs12}, gVs
endfor


for (i=0;i<nG;i+=1)
	for (j=1;j<nP;j+=1)
		if (spikes[j][i]==1)
			vv[j][i] = vThresh
		endif
		if (spikes2[j][i]==1)
			v2[j][i] = vThresh
		endif
	endfor
endfor








variable/g onROI = 0
variable/g nROI = nG


string/g vOuts = "vOut"
string/g v2Outs = "v2Out"
string/g qtsOuts = "qtsOut"
string/g qts2Outs = "qts2Out"
string/g spikesOuts = "spikesOut"
string/g spikes2Outs = "spikes2Out"
string/g gValOuts = "gValues"

string outV = eventName + "_rgcV"
string outV2 = eventName + "_rgcV2"
string outQTS = eventName + "_rgcQTS"
string outQTS2 = eventName + "_rgcQTS2"
string outSpikes = eventName + "_rgcSpikes"
string outSpikes2 = eventName + "_rgcSpikes2"
string gVSout = eventName + "_rgcvs"
string nSpikeOut = eventName + "_rgcNS"

duplicate/o nSpikes, $nSpikeOut
duplicate/o gvs, $gvsOut
duplicate/o vv,$vOuts, $outV
duplicate/o v2,$v2Outs, $outV2
duplicate/o qts,$qtsOuts, $outQTS
duplicate/o qts2,$qts2Outs, $outQTS2
duplicate/o spikes,$spikesOuts, $outSpikes
duplicate/o spikes2,$spikes2Outs, $outSpikes2
duplicate/o gVals, $gValOuts


killwaves/z vv,v2,qts,qts2,spikes,spikes2,gVals
plotShit()
end






function vecStren(spikeTime,spikeTime2)


wave spikeTime, spikeTime2

variable stimStart = 10
variable stimStop =40
variable freq = 5


duplicate/o spikeTime, ampsS
duplicate/o spikeTime2, ampsS2

ampsS = 1
ampsS2 = 1
variable vs1,vs2
pullEvents(stimStart,stimStop,spikeTime,ampsS)
wave stimEvents, stimAmps
if (dimsize(stimEvents,0)<2)
	 vs1 = 0
else
duplicate/o stimEvents, phaseTimes, sinPart, cosPart
phaseTimes = 0

variable i

for (i=0;i<dimsize(phaseTimes,0);i+=1)
	phaseTimes[i] = mod(stimEvents[i]-stimStart,1/freq)
	sinPart[i] = sin(phaseTimes[i] * 2 *pi * freq)
	cosPart[i] = cos(phaseTimes[i] * 2 * pi * freq)
endfor

 vs1 = 1/dimsize(phaseTimes,0) * sqrt(sum(cosPart)^2 + sum(sinPart)^2)
endif



pullEvents(stimStart,stimStop,spikeTime2,ampsS2)
wave stimEvents, stimAmps
if (dimsize(stimEvents,0)<2)
	vs2 = 0
else
duplicate/o stimEvents, phaseTimes, sinPart, cosPart
phaseTimes = 0
for (i=0;i<dimsize(phaseTimes,0);i+=1)
	phaseTimes[i] = mod(stimEvents[i]-stimStart,1/freq)
	sinPart[i] = sin(phaseTimes[i] * 2 *pi * freq)
	cosPart[i] = cos(phaseTimes[i] * 2 * pi * freq)
endfor
 vs2 = 1/dimsize(phaseTimes,0) * sqrt(sum(cosPart)^2 + sum(sinPart)^2)
endif
make/o/n=2 vs12
vs12 = {vs1,vs2}

end


////////////////////////////////////////
function plotShit()
string/g gValOuts
string/g vOuts
string/g v2Outs
string/g qtsOuts
string/g qts2Outs
string/g spikesOuts
string/g spikes2Outs


duplicate/o $vOuts,vv
duplicate/o $v2Outs, v2
duplicate/o $qtsOuts, qts
duplicate/o $qts2Outs, qts2
duplicate/o $spikesOuts,spikes
duplicate/o $spikes2Outs,spikes2


variable/g nROI
variable/g onROI

variable nTrace = dimsize(vv,0)
// Panel
duplicate/o/RMD=[][onROI] vv, tempV
duplicate/o/RMD=[][onROI] v2, tempV2
duplicate/o/RMD=[][onROI] spikes, tempSpikes
duplicate/o/RMD=[][onROI] spikes2, tempSpikes2
duplicate/o qts,tempqts
duplicate/o qts2,tempqts2
duplicate/o $gValOuts, gVals



Display/K=1/W=(100,0,800,400)/N=roiThingy tempV as "LIF Trace"
ModifyGraph margin(bottom)=100
ModifyGraph axisEnab(left)={0.4,0.6}
label left "Voltage (mV)"


AppendToGraph/L=qtAx tempqts
ModifyGraph axisEnab(qtAx)={0.7,0.9},freePos(qtAx)=0
ModifyGraph axisEnab(bottom)={0,0.45}
label qtAx "Estimated Event Quanta"

AppendToGraph/L=spikeAx tempSpikes
ModifyGraph axisEnab(spikeAx)={0,0.3},freePos(qtAx)=0
ModifyGraph freePos(spikeAx)=0


appendtograph/B=newBottom tempV2
AppendToGraph/B=newBottom/L=qtAx tempqts2
appendtograph/B=newBottom/L=spikeAx tempSpikes2
label spikeAx "RGC Spikes"



ModifyGraph axisEnab(bottom)={0,0.45}
ModifyGraph axisEnab(newBottom)={0.55,1}
ModifyGraph freePos(newBottom)=0
ModifyGraph freePos(spikeAx)=0



//Val Display
ValDisplay whichROI value = #"onROI"

// Good Button
Button prevButt, proc = prevButt, title = "Previous"
Button nextButt, proc =nextButt, title = "Next"
Button closeButt, proc = closeBu, title = "Close"


killwaves/z tempV, tempV2, tempSpikes, tempQts, tempQTS2, tempSpikes2
end

Function closeBu(ba) : buttonControl
struct WMButtonAction&ba
	switch(ba.eventcode)
		case 2:
			dowindow/k roiThingy
			killstrings/z vOuts, v2Outs, qtsOuts, qts2Outs, spikesOuts,spikes2Outs,gValOuts
			killwaves/z v,v2,qts,qts2,spikes,spikes2,gVals
	endswitch
end

Function nextButt(ba) : buttonControl
struct WMButtonAction&ba
	switch(ba.eventcode)
		case 2:
			dowindow/k roiThingy
			variable/g onROI
			variable/g nROI
			
			if (onROI < nROI-1)
				onROI+=1
				plotShit()
			else
				onROI = 0
				plotShit()
			endif
	endswitch
end


Function prevButt(ba) : buttonControl
struct WMButtonAction&ba
	switch(ba.eventcode)
		case 2:
			dowindow/k roiThingy
			variable/g onROI
			variable/g nROI
			
			if (onROI >0 )
				onROI-=1
				plotShit()
			else 
				onROI = nROI-1
				plotShit()
			endif
	endswitch
end