#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.



Function cycPer2(eventName,trig1,trig2,trig3,trig4,trig5,trig6,tstop)


string eventName
variable trig1, trig2, trig3, trig4, trig5, trig6, tStop
duplicate/o $eventName, events

string baseWindowName = eventName + "_EPC"
string wind1 = baseWindowName + "1"
string wind2 = baseWindowName + "2"
string wind3 = baseWindowName + "3"
string wind4 = baseWindowName + "4"
string wind5 = baseWindowName + "5"
string wind6 = baseWindowName + "6"


variable freq1 = 1// stimulus frequency
variable oneCyc1 = 1/freq1 // length of cycle
variable freq2 = 5 // stimulus frequency
variable oneCyc2 = 1/freq2 // length of cycle
variable freq3 = 10// stimulus frequency
variable oneCyc3 = 1/freq3 // length of cycle
variable freq4 = 20 // stimulus frequency
variable oneCyc4 = 1/freq4 // length of cycle
variable freq0 = 0.99// stimulus frequency
variable oneCyc0 = 1/freq0 // length of cycle

variable i 
variable j 

// Number of windows in each condition
variable window1 = (trig2-trig1) / oneCyc0
variable window2 = (trig3-trig2) / oneCyc1
variable window3 = (trig4-trig3) / oneCyc2
variable window4 = (trig5-trig4) / oneCyc3
variable window5 = (trig6-trig5) / oneCyc4
variable window6 = (tStop-trig6) / oneCyc0


variable windowBeg, windowEnd
variable nEvents = numpnts(events)-1
make/o/n=(window1) firstBlock
make/o/n=(window2) secondBlock
make/o/n=(window3) thirdBlock
make/o/n=(window4) fourthBlock
make/o/n=(window5) fifthBlock
make/o/n=(window6) sixthBlock
firstBlock = 0
secondBlock = 0
thirdBlock = 0
fourthBlock = 0
fifthBlock = 0
sixthBlock = 0


/// Calculate events per cycle of first window
for (i=0;i<window1;i+=1)
	windowBeg = trig1+ (i * oneCyc0)
	windowEnd = windowBeg+oneCyc0
	for (j=0;j<nEvents;j+=1)
		if (events[j] > windowBeg && events[j] < windowEnd)
			firstBlock[i] +=1
		endif
	endfor
endfor



/// Calculate events per cycle of second window
for (i=0;i<window2;i+=1)
	windowBeg = trig2+ (i * oneCyc1)
	windowEnd = windowBeg+oneCyc1
	for (j=0;j<nEvents;j+=1)
		if (events[j] > windowBeg && events[j] < windowEnd)
			secondBlock[i] +=1
		endif
	endfor
endfor


/// Calculate events per cycle of third window
for (i=0;i<window3;i+=1)
	windowBeg = trig3+ (i * oneCyc2)
	windowEnd = windowBeg+oneCyc2
	for (j=0;j<nEvents;j+=1)
		if (events[j] > windowBeg && events[j] < windowEnd)
			thirdBlock[i] +=1
		endif
	endfor
endfor


/// Calculate events per cycle of fourth window
for (i=0;i<window4;i+=1)
	windowBeg = trig4+ (i * oneCyc3)
	windowEnd = windowBeg+oneCyc3
	for (j=0;j<nEvents;j+=1)
		if (events[j] > windowBeg && events[j] < windowEnd)
			fourthBlock[i] +=1
		endif
	endfor
endfor



/// Calculate events per cycle of fifth window
for (i=0;i<window5;i+=1)
	windowBeg = trig5+ (i * oneCyc4)
	windowEnd = windowBeg+oneCyc4
	for (j=0;j<nEvents;j+=1)
		if (events[j] > windowBeg && events[j] < windowEnd)
			fifthBlock[i] +=1
		endif
	endfor
endfor


/// Calculate events per cycle of sixth window
for (i=0;i<window6;i+=1)
	windowBeg = trig6+ (i * oneCyc0)
	windowEnd = windowBeg+oneCyc0
	for (j=0;j<nEvents;j+=1)
		if (events[j] > windowBeg && events[j] < windowEnd)
			sixthBlock[i] +=1
		endif
	endfor
endfor

duplicate/o firstBlock, $wind1
duplicate/o secondBlock, $wind2
duplicate/o thirdBlock, $wind3
duplicate/o fourthBlock, $wind4
duplicate/o fifthBlock, $wind5
duplicate/o sixthBlock, $wind6


variable mean1 = mean($wind1)
variable mean2 = mean($wind2)
variable mean3 = mean($wind3)
variable mean4 = mean($wind4)
variable mean5 = mean($wind5)
variable mean6 = mean($wind6)

variable var1 = variance($wind1)
variable var2 = variance($wind2)
variable var3 = variance($wind3)
variable var4 = variance($wind4)
variable var5 = variance($wind5)
variable var6 = variance($wind6)

variable f1 = var1/mean1
variable f2 = var2/mean2
variable f3 = var3/mean3
variable f4 = var4/mean4
variable f5 = var5/mean5
variable f6 = var6/mean6

string m1 = "Mean Events per Cycle =" + num2str(mean1)
string m2 = "Mean Events per Cycle =" + num2str(mean2)
string m3 = "Mean Events per Cycle =" + num2str(mean3)
string m4 = "Mean Events per Cycle =" + num2str(mean4)
string m5 = "Mean Events per Cycle =" + num2str(mean5)
string m6 = "Mean Events per Cycle =" + num2str(mean6)

string s1 = "Variance Events per Cycle =" + num2str(var1)
string s2 = "Variance Events per Cycle =" + num2str(var2)
string s3 = "Variance Events per Cycle =" + num2str(var3)
string s4 = "Variance Events per Cycle =" + num2str(var4)
string s5 = "Variance Events per Cycle =" + num2str(var5)
string s6 = "Variance Events per Cycle =" + num2str(var6)

string f1a = "Fano factor =" + num2str(f1)
string f2a = "Fano factor =" + num2str(f2)
string f3a = "Fano factor =" + num2str(f3)
string f4a = "Fano factor =" + num2str(f4)
string f5a = "Fano factor =" + num2str(f5)
string f6a = "Fano factor =" + num2str(f6)


display/k=1 $wind1
Textbox/C/N=text1 m1
AppendText/N=text1 s1
AppendText/N=text1 f1a

display/k=1 $wind2
Textbox/C/N=text1 m2
AppendText/N=text1 s2
AppendText/N=text1 f2a


display/k=1 $wind3
Textbox/C/N=text1 m3
AppendText/N=text1 s3
AppendText/N=text1 f3a

display/k=1 $wind4
Textbox/C/N=text1 m4
AppendText/N=text1 s4
AppendText/N=text1 f4a

display/k=1 $wind5
Textbox/C/N=text1 m5
AppendText/N=text1 s5
AppendText/N=text1 f5a


display/k=1 $wind6
Textbox/C/N=text1 m6
AppendText/N=text1 s6
AppendText/N=text1 f6a




end


Function cycPer(eventName,trig1,trig2,trig3,tstop)


string eventName
variable trig1, trig2, trig3, tStop
duplicate/o $eventName, events

string baseWindowName = eventName + "_EPC"
string wind1 = baseWindowName + "1"
string wind2 = baseWindowName + "2"
string wind3 = baseWindowName + "3"
variable freq = 5 // stimulus frequency
variable oneCyc = 1/5 // length of cycle

variable i 
variable j 

// Number of windows in each condition
variable window1 = (trig2-trig1) / oneCyc
variable window2 = (trig3-trig2) / oneCyc
variable window3 = (tStop-trig3) / oneCyc


variable windowBeg, windowEnd
variable nEvents = numpnts(events)-1
make/o/n=(window1) firstBlock
make/o/n=(window2) secondBlock
make/o/n=(window3) thirdBlock
firstBlock = 0
secondBlock = 0
thirdBlock = 0


/// Calculate events per cycle of first window
for (i=0;i<window1;i+=1)
	windowBeg = trig1+ (i * oneCyc)
	windowEnd = windowBeg+oneCyc
	for (j=0;j<nEvents;j+=1)
		if (events[j] > windowBeg && events[j] < windowEnd)
			firstBlock[i] +=1
		endif
	endfor
endfor



/// Calculate events per cycle of second window
for (i=0;i<window2;i+=1)
	windowBeg = trig2+ (i * oneCyc)
	windowEnd = windowBeg+oneCyc
	for (j=0;j<nEvents;j+=1)
		if (events[j] > windowBeg && events[j] < windowEnd)
			secondBlock[i] +=1
		endif
	endfor
endfor


/// Calculate events per cycle of third window
for (i=0;i<window3;i+=1)
	windowBeg = trig3+ (i * oneCyc)
	windowEnd = windowBeg+oneCyc
	for (j=0;j<nEvents;j+=1)
		if (events[j] > windowBeg && events[j] < windowEnd)
			thirdBlock[i] +=1
		endif
	endfor
endfor


duplicate/o firstBlock, $wind1
duplicate/o secondBlock, $wind2
duplicate/o thirdBlock, $wind3

variable mean1 = mean($wind1)
variable mean2 = mean($wind2)
variable mean3 = mean($wind3)
variable var1 = variance($wind1)
variable var2 = variance($wind2)
variable var3 = variance($wind3)
variable f1 = var1/mean1
variable f2 = var2/mean2
variable f3 = var3/mean3

string m1 = "Mean Events per Cycle =" + num2str(mean1)
string m2 = "Mean Events per Cycle =" + num2str(mean2)
string m3 = "Mean Events per Cycle =" + num2str(mean3)

string s1 = "Variance Events per Cycle =" + num2str(var1)
string s2 = "Variance Events per Cycle =" + num2str(var2)
string s3 = "Variance Events per Cycle =" + num2str(var3)

string f1a = "Fano factor =" + num2str(f1)
string f2a = "Fano factor =" + num2str(f2)
string f3a = "Fano factor =" + num2str(f3)


display/k=1 $wind1
Textbox/C/N=text1 m1
AppendText/N=text1 s1
AppendText/N=text1 f1a

display/k=1 $wind2
Textbox/C/N=text1 m2
AppendText/N=text1 s2
AppendText/N=text1 f2a


display/k=1 $wind3
Textbox/C/N=text1 m3
AppendText/N=text1 s3
AppendText/N=text1 f3a

end

Function vesPerCyc(eventName,trig1,trig2,trig3,trig4,trig5,trig6,tstop,unitaryEvent)
string eventName
variable trig1, trig2, trig3, trig4, trig5, trig6, tStop, unitaryEvent





duplicate/o $eventName, events

string baseWindowName = eventName + "_EPC"
string wind1 = baseWindowName + "1"
string wind2 = baseWindowName + "2"
string wind3 = baseWindowName + "3"
string wind4 = baseWindowName + "4"
string wind5 = baseWindowName + "5"
string wind6 = baseWindowName + "6"


variable freq1 = 1// stimulus frequency
variable oneCyc1 = 1/freq1 // length of cycle
variable freq2 = 5 // stimulus frequency
variable oneCyc2 = 1/freq2 // length of cycle
variable freq3 = 10// stimulus frequency
variable oneCyc3 = 1/freq3 // length of cycle
variable freq4 = 20 // stimulus frequency
variable oneCyc4 = 1/freq4 // length of cycle
variable freq0 = 0.99// stimulus frequency
variable oneCyc0 = 1/freq0 // length of cycle


// Number of windows in each condition
variable window1 = (trig2-trig1) / oneCyc0
variable window2 = (trig3-trig2) / oneCyc1
variable window3 = (trig4-trig3) / oneCyc2
variable window4 = (trig5-trig4) / oneCyc3
variable window5 = (trig6-trig5) / oneCyc4
variable window6 = (tStop-trig6) / oneCyc0


variable windowBeg, windowEnd
variable nEvents = numpnts(events)-1
make/o/n=(window1+1) firstBlock
make/o/n=(window2) secondBlock
make/o/n=(window3) thirdBlock
make/o/n=(window4) fourthBlock
make/o/n=(window5) fifthBlock
make/o/n=(window6) sixthBlock
firstBlock = 0
secondBlock = 0
thirdBlock = 0
fourthBlock = 0
fifthBlock = 0
sixthBlock = 0


variable i,j

for (i=0;i<window1;i+=1)
	firstBlock[i] = sum(events,(trig1+(i*oneCyc0)),(trig1+(i*oneCyc0)+oneCyc0))
endfor
firstBlock /= unitaryEvent *oneCyc0

for (i=0;i<window2;i+=1)
	secondBlock[i] = sum(events,(trig2+(i*oneCyc1)),(trig2+(i*oneCyc1)+oneCyc1))
endfor
secondBlock /= unitaryEvent*oneCyc1
end

for (i=0;i<window3;i+=1)
thirdBlock[i] = sum(events,(trig3+(i*oneCyc2)),(trig3+(i*oneCyc2)+oneCyc2))
endfor
thirdBlock /= unitaryEvent*oneCyc2

for (i=0;i<window4;i+=1)
fourthBlock[i] = sum(events,(trig4+(i*oneCyc3)),(trig4+(i*oneCyc3)+oneCyc3))
endfor
fourthBlock /= unitaryEvent*oneCyc3

for (i=0;i<window5;i+=1)
fifthBlock[i] = sum(events,(trig5+(i*oneCyc4)),(trig5+(i*oneCyc4)+oneCyc4))
endfor
fifthBlock /= unitaryEvent*oneCyc4

for (i=0;i<window6;i+=1)
sixthBlock[i] = sum(events,(trig6+(i*oneCyc5)),(trig6+(i*oneCyc5)+oneCyc5))
endfor
fifthBlock /= unitaryEvent/oneCyc5


duplicate/o firstBlock, $wind1
duplicate/o secondBlock, $wind2
duplicate/o thirdBlock, $wind3
duplicate/o fourthBlock, $wind4
duplicate/o fifthBlock, $wind5
duplicate/o sixthBlock, $wind6


variable mean1 = mean($wind1)
variable mean2 = mean($wind2)
variable mean3 = mean($wind3)
variable mean4 = mean($wind4)
variable mean5 = mean($wind5)
variable mean6 = mean($wind6)

variable var1 = variance($wind1)
variable var2 = variance($wind2)
variable var3 = variance($wind3)
variable var4 = variance($wind4)
variable var5 = variance($wind5)
variable var6 = variance($wind6)

variable f1 = var1/mean1
variable f2 = var2/mean2
variable f3 = var3/mean3
variable f4 = var4/mean4
variable f5 = var5/mean5
variable f6 = var6/mean6

string m1 = "Mean Ves per Cycle =" + num2str(mean1)
string m2 = "Mean Ves per Cycle =" + num2str(mean2)
string m3 = "Mean Ves per Cycle =" + num2str(mean3)
string m4 = "Mean Ves per Cycle =" + num2str(mean4)
string m5 = "Mean Ves per Cycle =" + num2str(mean5)
string m6 = "Mean Ves per Cycle =" + num2str(mean6)

string s1 = "Variance Ves per Cycle =" + num2str(var1)
string s2 = "Variance Ves per Cycle =" + num2str(var2)
string s3 = "Variance Ves per Cycle =" + num2str(var3)
string s4 = "Variance Ves per Cycle =" + num2str(var4)
string s5 = "Variance Ves per Cycle =" + num2str(var5)
string s6 = "Variance Ves per Cycle =" + num2str(var6)

string f1a = "Fano factor =" + num2str(f1)
string f2a = "Fano factor =" + num2str(f2)
string f3a = "Fano factor =" + num2str(f3)
string f4a = "Fano factor =" + num2str(f4)
string f5a = "Fano factor =" + num2str(f5)
string f6a = "Fano factor =" + num2str(f6)


display/k=1 $wind1
Textbox/C/N=text1 m1
AppendText/N=text1 s1
AppendText/N=text1 f1a

display/k=1 $wind2
Textbox/C/N=text1 m2
AppendText/N=text1 s2
AppendText/N=text1 f2a


display/k=1 $wind3
Textbox/C/N=text1 m3
AppendText/N=text1 s3
AppendText/N=text1 f3a

display/k=1 $wind4
Textbox/C/N=text1 m4
AppendText/N=text1 s4
AppendText/N=text1 f4a

display/k=1 $wind5
Textbox/C/N=text1 m5
AppendText/N=text1 s5
AppendText/N=text1 f5a


display/k=1 $wind6
Textbox/C/N=text1 m6
AppendText/N=text1 s6
AppendText/N=text1 f6a




end



Function vesPerCyc2(eventName,trig1,trig2,trig3,tstop,unitaryEvent)
string eventName
variable trig1, trig2, trig3, tstop, unitaryEvent

duplicate/o $eventName, events

string baseWindowName = eventName + "_EPC"
string wind1 = baseWindowName + "1"
string wind2 = baseWindowName + "2"
string wind3 = baseWindowName + "3"
variable freq = 5 // stimulus frequency
variable oneCyc = 1/5 // length of cycle

// Number of windows in each condition

variable window1 = (trig2-trig1) / oneCyc
variable window2 = (trig3-trig2) / oneCyc
variable window3 = (tStop-trig3) / oneCyc


variable windowBeg, windowEnd
variable nEvents = numpnts(events)-1
make/o/n=(window1+1) firstBlock
make/o/n=(window2) secondBlock
make/o/n=(window3) thirdBlock

firstBlock = 0
secondBlock = 0
thirdBlock = 0

variable i,j

for (i=0;i<window1;i+=1)
	firstBlock[i] = sum(events,(trig1+(i*oneCyc)),(trig1+(i*oneCyc)+oneCyc))
endfor
firstBlock /= unitaryEvent*oneCyc

for (i=0;i<window2;i+=1)
	secondBlock[i] = sum(events,(trig2+(i*oneCyc)),(trig2+(i*oneCyc)+oneCyc))
endfor
secondBlock /= unitaryEvent*oneCyc
end

for (i=0;i<window3;i+=1)
thirdBlock[i] = sum(events,(trig3+(i*oneCyc)),(trig3+(i*oneCyc)+oneCyc))
endfor
thirdBlock /= unitaryEvent*oneCyc

duplicate/o firstBlock, $wind1
duplicate/o secondBlock, $wind2
duplicate/o thirdBlock, $wind3

variable mean1 = mean($wind1)
variable mean2 = mean($wind2)
variable mean3 = mean($wind3)
variable var1 = variance($wind1)
variable var2 = variance($wind2)
variable var3 = variance($wind3)
variable f1 = var1/mean1
variable f2 = var2/mean2
variable f3 = var3/mean3

string m1 = "Mean Ves per Cycle =" + num2str(mean1)
string m2 = "Mean Ves per Cycle =" + num2str(mean2)
string m3 = "Mean Ves per Cycle =" + num2str(mean3)

string s1 = "Variance Ves per Cycle =" + num2str(var1)
string s2 = "Variance Ves per Cycle =" + num2str(var2)
string s3 = "Variance Ves per Cycle =" + num2str(var3)

string f1a = "Fano factor =" + num2str(f1)
string f2a = "Fano factor =" + num2str(f2)
string f3a = "Fano factor =" + num2str(f3)


display/k=1 $winds1
Textbox/C/N=text1 m1
AppendText/N=text1 s1
AppendText/N=text1 f1a

display/k=1 $winds2
Textbox/C/N=text1 m2
AppendText/N=text1 s2
AppendText/N=text1 f2a


display/k=1 $winds3
Textbox/C/N=text1 m3
AppendText/N=text1 s3
AppendText/N=text1 f3a
