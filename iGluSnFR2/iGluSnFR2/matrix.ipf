#pragma rtGlobals=3		// Use modern global access method and strict wave access.

function ones_spikes(eventwave,w,sr)
wave eventwave
wave w
variable sr
variable i, nl=dimsize(w,0),nr=dimsize(eventwave,0)

	make /o/n=(dimsize(w,0),0) oneslocationspikes=0
	//oneslocationspikes[p][q]=0
		for (i=0;i<nr;i+=1)
			oneslocationspikes [eventwave[i]*sr][]=1
		endfor
end