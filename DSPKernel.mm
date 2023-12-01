/*
	Copyright (C) 2016 Apple Inc. All Rights Reserved.
	See LICENSE.txt for this sample’s licensing information
	
	Abstract:
	Utility code to manage scheduled parameters in an audio unit implementation.
*/

#import "../harmonizr-dsp/DSPKernel.hpp"

void DSPKernel::handleOneEvent(AURenderEvent const *event) {
    
    AUMIDIEvent e;
    
	switch (event->head.eventType) {
		case AURenderEventParameter:
		case AURenderEventParameterRamp: {
			AUParameterEvent const& paramEvent = event->parameter;
			
			startRamp(paramEvent.parameterAddress, paramEvent.value, paramEvent.rampDurationSampleFrames);
			break;
		}
			
		case AURenderEventMIDI:
            
            e = event->MIDI;
            if ((e.data[0] & 0xF0) == 0xC0)
            {
                pc_flag = 1;
                program_change = e.data[1];
                dispatch_semaphore_signal(sem);
            }
            else
            {
                handleMIDIEvent(event->MIDI);
                
                if ((e.data[0] & 0xF0) == 0xB0)
                {
                    cc_num = e.data[1];
                    cc_val = e.data[2];
                    cc_flag = 1;
                    dispatch_semaphore_signal(sem);
                }
            }
            
			break;
		
		default:
			break;
	}
}

void DSPKernel::sendMIDIOutput(AUEventSampleTime now, AUMIDIOutputEventBlock midiOut)
{
    for (int k = 0; k < n_output_events; k++)
    {
        //printf("%x %x %x\n", output_events[k].data[0],output_events[k].data[1],output_events[k].data[2]);
        midiOut(now, 0, output_events[k].length, output_events[k].data);
    }
    n_output_events = 0;
}

void DSPKernel::performAllSimultaneousEvents(AUEventSampleTime now, AURenderEvent const *&event, AUMIDIOutputEventBlock midiOut) {
	do {
		handleOneEvent(event);
        
        if (event->head.eventType == AURenderEventMIDI && midiOut)
        {
            midiOut(now, 0, event->MIDI.length, event->MIDI.data);
        }

		// Go to next event.
		event = event->head.next;
        
		// While event is not null and is simultaneous.
	} while (event && event->head.eventSampleTime == now);
}

/**
	This function handles the event list processing and rendering loop for you.
	Call it inside your internalRenderBlock.
*/
void DSPKernel::processWithEvents(AudioTimeStamp const *timestamp, AUAudioFrameCount frameCount, AURenderEvent const *events, AUMIDIOutputEventBlock midiOut) {

	AUEventSampleTime now = AUEventSampleTime(timestamp->mSampleTime);
    static AUEventSampleTime then = 0;
    //printf("%d frames\n", frameCount);
//    if ((now - then) != frameCount && then != 0)
//        printf("********** %lld\n", now - then);
    then = now;
    
	AUAudioFrameCount framesRemaining = frameCount;
	AURenderEvent const *event = events;
	
	while (framesRemaining > 0) {
		// If there are no more events, we can process the entire remaining segment and exit.
		if (event == nullptr) {
			AUAudioFrameCount const bufferOffset = frameCount - framesRemaining;
			process(framesRemaining, bufferOffset);
            sendMIDIOutput(now, midiOut);
			return;
		}

		AUAudioFrameCount const framesThisSegment = AUAudioFrameCount(event->head.eventSampleTime - now);
		
		// Compute everything before the next event.
		if (framesThisSegment > 0) {
			AUAudioFrameCount const bufferOffset = frameCount - framesRemaining;
			process(framesThisSegment, bufferOffset);
            sendMIDIOutput(now, midiOut);
							
			// Advance frames.
			framesRemaining -= framesThisSegment;

			// Advance time.
			now += AUEventSampleTime(framesThisSegment);
		}
		
		performAllSimultaneousEvents(now, event, midiOut);
	}
}

