//
//  HarmonizerDSPKermel.cpp
//  iOSHarmonizerFramework
//
//  Created by Matthew E Robbins on 12/1/23.
//

// NB This file must be set to type "Objective-C++ source in XCode"
#include <stdio.h>
#include <iostream>
#include "HarmonizerDSPKernel.hpp"
#include "Util.hpp"

void HarmonizerDSPKernel::init(int inChannels, int outChannels, double inSampleRate) {
    n_channels = outChannels;
    sampleRate = float(inSampleRate);
    fprintf(stderr,"**** init with %d channels! at %f Hz\n", n_channels, inSampleRate);
    
    in_buffers = (float **) calloc(inChannels, sizeof(float *));
    out_buffers = (float **) calloc(outChannels, sizeof(float *));
    
    filter = ButterworthFilter();
    looper.emplace(n_channels,60*sampleRate,(int)lrintf(0.05*sampleRate),(int)lrintf(0.05*sampleRate));
    
    noise_gate = NoiseGate(gate_thresh, 6.0f, sampleRate, 0.25);
    
    harmonizer.emplace(l2nfft, nvoices, maxT, sampleRate);
    
    freezers.reserve(n_channels);
    for (int k = 0; k < n_channels; k++) {
        freezers.emplace_back(13,5);
    }
    
    for (int k = 0; k < nvoices; k++)
    {
        if (k >= N_AUTO)
        {
            psolaVoices[k].setPan(((float)(k - 3) / (float)(nvoices - 3)) - 0.5);
        }
        psolaVoices[k].setFs(sampleRate);
    }
    
    // define equal tempered interval ratios
    
    intervals = interval_table + 24;
    
    for (int k = -23; k < 24; k++)
    {
        intervals[k] = powf(2.0, (float) (k) / 12);
    }
    
    memset(midinotes, 0, 128 * sizeof(float));
    memset(keys_down, 0, 128 * sizeof(int));
    
}

void HarmonizerDSPKernel::fini() {
    fprintf(stderr, "*** fini! ***\n");
    
    free(in_buffers);
    free(out_buffers);
    
    freezers.clear();
}

void HarmonizerDSPKernel::reset() {
    for (int k = 0; k < nvoices; k++)
    {
        if (k >= n_auto)
        {
            psolaVoices[k].setPan(((float)(k - 3) / (float)(nvoices - 3)) - 0.5);
            //voices[k].formant_ratio = ((float)(k - 3) / (float)(nvoices - 3)) + 0.5;
        }
    }
    
    rcnt = 256;
    T = 400;
}

void HarmonizerDSPKernel::setBuffers(float ** in, float ** out) {

    for (int k = 0; k < n_channels; k++)
    {
        in_buffers[k] = in[k];
        out_buffers[k] = out[k];
    }
}

void HarmonizerDSPKernel::setParameter(param_address_t address, param_value_t value) {
    switch (address) {
        case HarmParamKeycenter:
            root_key = (int) clamp(value,0.f,47.f);
            break;
        case HarmParamInversion:
            inversion = (int) clamp(value,0.f,3.f);
            harmonizer->setInversion(inversion);
            break;
        case HarmParamNvoices:
            n_auto = (int) clamp(value,1.f,4.f);
            harmonizer->setAutoVoiceCount(n_auto);
            for (int k = n_auto; k < N_AUTO; k++)
            {
                ui_voice_notes[k] = -1;
            }
            fprintf(stderr, "nvoices: %d\n", n_auto);
            break;
        case HarmParamAuto:
            autotune = (int) clamp(value,0.f,1.f);
            //fprintf(stderr, "autotune: %d\n", autotune);
            break;
        case HarmParamAutoStrength:
            corr_strength = clamp(value, 0.f, 1.f);
            break;
        case HarmParamMidi:
            midi_enable = (int) clamp(value,0.f,1.f);
            printf("set midi_enable to %d\n", midi_enable);
            break;
        case HarmParamMidiLink:
            midi_link = (int) clamp(value,0.f,1.f);
            break;
        case HarmParamMidiVelIgnore:
            midi_ignore_velocity = (int) clamp(value,0.f, 1.f);
            break;
        case HarmParamMidiKeyCC:
            midi_keycenter_cc = (int) clamp(value,0.f,127.f);
            break;
        case HarmParamMidiKeyCcOffset:
            midi_keycenter_cc_offset = (int) clamp(value, 0.f, 127.f);
            break;
        case HarmParamMidiQualCC:
            midi_keyquality_cc = (int) clamp(value,0.f,127.f);
            break;
        case HarmParamMidiQualCcOffset:
            midi_keyquality_cc_offset = (int) clamp(value, 0.f, 127.f);
            break;
        case HarmParamMidiNvoiceCC:
            midi_nvoices_cc = (int) clamp(value,0.f,127.f);
            break;
        case HarmParamMidiNvoiceCcRange:
            midi_nvoices_range = (int) clamp(value,0.f,127.f);
            break;
        case HarmParamMidiInvCC:
            midi_inversion_cc = (int) clamp(value,0.f,127.f);
            break;
        case HarmParamMidiInvCcRange:
            midi_inversion_range = (int) clamp(value,0.f,127.f);
            break;
        case HarmParamMidiFreezeCC:
            midi_freeze_cc = (int) clamp(value,0.f,127.f);
            break;
        case HarmParamMidiPC:
            midi_program_change_enable = (int) clamp(value,0.f,1.f);
            break;
        case HarmParamMidiMelOut:
            midi_transmit_melody = (int) clamp(value,0.f,1.f);
            break;
        case HarmParamMidiHarmOut:
            midi_transmit_harmony = (int) clamp(value,0.f,1.f);
            break;
        case HarmParamMidiLegato:
            //midi_legato = (int) clamp(value, 0.f, 1.f);
            break;
        case HarmParamMidiPedalFcn:
            midi_pedal_fcn = (int) clamp(value, 0.f, 2.f);
            break;
        case HarmParamMidiPedalInv:
            midi_pedal_inv = (int) clamp(value, 0.f, 1.f);
            break;
        case HarmParamTriad:
            triad = (int) clamp(value,-1.f,30.f);
            break;
        case HarmParamBypass:
            bypass = (int) clamp(value,0.f,1.f);
            printf("set bypass to %d\n", bypass);
            break;
        case HarmParamHgain:
            harmgain_target = clamp(value, 0.f, 2.f);
            break;
        case HarmParamVgain:
            voicegain_target = clamp(value, 0.f, 2.f);
            break;
        case HarmParamDryMix:
            dry_mix = clamp(value, 0.f, 1.f);
            break;
        case HarmParamSpeed:
            speed = clamp(value, 0.f, 1.f);
            break;
        case HarmParamTuning:
            baseTuning = value;
            break;
        case HarmParamThreshold:
            threshold = value;
            break;
        case HarmParamStereo:
            stereo_mode = value;
            break;
        case HarmParamSynth:
            synth_enable = value;
            fprintf(stderr, "synth_enable = %d\n", synth_enable);
            break;
        case HarmParamAlgorithm:
            algorithm = value;
            break;
        case HarmParamVibrato:
            vibrato = value;
            break;
        case HarmParamFreeze:
            for (int k = 0; k < n_channels; k++)
                freezers[k].freeze((bool) (value > 0));
            break;
        case HarmParamGateThresh:
            gate_thresh = value;
            noise_gate.set_thresh_db(gate_thresh);
            break;
        case HarmParamLoop:
            //int old_mode = loop_mode;
            looper->setMode(static_cast<Looper::loopMode>(value));
            break;
            
        case HarmParamInterval:
        default:
            int addr = (int) address - (int) HarmParamInterval;
            int scale_degree = addr / 4;
            
            chord_ratio_t * table = major_chord_table;
            if (scale_degree > 11)
                table = minor_chord_table;
            if (scale_degree > 23)
                table = blues_chord_table;
            
            float * ratios = (float *) &table[scale_degree%12];
            //fprintf(stderr, "addr = %d\n", addr);
            interval_offsets[addr] = (int) value;
            ratios[addr & 0x3] = intervals[(int) value];
            if (harmonizer) {
                harmonizer->setInterval(static_cast<std::size_t>(addr), static_cast<int>(value));
            }
            break;
    }
}

param_value_t HarmonizerDSPKernel::getParameter(param_address_t address) {
    switch (address) {
        case HarmParamKeycenter:
            return (float) root_key;
        case HarmParamInversion:
            return (float) inversion;
        case HarmParamNvoices:
            return (float) harmonizer->getAutoVoiceCount();
        case HarmParamAuto:
            return (float) autotune;
        case HarmParamAutoStrength:
            return (float) corr_strength;
        case HarmParamMidi:
            return (float) midi_enable;
        case HarmParamMidiLink:
            return (float) midi_link;
        case HarmParamMidiLegato:
            return (float) 0; //midi_legato;
        case HarmParamMidiVelIgnore:
            return (float) midi_ignore_velocity;
        case HarmParamMidiKeyCC:
            return (float) midi_keycenter_cc;
        case HarmParamMidiKeyCcOffset:
            return (float) midi_keycenter_cc_offset;
        case HarmParamMidiQualCC:
            return (float) midi_keyquality_cc;
        case HarmParamMidiQualCcOffset:
            return (float) midi_keyquality_cc_offset;
        case HarmParamMidiNvoiceCC:
            return (float) midi_nvoices_cc;
        case HarmParamMidiNvoiceCcRange:
            return (float) midi_nvoices_range;
        case HarmParamMidiInvCC:
            return (float) midi_inversion_cc;
        case HarmParamMidiInvCcRange:
            return (float) midi_inversion_range;
        case HarmParamMidiFreezeCC:
            return (float) midi_freeze_cc;
        case HarmParamMidiPC:
            return (float) midi_program_change_enable;
        case HarmParamMidiMelOut:
            return (float) midi_transmit_melody;
        case HarmParamMidiHarmOut:
            return (float) midi_transmit_harmony;
        case HarmParamMidiPedalFcn:
            return (float) midi_pedal_fcn;
        case HarmParamMidiPedalInv:
            return (float) midi_pedal_inv;
        case HarmParamTriad:
            return (float) triad;
        case HarmParamBypass:
            return (float) bypass;
        case HarmParamHgain:
            return harmgain_target;
        case HarmParamVgain:
            return voicegain_target;
        case HarmParamDryMix:
            return dry_mix;
        case HarmParamSpeed:
            return speed;
        case HarmParamTuning:
            return baseTuning;
        case HarmParamThreshold:
            return threshold;
        case HarmParamStereo:
            return stereo_mode;
        case HarmParamSynth:
            return synth_enable;
        case HarmParamAlgorithm:
            return algorithm;
        case HarmParamGateThresh:
            return noise_gate.get_thresh_db();
        case HarmParamVibrato:
            return vibrato;
        case HarmParamLoop:
            return (float) loop_mode;
        case HarmParamFreeze:
            return freezers[0].frozen();
        case HarmParamInterval:
        default:
            int addr = (int) address - (int) HarmParamInterval;
            int scale_degree = addr / 4;
            
            chord_ratio_t * table = major_chord_table;
            if (scale_degree > 11)
                table = minor_chord_table;
            if (scale_degree > 23)
                table = blues_chord_table;
            
            float * ratios = (float *) &table[scale_degree%12];
            return interval_offsets[addr];;
            return round(log2(ratios[addr & 0x3])*12);
    }
}

void HarmonizerDSPKernel::setPreset(int preset_ix_)
{
    int chords_intervals[] = {
            0,4,7,12, -1,3,6,11, 2,5,10,14, 1,4,9,13, 0,3,8,12, -1,2,7,11, 1,6,10,13, 0,5,9,12, -1,4,8,11, 0,3,7,10, 2,6,9,14, 1,5,8,13, // major
            0,3,7,12, -1,2,6,11, 1,5,10,13, 0,4,9,12, -1,3,8,11, -1,2,7,10, 1,6,9,13, 0,5,8,12, 0,4,7,11, 0,3,6,10, 0,5,9,14, 1,4,8,13, // minor
            0,4,10,12, -1,3,9,11, -2,2,8,10, 1,4,7,9, 0,3,6,8, 2,5,7,11, 1,4,6,10, 0,3,5,9, -1,2,4,8, 1,3,7,10, 0,2,6,9, -1,1,5,8, //dom
    };
    int diatonic_intervals[] = {
            0,4,7,12, -1,3,6,11, 0,5,10,12, 1,4,9,13, 0,3,8,12, 0,2,7,11, 1,6,10,13, 0,5,9,12, -1,4,8,11, 0,3,7,12, 1,2,6,13, 0,1,5,12, // major
            0,3,7,12, -1,2,6,11, 0,5,10,12, 0,4,9,12, -1,3,8,11, -2,2,7,10, 1,6,9,13, 0,5,8,12, -1,4,7,11, 0,3,6,10, 2,5,9,14, 1,4,8,13, // minor
            0,4,7,10, -1,3,9,11, -2,2,8,10, 1,4,7,9, 0,3,6,8, 2,5,7,11, 1,4,6,10, 0,3,5,9, -1,2,4,8, 1,3,7,10, 0,2,6,9, -1,1,5,8, //dom
    };
    int chromatic_intervals[] = {
            0,4,7,12, 0,3,6,12, 0,3,7,12, 0,3,9,12, 0,3,8,12, 0,4,7,12, 0,3,9,12, 0,5,9,12, 0,4,8,12, 0,5,8,12, 0,4,7,12, 0,3,6,12, // major
            0,3,7,12, 0,4,7,12, 0,3,9,12, 0,4,9,12, 0,3,8,12, 0,3,7,12, 0,6,9,12, 0,4,7,12, 0,4,7,12, 0,3,6,12, 0,5,9,12, 0,3,7,12, // minor
            0,4,7,12, 0,3,9,12, 0,2,8,12, 0,4,7,12, 0,3,6,12, 0,5,7,12, 0,4,6,12, 0,3,5,12, 0,2,4,12, 0,3,7,12, 0,2,6,12, 0,1,5,12, //dom
    };
    int barbershop_intervals[] = {
            0,4,7,12, 0,3,5,9, 0,3,5,9, 0,3,6,9, 0,3,8,12, 0,2,6,9, 0,3,5,9, 0,5,9,12, 0,3,6,9, 0,3,5,9, 0,3,6,9, 0,3,6,8, // major
            0,3,7,12, 0,4,7,10, 0,3,5,9, 0,4,9,12, 0,3,6,8, 0,3,7,9, 0,3,6,8, 0,5,8,12, 0,4,7,10, 0,3,6,10, 0,4,7,10, 0,3,6,8, // minor
            0,4,7,10, 0,3,6,9, 0,3,5,9, 0,3,6,9, 0,3,6,8, 0,2,6,9, 0,3,5,9, 0,3,5,9, 0,2,4,8, 0,3,6,9, 0,2,6,9, 0,4,7,10 //dom
    };
    int justmidi_intervals[] = {
        -12,-12,-12,-12, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, // "major"
        0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, // "minor"
        0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, // "dom"
    };
    int bohemian_intervals[] = {
        0,4,7,9, 0,3,6,8, 0,3,7,10, 0,3,6,9, 0,3,5,8, 0,4,7,9, 0,3,6,9, 0,2,5,9, 0,3,6,9, 0,3,5,8, 0,2,6,9, 0,1,5,8, // major
        0,3,7,10, 0,3,6,8, 0,3,6,9, 0,4,7,9, 0,3,5,8, 0,3,6,9, 0,3,6,9, 0,3,5,8, 0,3,6,9, 0,3,6,10, 0,4,7,10, 0,3,6,9, // minor
        0,4,7,10, 0,3,6,9, 0,3,5,8, 0,3,6,9, 0,3,6,8, 0,2,5,9, 0,3,5,9, 0,3,5,9, 0,2,6,9, 0,1,5,8, 0,2,6,9, 0,3,6,9 //dom
    };
    int bass_intervals[] = {
        -12,-12,-12,-12, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, // "major"
        0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, // "minor"
        0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, // "dom"
    };
    int fourths_intervals[] = {
        0,-5,7,12, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, // "major"
        0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, // "minor"
        0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, // "dom"
    };
    int modes_intervals[] = {
        0,4,7,11, 0,3,6,10, 0,3,7,10, 0,3,6,9, 0,3,7,10, 0,4,7,11, 0,3,6,10, 0,4,7,10, 0,4,8,11, 0,3,7,10, 0,4,7,11, 0,3,6,10, // major
        0,3,7,11, 0,3,7,10, 0,3,7,10, 0,4,8,11, 0,4,7,10, 0,4,7,10, 0,4,7,10, 0,4,7,10, 0,4,7,11, 0,3,6,10, 0,3,6,10, 0,3,6,10, // minor
        0,4,7,10, 0,3,7,10, 0,3,7,10, 0,3,6,10, 0,3,6,10, 0,4,7,11, 0,3,7,10, 0,3,7,10, 0,3,7,10, 0,3,7,10, 0,4,7,11, 0,4,7,10 //dom
    };

    int * intervals = chords_intervals;

    preset_ix = preset_ix_;
    switch (preset_ix_)
    {
        case HarmPresetChords:
            intervals = chords_intervals;
            setParameter(HarmParamNvoices,4);
            setParameter(HarmParamInversion,3);
            setParameter(HarmParamAuto,0);
            setParameter(HarmParamTriad,-1);
            break;
        case HarmPresetDiatonic:
            intervals = diatonic_intervals;
            setParameter(HarmParamNvoices,4);
            setParameter(HarmParamInversion,3);
            setParameter(HarmParamAuto,0);
            setParameter(HarmParamTriad,-1);
            break;
        case HarmPresetChromatic:
            intervals = chromatic_intervals;
            setParameter(HarmParamNvoices,4);
            setParameter(HarmParamInversion,3);
            setParameter(HarmParamAuto,0);
            setParameter(HarmParamTriad,-1);
            break;
        case HarmPresetBarbershop:
            intervals = barbershop_intervals;
            setParameter(HarmParamNvoices,4);
            setParameter(HarmParamInversion,2);
            setParameter(HarmParamAuto,0);
            setParameter(HarmParamTriad,-1);
            break;
        case HarmPresetMIDI:
            intervals = justmidi_intervals;
            setParameter(HarmParamNvoices,1);
            setParameter(HarmParamInversion,0);
            setParameter(HarmParamAuto,0);
            setParameter(HarmParamTriad,0);
            break;
        case HarmPresetBohemian:
            intervals = bohemian_intervals;
            setParameter(HarmParamNvoices,4);
            setParameter(HarmParamInversion,4);
            setParameter(HarmParamAuto,0);
            setParameter(HarmParamTriad,-1);
            break;
        case HarmPresetBass:
            intervals = bass_intervals;
            setParameter(HarmParamNvoices,1);
            setParameter(HarmParamInversion,1);
            setParameter(HarmParamAuto,0);
            setParameter(HarmParamTriad,0);
            break;
        case HarmPreset4ths:
            intervals = fourths_intervals;
            setParameter(HarmParamNvoices,1);
            setParameter(HarmParamInversion,1);
            setParameter(HarmParamAuto,0);
            setParameter(HarmParamTriad,0);
            break;
        case HarmPresetModes:
            intervals = modes_intervals;
            setParameter(HarmParamNvoices,4);
            setParameter(HarmParamInversion,3);
            setParameter(HarmParamAuto,0);
            setParameter(HarmParamTriad,-1);
            break;
        default:
            preset_ix = 0;
            return;
    }

    for (int i = 0; i < 144; i++)
    {
        setParameter(HarmParamInterval+i,(float) intervals[i]);
    }

}
int HarmonizerDSPKernel::getPreset()
{
    return preset_ix;
}

#ifdef __APPLE__
void HarmonizerDSPKernel::startRamp(AUParameterAddress address, AUValue value, AUAudioFrameCount duration) {
    // just set it now
    setParameter(address, value);
    return;
}

void HarmonizerDSPKernel::setBuffers(AudioBufferList* inBufferList, AudioBufferList* outBufferList) {

    for (int k = 0; k < n_channels; k++)
    {
        in_buffers[k] = (float*) inBufferList->mBuffers[k].mData;
        out_buffers[k] = (float*) outBufferList->mBuffers[k].mData;
    }
}

void HarmonizerDSPKernel::set_T(float t) {
    this->T = t;
}

void HarmonizerDSPKernel::handleMIDIEvent(midi_event_t const& midiEvent) {
    uint8_t status = midiEvent.data[0] & 0xF0;
    //uint8_t channel = midiEvent.data[0] & 0x0F; // works in omni mode.
    static int key_quality = 0;
    
//    if (channel != 0)
//        return;
    
    switch (status) {
        case 0x80 : { // note off
            uint8_t note = midiEvent.data[1];
            if (note > 127) break;
            remnote((int)note);
            break;
        }
        case 0x90 : { // note on
            uint8_t note = midiEvent.data[1];
            uint8_t veloc = midiEvent.data[2];
            D(std::cout << "note #" << (int) note << ": " << (int) veloc << "\n");
            if (note > 127 || veloc > 127) break;
            if (veloc == 0)
                remnote((int)note);
            else
                addnote((int)note,(int)veloc);
            break;
        }
        case 0xB0 : { // control
            uint8_t num = midiEvent.data[1];
            uint8_t val = midiEvent.data[2];
            D(std::cout << "cc #" << num << ": " << val << "\n");
            if (num == 11)
            {
                midigain = (float) val / 64.0;
            }
            if (num == 64)
            {
                if ((val > 0) ^ midi_pedal_inv)
                {
                    pedal_down();
                }
                else
                {
                    pedal_up();
                }
            }
            if (num == 123) { // all notes off

            }
            if (num == midi_keycenter_cc) // Keycenter change
            {
                if ((val >= midi_keycenter_cc_offset) && (val < (midi_keycenter_cc_offset + 12)))
                {
                    root_key = val - midi_keycenter_cc_offset + key_quality * 12;
                    key_quality = 0;
                    fprintf(stderr,"root key = %d\n", root_key);
                }
            }
            if (num == midi_keyquality_cc)
            {
                if (val > midi_keyquality_cc_offset)
                {
                    key_quality = val - midi_keyquality_cc_offset;
                    if (key_quality > 2)
                        key_quality = 2;
                }
            }
            if (num == midi_nvoices_cc)
            {
                setParameter(HarmParamNvoices, 1 + (val * 4) / 127);
            }
            if (num == midi_inversion_cc)
            {
                setParameter(HarmParamInversion, (val * n_auto) / 127);
                //inversion = 1 + (val * (n_auto-1)) / 127;
            }
            if (num == midi_freeze_cc)
            {
                setParameter(HarmParamFreeze, val == 0 ? 0.0 : 1.0);
            }
            if (num >= 20 && num <= 31)
            {
                
            }
            break;
        }
        case 0xC0 : {
            //uint8_t num = midiEvent.data[1];
            patch_number = midiEvent.data[1];
            
            break;
        }
    }
}

#endif

void HarmonizerDSPKernel::process(frame_count_t frameCount, frame_count_t bufferOffset) {
    //fprintf(stderr, "process!\n");

    int channelCount = n_channels;
    sample_count += frameCount;
    n_output_events = 0;
    
    float* in  = in_buffers[0] + bufferOffset;
    float* out[2] = {out_buffers[0] + bufferOffset, out_buffers[0] + bufferOffset};
    
    if (channelCount > 1)
    {
        out[1] = out_buffers[1] + bufferOffset;
    }
    
    //int n_computed = 0;
    
//    // For each sample.
//    for (int frameIndex = 0; frameIndex < frameCount; ++frameIndex)
//    {
//        raw_buffer.pushValue(in[frameIndex]);
//        filtered_buffer.pushValue(filter.compute_one(in[frameIndex]));
//        noise_gate.compute_one(in[frameIndex]);
//                
//        if (bypass)
//        {
//            out[0][frameIndex] = in[frameIndex] / 2;
//            out[1][frameIndex] = out[0][frameIndex];
//            continue;
//        }
//        else
//        {
//            out[0][frameIndex] = out[1][frameIndex] = 0;
//        }
//        
//        if (--rcnt == 0)
//        {
//            rcnt = 256;
//            float p = pitchEstimator.estimate(filtered_buffer);
//
//            if (p > 0.f)
//                T = p;
//            
//            voiced = (p != 0);
//            
//            if (synth_enable)
//            {
//                //get_model(cix - 3*T);
//                //get_minphase_pulse(-nfft);
//            }
//            
//            update_voices();
//        }
//        
//        if (pitchMarker.findMark(T,0.25) == true)
//        {
//            //std::cerr << "d=" << raw_buffer.getWriteIndex() - pitchMarker.mark << std::endl;
//            // set synth source
//            float * grainbuf = raw_buffer.getContiguous(pitchMarker.mark-T);
//            
//            for (int k = 0; k < nvoices; k++){
//                psolaVoices[k].setGrainSource(grainbuf, 0, 2*T);
//                //psolaVoices[k].win_enable = !synth_enable;
//            }
//                        
//            int nn = frameIndex - n_computed;
//            psola(out[0]+n_computed, out[1]+n_computed, nn);
//            n_computed += nn;
//        }
//                
//        float x = autotune ? simpleVoices[0].computeOne() : in[frameIndex];
//        voicegain += .001 * sgn(voicegain_target - voicegain);
//
//        
//        
//        out[0][frameIndex] = x * voicegain/2; //in[frameIndex] * voicegain/2;
//        if (stereo_mode != StereoModeSplit){
//            out[1][frameIndex] = out[0][frameIndex];
//        }
//        
//    }
//    
//    int nn = frameCount - n_computed;
//    psola(out[0]+n_computed, out[1]+n_computed, nn);
//    
//    filter.clear_bad();
    
    memcpy(out[0], in, frameCount * sizeof(float));
    memcpy(out[1], in, frameCount * sizeof(float));
    
    harmonizer->compute(in, out, channelCount, frameCount);
    
    midi_note_number = T_to_midi_note(harmonizer->getT(), sampleRate);
    harmonizer->getAutoNotes(ui_voice_notes, n_auto);
            
    // compute freezing stuff
    for (int ch = 0; ch < n_channels; ch++) {
        for (int frameIndex = bufferOffset; frameIndex < frameCount + bufferOffset; ++frameIndex)
        {
            float y = freezers[ch].process_frame(out_buffers[ch][frameIndex]);
            out_buffers[ch][frameIndex] += y;
        }
    }
    
    
    // compute looping stuff
    looper->compute(out_buffers, bufferOffset, frameCount);
}

void HarmonizerDSPKernel::psola(float *out, float *out2, int n)
{
    for (int sample_ix = 0; sample_ix < n; sample_ix++)
    {
        int first_psola_voice = 1;
        
        voicegain_target = dry_mix;
        if (triad >= 0)
        {
            first_psola_voice = 0;
            voicegain_target = 0;
        }
                
        // Ramp Gain
        harmgain += .001 * sgn(harmgain_target - harmgain);
        
        for (int vix = first_psola_voice; vix < nvoices; vix++)
        {
            float target_gain = noise_gate.get_gain() *
                ((voiced && psolaVoices[vix].midinote > 0) ? 1.0 : 0.0);
            
            psolaVoices[vix].setGain(target_gain);
            
            if (vix >= n_auto && !midi_enable)
                continue;
            
//            float u = (algorithm == AlgorithmPSOLA) ? psolaVoices[vix].synthesizeOne() :
//                simpleVoices[vix].computeOne();
            
            //volatile float uu = psolaVoices[vix].synthesizeOne();
            volatile stereoSample u = psolaVoices[vix].render();
            
            u.l *= harmgain;
            u.r *= harmgain;
            switch (stereo_mode)
            {
                case StereoModeNormal:
                    out[sample_ix] += u.l;
                    out2[sample_ix] += u.r;
                    break;
                case StereoModeMono:
                    out[sample_ix] += u.l+u.r;
                    out2[sample_ix] = out[sample_ix];
                    break;
                case StereoModeSplit:
                    out2[sample_ix] += u.l+u.r;
                    break;
            }
        }
    }
}

void HarmonizerDSPKernel::addnote(int note, int vel)
{
    int min_dist = 129;
    int dist;
    int min_ix = -1;
    int vel_ = midi_ignore_velocity ? 100 : vel;
    
    keys_down[note] = 1;
    
//        if (!midi_enable)
//        {
//            return;
//        }
    
    // look for one with the same note and take that if we can, or the empty one with
    // the closest last note to the one we want
    for (int k = n_auto; k < nvoices; k++)
    {
        dist = (psolaVoices[k].midinote == note) ? 0 : abs(psolaVoices[k].midinote.prev() - note);
        if ((dist < min_dist && psolaVoices[k].midinote < 0))
        {
            min_ix = k;
            min_dist = dist;
        }
    }

    psolaVoices[min_ix].setMidiNote(note, vel_);
    
    harmonizer->addMidiNote(note, vel_);
    update_midi();
}

void HarmonizerDSPKernel::remnote(int note)
{
    keys_down[note] = 0;
    
    if (midi_pedal)
    {
        return;
    }
    
    for (int k = n_auto; k < nvoices; k++)
    {
        if (psolaVoices[k].midinote == note)
        {
            psolaVoices[k].midinote = -1;
        }
    }
    
    update_midi();
}

void HarmonizerDSPKernel::pedal_down()
{
    switch (midi_pedal_fcn) {
        case PedalFreeze:
            for (SpectralProcessor& freezer : freezers)
                freezer.freeze(true);
            break;
        case PedalNotes:
            midi_pedal = 1;
    }
}
    
void HarmonizerDSPKernel::pedal_up()
{
    switch (midi_pedal_fcn) {
        case PedalFreeze:
            for (SpectralProcessor& freezer : freezers)
                freezer.freeze(false);
            break;
        case PedalNotes:
            midi_pedal = 0;

            for (int k = n_auto; k < nvoices; k++)
            {
                if (!keys_down[psolaVoices[k].midinote])
                {
                    psolaVoices[k].midinote = -1;
                }
            }
            break;
    }
}

void HarmonizerDSPKernel::analyze_harmony(void)
{
    //float intervals[127];
    int octave[12];
    
    memset(octave, 0, 12 * sizeof(int));
    
    int n = 0;
    for (int j = n_auto; j < nvoices; j++)
    {
        if (psolaVoices[j].midinote < 0)
            continue;
        
        midinotes[n++] = (float) psolaVoices[j].midinote;
        
        octave[psolaVoices[j].midinote % 12] = 1;
    }
    if (n == 0)
        return;

#ifdef __APPLE__
    vDSP_vsort(midinotes, (vDSP_Length) n, 1);
#else
    std::sort(midinotes, midinotes + n);
#endif
    for (int j = 0; j < n; j++)
    {
        // ignore doubles
        int dbl = 0;
        for (int k = 0; k < j; k++)
        {
            if (((int)(midinotes[j] - midinotes[k]) % 12) == 0)
            {
                dbl = 1; break;
            }
        }
        
        if (dbl)
            continue;
        
        int ix = (int) midinotes[j] % 12;
        
        if (octave[(ix+4)%12] && octave[(ix+10)%12])
        {
            // 7th
            root_key = ix+ 24;
            break;
        }
        
        if (octave[(ix+3)%12] && octave[(ix+6)%12])
        {
            // dim
            root_key = ((ix+2)%12) + 24;
            break;
        }
        
        if (octave[(ix+3)%12] && octave[(ix+7)%12])
        {
            //min
            root_key = ix + 12;
            break;
        }
        
        if (octave[(ix+4)%12] && octave[(ix+7)%12] && octave[(ix+11)%12])
        {
            root_key = ((ix+7)%12);
            break;
        }
        
        if (octave[(ix+4)%12] && octave[(ix+7)%12])
        {
            //maj
            root_key = ix;
            break;
        }
        if (octave[(ix+4)%12])
        {
            //maj
            root_key = ix;
            break;
        }
        if (octave[(ix+3)%12])
        {
            //min
            root_key = ix+12;
            break;
        }
        
        if (octave[(ix+4)%12] && octave[(ix+8)%12])
        {
            //aug
        }
   
    }
}

void HarmonizerDSPKernel::send_note_on(int nn, int vel)
{
    if (n_output_events > max_output_events -1)
        return;
    // queue MIDI note on messages
    output_events[n_output_events].length = 3;
    output_events[n_output_events].data[0] = 0x90;
    output_events[n_output_events].data[1] = nn;
    output_events[n_output_events].data[2] = vel;
    n_output_events++;
}
    
void HarmonizerDSPKernel::send_note_off(int nn, int vel)
{
    if (n_output_events > max_output_events -1)
        return;
    
    // queue MIDI note off
    output_events[n_output_events].length = 3;
    output_events[n_output_events].data[0] = 0x80;
    output_events[n_output_events].data[1] = nn;
    output_events[n_output_events].data[2] = vel;
    n_output_events++;
}

void HarmonizerDSPKernel::update_midi() {
    if (midi_link && (sample_count - midi_changed_sample_num) > (int) sampleRate / 50)
    {
        analyze_harmony();
    }
    midi_changed_sample_num = sample_count;
}

void HarmonizerDSPKernel::calculate_voices()
{
    
}

void HarmonizerDSPKernel::list_intervals()
{
    for (int k = 0; k < n_auto; k++) {
        std::cout << "voice " << k << ": \n"
        "\tmidi note: " << psolaVoices[k].midinote << "\n";
    }
}

void HarmonizerDSPKernel::set_voiced(bool voiced) {
    this->voiced = (int) voiced;
}

void HarmonizerDSPKernel::update_voices (void)
{
    static int was_voiced = 0;
    
    if (!voiced)
    {
        if (was_voiced)
        {
            //send_note_off(midi_note_number, 100);
            for (int k = 0; k < n_auto; k++)
            {
                send_note_off(ui_voice_notes[k], 100);
            }
        }
        
        midi_note_number = -1;
        
        for (int k = 0; k < n_auto; k++)
        {
            psolaVoices[k].midinote = ui_voice_notes[k] = -1;
        }

        was_voiced = 0;
        return;
    }
    
    // half steps
    float note_f = 12.0 * log2f (sampleRate / (T * baseTuning));
    // integer half steps from a440 (baseTuning)
    int nn = (int) round(note_f);
    // microtonal error
    float err = note_f - (float) nn;
    
    // convert to midi note number, A440 = midi note 69
    midi_note_number = nn + 69;
    
    // hysteresis:
    // we want to avoid flipping back and forth
        
    int root = root_key % 12;
    int quality = root_key / 12;
    
    // this is just (old_midi_note_number - root) % 12 -- why not current?
    //int interval = (last_nn + 69 - root) % 12;
    int interval = (midi_note_number - root) % 12;
    
    // convert interval table to midi notes for auto voices. Send MIDI off/on messages if note changes
    for (int k = 0; k < n_auto; k++)
    {
        int midinote = midi_note_number + interval_offsets[k + (interval*4) + (quality*48)];
        if (k > inversion){
            midinote -= 12;
        }
        
        psolaVoices[k].setMidiNote(midinote, 65);
        
        if (psolaVoices[k].midinote.diff() != 0)
        {
            send_note_on(ui_voice_notes[k], 100);
            if (psolaVoices[k].midinote.prev() >= 0)
                send_note_off(psolaVoices[k].midinote.prev(), 100);
        }
                
        // let the UI know
        ui_voice_notes[k] = psolaVoices[k].midinote;
    }
    
    return;
}


float HarmonizerDSPKernel::loopPosition() {
    return looper->position();
}

int HarmonizerDSPKernel::setLoopMode(int mode) {
    looper->setMode(static_cast<Looper::loopMode>(mode));
    return mode;
}
int HarmonizerDSPKernel::getLoopMode() {
    return looper->loop_mode;
}

