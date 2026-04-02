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


void HarmonizerDSPKernel::init(int inChannels, int outChannels, double inSampleRate) {
    n_channels = outChannels;
    fprintf(stderr,"**** init with %d channels! at %f Hz\n", n_channels, inSampleRate);
    int log2nfft = 11;
    int nfft = 1 << log2nfft;
    sampleRate = float(inSampleRate);
#ifdef __APPLE__
    fft_s = vDSP_create_fftsetup(log2nfft, 2);
    
    fft_in.realp = (float *) calloc(nfft, sizeof(float));
    fft_in.imagp = (float *) calloc(nfft, sizeof(float));
    
    fft_out.realp = (float *) calloc(nfft, sizeof(float));
    fft_out.imagp = (float *) calloc(nfft, sizeof(float)); 
    
    fft_out2.realp = (float *) calloc(nfft, sizeof(float));
    fft_out2.imagp = (float *) calloc(nfft, sizeof(float));
    
    fft_buf.realp = (float *) calloc(nfft, sizeof(float));
    fft_buf.imagp = (float *) calloc(nfft, sizeof(float));
    
    A_model.realp = (float *) calloc(100, sizeof(float));
    A_model.imagp = (float *) calloc(100, sizeof(float));
    
    Hann.realp = (float *) calloc(nfft, sizeof(float));
    Hann.imagp = (float *) calloc(nfft, sizeof(float));
    
    spec_env.realp = (float *) calloc(nfft, sizeof(float));
    spec_env.imagp = (float *) calloc(nfft, sizeof(float));
    fft_in.realp[0] = 1.0;
    
#else
    fft_s = kiss_fft_alloc(nfft,0,NULL,0);
    ifft_s = kiss_fft_alloc(nfft,1,NULL,0);
    fft_in = (kiss_fft_cpx *) calloc(nfft, sizeof(kiss_fft_cpx));
    fft_out = (kiss_fft_cpx *) calloc(nfft, sizeof(kiss_fft_cpx));
    fft_out2 = (kiss_fft_cpx *) calloc(nfft, sizeof(kiss_fft_cpx));
    Hann = (kiss_fft_cpx *) calloc(nfft, sizeof(kiss_fft_cpx));
    spec_env = (kiss_fft_cpx *) calloc(nfft, sizeof(kiss_fft_cpx));
#endif
    
    ncbuf = 2*nfft;
    cbuf = (float *) calloc(ncbuf + 3, sizeof(float));
    fbuf = (float *) calloc(ncbuf + 3, sizeof(float));
    
    filter = ButterworthFilter();
    
    cmdf = (float *) calloc(maxT, sizeof(float));
    n_snr = sampleRate/100;
    snr_buf = (float *) calloc(n_snr, sizeof(float));
    
    synth_pulse = (float **) calloc(n_synth_pulse, sizeof(float *));
    
    for (int k = 0; k < n_synth_pulse; k++)
        synth_pulse[k] = (float *) calloc(2*maxT, sizeof(float));
    
    lms_h = (float *) calloc(lms_n, sizeof(float));
    
    nvoices = 16;
    voices = (voice_t *) calloc(nvoices, sizeof(voice_t));
    voice_ix = 1;
    
    in_buffers = (float **) calloc(inChannels, sizeof(float *));
    out_buffers = (float **) calloc(outChannels, sizeof(float *));
    
    fd_lpf = (float *) calloc(nfft, sizeof(float));
    
    int fc = (int) ((4000. / sampleRate) * nfft); // cutoff frequency chosen to preserve formants
    int bw = fc / 2;
    //fprintf(stderr, "fc = %d\n", fc);
    for (int k = 1; k < nfft/2; k++)
    {
        if (k <= fc-bw)
            fd_lpf[k] = fd_lpf[nfft-k] = 1.0;
        else if (k < fc+bw)
            fd_lpf[k] = fd_lpf[nfft-k] = 0.5 - 0.5 * cos(M_PI * (k - (fc - bw)) / (2*bw + 1));
        //fprintf(stderr, "%f\n", fd_lpf[k]);
    }
    
//    loop_buf = (float **) calloc(2, sizeof(float *));
//    loop_max = (int) (sampleRate * 60);
//    for (int k = 0; k < 2; k++)
//    {
//        loop_buf[k] = (float *) calloc(loop_max, sizeof(float));
//    }
    
    looper = Looper(n_channels,60*sampleRate,(int)lrintf(0.05*sampleRate),(int)lrintf(0.05*sampleRate));
    
    for (int k = 0; k < nvoices; k++)
    {
        voices[k].midinote = -1;
        voices[k].lastnote = 0;
        voices[k].midinote_ = 69.0;
        voices[k].error = 0;
        voices[k].ratio = 1;
        voices[k].target_ratio = 1.0;
        voices[k].formant_ratio = 1;
        voices[k].nextgrain = 250;
        voices[k].ix1 = 0;
        voices[k].ix2 = 0;
        voices[k].xfade_ix = 0;
        voices[k].xfade_dur = 0;
        voices[k].vibrato_rate = 4.0 + ( (float)rand( ) / (float)RAND_MAX ) * 1.0;
        voices[k].vibrato_amp = 0.01;
        voices[k].vib_phase = 0.0;
        
        voices[k].gain = 1.0;
        voices[k].target_gain = 1.0;
        
        if (k >= N_AUTO)
        {
            voices[k].pan = ((float)(k - 3) / (float)(nvoices - 3)) - 0.5;
            //voices[k].formant_ratio = ((float)(k - 3) / (float)(nvoices - 3)) + 0.5;
            //psolaVoices[k].enable = false;
        }
        
        psolaVoices[k].setGrainSource(fd_lpf, fc, 0);
        psolaVoices[k].T = 200;
        psolaVoices[k].ratio = 1.0;
    }
    
    voices[1].formant_ratio = 0.99;
    voices[2].formant_ratio = 1.01;
    
    voices[0].midinote = 0;
    
    ngrains = 5 * nvoices;
    grains = new grain_t[ngrains];
    
    for (int k = 0; k < ngrains; k++)
    {
        grains[k].size = -1;
        grains[k].start = -1;
        grains[k].ix = 0;
        grains[k].gain = 1;
    }
        
    // create grain window with zero-padded shoulders
    grain_window = new float[graintablesize + 6];
    memset(grain_window, 0, (graintablesize + 6) * sizeof(float));
    
    for (int k = 0; k < graintablesize; k++)
    {
        grain_window [3 + k] = 0.5 * (1 - cosf (2 * M_PI * k / graintablesize));
    }
    
    #ifdef __APPLE__
    int pulsewin = 20;
    for (int k = 0; k < pulsewin; k++)
    {
        fft_in.realp[k] = window_value((float)k/(pulsewin-1));
//            fprintf(stderr,"%f\n", fft_in.realp[k]);
    }
    
    vDSP_fft_zopt(fft_s, &fft_in, 1, &Hann, 1, &fft_buf, 11, 1);
    #else
    #endif
//        memcpy(cbuf, TestAudioData, 2048*sizeof(float));
//        T = 200;
//        get_minphase_pulse(0);
//        exit(127);
    
    fft_mag = new float[nfft];
    fft_mag_db = new float[nfft];
    
    // define equal tempered interval ratios
    
    intervals = interval_table + 24;
    
    for (int k = -23; k < 24; k++)
    {
        intervals[k] = powf(2.0, (float) (k) / 12);
    }
    
    memset(midinotes, 0, 128 * sizeof(float));
    memset(keys_down, 0, 128 * sizeof(int));

//        int chords_intervals[] = {0,4,7,12, -1,3,6,11, 2,5,10,14, 1,4,9,13, 0,3,8,12, -1,2,7,11, 1,6,10,13, 0,5,9,12, -1,4,8,11, 0,3,7,10, 2,6,9,14, 1,5,8,13, // major
//                           0,3,7,12, -1,2,6,11, 1,5,10,13, 0,4,9,12, -1,3,8,11, -1,2,7,10, 1,6,9,13, 0,5,8,12, 0,4,7,11, 0,3,6,10, 0,5,9,14, 1,4,8,13, // minor
//                           0,4,10,12, -1,3,9,11, -2,2,8,10, 1,4,7,9, 0,3,6,8, 2,5,7,11, 1,4,6,10, 0,3,5,9, -1,2,4,8, 1,3,7,10, 0,2,6,9, -1,1,5,8, //dom
//        };
//        for (int i = 0; i < 144; i++)
//        {
//            setParameter(HarmParamInterval+i,(float) chords_intervals[i]);
//        }
    
}

void HarmonizerDSPKernel::fini() {
#ifdef __APPLE__
    vDSP_destroy_fftsetup(fft_s);
    free(fft_in.realp);
    free(fft_in.imagp);
    free(fft_out.realp);
    free(fft_out.imagp);
    free(fft_out2.realp);
    free(fft_out2.imagp);
    free(fft_buf.realp);
    free(fft_buf.imagp);
    free(A_model.realp);
    free(A_model.imagp);
#else
    kiss_fft_free(fft_s);
    free(fft_in);
    free(fft_out);
    free(fft_out2);
    free(fft_buf);
#endif

    fprintf(stderr, "*** fini! ***\n");
    delete grain_window;
    delete fft_mag;
    delete fft_mag_db;
    delete grains;

    free(cbuf);
    free(voices);
    free(fd_lpf);
    
    free(in_buffers);
    free(out_buffers);
}

void HarmonizerDSPKernel::reset() {
    for (int k = 0; k < nvoices; k++)
    {
        voices[k].midinote = -1;
        voices[k].error = 0;
        voices[k].ratio = 1;
        voices[k].target_ratio = 1.0;
        voices[k].formant_ratio = 1;
        voices[k].nextgrain = 250;
        
        if (k >= 3)
        {
            voices[k].pan = ((float)(k - 3) / (float)(nvoices - 3)) - 0.5;
            //voices[k].formant_ratio = ((float)(k - 3) / (float)(nvoices - 3)) + 0.5;
        }
    }
    voice_ix = 1;
    
    cix = 0;
    rix = 0;
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
            break;
        case HarmParamNvoices:
            n_auto = (int) clamp(value,1.f,4.f);
            for (int k = n_auto; k < N_AUTO; k++)
            {
                voice_notes[k] = -1;
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
            midi_legato = (int) clamp(value, 0.f, 1.f);
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
        case HarmParamVibrato:
            vibrato = value;
            break;
        case HarmParamLoop:
            //int old_mode = loop_mode;
            looper.setMode(static_cast<Looper::loopMode>(value));
           
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
            return (float) n_auto;
        case HarmParamAuto:
            return (float) autotune;
        case HarmParamAutoStrength:
            return (float) corr_strength;
        case HarmParamMidi:
            return (float) midi_enable;
        case HarmParamMidiLink:
            return (float) midi_link;
        case HarmParamMidiLegato:
            return (float) midi_legato;
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
        case HarmParamMidiPC:
            return (float) midi_program_change_enable;
        case HarmParamMidiMelOut:
            return (float) midi_transmit_melody;
        case HarmParamMidiHarmOut:
            return (float) midi_transmit_harmony;
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
        case HarmParamVibrato:
            return vibrato;
        case HarmParamLoop:
            return (float) loop_mode;
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

void HarmonizerDSPKernel::handleMIDIEvent(midi_event_t const& midiEvent) {
    uint8_t status = midiEvent.data[0] & 0xF0;
    uint8_t channel = midiEvent.data[0] & 0x0F; // works in omni mode.
    static int key_quality = 0;
    
    if (channel != 0)
        return;
    
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
            //fprintf(stderr, "cc #%d: %d\n", num, val);
            if (num == 11)
            {
                midigain = (float) val / 64.0;
            }
            if (num == 64)
            {
                if ((float) (val > 0))
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
    float* out = out_buffers[0] + bufferOffset;
    float* out2 = out;
    
    if (channelCount > 1)
    {
        out2 = out_buffers[1] + bufferOffset;
    }
    
    int n_computed = 0;
    
    // For each sample.
    for (int frameIndex = 0; frameIndex < frameCount; ++frameIndex)
    {
        raw_buffer.pushValue(in[frameIndex]);
        filtered_buffer.pushValue(filter.compute_one(in[frameIndex]));
        noise_gate.compute_one(in[frameIndex]);
                
        if (bypass)
        {
            out[frameIndex] = in[frameIndex] / 2;
            out2[frameIndex] = out[frameIndex];
            continue;
        }
        else
        {
            out[frameIndex] = out2[frameIndex] = 0;
        }
        
        if (--rcnt == 0)
        {
            rcnt = 256;
            int oldT = T;
            float p = pitchEstimator.estimate(filtered_buffer);
            if (p == 0)
                T = oldT + 0.1 * (400 - oldT);
            else if (p < minT)
                T = minT;
            else
                T = p;
            
            voiced = (p != 0);
            
            if (synth_enable)
            {
                //get_model(cix - 3*T);
                get_minphase_pulse(-nfft);
            }
            
            update_voices();
        }
        
        if (pitchMarker.findMark(T,0.25))
        {
            //std::cerr << "d=" << raw_buffer.getWriteIndex() - pitchMarker.mark << std::endl;
            // set synth source
            for (int k = 0; k < nvoices; k++){
                float * buf = synth_enable ? synth_pulse[synth_pulse_ix] : raw_buffer.getContiguous(pitchMarker.mark-T);
                
                psolaVoices[k].setGrainSource(buf, 0, 2*T);
                //psolaVoices[k].win_enable = !synth_enable;
            }
                        
            int nn = frameIndex - n_computed;
            psola(out+n_computed, out2+n_computed, nn);
            n_computed += nn;
        }
        
        float at = simpleVoices[0].computeOne();
        
        float x = autotune ? at : in[frameIndex];
        out[frameIndex] = x * voicegain/2; //in[frameIndex] * voicegain/2;
        if (stereo_mode != StereoModeSplit){
            out2[frameIndex] = out[frameIndex];
        }
    }
    
    int nn = frameCount - n_computed;
    psola(out+n_computed, out2+n_computed, nn);
    
    filter.clear_bad();
            
    // compute looping stuff
    
    looper.compute(out_buffers, bufferOffset, frameCount);
}

inline float HarmonizerDSPKernel::measure_snr(float in)
{
    int ixp = snr_ix + 1;
    if (ixp >= n_snr) ixp = 0;
    
    snr_buf[snr_ix] = in*in;
    mean_sq += (snr_buf[snr_ix] - snr_buf[ixp])/n_snr;
    //fprintf(stderr, "%f\n", sumsq);
    snr_ix = ixp;
    
    rms = sqrtf(mean_sq);
    if ((nse_floor < rms) && !voiced)
    {
        nse_floor += (rms - nse_floor) * 0.000001;
    }
    else if (nse_floor > rms)
    {
        nse_floor -= (nse_floor - rms) * 0.001;
    }
    
    return rms/nse_floor;
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
        
        int nonvoiced_count = (int) (voicegain > 0);
        
        // Ramp Gain
        voicegain += .001 * sgn(voicegain_target - voicegain);
        harmgain += .001 * sgn(harmgain_target - harmgain);
        
        for (int vix = first_psola_voice; vix < nvoices; vix++)
        {
            voices[vix].target_gain =  noise_gate.get_gain() *
                (((voiced || nonvoiced_count == 0) && voices[vix].midinote > 0) ? (float) voices[vix].midivel / 127.f : 0.0);
            
            if (voices[vix].target_gain > 0)
                nonvoiced_count++;
            
            if (vix == 0 && (autotune || triad >= 0)) voices[vix].target_gain = dry_mix;
            //voices[vix].gain += .001 * sgn(voices[vix].target_gain - voices[vix].gain);
            voices[vix].gain = inc_to_target(voices[vix].gain, voices[vix].target_gain, 0.9, 0.001, -0.0004);
            
            if (voices[vix].gain < 0.0001){
                continue;
            }
            
            if (vix >= n_auto && !midi_enable)
                continue;
            
            float u = psolaVoices[vix].synthesizeOne() * harmgain * voices[vix].gain;
            //float u = simpleVoices[vix].computeOne() * harmgain * voices[vix].gain;
            
            switch (stereo_mode)
            {
                case StereoModeNormal:
                    out[sample_ix] += u * (voices[vix].pan + 1.0)/2.f;
                    out2[sample_ix] += u * (-voices[vix].pan + 1.0)/2.f;
                    break;
                case StereoModeMono:
                    out[sample_ix] += u/2.f;
                    out2[sample_ix] = out[sample_ix];
                    break;
                case StereoModeSplit:
                    out2[sample_ix] += u;
                    break;
            }
        }
    }
}

#ifdef __APPLE__
    
float HarmonizerDSPKernel::get_minphase_pulse(int start_ix)
{
    memset(fft_in.realp, 0, nfft * sizeof(float));
    memset(fft_in.imagp, 0, nfft * sizeof(float));
    static int count = 0;
    float * input = raw_buffer.getContiguousRelative(start_ix);
    for (int k = 0; k < nfft; k++)
    {
        //int ix = (start_ix + k) & cmask;
        fft_in.realp[k] = input[k] * window_value((float)k/(nfft));
    }
    
    vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out, 1, &fft_buf, 11, 1);
    // log(abs(fft))
    for (int k = 0; k < nfft; k++)
    {
        float r1,c1;
        r1 = fft_out.realp[k]; c1 = fft_out.imagp[k];
        float mag = (r1*r1 + c1*c1);
        
        fft_in.realp[k] = logf(mag+0.00001)/2; // factor of 2 is square root
        fft_in.imagp[k] = 0;
    }

    // cepstrum
    vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out, 1, &fft_buf, 11, -1);
    
    // fold and window
    fft_in.realp[0] = fft_out.realp[0]/nfft;

    int cutoff = 300;
    if (voiced)
    {
        cutoff = roundf(T * 0.7);
    }
    for (int k = 1; k < cutoff; k++)
    {
        fft_in.realp[k] = 0.f;
        fft_in.imagp[k] = 0;
    }
    for (int k = cutoff; k < nfft; k++)
    {
        fft_in.realp[k] = (fft_out.realp[k]*2)/nfft;
        fft_in.imagp[k] = 0;
    }
    
    // fft_out contains smooth envelope (fft)
    vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out, 1, &fft_buf, 11, 1);
    
    std::default_random_engine de(sample_count); //seed
    std::normal_distribution<float> nd(0, 1); //zero mean, 1 std
    
    for (int k = 0; k < nfft; k++)
    {
        fft_in.realp[k] = nd(de); fft_in.imagp[k] = 0;
    }
    
    vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out2, 1, &fft_buf, 11, 1);
    // fft_out2 contains random spectrum, convolve with envelope
    for (int k = 0; k < nfft; k++)
    {
        float ex = expf(fft_out.realp[k]);
        float flt = 1; //(1 + 9*(1 - fd_lpf[k]))/10;
        if (voiced)
        {
            flt = 1 - fd_lpf[k];
        }
        fft_in.realp[k] = flt * ex * fft_out2.realp[k]/nfft;
        fft_in.imagp[k] = flt * ex * fft_out2.imagp[k]/nfft;
    }
    
    vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out2, 1, &fft_buf, 11, -1);
    // fft_out2 now contains noise matching original sound spectrum
    
    for (int k = 0; k < nfft; k++)
    {
        spec_env.realp[k] = 0.9*spec_env.realp[k] + 0.1*fft_out.realp[k];
        spec_env.imagp[k] = 0.9*spec_env.imagp[k] + 0.1*fft_out.realp[k];

        float ex = expf(spec_env.realp[k]);
        //float rnd = nd(de);
        
        if (voiced)
        {
            ex *= fd_lpf[k];
        }
        
        float fr = fabs(1. - ((float) k / (float) (1 + nfft/2)));

        if (!voiced)
        {
            fr = 0.;
        }
        
        fft_in.realp[k] = ex * cos(spec_env.realp[k]);
        fft_in.imagp[k] = ex * sin(spec_env.realp[k]);
    }
    
    // get impulse response (ifft)
    vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out, 1, &fft_buf, 11, -1); //inverse
    
    if (++synth_pulse_ix >= n_synth_pulse)
        synth_pulse_ix = 0;
    
    if (voiced && count < 10)
    {
        count++;
    }
    
    if (!voiced && count > 0)
    {
        count--;
    }
    
    float mix = (float) count / 100.0;
    //fprintf(stderr, "mix = %f\n", mix);
    
    for (int k = 0; k < 2*maxT; k++)
    {
        synth_pulse[synth_pulse_ix][k] = mix*fft_out.realp[k]/sqrtf(nfft)/16;
        //synth_pulse[synth_pulse_ix][k] += fft_out2.realp[k]/16;
    }
            
    return 0.0;
}
    
float HarmonizerDSPKernel::get_model(int start_ix)
{
    memset(fft_in.realp, 0, nfft * sizeof(float));
    memset(fft_in.imagp, 0, nfft * sizeof(float));
    memset(A_model.realp, 0, 100 * sizeof(float));
    memset(A_model.imagp, 0, 100 * sizeof(float));
    
    float f0 = (float)nfft/(float)T;
    fprintf(stderr, "f0 = %f\n", f0);
    int h_ix = 0;
    float ratio = 1.;
    int miss_cnt = 0;
    float h_max = -500.;
    
    for (int k = 0; k < 3*T; k++)
    {
        int ix = (start_ix + k) & cmask;
        fft_in.realp[k] = cbuf[ix] * window_value((float)k/(3*T-1));
        //fprintf(stderr, "%f\n",fft_in.realp[k]);
    }
    
    vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out, 1, &fft_buf, 11, 1);
    
    for (int k = 0; k < nfft; k++)
    {
        fft_mag[k] = fft_out.realp[k]*fft_out.realp[k] + fft_out.imagp[k]*fft_out.imagp[k];
        fft_mag_db[k] = 10.*log10(fft_mag[k]);
    }
    
    float fix = 0;
    
    while (fix < nfft/8 && h_ix < 30)
    {
        fix += f0;
        h_ix++;
        //fprintf(stderr, "fix=%f\n",fix);
        float max = -500;;
        int max_ix = 0;
        
        for (int k = floorf(fix - f0/4); k < ceilf(fix + f0/4); k++)
        {
            if (fft_mag_db[k] > max)
            {
                max = fft_mag_db[k]; max_ix = k;
            }
        }
        
        if (max > fft_mag_db[max_ix+1] && max > fft_mag_db[max_ix-1])
        {
            float pv;
            float pix = max_ix + quadratic_peak(&pv, fft_mag_db, max_ix);
            
            if (pv > h_max) h_max = pv;
            
            //fprintf(stderr, "px = %f\n",pix);
            //fprintf(stderr, "pv = %f\n",pv);
            
            f0 = f0 * (1-ratio) + ratio * (float)(pix)/h_ix;
            A_model.realp[h_ix] = linear_interp(fft_out.realp, pix);
            A_model.imagp[h_ix] = linear_interp(fft_out.imagp, pix);
            fprintf(stderr, "A[%d] = %f\n", h_ix, sqrtf(A_model.realp[h_ix]*A_model.realp[h_ix] + A_model.imagp[h_ix]*A_model.imagp[h_ix]));
        }
        else
        {
            if(++miss_cnt > 3)
            {
                break;
            }
        }
        
        ratio /= 1.2;
    }
    
    fprintf(stderr,"f0_hat = %f\n\n", f0);
    
    return 0.0;
}

#else

#endif

void HarmonizerDSPKernel::addnote(int note, int vel)
{
    int min_dist = 129;
    int dist;
    int min_ix = -1;
    midi_changed_sample_num = sample_count;
    midi_changed = 1;
    
    int vel_ = midi_ignore_velocity ? 100 : vel;
    
    keys_down[note] = 1;
    
//        if (!midi_enable)
//        {
//            return;
//        }
    
    if (!midi_legato)
    {
        // look for one with the same note and take that if we can, or the empty one with
        // the closest last note to the one we want
        for (int k = n_auto; k < nvoices; k++)
        {
            dist = (voices[k].midinote == note) ? 0 : abs(voices[k].lastnote - note);
            if ((dist < min_dist && voices[k].midinote < 0))
            {
                min_ix = k;
                min_dist = dist;
            }
        }
        //fprintf(stderr, "%d,%d,%d\n", min_dist, min_ix,n_auto);

        if (min_dist >= 0)
        {
            voices[min_ix].lastnote = voices[min_ix].midinote;
            voices[min_ix].midinote = note;
            voices[min_ix].midivel = vel_;
            voices[min_ix].sample_num = sample_count;
            //voices[min_ix].ratio = 1.0;
            //voices[min_ix].nextgrain = 0;
            return;
        }
    }
    
    // otherwise, we are stealing
    min_dist = 129;
    min_ix = -1;
    for (int k = n_auto; k < nvoices; k++)
    {
        dist = abs(voices[k].midinote - note);
        if (dist < min_dist)
        {
            min_dist = dist;
            min_ix = k;
        }
    }
    
    voices[min_ix].lastnote = voices[min_ix].midinote;
    voices[min_ix].midinote = note;
    voices[min_ix].midivel = vel_;
    voices[min_ix].sample_num = sample_count;
    //voices[min_ix].nextgrain = 0;
    
    if (++voice_ix > nvoices)
        voice_ix = n_auto;
    
    return;
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
        if (voices[k].midinote == note)
        {
            voices[k].lastnote = note;
            voices[k].midinote = -1;
        }
    }
    
    midi_changed_sample_num = sample_count;
    midi_changed = 1;
}

void HarmonizerDSPKernel::pedal_down()
{
    midi_pedal = 1;
}
    
void HarmonizerDSPKernel::pedal_up()
{
    midi_pedal = 0;
    for (int k = n_auto; k < nvoices; k++)
    {
        if (!keys_down[voices[k].midinote])
        {
            voices[k].lastnote = voices[k].midinote;
            voices[k].midinote = -1;
        }
    }
}

void HarmonizerDSPKernel::analyze_harmony(void)
{
    //float intervals[127];
    int octave[12];
    
    memset(octave, 0, 12 * sizeof(int));
    
    int n = 0;
    for (int j = 3; j < nvoices; j++)
    {
        if (voices[j].midinote < 0)
            continue;
        
        midinotes[n++] = (float) voices[j].midinote;
        
        octave[voices[j].midinote % 12] = 1;
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

void HarmonizerDSPKernel::update_voices (void)
{
    static int last_nn = 0;
    
    voices[0].error = 0;
    //voices[0].ratio = 1;
    voices[0].target_ratio = 1;
    voices[0].formant_ratio = 1.0;
    voices[0].midivel = 127;
    voices[0].midinote = 0;
    
    if (n_auto > 1)
    {
        voices[1].midinote = -1;
        voices[1].midivel = 65;
        voices[1].pan = 0.5;
    }
    if (n_auto > 2)
    {
        voices[2].midinote = -1;
        voices[2].midivel = 65;
        voices[2].pan = -0.5;
    }
    if (n_auto > 3)
    {
        voices[3].midinote = -1;
        voices[3].midivel = 65;
        voices[3].pan = -0.5;
    }
    
    if (midi_changed && (sample_count - midi_changed_sample_num) > (int) sampleRate / 50)
    {
        midi_changed = 0;
        //printf("%d samples since last midi note\n", sample_count - midi_changed_sample_num);
        if (midi_link) { analyze_harmony(); }
    }
    
    static int was_voiced = 0;
    
    if (!voiced)
    {
        if (was_voiced)
        {
            //send_note_off(midi_note_number, 100);
            for (int k = 0; k < n_auto; k++)
            {
                send_note_off(voice_notes[k], 100);
            }
        }
        
        note_number = -1.0;
        midi_note_number = -1.0;
        last_nn = -1;
        
        for (int k = 0; k < n_auto; k++)
        {
            voice_notes[k] = -1;
        }
        voices[0].ratio = 1;
        was_voiced = 0;
        
        return;
    }
    
    float f = log2f (sampleRate / (T * baseTuning));
    
    float note_f = f * 12.0;
    int nn = (int) round(note_f);
    
    if (nn != last_nn && fabs(note_f - (float) last_nn) < .65)
    {
        nn = last_nn;
    }
    
    float err = note_f - (float) nn;
    
    //int old_midi_note_number = midi_note_number;
    midi_note_number = nn + 69;
    
    //fprintf(stderr, "%f,%d\n",note_f,last_nn);
    
    if (!was_voiced)
    {
        //send_note_on(midi_note_number, 100);
    }
    else
    {
        // note didn't change
        if (fabs(note_f - (float) last_nn) < 0.6)
        {
            midi_note_number = last_nn + 69;
        }
        else // note changed
        {
            last_nn = nn;
//                send_note_on(midi_note_number, 100);
//                send_note_off(old_midi_note_number, 100);
        }
    }
    
    int root = root_key % 12;
    int quality = root_key / 12;
    
    int interval = (last_nn + 69 - root) % 12;
    note_number = (nn + 69) % 12 + (note_f - nn);
    
    //last_nn = nn;
    
    // for autotune, find nearest 0-interval for voice 0
    
    voice_notes_old[0] = voice_notes[0];
    int inc = err > 0 ? 1 : -1;
    int ix = interval;
    int ix2 = interval;
    for (int k = 0; k < 12; k++)
    {
        if (interval_offsets[((ix+12)%12)*4 + quality*48] == 0)
        {
            voice_notes[0] = midi_note_number + ix - interval;
            voices[0].midinote = voice_notes[0];
            break;
        }
        
        if (interval_offsets[((ix2+12)%12)*4 + quality*48] == 0)
        {
            voice_notes[0] = midi_note_number + ix2 - interval;
            voices[0].midinote = voice_notes[0];
            break;
        }
        
        ix += inc; ix2 -= inc;
    }
    
    // convert interval table to midi notes for auto voices.
    for (int k = 1; k < n_auto; k++)
    {
        voice_notes_old[k] = voice_notes[k];
        voice_notes[k] = midi_note_number + interval_offsets[k + (interval*4) + (quality*48)];
        
        if (k > inversion){
            voice_notes[k] -= 12;
        }
        
        if (voice_notes[k] != voice_notes_old[k])
        {
            send_note_on(voice_notes[k], 100);
            if (voice_notes_old[k] >= 0)
                send_note_off(voice_notes_old[k], 100);
        }
        
        voices[k].midinote = voice_notes[k];
        voices[k].midivel = 65;
    }
    
    if (!auto_enable)
    {
        for (int k = 0; k < n_auto; k++)
        {
            simpleVoices[k].ratio = voices[k].ratio = voices[k].target_ratio = 1.0;
        }
    }
    
    // compute target resampling ratios for all voices.
    for (int k = 0; k < nvoices; k++)
    {
        //fprintf(stderr, "%d\t", voices[k].midinote);
        if (voices[k].midinote < 0)
            continue;
        
        // clear out lingering notes
        if (k >= n_auto && voices[k].midinote >= 0 && !keys_down[voices[k].midinote] && !midi_pedal)
        {
            voices[k].midinote = -1;
        }
        
        if (triad >= 0 && k < n_auto)
        {
            voices[k].target_ratio = major_chord_table[0][k];
            if (autotune)
                voices[k].target_ratio *= powf(2.0, -err/12.0);
            
            if (k > inversion)
            {
                voices[k].target_ratio /= 2;
            }
            
            if (k == 0 && voices[k].target_ratio != 1)
            {
                voices[0].midinote = 0; // hack to turn on voice
            }
            
            voice_notes[k] = midi_note_number + interval_offsets[k];
            
            if (k > inversion)
            {
                voice_notes[k] -= 12;
            }
            
            //fprintf(stderr,"%d: %f\n", k, major_chord_table[0][k]);
            continue;
        }
        
        if ((voiced && !was_voiced) || voices[k].lastnote < 0)
            voices[k].midinote_ = voices[k].midinote;
        else
        {
            float diff = 1 * (voices[k].midinote - voices[k].midinote_);
            if (fabs(diff) > speed)
                diff = sgn(diff) * speed;
            
            voices[k].midinote_ += diff;
        }
        
        float error_hsteps = (voices[k].midinote_ - 69) - note_f;
        
        if (k == 0)
        {
            error_hsteps *= corr_strength;
        }
        
        voices[k].target_ratio = powf(2.0, error_hsteps/12);
    }
    
    //fprintf(stderr, "\n");
    
    was_voiced = voiced;

    //float frac = 0.99 - (speed * 0.29);
    float v1frac = 0.1;
    
    for (int k = 0; k < nvoices; k++)
    {
//            if (k < n_auto && k >= start)
//                voices[k].target_ratio /= error_ratio; // for autoharm
        
        if (k == 0)
        {
            simpleVoices[k].ratio = voices[k].ratio = v1frac * voices[k].ratio + (1-v1frac) * voices[k].target_ratio;
            psolaVoices[k].T = T / simpleVoices[k].ratio;
        }
        else
        {
            simpleVoices[k].ratio = voices[k].ratio = v1frac * voices[k].ratio + (1-v1frac) * voices[k].target_ratio;
            psolaVoices[k].T = T / simpleVoices[k].ratio;
        }
    }
}


float HarmonizerDSPKernel::loopPosition() {
    return looper.position();
}

int HarmonizerDSPKernel::setLoopMode(int mode) {
    looper.setMode(static_cast<Looper::loopMode>(mode));
    return mode;
}
int HarmonizerDSPKernel::getLoopMode() {
    return looper.loop_mode;
}

inline float HarmonizerDSPKernel::window_value(float f)
{
    if (f <= 0 || f >= 1)
    {
        return 0.0;
    }
    
    float wi = 3 + graintablesize * f;
    int i = (int) wi;
    
    float w = cubic (grain_window + i, wi - i);
//        fprintf(stderr, "index: %d: %f\n", i, w);
//        fprintf(stderr, "done\n");
    
    return w;
}



// :MARK utility Functions
inline float cubic (float *v, float a)
{
    float b, c;
    
    b = 1 - a;
    c = a * b;
    return (1.0f + 1.5f * c) * (v[1] * b + v[2] * a)
    - 0.5f * c * (v[0] * b + v[1] + v[2] + v[3] * a);
}

float quadratic_peak(float *pv, float *v, int ix)
{
    float pix = 0.5*(v[ix-1] - v[ix+1])/(v[ix+1] + v[ix-1] - 2*v[ix]);
    *pv = v[ix] - .25*(v[ix-1]-v[ix+1])*pix;
    return pix;
}

float linear (float *v, float a)
{
    return v[0] * (1 - a) + v[1] * a;
}

float linear_interp(float *v, float ix)
{
    int iix = (int) floorf(ix);
    return linear(v + iix, ix - iix);
}

float cubic_interp(float *v, float ix)
{
    int iix = (int) floorf(ix);
    return cubic(v + iix - 1, ix - iix);
}
