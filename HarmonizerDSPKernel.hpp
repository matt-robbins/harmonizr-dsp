/*
	<samplecode>
		<abstract>
			A DSPKernel subclass implementing the realtime signal processing algorithm.
		</abstract>
	</samplecode>
*/

#ifndef HarmonizerDSPKernel_hpp
#define HarmonizerDSPKernel_hpp


#import <vector>
#import <cmath>
#import <sys/time.h>

#ifdef __APPLE__
#import "../harmonizr-dsp/DSPKernel.hpp"
#import "ParameterRamper.hpp"
#include <Accelerate/Accelerate.h>
#include <dispatch/dispatch.h>

#include "TestAudioData.h"

typedef AUAudioFrameCount frame_count_t;
typedef AUParameterAddress param_address_t;
typedef AUValue param_value_t;
typedef AudioTimeStamp timestamp_t;

typedef AUMIDIEvent midi_event_t;

#else
#include "kiss_fft.h"
#import <algorithm>
#include <android/log.h>
typedef int32_t frame_count_t;
typedef int32_t param_address_t;
typedef float param_value_t;
typedef float timestamp_t;

template <typename T>
T clamp(T input, T low, T high) {
    return std::min(std::max(input, low), high);
}

#endif

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

typedef struct grain_s
{
    float size;
    float start;
    float ix;
    float ratio;
    float gain;
    float pan;
    int vix;
    int synth_ix;
} grain_t;

typedef struct voice_s
{
    float error;
    float ratio;
    float target_ratio;
    float formant_ratio;
    float nextgrain;
    float pan;
    float ix1;
    float ix2;
    float gain;
    float target_gain;
    float vibrato_rate;
    float vibrato_amp;
    float vib_phase;
    int xfade_ix;
    int xfade_dur;
    int midinote;
    float midinote_;
    int midivel;
    int lastnote;
    unsigned int sample_num;
} voice_t;

typedef struct triad_ratio_s
{
    float r1;
    float r2;
} triad_ratio_t;

typedef float chord_ratio_t[4];

typedef enum triads
{
    TRIAD_MAJOR_R = 0,
    TRIAD_MAJOR_1,
    TRIAD_MAJOR_2,
    TRIAD_MINOR_R,
    TRIAD_MINOR_1,
    TRIAD_MINOR_2,
    TRIAD_DIM_R,
    TRIAD_DIM_1,
    TRIAD_DIM_2,
    TRIAD_SUS4_R,
    TRIAD_SUS4_1,
    TRIAD_SUS4_2,
    TRIAD_SUS2_R,
    TRIAD_SUS2_1,
    TRIAD_SUS2_2,
    TRIAD_AUG,
    TRIAD_5TH,
    TRIAD_OCT
} triad_t;

enum {
    HarmParamKeycenter = 0,
    HarmParamInversion,
    HarmParamNvoices,
    HarmParamAuto,
    HarmParamAutoStrength,
    HarmParamMidi,
    HarmParamMidiLink,
    HarmParamMidiLegato,
    HarmParamMidiKeyCC,
    HarmParamMidiKeyCcOffset,
    HarmParamMidiQualCC,
    HarmParamMidiQualCcOffset,
    HarmParamMidiNvoiceCC,
    HarmParamMidiNvoiceCcRange,
    HarmParamMidiInvCC,
    HarmParamMidiInvCcRange,
    HarmParamMidiPC,
    HarmParamMidiHarmOut,
    HarmParamMidiMelOut,
    HarmParamTriad,
    HarmParamBypass,
    HarmParamDouble,
    HarmParamHgain,
    HarmParamVgain,
    HarmParamDryMix,
    HarmParamSpeed,
    HarmParamTuning,
    HarmParamThreshold,
    HarmParamStereo,
    HarmParamSynth,
    HarmParamInterval
};

enum {
    HarmPresetChords=0,
    HarmPresetDiatonic,
    HarmPresetChromatic,
    HarmPresetBarbershop,
    HarmPresetMIDI,
    HarmPresetBohemian,
    HarmPresetBass,
    HarmPreset4ths,
    HarmPresetModes
};

enum {
    StereoModeNormal=0,
    StereoModeMono,
    StereoModeSplit
};

#define N_AUTO 4

static inline double squared(double x) {
    return x * x;
}

static inline float inc_to_target(float value, float target, float c, float maxrate_up, float maxrate_down)
{
    float diff = c * (target - value);
    int sign = sgn(diff);
    if (sign > 0 && diff > maxrate_up)
        diff = maxrate_up;
    
    else if (sign < 0 && diff < maxrate_down)
        diff = maxrate_down;
    
    return value + diff;
}

/*
	HarmonizerDSPKernel
	Performs our filter signal processing.
	As a non-ObjC class, this is safe to use from render thread.
*/

#ifdef __APPLE__
class HarmonizerDSPKernel : public DSPKernel {
#else
class HarmonizerDSPKernel {
#endif
public:
    int n_channels = 0;
    // MARK: Member Functions

    HarmonizerDSPKernel() {}
	
	void init(int channelCount, double inSampleRate) {
		n_channels = channelCount;
        fprintf(stderr,"**** init with %d channels! at %f Hz\n", n_channels, inSampleRate);
		
		sampleRate = float(inSampleRate);
#ifdef __APPLE__
        fft_s = vDSP_create_fftsetup(11, 2);
        
        fft_in.realp = (float *) calloc(2048, sizeof(float));
        fft_in.imagp = (float *) calloc(2048, sizeof(float));
        
        fft_out.realp = (float *) calloc(2048, sizeof(float));
        fft_out.imagp = (float *) calloc(2048, sizeof(float));
        
        fft_out2.realp = (float *) calloc(2048, sizeof(float));
        fft_out2.imagp = (float *) calloc(2048, sizeof(float));
        
        fft_buf.realp = (float *) calloc(2048, sizeof(float));
        fft_buf.imagp = (float *) calloc(2048, sizeof(float));
        
        A_model.realp = (float *) calloc(100, sizeof(float));
        A_model.imagp = (float *) calloc(100, sizeof(float));
        
        Hann.realp = (float *) calloc(2048, sizeof(float));
        Hann.imagp = (float *) calloc(2048, sizeof(float));
        
        fft_in.realp[0] = 1.0;
        
        fprintf(stderr, "imagp = %p\n", fft_in.imagp);
        

#else
        fft_s = kiss_fft_alloc(2048,0,NULL,0);
        ifft_s = kiss_fft_alloc(2048,1,NULL,0);
        fft_in = (kiss_fft_cpx *) calloc(2048, sizeof(kiss_fft_cpx));
        fft_out = (kiss_fft_cpx *) calloc(2048, sizeof(kiss_fft_cpx));
        fft_out2 = (kiss_fft_cpx *) calloc(2048, sizeof(kiss_fft_cpx));
        Hann = (kiss_fft_cpx *) calloc(2048, sizeof(kiss_fft_cpx));
#endif
        
        ncbuf = 4096;
        cbuf = (float *) calloc(ncbuf + 3, sizeof(float));
        
        spec_env = (float *) calloc(2048, sizeof(float));
        synth_pulse = (float **) calloc(n_synth_pulse, sizeof(float *));
        
        for (int k = 0; k < n_synth_pulse; k++)
            synth_pulse[k] = (float *) calloc(2*maxT, sizeof(float));
        
        lms_h = (float *) calloc(lms_n, sizeof(float));
        
        nvoices = 16;
        voices = (voice_t *) calloc(nvoices, sizeof(voice_t));
        voice_ix = 1;
        
        in_buffers = (float **) calloc(channelCount, sizeof(float *));
        out_buffers = (float **) calloc(channelCount, sizeof(float *));
        
        for (int k = 0; k < nvoices; k++)
        {
            voices[k].midinote = -1;
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
            
            if (k >= 3)
            {
                voices[k].pan = ((float)(k - 3) / (float)(nvoices - 3)) - 0.5;
                //voices[k].formant_ratio = ((float)(k - 3) / (float)(nvoices - 3)) + 0.5;
            }
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
        
        memset(Tbuf, 0, 5*sizeof(float));
        Tix = 0;
        
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
        
        // define triad ratios
        
        triads[TRIAD_MAJOR_R].r1 = intervals[4];
        triads[TRIAD_MAJOR_R].r2 = intervals[7];
        triads[TRIAD_MAJOR_1].r1 = intervals[5];
        triads[TRIAD_MAJOR_1].r2 = intervals[9];
        triads[TRIAD_MAJOR_2].r1 = intervals[3];
        triads[TRIAD_MAJOR_2].r2 = intervals[8];
        
        triads[TRIAD_MINOR_R].r1 = intervals[3];
        triads[TRIAD_MINOR_R].r2 = intervals[7];
        triads[TRIAD_MINOR_1].r1 = intervals[5];
        triads[TRIAD_MINOR_1].r2 = intervals[8];
        triads[TRIAD_MINOR_2].r1 = intervals[4];
        triads[TRIAD_MINOR_2].r2 = intervals[9];
        
        triads[TRIAD_DIM_R].r1 = intervals[3];
        triads[TRIAD_DIM_R].r2 = intervals[6];
        triads[TRIAD_DIM_1].r1 = intervals[6];
        triads[TRIAD_DIM_1].r2 = intervals[9];
        triads[TRIAD_DIM_2].r1 = intervals[3];
        triads[TRIAD_DIM_2].r2 = intervals[9];
        
        triads[TRIAD_SUS4_R].r1 = intervals[5];
        triads[TRIAD_SUS4_R].r2 = intervals[7];
        triads[TRIAD_SUS4_1].r1 = intervals[2];
        triads[TRIAD_SUS4_1].r2 = intervals[7];
        triads[TRIAD_SUS4_2].r1 = intervals[5];
        triads[TRIAD_SUS4_2].r2 = intervals[10];
        
        triads[TRIAD_SUS2_R].r1 = intervals[2];
        triads[TRIAD_SUS2_R].r2 = intervals[7];
        
        triads[TRIAD_SUS2_1].r1 = intervals[5];
        triads[TRIAD_SUS2_1].r2 = intervals[10];
        
        triads[TRIAD_SUS2_2].r1 = intervals[5];
        triads[TRIAD_SUS2_2].r2 = intervals[7];
        
        triads[TRIAD_AUG].r1 = intervals[4];
        triads[TRIAD_AUG].r2 = intervals[8];
        
        triads[TRIAD_5TH].r1 = 3.0/2.0;
        triads[TRIAD_5TH].r2 = 2.0;
        
        triads[TRIAD_OCT].r1 = 0.5;
        triads[TRIAD_OCT].r2 = 2.0;
        
        memset(midinotes, 0, 128 * sizeof(float));
        memset(keys_down, 0, 128 * sizeof(int));

        int chords_intervals[] = {0,4,7,12, -1,3,6,11, 2,5,10,14, 1,4,9,13, 0,3,8,12, -1,2,7,11, 1,6,10,13, 0,5,9,12, -1,4,8,11, 0,3,7,10, 2,6,9,14, 1,5,8,13, // major
                           0,3,7,12, -1,2,6,11, 1,5,10,13, 0,4,9,12, -1,3,8,11, -1,2,7,10, 1,6,9,13, 0,5,8,12, 0,4,7,11, 0,3,6,10, 0,5,9,14, 1,4,8,13, // minor
                           0,4,10,12, -1,3,9,11, -2,2,8,10, 1,4,7,9, 0,3,6,8, 2,5,7,11, 1,4,6,10, 0,3,5,9, -1,2,4,8, 1,3,7,10, 0,2,6,9, -1,1,5,8, //dom
        };
        for (int i = 0; i < 144; i++)
        {
            setParameter(HarmParamInterval+i,(float) chords_intervals[i]);
        }
        
	}
    
    void fini() {
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

        delete grain_window;
        delete fft_mag;
        delete fft_mag_db;
        delete grains;

        free(cbuf);
        free(voices);
        
        free(in_buffers);
        free(out_buffers);
    }
	
	void reset() {
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
        
        pitchmark[0] = 0;
        pitchmark[1] = -1;
        pitchmark[2] = -1;
	}
    
    float cubic (float *v, float a)
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
    
    inline float window_value(float f)
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
	
	void setParameter(param_address_t address, param_value_t value) {
        switch (address) {
            case HarmParamKeycenter:
                root_key = (int) clamp(value,0.f,47.f);
                break;
            case HarmParamInversion:
                inversion = (int) clamp(value,0.f,3.f);
                break;
            case HarmParamNvoices:
                n_auto = (int) clamp(value,1.f,4.f);
                //fprintf(stderr, "nvoices: %d\n", n_auto);
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

	param_value_t getParameter(param_address_t address) {
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

    void setBuffers(float ** in, float ** out) {

        for (int k = 0; k < n_channels; k++)
        {
            in_buffers[k] = in[k];
            out_buffers[k] = out[k];
        }
    }


    void setPreset(int preset_ix_)
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
    int getPreset()
    {
        return preset_ix;
    }

#ifdef __APPLE__

	void startRamp(AUParameterAddress address, AUValue value, AUAudioFrameCount duration) override {
        return;
	}

	void setBuffers(AudioBufferList* inBufferList, AudioBufferList* outBufferList) {

        for (int k = 0; k < n_channels; k++)
        {
            in_buffers[k] = (float*) inBufferList->mBuffers[k].mData;
            out_buffers[k] = (float*) outBufferList->mBuffers[k].mData;
        }

    }
    
    virtual void handleMIDIEvent(midi_event_t const& midiEvent) override {
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

#else

#endif
    
    void psola(float *out, float *out2)
    {
        int first_psola_voice = 0;
        
        voicegain_target = 0;
        if (!autotune && triad < 0)
        {
            first_psola_voice = 1;
            voicegain_target = dry_mix;
        }
        
        int nonvoiced_count = (int) (voicegain > 0);
        
        // Ramp Gain
        voicegain += .001 * sgn(voicegain_target - voicegain);
        harmgain += .001 * sgn(harmgain_target - harmgain);
        
        for (int vix = first_psola_voice; vix < nvoices; vix++)
        {
            voices[vix].target_gain = ((voiced || nonvoiced_count == 0) && voices[vix].midinote > 0) ? 1.0 : 0.0;
            
            if (voices[vix].target_gain > 0)
                nonvoiced_count++;
            
            if (vix == 0 && (autotune || triad >= 0)) voices[vix].target_gain = dry_mix;
            //voices[vix].gain += .001 * sgn(voices[vix].target_gain - voices[vix].gain);
            voices[vix].gain = inc_to_target(voices[vix].gain, voices[vix].target_gain, 0.9, 0.001, -0.0004);
            
            if (voices[vix].gain < 0.001)
            {
                continue;
            }
            
            if (vix >= n_auto && !midi_enable)
                continue;
            
            float unvoiced_offset = 0;
            //                if (!voiced && vix > 5)
            //                {
            //                    continue;
            //                    //unvoiced_offset = T * (0.5 - (float) rand() / RAND_MAX);
            //                }
            
            float midigain_local = 1.0;
            if (vix >= 3)
                midigain_local = midigain;
            
            if (--voices[vix].nextgrain < 0)
            {
                // search for the first open grain
                int found_grain = 0;
                for (int k = 0; k < ngrains; k++)
                {
                    if (grains[k].size < 0)
                    {
                        grains[k].size = 2 * T;
                        grains[k].start = pitchmark[0] - voices[vix].nextgrain - T + unvoiced_offset;
                        grains[k].ratio = voices[vix].formant_ratio;
                        grains[k].synth_ix = synth_pulse_ix;
                        //memcpy(grains[k].data,synth_pulse,roundf(2*T)*sizeof(float));
                        
                        grains[k].ix = 0;
                        grains[k].gain = midigain_local * (float) voices[vix].midivel / 127.0;
                        grains[k].pan = voices[vix].pan;
                        grains[k].vix = vix;
                        
                        if (vix == 0)
                        {
                            grains[k].gain = 1.0;
                        }
                        
                        if (!voiced)
                        {
                            grains[k].ratio = 1.0; //voices[vix].ratio;
                        }
                        else
                        {
                            // for low transpositions, increase gain
                            if (voices[vix].ratio < 1)
                                grains[k].gain *= powf(1/voices[vix].ratio,0.5);
                            
                            // for high transpositions, start shortening the blips.
                            //                                if (voices[vix].ratio > 1.7)
                            //                                {
                            //                                    //grains[k].size = T/2;
                            //                                    grains[k].ratio *= (1 + (voices[vix].ratio - 1.7)/2);
                            //                                    //grains[k].gain *= powf(voices[vix].ratio,0.5);
                            //                                }
                        }
                        
                        //voices[vix].nextgrain += voiced ? (T / voices[vix].ratio) : T;
                        
                        float v_per = (T / voices[vix].ratio);
                        
                        while (voices[vix].nextgrain < 0)
                            voices[vix].nextgrain += v_per;
                        
                        voices[vix].vib_phase += 2*M_PI * (v_per / (sampleRate/voices[vix].vibrato_rate));
                        if (voices[vix].vib_phase > 2*M_PI)
                            voices[vix].vib_phase -= 2*M_PI;
                        
                        float vib_delay = voices[vix].vibrato_amp * v_per * sinf(voices[vix].vib_phase);
                        voices[vix].nextgrain += vib_delay;
                        
                        //printf("phase = %f\n", voices[vix].vib_phase);
                        if (k > maxgrain)
                            maxgrain = k;
                        
                        found_grain = true;
                        break;
                    }
                }
                
                if (found_grain == false)
                {
                    fprintf(stderr, "couldn't find grain!\n");
                }
                for (int k = maxgrain; k > 0; k--)
                {
                    if (grains[k].size >= 0)
                        break;
                    
                    maxgrain = k;
                }
            }
        }
        
        for (int ix = 0; ix <= maxgrain; ix++)
        {
            grain_t g = grains[ix];
            
            // if this grain has been "triggered", it's size is > 0
            if (g.size > 0)
            {
                float fi = g.start + g.ix;
                
                if (fi >= ncbuf)
                    fi -= ncbuf;
                else if (fi < 0)
                    fi += ncbuf;
                
                int i = (int) fi;
                //float u = cbuf[i];
                float u = cubic (cbuf + i, fi - i);
                
                float wi = 2 + graintablesize * (g.ix / g.size);
                i = (int) wi;
                float w = cubic (grain_window + i, wi - i);
                
                //w = 1;
                
                // shrink left side of the window to create a sharper attack,
                // which should translate to a cleaner sound for high transpositions.
                float mix = 1.0;
                if (voices[g.vix].ratio > 1.8)
                {
                    mix = fmax(0.0, 1.0 - (voices[g.vix].ratio - 1.8));
                }
//                if (voices[g.vix].ratio > 2)
//                {
//                    mix = 0.0;
//                }
                
                if (synth_enable)
                {
                    u = (mix * u) + (1 - mix) * cubic(synth_pulse[g.synth_ix] + (int)(g.ix)+3, g.ix - floorf(g.ix));
                    //u = synth_pulse[(int) g.ix];
                }
                
                
                float f = g.ix/g.size;
                
                if (f < 0.5 && voices[g.vix].ratio > 1)
                {
                    f = 0.5 - (0.5 - f)*voices[g.vix].ratio;
                }
                
                w = window_value(f);
                
                float hgain = g.vix ? harmgain : 1.0;
                
                switch (stereo_mode)
                {
                    case StereoModeNormal:
                        *out += u * hgain * w * g.gain * voices[g.vix].gain * (g.pan + 1.0)/2;
                        *out2 += u * hgain * w * g.gain * voices[g.vix].gain * (-g.pan + 1)/2;
                        break;
                    case StereoModeMono:
                        *out += u * hgain * w * g.gain * voices[g.vix].gain;
                        *out2 = *out;
                        break;
                    case StereoModeSplit:
                        *out2 += u * hgain * w * g.gain * voices[g.vix].gain;
                        break;
                }
                
                g.ix += g.ratio;
                
                if (g.ix > g.size)
                {
                    g.size = -1;
                    //printf("ending grain %d\n", ix);
                }
                
                grains[ix] = g;
            }
        }
    }

#ifdef __APPLE__
	void process(frame_count_t frameCount, frame_count_t bufferOffset) override {
#else
    void process(frame_count_t frameCount, frame_count_t bufferOffset) {
#endif

		int channelCount = n_channels;
        sample_count += frameCount;
        n_output_events = 0;
        
        // 1st order filters
        static float inbuf = 0;
        static float outbuf1 = 0;
        static float outbuf2 = 0;
        
        if (bufferOffset != 0)
            fprintf(stderr, "buffer_offset = %d\n", bufferOffset);
        
        // For each sample.
		for (int frameIndex = 0; frameIndex < frameCount; ++frameIndex)
        {
            int frameOffset = int(frameIndex + bufferOffset);
						
            float* in  = in_buffers[0] + frameOffset;
            float* out = out_buffers[0] + frameOffset;
            float* out2 = out;
            
            if (channelCount > 1)
            {
                out2 = out_buffers[1] + frameOffset;
            }
            
            // preemphasis
            cbuf[cix] = (*in - inbuf * 0.9);
            inbuf = *in;
            *in = cbuf[cix];
            
            // TODO: lms equalizer
            
            if (bypass)
            {
                *out = *in / 2;
                *out2 = *out;
                continue;
            }
            else
            {
                *out = *out2 = 0;
            }
            
            if (cix < 3)
            {
                cbuf[ncbuf+cix] = cbuf[cix];
            }
            if (++cix >= ncbuf)
                cix = 0;
            
            if (--rcnt == 0)
            {
                rcnt = 256;
                int oldT = T;
                float p = estimate_pitch(cix - 2*maxT);
                if (p > 0)
                    T = p;
                else
                    T = oldT;
                
                voiced = (p != 0);
                
                if (voiced && synth_enable)
                {
                    //get_model(cix - 3*T);
                    get_minphase_pulse(cix - nfft);
                }
                
                update_voices();
            }
            
            float dp = cix - 2*maxT - pitchmark[0];
            if (dp < 0)
                dp += ncbuf;
            
            if (dp > (T + T/4))
            {
                findmark();
            
                //printf("pitchmark[0,1,2] = %.2f,%.2f,%.2f\ninput = %d\n", pitchmark[0],pitchmark[1],pitchmark[2],cix);
            }
            
            *out = *in * voicegain / 2;
            if (stereo_mode != StereoModeSplit)
            {
                *out2 = *out;
            }
            
            psola(out, out2);
            
            float a = .9;

            *out = (*out + outbuf1 * a);
            *out2 = (*out2 + outbuf2 * a);

            outbuf1 = *out;
            outbuf2 = *out2;
		}
	}

#ifdef __APPLE__
    float estimate_pitch(int start_ix)
    {
        if (fft_in.imagp <= (float *) 0)
        {
            return 0.0;
        }
        //fprintf(stderr, "estimate pitch : imagp = %p\n", fft_in.imagp);
        memset(fft_in.realp, 0, nfft * sizeof(float));
        memset(fft_in.imagp, 0, nfft * sizeof(float));
        
        for (int k = 0; k < maxT; k++)
        {
            int ix = (start_ix + k) & cmask;
            fft_in.realp[k] = cbuf[ix];
        }
        
        vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out, 1, &fft_buf, 11, 1);
        
        for (int k = maxT; k < 2*maxT; k++)
        {
            int ix = (start_ix + k) & cmask;
            fft_in.realp[k] = cbuf[ix];
        }
        vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out2, 1, &fft_buf, 11, 1);
        
        // conjugate small window and correlate with large window
        for (int k = 0; k < nfft; k++)
        {
            float r1,c1,r2,c2;
            r1 = fft_out.realp[k]; c1 = -fft_out.imagp[k];
            r2 = fft_out2.realp[k]; c2 = fft_out2.imagp[k];

            float factor = (float) k / (float) nfft;
            if (k < nfft)
            {
                fft_in.realp[k] = (r1*r2 - c1*c2);// * (1 - factor);
                fft_in.imagp[k] = (r1*c2 + r2*c1);// * (1 - factor);
            }
            else
            {
                fft_in.realp[k] = 0;
                fft_in.imagp[k] = 0;
            }
        }
        // inverse transform
        vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out, 1, &fft_buf, 11, -1);
        
        float sumsq_ = fft_out.realp[0]/nfft;
        float sumsq = sumsq_;
        
        float df,cmdf,cmdf1,cmdf2, sum = 0;
        
        float period = 0.0;
        
        cmdf2 = cmdf1 = cmdf = 1;
        for (int k = 1; k < maxT; k++)
        {
            int ix1 = (start_ix + k) & cmask;
            int ix2 = (start_ix + k + maxT) & cmask;
            
            sumsq -= cbuf[ix1]*cbuf[ix1];
            sumsq += cbuf[ix2]*cbuf[ix2];
            
            df = sumsq + sumsq_ - 2 * fft_out.realp[k]/nfft;
            sum += df;
            cmdf2 = cmdf1; cmdf1 = cmdf;
            cmdf = (df * k) / sum;
            
            if (k > 0 && cmdf2 > cmdf1 && cmdf1 < cmdf && cmdf1 < threshold && k > 20)
            {
                period = (float) (k-1) + 0.5*(cmdf2 - cmdf)/(cmdf2 + cmdf - 2*cmdf1); break;
            }
        }
        
//        if (sumsq < 0.01)
//        {
//            period = 0;
//        }
        
        Tbuf[Tix++] = period;
        
        //fprintf(stderr, "%f\n", sumsq);
        
        if (Tix >= nmed)
            Tix = 0;
        
        memcpy(Tsrt, Tbuf, nmed * sizeof(float));
        vDSP_vsort(Tsrt, (vDSP_Length) nmed, 1);
        
        return Tsrt[nmed/2];
    }
        
    float get_minphase_pulse(int start_ix)
    {
        memset(fft_in.realp, 0, nfft * sizeof(float));
        memset(fft_in.imagp, 0, nfft * sizeof(float));
        for (int k = 0; k < nfft; k++)
        {
            int ix = (start_ix + k) & cmask;
            fft_in.realp[k] = cbuf[ix] * window_value((float)k/(nfft));
        }
        
        vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out, 1, &fft_buf, 11, 1);
        // log(abs(fft))
        for (int k = 0; k < nfft; k++)
        {
            float r1,c1;
            r1 = fft_out.realp[k]; c1 = fft_out.imagp[k];
            
            fft_in.realp[k] = logf((r1*r1 + c1*c1)+0.00001)/2; // factor of 2 is square root
            fft_in.imagp[k] = 0;
        }

        // cepstrum
        vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out, 1, &fft_buf, 11, -1);
        
        // fold and window
        fft_in.realp[0] = fft_out.realp[0]/nfft;

        for (int k = 1; k < nfft; k++)
        {
            if (k >= roundf(T*0.4)-1)
                fft_in.realp[k] = 0.f;
            else
            {
                fft_in.realp[k] = (fft_out.realp[k]*2)/nfft;
            }
            
            fft_in.imagp[k] = 0;
        }
        
        // smooth envelope (fft)
        vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out, 1, &fft_buf, 11, 1);
        
        for (int k = 0; k < nfft; k++)
        {
            spec_env[k] = 0.9*spec_env[k] + .1*fft_out.realp[k];

            float ex = expf(spec_env[k]);
            
            fft_in.realp[k] = ex * Hann.realp[k]/2;
            fft_in.imagp[k] = ex * Hann.imagp[k]/2;
        }
        
        // get impulse response (ifft)
        vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out, 1, &fft_buf, 11, -1); //inverse
        int ix = nfft-round(T);
        
        if (++synth_pulse_ix >= n_synth_pulse)
            synth_pulse_ix = 0;
        
        for (int k = 0; k < 2*maxT; k++)
        {
            if (ix >= nfft)
                ix = 0;
            synth_pulse[synth_pulse_ix][k] = fft_out.realp[ix]/nfft;
            ix++;
            //fprintf(stderr,"%f\n",synth_pulse[k]);
        }
        
        return 0.0;
    }
        
    float get_model(int start_ix)
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
    float estimate_pitch(int start_ix) {
        memset(fft_in, 0, nfft * sizeof(kiss_fft_cpx));

        for (int k = 0; k < maxT; k++) {
            int ix = (start_ix + k) & cmask;
            fft_in[k].r = cbuf[ix];
        }

        kiss_fft(fft_s, fft_in, fft_out);

        //memset(fft_in, 0, nfft * sizeof(kiss_fft_cpx));

        for (int k = maxT; k < 2 * maxT; k++) {
            int ix = (start_ix + k) & cmask;
            fft_in[k].r = cbuf[ix];
        }

        kiss_fft(fft_s, fft_in, fft_out2);

        // conjugate small window and correlate with large window
        for (int k = 0; k < nfft; k++) {
            float r1, c1, r2, c2;
            r1 = fft_out[k].r;
            c1 = -fft_out[k].i;
            r2 = fft_out2[k].r;
            c2 = fft_out2[k].i;

            fft_in[k].r = r1 * r2 - c1 * c2;
            fft_in[k].i = r1 * c2 + r2 * c1;
        }
        // inverse transform
        kiss_fft(ifft_s, fft_in, fft_out);

        float sumsq_ = fft_out[0].r / nfft;
        float sumsq = sumsq_;

        float df, cmdf, cmdf1, cmdf2, sum = 0;

        float period = 0.0;

        cmdf2 = cmdf1 = cmdf = 1;
        for (int k = 1; k < maxT; k++) {
            int ix1 = (start_ix + k) & cmask;
            int ix2 = (start_ix + k + maxT) & cmask;

            sumsq -= cbuf[ix1] * cbuf[ix1];
            sumsq += cbuf[ix2] * cbuf[ix2];

            df = sumsq + sumsq_ - 2 * fft_out[k].r / nfft;
            sum += df;
            cmdf2 = cmdf1;
            cmdf1 = cmdf;
            cmdf = (df * k) / sum;

            if (k > 0 && cmdf2 > cmdf1 && cmdf1 < cmdf && cmdf1 < threshold && k > 20) {
                period = (float) (k - 1) + 0.5 * (cmdf2 - cmdf) / (cmdf2 + cmdf - 2 * cmdf1);
                break;
            }
        }

        Tbuf[Tix++] = period;

        if (Tix >= nmed)
            Tix = 0;

        memcpy(Tsrt, Tbuf, nmed * sizeof(float));

        std::sort(Tsrt, Tsrt+nmed);

        return Tsrt[nmed / 2];
    }

#endif
    
    void findmark (void)
    {
        int mask = ncbuf - 1;
        
        memmove(pitchmark + 1, pitchmark, 2 * sizeof(float));
        pitchmark[0] = pitchmark[1] + T;
        
        if (pitchmark[0] >= ncbuf - 1)
            pitchmark[0] -= ncbuf;
        
        if (!voiced)
            return;
        
        // find center of mass around next pitch mark
        float mean = 0.0;
        float min = HUGE_VALF;
        
        int srch_n = (int)T/4;
        
        for (int k = -srch_n; k < srch_n; k++)
        {
            int ix = ((int) pitchmark[0] + k) & mask;
            if (cbuf[ix] < min)
                min = cbuf[ix];
        }
        
        mean = 0;
        float sum = 0;
        for (int k = -srch_n; k < srch_n; k++)
        {
            int ix = ((int) pitchmark[0] + k) & mask;
            mean += (float) k * (cbuf[ix] - min);
            sum += (cbuf[ix] - min);
        }
        
        if (sum == 0)
            pitchmark[0] += T;
        else
        {
            pitchmark[0] += (mean/sum);
        }
        
        if (pitchmark[0] < 0)
            pitchmark[0] += ncbuf;
        else if (pitchmark[0] > ncbuf - 1)
            pitchmark[0] -= ncbuf;
    }
    
    void addnote(int note, int vel)
    {
        int min_dist = 129;
        int dist;
        int min_ix = -1;
        midi_changed_sample_num = sample_count;
        midi_changed = 1;
        
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
                dist = abs(voices[k].lastnote - note);
                if ((dist < min_dist && voices[k].midinote < 0) || voices[k].midinote == note)
                {
                    min_ix = k;
                    min_dist = dist;
                }
            }
            
            if (min_dist >= 0)
            {
                voices[min_ix].lastnote = voices[min_ix].midinote;
                voices[min_ix].midinote = note;
                voices[min_ix].midivel = vel;
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
        voices[min_ix].midivel = vel;
        voices[min_ix].sample_num = sample_count;
        //voices[min_ix].nextgrain = 0;
        
        if (++voice_ix > nvoices)
            voice_ix = n_auto;
        
        return;
    }
    
    void remnote(int note)
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
        
    void pedal_down()
    {
        midi_pedal = 1;
    }
        
    void pedal_up()
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
    
    void analyze_harmony(void)
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
        
    void send_note_on(int nn, int vel)
    {
        // queue MIDI note on messages
        output_events[n_output_events].length = 3;
        output_events[n_output_events].data[0] = 0x90;
        output_events[n_output_events].data[1] = nn;
        output_events[n_output_events].data[2] = vel;
        n_output_events++;
    }
        
    void send_note_off(int nn, int vel)
    {
        // queue MIDI note off
        output_events[n_output_events].length = 3;
        output_events[n_output_events].data[0] = 0x80;
        output_events[n_output_events].data[1] = nn;
        output_events[n_output_events].data[2] = vel;
        n_output_events++;
    }
    
    void update_voices (void)
    {
        static int last_nn = 0;
        
        voices[0].error = 0;
        //voices[0].ratio = 1;
        voices[0].target_ratio = 1;
        voices[0].formant_ratio = 1.0;
        voices[0].midivel = 127;
        voices[0].midinote = 0;
        
        voices[1].midinote = -1;
        voices[1].midivel = 65;
        voices[1].pan = 0.5;
        
        voices[2].midinote = -1;
        voices[2].midivel = 65;
        voices[2].pan = -0.5;
        
        voices[3].midinote = -1;
        voices[3].midivel = 65;
        voices[3].pan = -0.5;
        
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
            
            was_voiced = 0;
            return;
        }
        
        float f = log2f (sampleRate / (T * baseTuning));
        
        float note_f = f * 12.0;
        int nn = (int) round(note_f);
        
        int old_midi_note_number = midi_note_number;
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
        
        // convert interval table to midi notes for auto voices.
        for (int k = 0; k < n_auto; k++)
        {
            voice_notes_old[k] = voice_notes[k];
            voice_notes[k] = midi_note_number + interval_offsets[k + (interval*4) + (quality*48)];
            
            if (k > inversion)
            {
                voice_notes[k] -= 12;
            }
            
            if (voice_notes[k] != voice_notes_old[k])
            {
                send_note_on(voice_notes[k], 100);
                if (voice_notes_old[k] >= 0)
                    send_note_off(voice_notes_old[k], 100);
            }
            
            voices[k].midinote = voice_notes[k];
        }
        
        if (!auto_enable)
        {
            for (int k = 0; k < n_auto; k++)
            {
                voices[k].ratio = voices[k].target_ratio = 1.0;
            }
        }
        
        // compute target resampling ratios for all voices.
        for (int k = 0; k < nvoices; k++)
        {
            if (voices[k].midinote < 0)
                continue;
            
            if (triad >= 0 && k < n_auto)
            {
                voices[k].target_ratio = major_chord_table[0][k];
                
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
                float diff = .9 * (voices[k].midinote - voices[k].midinote_);
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
        
        was_voiced = voiced;
    
        //float frac = 0.99 - (speed * 0.29);
        float v1frac = 0.5;
        
        for (int k = 0; k < nvoices; k++)
        {
//            if (k < n_auto && k >= start)
//                voices[k].target_ratio /= error_ratio; // for autoharm
            
            if (k == 0)
            {
                voices[k].ratio = v1frac * voices[k].ratio + (1-v1frac) * voices[k].target_ratio;
            }
            else
            {
                voices[k].ratio = v1frac * voices[k].ratio + (1-v1frac) * voices[k].target_ratio;
            }
        }
    }
	
    // MARK: Member Variables

private:
	//std::vector<FilterState> channelStates;
    int nfft = 2048;
    int l2nfft = 11;
#ifdef __APPLE__
    FFTSetup fft_s;
    DSPSplitComplex fft_in, fft_out, fft_out2, fft_buf, A_model, Hann;
#else
    kiss_fft_cfg fft_s, ifft_s;
    kiss_fft_cpx *fft_in, *fft_out, *fft_out2, *fft_buf, *A_model, *Hann;
#endif
    float * cbuf;
    float * spec_env;
    float * fft_mag;
    float * fft_mag_db;
    float * lms_h;
    int lms_n = 10;
    float ** synth_pulse;
    int n_synth_pulse = 10;
    int synth_pulse_ix = 0;
    int ncbuf = 4096;
    int cix = 0;
    int rix = 0;
    int maxgrain = 0;
    float rcnt = 256;
    float T = 400;
    int nmed = 5;
    float Tbuf[5];
    float Tsrt[5];
    int nAmpl = 50;
    float Abuf[50];
    int Tix;
    float pitchmark[3] = {0,-1,-1};
    int maxT = 600; // note nfft should be bigger than 3*maxT
    int cmask = ncbuf - 1;
    int voiced = 0;
	float sampleRate = 44100.0;
    float baseTuning = 440.0;
    int keycenter = 0;
    float midinotes[128];
    float midigain = 1.0;
    float harmgain = 0.0;
    float harmgain_target = 1.0;
    float voicegain = 0.0;
    float voicegain_target = 1.0;
    float dry_mix = 1.0;
    float speed = 1.0;
    float corr_strength = 0.5;
    float threshold = 0.2;
    int autotune = 1;
    int bypass = 0;
    
    int nvoices = 16;
    int voice_ix = 3;
    
    int chord_quality = 0;
    voice_t * voices;
    int inversion = 2;
    int midi_enable = 1;
    int synth_enable = 0;
    
    int midi_keycenter_cc = 16;
    int midi_keycenter_cc_offset = 1;
    int midi_keyquality_cc = 17;
    int midi_keyquality_cc_offset = 1;
    int midi_nvoices_cc = 18;
    int midi_nvoices_range = 0;
    int midi_inversion_cc = 19;
    int midi_inversion_range = 0;
    int midi_program_change_enable = 1;
    int midi_transmit_harmony = 0;
    int midi_transmit_melody = 0;
        
    int midi_legato = 0;
    int midi_pedal = 0;
    int auto_enable = 1;
    int midi_link = 1;
    int stereo_mode = StereoModeSplit;
    int n_auto = 4;
    int triad = -1;
    float interval_table[48];
    int interval_offsets[144];
    float * intervals;
    triad_ratio_t triads[18];
    chord_ratio_t major_chord_table[12];
    chord_ratio_t minor_chord_table[12];
    chord_ratio_t blues_chord_table[12];
    
    int ngrains;
    int grain_ix = 0;
    grain_t * grains;
    
    int graintablesize = maxT;
    float * grain_window;
    
    unsigned int sample_count = 0;
    unsigned int midi_changed_sample_num = 0;
    unsigned int midi_changed = 1;

    int preset_ix = 0;
        
    int voice_notes_old[N_AUTO];

//    AudioBufferList* inBufferListPtr = nullptr;
//    AudioBufferList* outBufferListPtr = nullptr;

public:

    float ** in_buffers;
    float ** out_buffers;

    float note_number = -1.0;
    float midi_note_number;
    int voice_notes[N_AUTO];
    int keys_down[128];
    int root_key = 0;
        
    int patch_number = 3;

    std::string preset_names[9] = {"Chords","Diatonic","Chromatic","Barbershop","JustMidi","Bohemian?","Bass!","4ths","Modes"};
    midi_event_t output_events[10];
    int n_output_events;
};

#endif /* FilterDSPKernel_hpp */
