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
#import <complex>
#import <random>
#import <sys/time.h>
#include "ButterworthFilter.hpp"
#include "CircularAudioBuffer.hpp"
#include "PitchEstimator.hpp"

#ifdef __APPLE__
#import "DSPKernel.hpp"
#import "ParameterRamper.hpp"
#include <Accelerate/Accelerate.h>
#include <dispatch/dispatch.h>

//#include "TestAudioData.h"

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

#endif //#ifdef APPLE

// Utility functions
float linear (float *v, float a);
float linear_interp(float *v, float ix);
inline float cubic (float *v, float a);
float cubic_interp(float *v, float ix);
float quadratic_peak(float *pv, float *v, int ix);

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

#define N_AUTO 4
#define N_MIDI_BUF 100

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

typedef float chord_ratio_t[4];

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
    HarmParamVibrato,
    HarmParamLoop,
    HarmParamInterval
};

enum loopMode {
    LoopStopped=0,
    LoopRec,
    LoopPlay,
    LoopPlayRec,
    LoopPause
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
    
// MARK: Member Functions
private:
    inline float window_value(float f);
    void psola(float *out, float *out2, int n);
    void findmark (void);
    
    void update_voices (void);
    
    inline float measure_snr(float in);
    float pitch_resample();
    
   // float estimate_pitch(int start_ix);
   // float estimate_pitch2(int start_ix);
    float get_minphase_pulse(int start_ix);
    float get_model(int start_ix);
    
    void analyze_harmony(void);
    
    void send_note_on(int nn, int vel);
    void send_note_off(int nn, int vel);
    
public:
    int n_channels = 0;

    HarmonizerDSPKernel() : 
        raw_buffer { CircularAudioBuffer(nfft*2)},
        filtered_buffer { CircularAudioBuffer(nfft*2)},
        pitchEstimator { PitchEstimatorYIN(maxT,l2nfft,threshold,nmed) }
    {
        
        fprintf(stderr, "bufsize = %d\n", raw_buffer.getSize());
    }
	
    void init(int inChannels, int outChannels, double inSampleRate);
    void fini();
    void reset();
    
    void setParameter(param_address_t address, param_value_t value);
    param_value_t getParameter(param_address_t address);
    float loopPosition();

    void setBuffers(float ** in, float ** out);

    void setPreset(int preset_ix_);
    int getPreset();
    
    void addnote(int note, int vel);
    void remnote(int note);
    void pedal_down();
    void pedal_up();
    
#ifdef __APPLE__

    void startRamp(AUParameterAddress address, AUValue value, AUAudioFrameCount duration) override;
    void setBuffers(AudioBufferList* inBufferList, AudioBufferList* outBufferList);
    virtual void handleMIDIEvent(midi_event_t const& midiEvent) override;
    
    void process(frame_count_t frameCount, frame_count_t bufferOffset) override;

#else
    void process(frame_count_t frameCount, frame_count_t bufferOffset);
#endif
        
    
    // MARK: Member Variables
private:
	//std::vector<FilterState> channelStates;
    int l2nfft = 11;
    int nfft = 1 << l2nfft;
    int maxT = 600; // note nfft should be bigger than 3*maxT
    float threshold = 0.2; // for YIN pitch estimator
    int nmed = 7; // length of median filter for YIN estimator
    
#ifdef __APPLE__
    FFTSetup fft_s;
    DSPSplitComplex fft_in, fft_out, fft_out2, fft_buf, A_model, Hann, spec_env;
#else
    kiss_fft_cfg fft_s, ifft_s;
    kiss_fft_cpx *fft_in, *fft_out, *fft_out2, *fft_buf, *A_model, *Hann, *spec_env;
#endif
        
    ButterworthFilter filter;
    CircularAudioBuffer raw_buffer;
    CircularAudioBuffer filtered_buffer;
    PitchEstimatorYIN pitchEstimator;
    
    float * in_filt;
    float * cbuf;
    float * fbuf;
    float * cmdf;
    float * snr_buf;
    float * fd_lpf;
    float * fft_mag;
    float * fft_mag_db;
    float * lms_h;
    int lms_n = 10;
    float ** synth_pulse;
    float ** loop_buf;
    int loop_max;
    int loop_n = 0;
    int loop_ix = 0;
    int loop_xf = 0;
    int loop_xfn = 500;
    int n_synth_pulse = 20;
    int synth_pulse_ix = 0;
    int ncbuf = 4096;
    int n_snr = 400;
    int snr_ix = 0;
    int cix = 0;
    int rix = 0;
    int maxgrain = 0;
    
    float mean_sq = 0;
    float nse_floor = 1.0;
    float rcnt = 256;
    float T = 400;
    float Tbuf[7];
    float Tsrt[7];
    int nAmpl = 50;
    float Abuf[50];
    int Tix;
    float pitchmark[3] = {0,-1,-1};
    
    int cmask = ncbuf - 1;
    int voiced = 0;
    float noise_pct = 0.0;
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
    float vibrato = 0.;
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
    
public:
    int voice_notes[N_AUTO];
    float note_number = -1.0;
    float midi_note_number;
    int keys_down[128];
    int root_key = 0;
    int patch_number = 0;
    float rms = 0;
    
    float ** in_buffers;
    float ** out_buffers;

    int loop_mode = LoopStopped;
        
    std::string preset_names[9] = {"Chords","Diatonic","Chromatic","Barbershop","JustMidi","Bohemian?","Bass!","4ths","Modes"};
    midi_event_t output_events[N_MIDI_BUF];
    int max_output_events = N_MIDI_BUF;
    int n_output_events;
};

#endif /* FilterDSPKernel_hpp */
