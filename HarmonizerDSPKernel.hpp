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
#include "Window.hpp"
#include "PitchEstimator.hpp"
#include "PitchMarker.hpp"
#include "SimplePitchShifter.hpp"
#include "GranularSynth.hpp"
#include "SpectralProcessor.hpp"
#include "Looper.h"
#include "NoiseGate.hpp"
#include "Harmonizer.hpp"

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

#define N_AUTO 4
#define N_MIDI_BUF 100

typedef struct voice_s
{
    float ratio;
    float target_ratio;
    float pan;
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
    HarmParamGateThresh,
    HarmParamAlgorithm,
    HarmParamMidi,
    HarmParamMidiLink,
    HarmParamMidiLegato,
    HarmParamMidiVelIgnore,
    HarmParamMidiKeyCC,
    HarmParamMidiKeyCcOffset,
    HarmParamMidiQualCC,
    HarmParamMidiQualCcOffset,
    HarmParamMidiNvoiceCC,
    HarmParamMidiNvoiceCcRange,
    HarmParamMidiInvCC,
    HarmParamMidiInvCcRange,
    HarmParamMidiFreezeCC,
    HarmParamMidiPC,
    HarmParamMidiHarmOut,
    HarmParamMidiMelOut,
    HarmParamMidiPedalFcn,
    HarmParamMidiPedalInv,
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
    HarmParamInputChannel,
    HarmParamSynth,
    HarmParamVibrato,
    HarmParamFreeze,
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

enum pedalMode {
    PedalFreeze=0,
    PedalNotes,
    PedalBoth
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

enum {
    AlgorithmPSOLA=0,
    AlgorithmResample
};

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
    
    void analyze_harmony(void);
    void send_note_on(int nn, int vel);
    void send_note_off(int nn, int vel);
    
public:
    int n_channels = 0;

    HarmonizerDSPKernel() : 
        raw_buffer { CircularAudioBuffer(l2nfft+1) },
        filtered_buffer { CircularAudioBuffer(l2nfft+1) },
        noise_gate { NoiseGate()},
        pitchEstimator { PitchEstimatorYIN(maxT,l2nfft,threshold,nmed) },
        pitchMarker { PitchMarker(raw_buffer,maxT) },
        window { Window(Window::Hann, nfft)},
        looper { Looper()}
    {
        for (int i = 0; i < nvoices; i++){
            simpleVoices.push_back(SimplePitchShifter(raw_buffer,window,maxT));
            psolaVoices.push_back(PsolaVoice());
        }
        fprintf(stderr, "bufsize = %d\n", raw_buffer.getSize());
    }
	
    void init(int inChannels, int outChannels, double inSampleRate);
    void fini();
    void reset();
    
    void setParameter(param_address_t address, param_value_t value);
    param_value_t getParameter(param_address_t address);
    float loopPosition();
    int getLoopMode();
    int setLoopMode(int mode);

    void setBuffers(float ** in, float ** out);

    void setPreset(int preset_ix_);
    int getPreset();
    
    void addnote(int note, int vel);
    void remnote(int note);
    void update_midi();
    void calculate_voices();
    void update_voices (void);
    void set_voiced(bool voiced);
    void list_intervals();

    void pedal_down();
    void pedal_up();
    
    void set_T(float t);
    
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
    int l2nfft = 12;
    int nfft = 1 << l2nfft;
    int maxT = 1000; // note nfft should be bigger than 3*maxT
    int minT = 25; // corresponds to A6 (basically impossible)
    float threshold = 0.2; // for YIN pitch estimator
    int nmed = 7; // length of median filter for YIN estimator
        
    ButterworthFilter filter;
    CircularAudioBuffer raw_buffer;
    CircularAudioBuffer filtered_buffer;
    NoiseGate noise_gate;
    Window window;
    PitchEstimatorYIN pitchEstimator;
    PitchMarker pitchMarker;
    std::vector<SimplePitchShifter> simpleVoices;
    std::vector<PsolaVoice> psolaVoices;
    
    std::vector<SpectralProcessor> freezers;
    
    Looper looper;
    
    int nvoices = 16;
    
    float rcnt = 256;
    float T = 400;
    
    int voiced = 0;
    float noise_pct = 0.0;
    float sampleRate = 44100.0;
    float baseTuning = 440.0;
    int keycenter = 0;
    float midinotes[128];
    float midigain = 1.0;
    float harmgain = 1.0;
    float harmgain_target = 1.0;
    float voicegain = 1.0;
    float voicegain_target = 1.0;
    float dry_mix = 1.0;
    float speed = 1.0;
    float corr_strength = 0.5;
    float vibrato = 0.;
    
    float gate_thresh = -40.f;
    
    int autotune = 1;
    int bypass = 0;
    
    int chord_quality = 0;
    int inversion = 2;
    int midi_enable = 1;
    int synth_enable = 0;
    
    int algorithm = AlgorithmPSOLA;
    
    int midi_keycenter_cc = 16;
    int midi_keycenter_cc_offset = 1;
    int midi_keyquality_cc = 17;
    int midi_keyquality_cc_offset = 1;
    int midi_nvoices_cc = 18;
    int midi_nvoices_range = 0;
    int midi_inversion_cc = 19;
    int midi_inversion_range = 0;
    int midi_freeze_cc = 20;
    int midi_program_change_enable = 1;
    int midi_transmit_harmony = 0;
    int midi_transmit_melody = 0;
        
    int midi_ignore_velocity = 0;
    int midi_pedal = 0;
    int midi_link = 1;
    int midi_pedal_fcn = 0;
    int midi_pedal_inv = 0;
    int stereo_mode = StereoModeSplit;
    int in_channel = 0;
    int n_auto = 4;
    int triad = -1;
    float interval_table[48];
    int interval_offsets[144];
    float * intervals;
    chord_ratio_t major_chord_table[12];
    chord_ratio_t minor_chord_table[12];
    chord_ratio_t blues_chord_table[12];
    
    unsigned int sample_count = 0;
    unsigned int midi_changed_sample_num = 0;

    int preset_ix = 0;
    
public:
    int ui_voice_notes[N_AUTO];
    int midi_note_number = -1;
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
