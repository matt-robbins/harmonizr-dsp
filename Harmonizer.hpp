#ifndef _harmonizer_h
#define _harmonizer_h
#include "CircularAudioBuffer.hpp"
#include "GranularSynth.hpp"
#include "PitchEstimator.hpp"
#include "PitchMarker.hpp"
#include "IntervalTable.hpp"
#include "ButterworthFilter.hpp"
#include "Util.hpp"
#include <vector>
#include <iostream>

class PsolaVoice : public GranularSynth {
public:
    PsolaVoice(int maxOctave=4) : note_track{ IirTracker(1000) }, vel_track{ IirTracker(1000)} {}
    ~PsolaVoice() {}
    
    void setMidiNote(int nn, int vel) {
        midinote = nn;
        midivel = vel;
        T = midi_note_to_T(static_cast<float>(midinote), fs);
    }
    
    int getMidiNote() {
        if (midinote < 0)
            return midinote.prev();
        
        return midinote;
    }
    bool isOn() {
        return midinote >= 0;
    }
    
    void setGain(float gain) {
        GranularSynth::setGain(gain * static_cast<float>(midivel)/100.f);
    }
    
    void setPan(float pan) {
        this->pan = pan;
    }
    
    struct stereoSample {
        float l = 0.f;
        float r = 0.f;
    };
    
    stereoSample render() {
        if (midinote > 0)
        {
            T = midi_note_to_T(static_cast<float>(midinote), fs);
        }
        
        float u = GranularSynth::synthesizeOne();
        return stereoSample {u * (pan + 1.f)/2.f,u * (pan - 1.f)/2.f};
    }
    
    StateMemVariable<int> midinote = -1;
private:
    
    IirTracker note_track;
    IirTracker vel_track;
    int midivel = 0;
    
    float pan = 0.0;
};

class InputVoice {
public:
    InputVoice(float sampleRate, float baseTuning) : m_sampleRate{sampleRate}, m_baseTuning{baseTuning} {}
    void calculate(float T) {
        m_midiF = m_midiA4 + log2f (m_sampleRate / (T * m_baseTuning)) * m_nTet;
        m_midiN = (int) roundf(m_midiF);
    }
    float getMidiF() {
        return m_midiF;
    }
    int getMidiN() {
        return m_midiN.value();
    }
private:
    float m_midiF;
    StateMemVariable<int> m_midiN = {-1};
    StateMemVariable<bool> m_voiced = {false};
    float m_T;
    float m_sampleRate;
    float m_baseTuning;
    float m_midiA4 = 69.0;
    float m_nTet = 12.0f;
    float m_hystThresh = 0.65;
    
};

class Harmonizer {
public:
    Harmonizer(int l2bufsize, int nvoices, int maxT, float sampleRate) :
        maxT{maxT},
        buffer{l2bufsize+1},
        fbuffer{l2bufsize+1},
        pEst{maxT, l2bufsize - 1, 0.2, 7},
        pMark{buffer,(float)maxT},
        table{n_auto, n_qual, n_tet},
        sampleRate{sampleRate},
        input_voice{sampleRate,baseTuning}
    {
        for (int k = 0; k < nvoices; k++) {
            voices.push_back(PsolaVoice());
        }
    }
    ~Harmonizer() {}

    void compute(float *in[], float *out[], int nch, int N);
    void setVoiceT(int voice_n, float T);
    void setPitchEstPeriod(int per);
    void updateVoices();
    void addMidiNote(int nn, int vel);
    void sendMidiNoteOn(int nn, int vel);
    void sendMidiNoteOff(int nn, int vel);
    
    void setRootKey(int rk) {
        root_key = rk;
    }
    void setKeyQuality(int kq) {
        key_quality = kq;
    }

private:
    int maxT = 400;
    float T = 100;
    StateMemVariable<bool> m_voiced = {false};
    Counter m_est_counter = {256};
    float sampleRate = 44100;
    float baseTuning = 440.0;
    int n_auto = 4;
    int n_tet = 12;
    int n_qual = 3;
    
    int stereo_mode = 0;
    unsigned int sample_count = 0;
    
    bool midi_ignore_velocity = false;
    
    int root_key = 0;
    int key_quality = 0;
    int inversion = 0;
    
    CircularAudioBuffer buffer, fbuffer;
    ButterworthFilter filter = {};
    PitchEstimatorYIN pEst;
    PitchMarker pMark;
    IntervalTable table;
    std::vector<PsolaVoice> voices;
    std::vector<PsolaVoice> midi_voices;
    InputVoice input_voice;
};


#endif
