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

struct stereoSample {
    float l = 0.f;
    float r = 0.f;
    stereoSample& operator+=(const stereoSample& rhs) {
        l += rhs.l;
        r += rhs.r;
        return *this;
    }
};

class PsolaVoice : public GranularSynth {
public:
    PsolaVoice(int fs) : note_track{ IirTracker(1000) }, vel_track{ IirTracker(1000)} {
        this->fs = fs;
    }
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
    
    stereoSample render() {
        if (midinote > 0)
        {
            T = midi_note_to_T(static_cast<float>(midinote), fs);
        }
        if (!isOn()) {
            return stereoSample {};
        }
        
        float u = GranularSynth::synthesizeOne();
        float pan_ = (pan + 1.f)/2.f;
        return stereoSample {u * (1.f - pan_),u * pan_};
    }
    
    StateMemVariable<int> midinote = -1;
private:
    
    IirTracker note_track;
    IirTracker vel_track;
    int midivel = 0;
    
    float pan = 0.0;
};

class Harmonizer {
public:
    Harmonizer(int l2bufsize, int nvoices, int maxT, float sampleRate) :
        maxT{maxT},
        buffer{l2bufsize+1},
        fbuffer{l2bufsize+1},
        pEst{maxT, l2bufsize - 1, 0.2, 7},
        pMark{buffer,(float)maxT},
        table{N_AUTO, N_QUAL, N_TET},
        sampleRate{sampleRate}
    {
        for (int k = 0; k < nvoices; k++) {
            voices.push_back(PsolaVoice(sampleRate));
        }
    }
    ~Harmonizer() {}

    stereoSample compute_one(float in);
    void compute(float *in, float *out[], int nch, int N);
    void setVoiceT(int voice_n, float T);
    void setPitchEstPeriod(int per);
    void updateVoices();
    void addMidiNote(int nn, int vel);
    void sendMidiNoteOn(int nn, int vel);
    void sendMidiNoteOff(int nn, int vel);
    
    void setRootKey(int rk) {
        if (rk < 0 || rk >= N_TET) {
            throw std::runtime_error("Key root must be between 0 and " + std::to_string(N_TET-1) + "!");
        }
        root_key = rk;
    }
    void setKeyQuality(int kq) {
        if (kq < 0 || kq >= N_QUAL) {
            throw std::runtime_error("Key quality must be between 0 and " + std::to_string(N_QUAL-1) + "!");
        }
        key_quality = kq;
    }
    
    void setAutoIntervals(int * intervals, int N) {
        if (N_AUTO * N_QUAL * N_TET != N) {
            throw std::runtime_error("Interval Table size must match!");
        }
        
        table.setAll(intervals, N);
    }
    
    void setInterval(std::size_t interval_n, int value) {
        table.setRaw(interval_n, value);
    }
    
    void setAutoVoiceCount(int n) {
        n_auto = n;
    }
    
    void setInversion(int n) {
        inversion = n;
    }
    
    int getAutoVoiceCount() {
        return n_auto;
    }
    
    float getT() {
        return T;
    }
    
    void getAutoNotes(int * notes, int N)
    {
        for (int k = 0; k < n_auto; k++) {
            notes[k] = (int) lrintf(voices[k].midinote);
        }
    }

private:
    static constexpr int N_AUTO = 4;
    static constexpr int N_QUAL = 3;
    static constexpr int N_TET = 12;
    int maxT = 1000;
    float T = 400;
    StateMemVariable<bool> m_voiced = {false};
    Counter m_est_counter = {256};
    float sampleRate = 44100;
    float baseTuning = 440.0;
    int n_auto = N_AUTO;
    
    int stereo_mode = 0;
    unsigned int sample_count = 0;
    
    bool midi_ignore_velocity = false;
    
    int root_key = 0;
    int key_quality = 0;
    int inversion = 2;
    
    CircularAudioBuffer buffer, fbuffer;
    ButterworthFilter filter = {};
    PitchEstimatorYIN pEst;
    PitchMarker pMark;
    IntervalTable table;
    std::vector<PsolaVoice> voices;
};


#endif
