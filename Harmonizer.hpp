#ifndef _harmonizer_h
#define _harmonizer_h
#include "CircularAudioBuffer.hpp"
#include "GranularSynth.hpp"
#include "PitchEstimator.hpp"
#include "PitchMarker.hpp"
#include <vector>
#include <iostream>

class PsolaVoice {
    public:
        PsolaVoice(int maxOctave=4)
            : g(20) {
        }
        ~PsolaVoice() {}
        void updateT(float newT);
        GranularSynth g;

    private:
        int graintable_size = 20;
        float gain = 1.0;
        float gain_target = 1.0;
        float pitch_ratio = 1.0;
        float pitch_ratio_target = 1.0;
        float formant_ratio = 1.0;
        float formant_ratio_target = 1.0;
        float nextgrain = 200.0;
        float pan = 0.0;

        float vibrato_rate = 5.0; //Hz
        float vibrato_amp = 0.0;
        float vib_phase = 0.0;
        int midinote = -1; // -1 means note is off;
        int prev_midinote = -1;
        float midinote_ = 69.0; // midi note 69 is the pitch A(440)
        int midivel;

        float T = 200.0; // T<=0 means we're unvoiced
};


class Harmonizer {
    public:
        Harmonizer(int l2bufsize, int nvoices, int maxT) : 
            maxT(maxT),
            buffer(l2bufsize),
            pEst(maxT, l2bufsize - 1, 0.2, 7),
            pMark(buffer,maxT)
        {
            for (int k = 0; k < nvoices; k++) {
                voices.push_back(PsolaVoice());
            }
        }
        ~Harmonizer() {}

        void compute(float *in[], float *out[], int nch, int N);
        void setVoiceT(int voice_n, float T);
        void setPitchEstPeriod(int per);

    private:
        int nvoices = 0;
        int maxT = 20;
        float T = 10;
        int p_est_period = 256;
        int p_est_ix = 0;
        CircularAudioBuffer buffer;
        PitchEstimatorYIN pEst;
        PitchMarker pMark;
        std::vector<PsolaVoice> voices;
};


#endif