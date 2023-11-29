#ifndef _harmonizer_h
#define _harmonizer_h
#include "CircularAudioBuffer.hpp"
#include "GranularSynth.hpp"
#include "PitchEstimator.hpp"
#include "PitchMarker.hpp"
#include <vector>

class PsolaVoice {
    public:
        PsolaVoice(int maxOctave=4);
        ~PsolaVoice();
        void updateT(float newT);
        void compute(float *out[], int nch, int n);
    private:
        GranularSynth &g;
        CircularAudioBuffer &buffer;
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
        Harmonizer(int Bufsize, int nvoices);
        ~Harmonizer();

        void compute(float *out[], int nch, int N);
    private:
        int nvoices = 0;
        CircularAudioBuffer &buffer;
        PitchEstimatorYIN &pEst;
        PitchMarker &pMark;
        std::vector<PsolaVoice> voices;
};


#endif