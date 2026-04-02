#include "Harmonizer.hpp"

void Harmonizer::compute(float *in[], float *out[], int nch, int N){
    for (int ix = 0; ix < N; ix++) {
        buffer.pushValue(in[0][ix]);
        float p = 0;
        if (++p_est_ix > p_est_period) {
            p_est_ix = 0;

            p = pEst.estimate(buffer);
            if (p > 0) T = p;

            //std::cout << "p = " << T << "\n";
        }
        
        if (pMark.findMark(T,0.35)) {
            //std::cerr << "mark =" << pMark.mark << std::endl;
            //std::cout << "mark=" << pMark.mark << " wp=" << buffer.getWriteIndex() - pMark.mark <<"\n";
            //float *data = buffer.getContiguous(0);
            // for (int k = 0; k < 2*T+1; k++) {
            //     std::cout << data[k + (int) (pMark.mark - T)] << ",";
            // }
            // std::cout << std::endl;

            for (PsolaVoice &voice : voices) {
                voice.g.setGrainSource(buffer.getContiguous(0), pMark.mark-T,2*T);
                //std::cerr << "set grain source " << k << " to " << voices[k].g.getGrainSource() << std::endl;
            }
        }
        float u = 0;
        for (PsolaVoice &voice : voices) {
            //std::cerr << "synthesizing on " << k << std::endl;
            u += voice.g.synthesizeOne();
        }
        out[0][ix] = u;
    }
}

void Harmonizer::setVoiceT(int voice_n, float T) {
    //std::cerr << "voices:" << voices.capacity() << std::endl;
    //std::cerr << "setting period to " << T << std::endl;
    voices[voice_n].g.T = T;
}

void Harmonizer::setPitchEstPeriod(int per) {
    p_est_period = per;
}
