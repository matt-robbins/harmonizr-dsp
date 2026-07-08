#include "Harmonizer.hpp"

void Harmonizer::compute(float *in[], float *out[], int nch, int N){
    filter.clear_bad();
    
    for (int ix = 0; ix < N; ix++) {
        buffer.pushValue(in[0][ix]);
        fbuffer.pushValue(filter.compute_one(in[0][ix]));
        float p = 0;
        if (m_est_counter.update())
        {
            p = pEst.estimate(fbuffer);
            m_voiced = p > 0.0f;
            if (m_voiced == true) {
                T = p;
            }
            
            updateVoices();
        }
        
        if (pMark.findMark(T,0.35)) {
            for (PsolaVoice &voice : voices) {
                voice.setGrainSource(buffer.getContiguous(0), pMark.mark-T,2*T);
            }
        }
        float u = 0;
        for (PsolaVoice &voice : voices) {
            u += voice.synthesizeOne();
        }
        out[0][ix] = u;
    }
}

void Harmonizer::setVoiceT(int voice_n, float T) {
    voices[voice_n].T = T;
}

void Harmonizer::setPitchEstPeriod(int per) {
    m_est_counter.setPeriod(per);
}

// NOTE: this never happens *during* processing. Always between blocks. Block lengths change.
void Harmonizer::addMidiNote(int note, int vel)
{
    int vel_ = midi_ignore_velocity ? 100 : vel;
    int min_ix = -1;
    int min_dist = 129;
    
    
    // look for one with the same note and take that if we can, or the empty one with
    // the closest last note to the one we want
    for (int k = n_auto; k < voices.size(); k++)
    {
        int dist = abs(voices[k].midinote - note);
        
        if (dist == 0 && voices[k].isOn()) {
            // this voice is already what we're looking for
            voices[k].setMidiNote(note,vel_);
            return;
        }
        
        if (dist < min_dist)
        {
            min_ix = k; min_dist = dist;
        }
    }

    if (min_ix > 0)
    {
        voices[min_ix].setMidiNote(note,vel_);
        return;
    }

}

// TODO: implement these
void Harmonizer::sendMidiNoteOn(int nn, int vel) {
    return;
}
void Harmonizer::sendMidiNoteOff(int nn, int vel) {
    return;
}

void Harmonizer::updateVoices ()
{
    if (m_voiced == false)
    {
        if (m_voiced.prev() == true)
        {
            for (int k = 0; k < n_auto; k++)
            {
                sendMidiNoteOff(voices[k].getMidiNote(), 100);
            }
        }
        
        return;
    }
    
    input_voice.calculate(T);
    
    for (int k = 0; k < n_auto; k++)
    {
        int midinote = input_voice.getMidiN() + table.getInterval(k,key_quality,root_key);
        if (k > inversion){
            midinote -= 12;
        }
        
        voices[k].setMidiNote(midinote, 65);
        
        if (voices[k].midinote.diff() != 0)
        {
            sendMidiNoteOn(midinote, 100);
            if (voices[k].midinote.prev() >= 0)
                sendMidiNoteOff(voices[k].midinote.prev(), 100);
        }
                
        // let the UI know
        //ui_voice_notes[k] = voices[k].midinote;
    }
    
    
}
