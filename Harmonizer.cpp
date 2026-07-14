#include "Harmonizer.hpp"
#include <span>

stereoSample Harmonizer::compute_one(float in) {
    buffer.pushValue(in);
    fbuffer.pushValue(filter.compute_one(in));
    
    if (m_est_counter.update())
    {
        float p = pEst.estimate(fbuffer);
        m_voiced = p > 0.0f;
        if (m_voiced == true) {
            T = p;
        }
        
        updateVoices();
    }
    
    if (pMark.findMark(T,m_voiced)) {
        for (PsolaVoice &voice : voices) {
            voice.setGrainSource(buffer.getContiguous(0), pMark.mark-T,2*T);
            voice.setGain(m_voiced ? 1.0 : 0.0);
        }
    }
    stereoSample u = {0.f, 0.f};
    for (PsolaVoice &voice : voices) {
         u += voice.render();
    }
    
    return u;
}

void Harmonizer::compute(float *in, float *out[], int nch, int N){
    
    if (nch != 2) {
        throw std::runtime_error("Harmonizer needs 2 channels of output!");
    }
    
    filter.clear_bad();
    
    for (int ix = 0; ix < N; ix++) {
        stereoSample u = compute_one(in[ix]);
        out[0][ix] = u.l;
        out[1][ix] = u.r;
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
        int dist = (voices[k].midinote == note) ? 0 : abs(voices[k].midinote.prev() - note);
        if (dist == 0 && voices[k].isOn()) {
            // this voice is already what we're looking for
            voices[k].setMidiNote(note,vel_);
            return;
        }
        
        if (dist < min_dist) {
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
        if (m_voiced.prev() == true) {
            for (int k = 0; k < n_auto; k++) {
                sendMidiNoteOff(voices[k].getMidiNote(), 100);
            }
        }
        
        return;
    }
    
    int input_note = static_cast<int>(std::lrint(T_to_midi_note(T, sampleRate)));
    
    for (int k = 0; k < n_auto; k++)
    {
        int midinote = input_note + table.getInterval(k,key_quality,(input_note + N_TET - root_key) % N_TET);
        if (k > inversion){
            midinote -= N_TET;
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
