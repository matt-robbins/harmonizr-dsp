//
//  NoiseGate.hpp
//  Harmonizer
//
//  Created by Matthew E Robbins on 4/1/26.
//

class NoiseGate {
public:
    enum gateState {
        Open = 0,
        Closed
    };
    
    NoiseGate(float thresh_db, float hyst_db, int fs, int hold_t);
    
    void compute(float * out, float * in, int n);
    float compute_one(float in);
    float get_gain();
    
private:
    float env = 0.0;
    float gain = 0.0;
    float thresh = 0.0;
    float thresh_db = -100.0;
    float hyst_db = 10.0;
    float hyst = 10;
    float target_gain = 0.0;
    int t = 0;
    int hold_t = 48000;
    int fs = 48000;
    gateState state = Open;
};
