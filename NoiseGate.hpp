//
//  NoiseGate.hpp
//  Harmonizer
//
//  Created by Matthew E Robbins on 4/1/26.
//

class IirTracker {
public:
    IirTracker(float p=0.999f) : p{p} { };
    float compute(float in);
private:
    float y[2] = {0.f, 0.f};
    float p;
};

class NoiseGate {
public:
    enum gateState {
        Open = 0,
        Closed
    };
    
    NoiseGate(float thresh_db=-36.f, float hyst_db=3.f, float fs=441000.f, float hold_s=0.25f);
    
    void compute(float * out, float * in, int n);
    float compute_one(float in);
    float get_gain();
    void set_thresh_db(float thresh_db);
    float get_thresh_db();
    
private:
    float env = 0.0;
    float gain = 0.0;
    float thresh = 0.0;
    float thresh_db = -100.0;
    float hyst_db = 10.0;
    float hyst = 10;
    float target_gain = 0.0;
    int t = 0;
    float hold_s = 1.f;
    float fs = 48000;
    int hold_t = 0;
    gateState state = Open;
    
    IirTracker env_tracker = IirTracker(0.999);
    IirTracker gain_tracker = IirTracker(0.999);
    
};


