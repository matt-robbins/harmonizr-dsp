#ifndef _window_
#define _window_

class Window {
public:
    enum winType { Hann=0 };
    Window(winType type, int table_size);
    ~Window();

    float value(float position); // compute window value
    
private:
    winType type;
    int N;
    float * table;
};

#endif