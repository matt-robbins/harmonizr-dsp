//
//  IntervalTable.hpp
//  Harmonizer
//
//  Created by Matthew E Robbins on 6/24/26.
//
#include <vector>

class IntervalTable
{
public:
    IntervalTable(int n_voices, int n_types, int n_tet) : m_types{n_types}, m_tet{n_tet}, m_voices{n_voices} {
        intervals.resize(n_voices * n_types * n_tet);
    }

    void setInterval(int voice, int type, int root, int interval) {
        intervals[addr(voice, type, root)] = interval;
    }
    
    int getInterval(int voice, int type, int root) {
        return intervals[addr(voice, type, root)];
    }
    
    void setAll(int * source, int N) {
        memcpy(intervals.data(), source, N * sizeof(int));
    }
    
    int getRoot(int raw_kcenter) {
        return raw_kcenter % m_tet;
    }
    
    int getQuality(int raw_kcenter) {
        return raw_kcenter / m_tet;
    }
    
    int getRawKcenter(int root, int quality) {
        return quality * m_tet + root;
    }
    
private:
    int addr(int voice, int type, int root) {
        return voice + root * m_voices + type * m_voices * m_tet;
    }
    std::vector<int> intervals;
    int m_types = 3;
    int m_tet = 12;
    int m_voices = 4;
};
