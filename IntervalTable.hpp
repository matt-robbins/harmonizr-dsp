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
    IntervalTable(int n_voices, int n_quals, int n_tet) : m_quals{n_quals}, m_tet{n_tet}, m_voices{n_voices} {
        intervals.resize(n_voices * n_quals * n_tet);
    }

    void setInterval(int voice, int qual, int root, int interval) {
        intervals[addr(voice, qual, root)] = interval;
    }
    
    void setRaw(std::size_t address, int interval) {
        if (address >= intervals.size()) {
            std::cerr << "not writing interval with address " << address << "\n";
            return;
        }
        intervals[address] = interval;
    }
    
    int getInterval(int voice, int qual, int root) {
        return intervals[addr(voice, qual, root)];
    }
    
    void setAll(int * source, int N) {
        memcpy(intervals.data(), source, N * sizeof(int));
    }
    
    int getVoice(int raw_kcenter) {
        return raw_kcenter % m_voices;
    }
    
    int getRoot(int raw_kcenter) {
        return (raw_kcenter / m_voices) % m_tet;
    }
    
    int getQuality(int raw_kcenter) {
        return (raw_kcenter / (m_tet * m_voices)) % m_quals;
    }
    
    int getRawKcenter(int voice, int qual, int root) {
        return voice + root * m_voices + qual * m_voices * m_tet;
    }
    
private:
    std::size_t addr(int voice, int qual, int root) {
        // error checking here is super important
        
        if (voice >= m_voices || voice < 0) {
            throw std::runtime_error("voice # out of range.");
        }
        
        if (qual >= m_quals || qual < 0) {
            throw std::runtime_error("quality # out of range.");
        }
        
        if (root >= m_tet || root < 0) {
            throw std::runtime_error("root key # is out of range.");
        }
        
        return static_cast<std::size_t>(voice + root * m_voices + qual * m_voices * m_tet);
    }
    std::vector<int> intervals;
    int m_quals = 3;
    int m_tet = 12;
    int m_voices = 4;
};
