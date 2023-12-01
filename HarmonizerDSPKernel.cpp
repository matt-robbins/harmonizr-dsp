//
//  HarmonizerDSPKermel.cpp
//  iOSHarmonizerFramework
//
//  Created by Matthew E Robbins on 12/1/23.
//

// NB This file must be set to type "Objective-C++ source in XCode"
#include <stdio.h>
#include "HarmonizerDSPKernel.hpp"


float HarmonizerDSPKernel::loopPosition() {
    if (loop_n <= 0)
    {
        return (float) loop_ix / (float) loop_max;
    }
    return (float) loop_ix / (float) loop_n;
}
