//
//  main.h
//  speechProcess
//
//  Created by Collin Burger on 7/8/14.
//
//

#ifndef speechProcess_main_h
#define speechProcess_main_h

#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv/cvaux.hpp>
#include <iostream>
#include <fstream>
#include <sndfile.h>

struct snippet
{
    unsigned int start,stop; //start and stop of snippet in samples
    unsigned int startSegment,endSegment,length; //which segments of the audio the snippet belongs to
};

void findStartStopSamples(std::vector<snippet *> *channelSnips, int samplesPerSeg, float *channelData, int sampleRate, int numSegs, int numFrames);

#endif
