//
//  main.cpp
//  speechProcess
//
//  Created by Collin Burger on 7/8/14.
//
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include <cmath>
#include <sndfile.h>

using namespace std;

struct snippet
{
    unsigned int start,stop; //start and stop of snippet in samples
    unsigned int startSegment,endSegment,length; //which segments of the audio the snippet belongs to
};

void findStartStopSamples(vector<snippet *> *channelSnips, int samplesPerSeg, float *channelData, int sampleRate, int numSegs, int numFrames);

int main(int argc, char** argv)
{
    if (argc < 2) {
        printf("Missing argument. Must include path of wav file to be processed.\n");
        return -1;
    }
    else if (argc > 2){
        printf("Too many argruments. Only indicate the path of the .wav file.\n");
        return -1;
    }
    
    
    SF_INFO soundInfo;
    char *wavFileName = argv[1];
    SNDFILE *soundFile = sf_open(wavFileName, SFM_READ, &soundInfo); //open wav file
    if(NULL == soundFile){
        printf("Error opening wav file.\n");
        return -1;
    }
    else if(sf_error(soundFile) != 0){
        string err = sf_strerror(soundFile);
        cout << err << endl;
        return -1;
    }
    
    if(soundInfo.channels != 2){
        printf("Audio file must have 2 channels.\n");
        return -1;
    }
    
    long dataLength = soundInfo.frames * soundInfo.channels; //total amount of samples
    float *audioData = new float[dataLength]; //allocate space for the audio data
    sf_readf_float(soundFile, audioData, soundInfo.frames); //read samples into data array as floats
    
    //create an array for each channel and copy the samples into them
    float *channelOneData = new float[dataLength/2];
    float *channelTwoData = new float[dataLength/2];
    //the samples for each channel alternate
    for (int i = 0; i < dataLength-1; i+=2){
        int indexTwo = ceil((float)(i+1)/2);
        channelOneData[i/2] = audioData[i];
        channelTwoData[indexTwo] = audioData[i+1];
    }
    
    //Section for writing files to check the channels
    /*SF_INFO soundWrite;
    soundWrite.samplerate = soundInfo.samplerate;
    soundWrite.channels = 1;
    soundWrite.frames = soundInfo.frames;
    soundWrite.format = soundInfo.format;
    soundWrite.sections = soundInfo.sections;
    soundWrite.seekable = soundInfo.seekable;
    SNDFILE *channelOneWrite = sf_open("channelOne.wav", SFM_WRITE, &soundWrite);
    sf_write_float(channelOneWrite, channelOneData, dataLength/2);
    SNDFILE *channelTwoWrite = sf_open("channelTwo.wav", SFM_WRITE, &soundWrite);
    sf_write_float(channelTwoWrite,channelTwoData,dataLength/2);
    */
    
    //calculate the power of the each segment of signal
    int numSegs = ceil((float)soundInfo.frames/soundInfo.samplerate) * 2; //amount of 1/2 second segments
    int samplesPerSeg = (int)soundInfo.samplerate/2;
    float *sumSegsChannelOne = new float[numSegs];
    float *sumSegsChannelTwo = new float[numSegs];
    int segment = 0;
    for (int i = 0; i < soundInfo.frames; i+= soundInfo.samplerate/2) {
        int j;
        for (j = 0; j < soundInfo.samplerate/2; j++){
            int index = i+j;
            if (index > soundInfo.frames || index > soundInfo.frames) //test for end of array
                break;
            
            //power calculations
            sumSegsChannelOne[segment] += pow(abs(channelOneData[index]), 2.0);
            sumSegsChannelTwo[segment] += pow(abs(channelTwoData[index]), 2.0);
        }
        segment++;
    }
    
    
    //this section analyzes the adjacent segments to identify a section of audio with activity
    bool *maskChannelOne = new bool[numSegs];
    bool *maskChannelTwo = new bool[numSegs];
    int maskCountChannelOne = 0;
    int maskCountChannelTwo = 0;
    vector<struct snippet*> snipsChannelOne;
    vector<struct snippet*> snipsChannelTwo;
    for (int i = 0; i < numSegs; i++) {
        if (sumSegsChannelOne[i] > 1.0){ //test for audio activity in segment of channel one
            maskChannelOne[i] = true;
            maskCountChannelOne++;
            if (maskCountChannelOne > 60)
                maskCountChannelOne = 0;
        }
        else{
            maskChannelOne[i] = false;
            if (maskCountChannelOne >= 3) {
                snippet *snippy = new snippet;
                snippy->length = maskCountChannelOne;
                snippy->startSegment = i - maskCountChannelOne;
                snippy->endSegment = i - 1;
                snippy->start = 0;
                snippy->stop = 0;
                snipsChannelOne.push_back(snippy);
            }
            maskCountChannelOne = 0;
        }
        
        
        if (sumSegsChannelTwo[i] > 1.0){ //test for audio activity in segment of channel two
            maskChannelTwo[i] = true;
            maskCountChannelTwo++;
            if (maskCountChannelTwo > 60)
                maskCountChannelTwo = 0;
        }
        else{
            maskChannelTwo[i] = false;
            if (maskCountChannelTwo >= 3) {
                snippet *snippy = new snippet;
                snippy->length = maskCountChannelTwo;
                snippy->startSegment = i - maskCountChannelTwo;
                snippy->endSegment = i - 1;
                snippy->start = 0;
                snippy->stop = 0;
                snipsChannelTwo.push_back(snippy);
            }
            maskCountChannelTwo = 0;
        }
        //cout << "sec " << (double)i/2 << ": " << maskChannelTwo[i] << endl;
    }

    findStartStopSamples(&snipsChannelOne,samplesPerSeg,channelOneData,soundInfo.samplerate,numSegs,(int)soundInfo.frames);
    findStartStopSamples(&snipsChannelTwo, samplesPerSeg, channelTwoData, soundInfo.samplerate, numSegs, (int)soundInfo.frames);
    
    //csv output
    stringstream fileOut;
    string wavFileString = wavFileName;
    int dot = (int)wavFileString.find_last_of(".");
    string rawName = wavFileString.substr(0,dot);
    
    fileOut << rawName << "_channel_1.csv";
    string fileOutOne = fileOut.str();
    ofstream outputOne (fileOutOne,ofstream::out);
    for (int i = 0; i < snipsChannelOne.size(); i++) { //change output to milliseconds
        outputOne << floor((float)(snipsChannelOne.at(i)->start)/soundInfo.samplerate*1000) << "," <<  ceil((float)snipsChannelOne.at(i)->stop/soundInfo.samplerate*1000) << endl;
    }
    outputOne.close();
    
    fileOut.str("");
    fileOut << rawName << "_channel_2.csv";
    string fileOutTwo = fileOut.str();
    ofstream outputTwo (fileOutTwo,ofstream::out);
    for (int i = 0; i < snipsChannelTwo.size(); i++) {
        outputTwo << floor((float)snipsChannelTwo.at(i)->start/soundInfo.samplerate*1000) << "," << ceil((float)snipsChannelTwo.at(i)->stop/soundInfo.samplerate*1000) << endl;
    }
    outputTwo.close();
    
    
    //close and free everything;
    sf_close(soundFile);
    //sf_close(channelOneWrite);
    //sf_close(channelTwoWrite);
    delete audioData;
    delete sumSegsChannelOne;
    delete sumSegsChannelTwo;
    delete channelOneData;
    delete channelTwoData;
    delete maskChannelOne;
    delete maskChannelTwo;
    for(std::vector<struct snippet*>::iterator it = snipsChannelOne.begin(); it != snipsChannelOne.end(); ++it) {
        free(*it);
    }
    snipsChannelOne.clear();
    for(std::vector<struct snippet*>::iterator it = snipsChannelTwo.begin(); it != snipsChannelTwo.end(); ++it) {
        free(*it);
    }
    snipsChannelTwo.clear();
    return 0;
}


//--------------------------------------------------------
//This function makes sure that the audio snippets are not truncated at the beginning or end
void findStartStopSamples(vector<snippet *> *channelSnips, int samplesPerSeg, float *channelData, int sampleRate, int numSegs, int numSamples){
    float *dataTemp = new float[samplesPerSeg];
    for (int i = 0; i < channelSnips->size(); i++) {
        snippet *current = channelSnips->at(i);
        if (current->startSegment == 0)
            current->start = 0;
        else{
            memcpy(dataTemp,channelData + (current->startSegment - 1)*samplesPerSeg,samplesPerSeg*sizeof(float));
            for (int j = 0; j < samplesPerSeg; j+= sampleRate/1000){ //skip by millisecond
                float milliPower = 0;
                for (int k = 0; k < 8; k++)
                    milliPower += pow(abs(dataTemp[j+k]),2.0);
                if (milliPower > 0.001) {
                    current->start = (current->startSegment - 1)*samplesPerSeg + j; //store in samples
                    break;
                }
            }
            if (current->start == 0)
                current->start = (current->startSegment)*samplesPerSeg; //audio is not truncated at beginning
        }
        if (current->endSegment == numSegs - 1)
            current->stop = numSamples;
        else{
            memcpy(dataTemp,channelData + (current->endSegment + 1)*samplesPerSeg,samplesPerSeg*sizeof(float));
            for (int j = samplesPerSeg - 1; j >= 0 ; j-= sampleRate/1000){ //go by millisecond
                float milliPower = 0;
                for (int k = 7; k >= 0; k--)
                    milliPower += pow(abs(dataTemp[j+k]),2.0);
                if (milliPower > 0.001) {
                    current->stop = (current->endSegment + 1)*samplesPerSeg + j; //in samples
                    break;
                }
            }
            if (current->stop == 0)
                current->stop = (current->endSegment + 1)*samplesPerSeg; //audio is not truncated at end
        }
    }
    delete dataTemp;
}


