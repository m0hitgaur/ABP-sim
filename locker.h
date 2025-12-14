#pragma once
#include <iostream>
#include <time.h>
using namespace std;

void print_progress(float progress,float timestarted) 
{   float elapsed_time=time(NULL) - timestarted;
    float avg_time_remaining=0;
    if(progress>0.0)
    {avg_time_remaining= elapsed_time*(1.0/progress -1.0);}
    const int bar_width = 50;
    cout << "\rProgress: [";
    int pos = bar_width * progress;
    for (int i = 0; i < bar_width; i++) {
        if (i < pos) cout << "=";
        else if (i == pos) cout << ">";
        else cout << " ";
    }
    if(avg_time_remaining<1e-1)cout << "] " << int(progress * 100.0) << " %"<<"  Estimated time remaining: N/A secs"<<flush;
    else cout << "] " << int(progress * 100.0) << " %"<<"  Estimated time remaining: "<< avg_time_remaining<<" secs ("<<avg_time_remaining/60<<" mins) ("<<avg_time_remaining/3600<<" hrs)"<<flush;
    
}