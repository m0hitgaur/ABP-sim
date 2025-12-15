#pragma once
#include <iostream>
#include <time.h>
using namespace std;

void print_progress(double progress,double timestarted) 

{   double current_time=static_cast <double>(time(NULL));
    double elapsed_time=current_time - timestarted;
    if(progress>0.0)
        {double avg_time_remaining= elapsed_time*(1.0/progress -1.0);
        const int bar_width = 50;
        cout << "\rProgress: [";
        int pos = bar_width * progress;
        for (int i = 0; i < bar_width; i++) {
            if (i < pos) cout << "=";
            else if (i == pos) cout << ">";
            else cout << " ";
        }
        cout << "] " << int(progress * 100.0) << " %"<<"  Estimated time remaining: "<< avg_time_remaining<<" secs ("<<avg_time_remaining/60<<" mins) ("<<avg_time_remaining/3600<<" hrs)"<<flush;
         }

    else  cout << "] " << int(progress * 100.0) << " %"<<"  Estimated time remaining: N/A secs"<<flush;
     
}