#ifndef __SELECTION_TRAJECTORY_H
#define __SELECTION_TRAJECTORY_H

void selection_trajectory(vector<double> &freq, int &tt, double m, int generations, int n) 
{
    //double m = 0.1;
    double s = 0.01;
    //int n = 1000;
    //int generations = 200;
    double r = 0.0001;

    int t0 = 0;
    //int tt = 0;
    int t0min = 0;
    bool found = false;
    //int mind = 2*m;

    double f;
    //vector<double> freq;

    while (found == false) {
        f = 1 / (1 + 2 * n * s * exp(-s*t0));
        if (f > m) {
            if (tt == 0) {
                tt = t0;
                // cout << f << " " << tt << endl;
            }
            else if (t0 == (tt + generations)) {
                found = true;
            }
        }
        freq.push_back(f);
        t0++;
    }

    //return freq;
}   

#endif
