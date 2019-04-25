#ifndef __SELECTION_TRAJECTORY_H
#define __SELECTION_TRAJECTORY_H

void selection_trajectory(vector<double> &freq, double s, int tt, double m, int generations, int n) 
{
    int t0 = 0;
    int t0min = 0;
    bool found = false;
    double f;

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
