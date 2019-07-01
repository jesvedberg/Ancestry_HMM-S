#ifndef __SELECTION_TRAJECTORY_H
#define __SELECTION_TRAJECTORY_H


// generates vector with allele frequency of selected allele over time/generations
void selection_trajectory(vector<double> &freq, double s, int tt, double m, int generations, int n) 
{
    // returns flat vector if selection is 0 (ie, no change in ellele frequency over time)
    if ( s == 0) {
        freq.assign(generations,0.5); // maybe use a different value than 0.5???
        return;
    }

    int t0 = 0;
    int t0min = 0;
    bool found = false;
    double f;

    // loops over generations until initial frequency (m) is reached + number of generations (gen) has passed
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
            freq.push_back(f);
        }
        t0++;
    }

    //return freq;
}   

#endif
