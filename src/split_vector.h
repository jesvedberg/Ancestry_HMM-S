#ifndef __SPLIT_VECTOR_H
#define __SPLIT_VECTOR_H

/// splits vector (of recombination rates) into two vectors at the selected site.
/// the back vector is generated in reverse order

// min function because of namespace collision between std and arma
int int_min(int a, int b){
    int minout;
    if (a < b) {
        minout = a;
    }
    else {
        minout = b;
    }
    return minout;
}

void split_vector(int sel_site, vector<double> &whole_vec, vector<double> &back_vec, vector<double> &fwd_vec, cmd_line &options) 
{

    // trim vector if size is specified in morgans
    if (options.win_unit == "m") {
        double sum_morgans_fwd;
        double sum_morgans_back;
        for (int i = sel_site; i < whole_vec.size(); i++) {
            sum_morgans_fwd += whole_vec[i];
            if (sum_morgans_fwd <= options.win_morgan) {
                fwd_vec.push_back(whole_vec[i]) ;
            }
            else {
                break;
            }
        }
        for (int i = sel_site; i > 0; i--) {
            sum_morgans_back += whole_vec[i];
            if (sum_morgans_back <= options.win_morgan) {
                back_vec.push_back(whole_vec[i]) ;
            }
            else {
                break;
            }
        }
    }

    // trim vector if size is specified in percent
    else if (options.win_unit == "p") {
        int percent_size = whole_vec.size() * (options.win_percent/100);
        int trim_size_fwd = int_min((whole_vec.size() - sel_site) , percent_size); // Check for off by 1 error
        int trim_size_back = int_min(sel_site, percent_size);

        for (int i = sel_site; i < (sel_site + trim_size_fwd); i++) {
            fwd_vec.push_back(whole_vec[i]) ;
        }

        for (int i = sel_site; i > (sel_site - trim_size_back); i--) {
            back_vec.push_back(whole_vec[i]) ;
        }

    }

    //cout << "Split vector lengths: " << fwd_vec.size() << ", " << back_vec.size() << endl;
}

/*
// test. remove
void split_vector_int(int sel_site, vector<int> &whole_vec, vector<int> &back_vec, vector<int> &fwd_vec) 
{
    for (int i = sel_site; i < whole_vec.size(); i++) {
        fwd_vec.push_back(whole_vec[i]) ;
    }

    // Inverts vector from the selected site to the start. Removes the first site, as it has recombination rate 0.5
    for (int i = sel_site; i > 0; i--) {
        back_vec.push_back(whole_vec[i]) ;
    }
}
 */

#endif