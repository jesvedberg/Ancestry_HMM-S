#ifndef __SELECTION_FORWARD_H
#define __SELECTION_FORWARD_H

// Forward algoritm modified for selection inferrence.
double markov_chain::selection_forward_probabilities( map<int, vector<mat> > &transition_probabilites, vector<mat> &interploidy_transitions, selection &point, bool go_downstream  ) {
        
    /// return log likelihood which is sum of cts
    double lnl = 0 ;
    
    /// clear the fw probs matrix
    alphas.resize( transition_probabilites[ploidy_switch[0]].size() ) ;

    /// ploidy index to tract where in path we are
    int ploidy_index = 0 ;
    
    //// set all values to zero, but mostly just reize
    alphas[0].resize( transition_probabilites[ploidy_switch[0]][1].n_cols ) ;
    
    /// get initial state set
    alphas[0] = emission_probabilities[point.pos] * start_prob ;
    lnl += normalize( alphas[0] ) ;
    //cerr << "BEFORE: lnL: " << lnl <<  "  " << alphas[0] <<  "  " << alphas.size() <<  "  " << alphas[0].size() << endl;
    
    /// do all other sites
    /// Checks if going upstream or downstream from the selected site.
    if (go_downstream == true) {
        //cerr << "revloop1" << endl;
        selection_forward_loop_reverse(transition_probabilites, interploidy_transitions, point, lnl, ploidy_index) ;
    }
    else {
        //cerr << "fwdloop" << endl;
        selection_forward_loop(transition_probabilites, interploidy_transitions, point, lnl, ploidy_index) ;
    }
    //cerr << "BP4.2: selection forward 2." << endl;
    return lnl ;
}

// Loop in forward algorithm going downstream from the selected site
void markov_chain::selection_forward_loop( map<int, vector<mat> > &transition_probabilites, vector<mat> &interploidy_transitions, selection &point, double &lnl, int ploidy_index) {
    // WARNING: Check what +1 index does. May be unnecessary.
    //cerr << "Emission probabilities: " << emission_probabilities.size() << " " << point.pos << endl;

    int j = 1;

    for ( int i = point.pos+1 ; i < emission_probabilities.size() ; i ++ ) {
        //cerr << "fwdloop1" << endl;
        /// if we're at or past the next switch position
        bool ploidy_change = false ;
        /*if ( i >= ploidy_switch_position[ploidy_index+1] ) {
            ploidy_index ++ ;
            if ( ploidy_switch[ploidy_index] != ploidy_switch[ploidy_index-1] ) {
                ploidy_change = true ;
            }
        }*/
        /// resize matrix
        alphas[j].resize( transition_probabilites[ploidy_switch[ploidy_index]][1].n_cols ) ;
        //cerr << "n_cols: " << transition_probabilites[ploidy_switch[ploidy_index]].size() << endl;
        
        /// requires slightly different math if we are transitioning in ploidy between two adjacent sites
        if ( ploidy_change == true ) {
            /// transitions across a chromosome boundary will have low self-self rates
            if ( transition_probabilites[ploidy_switch[ploidy_index]][i](0,0) < 0.75 ) {
                alphas[j].fill( 1 ) ;
            }
            
            //// otherwise, this is a transition across ploidy types on the same chromosome use the interploidy transition rates
            else {
                alphas[j] = interploidy_transitions[ploidy_switch[ploidy_index-1]-1] * alphas[j-1] % emission_probabilities[i] ;
            }
        }
        
        /// otehrwise business as ususal
        else {
            alphas[j] = transition_probabilites[ploidy_switch[ploidy_index]][j] * alphas[j-1] % emission_probabilities[i] ;
        }

        /// normalize and updated likelihood
        lnl += normalize( alphas[j] ) ;
        //cerr << transition_probabilites[ploidy_switch[ploidy_index]][j] << endl;
        //cerr << i << "  " << j << "  " << "lnL: " << lnl <<  endl << alphas[j] << endl;
        j++;
    }
}

// Loop in forward algorithm going upstream from the selected site
void markov_chain::selection_forward_loop_reverse( map<int, vector<mat> > &transition_probabilites, vector<mat> &interploidy_transitions, selection &point, double &lnl, int ploidy_index) {
    // WARNING: Check what +1 index does. May be unnecessary.
    //cerr << "Emission probabilities: " << emission_probabilities.size() << " " << point.pos << endl;

    int j = 1;

    for ( int i = point.pos-1 ; i > 0 ; i -- ) {
        //cerr << "revloop1" << endl;
        /// if we're at or past the next switch position
        bool ploidy_change = false ;
        /*if ( i >= ploidy_switch_position[ploidy_index+1] ) {
            ploidy_index ++ ;
            if ( ploidy_switch[ploidy_index] != ploidy_switch[ploidy_index-1] ) {
                ploidy_change = true ;
            }
        }*/
        /// resize matrix
        alphas[j].resize( transition_probabilites[ploidy_switch[ploidy_index]][1].n_cols ) ;
        //cerr << "n_cols: " << transition_probabilites[ploidy_switch[ploidy_index]].size() << endl;
        
        /// requires slightly different math if we are transitioning in ploidy between two adjacent sites
        if ( ploidy_change == true ) {
            /// transitions across a chromosome boundary will have low self-self rates
            if ( transition_probabilites[ploidy_switch[ploidy_index]][i](0,0) < 0.75 ) {
                alphas[j].fill( 1 ) ;
            }
            
            //// otherwise, this is a transition across ploidy types on the same chromosome use the interploidy transition rates
            else {
                alphas[j] = interploidy_transitions[ploidy_switch[ploidy_index-1]-1] * alphas[j-1] % emission_probabilities[i] ;
            }
        }
        
        /// otehrwise business as ususal
        else {
            alphas[j] = transition_probabilites[ploidy_switch[ploidy_index]][j] * alphas[j-1] % emission_probabilities[i] ;
        }

        /// normalize and updated likelihood
        lnl += normalize( alphas[j] ) ;
        //cerr << transition_probabilites[ploidy_switch[ploidy_index]][j] << endl;
        //cerr << i << "  " << j << "  " << "lnL: " << lnl <<  endl << alphas[j] << endl;
        j++;
    }
}

#endif
