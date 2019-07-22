#ifndef __PRINT_USAGE_H
#define __PRINT_USAGE_H

void print_usage() {
    
    cerr << endl << endl << "ancestry_hmm usage:" << endl << endl ;
    cerr << "\trequired:" << endl ;
    cerr << "\t\t-i [string]\t\tinput file name" << endl ;
    cerr << "\t\t-s [string]\t\tsample id and ploidy file" << endl ;
    cerr << "\t\t-a [int] [float] [float] ..." << endl ;
    cerr << "\t\t\tnumber of ancestral populations and ancestry proportion attributable to each" << endl ;
    cerr << "\t\t-p [int] [int] [float]" << endl ;
    cerr << "\t\t\tancestry pulse with format, ancestral population, time," << endl ;
    cerr << "\t\t\tand proportion of final ancestry from this pulse" << endl ;
    cerr << "\t\t\tnegative time or proportions indicate that parameters are to be estimated" << endl << endl ;
    
    cerr << "\toptional:" << endl ;
    cerr << "\t\t--help\t\t\tprint this help statement" << endl ;
    cerr << "\t\t--ne [int]\t\teffective population size of the admixed population" << endl ;
    cerr << "\t\t-g\t\t\tsamples are specified with genotypes rather than read counts" << endl ;
    cerr << "\t\t--precision [int]\tmodify float and double precision to int" << endl ;
    cerr << "\t\t-v\t\t\tviterbi decoding" << endl ;
    cerr << "\t\t-b [int] [int]\t\tnumber of bootstraps and bootstrap block size in number of SNPs" << endl ;
    cerr << "\t\t--tmax [int]\t\tmaximum time of an admixture pulse" << endl ;
    cerr << "\t\t--tmin [int]\t\tminimum time of an admixture pulse" << endl ;
    cerr << "\t\t--tolerance [float]\tdistance in lnL units to just convergence" << endl ;
    cerr << "\t\t-e [float]\t\terror rates" << endl ;
    cerr << "\t\t-E\t\t\tsite specific error rates are included" << endl ;
    cerr << "\t\t--fix\t\t\tancestral allele frequencies are certain" << endl << endl ;
    
    cerr << "\toptional and relevant only for multiple pulse models:" << endl ;
    cerr << "\t\t--output-ancestry\toutput ancestry posteriors rather than pulses" << endl ;
    cerr << "\t\t-r [int]\t\tnumber of random restarts during nelder-mead optimization" << endl ;
    cerr << "\t\t--pmax [int]\t\tmaximum proportion ancestry in an admixture pulse" << endl ;
    cerr << "\t\t--pmin [int]\t\tminimum proportion ancestry in an admixture pulse" << endl << endl;

    /// Parameters for adaptive introgression
    cerr << "\toptional and relevant only for detecting and quantifying selection:" << endl ;
    cerr << "\t\t--chr [string]" << endl ;
    cerr << "\t\t\t specify chromosome that will be analyzed" << endl ;
    cerr << "\t\t--chr_win [int] [int]" << endl ;
    cerr << "\t\t\t limit region on chromosome that will be analyzed (optional)" << endl ;
    cerr << "\t\t--grid [int] [int] [int] [float] [float] [float]" << endl ;
    cerr << "\t\t\t calculate likelihood rations in a grid. parameters: chromosomal position start, stop, step, selection coefficient start, stop, step." << endl ;
    cerr << "\t\t--gss [int] [int] [int] [float] [float]" << endl ;
    cerr << "\t\t\t golden section search for optimal selection coeffient at each site. parameters: chromosomal position start, stop, step, selection coefficient start, stop" << endl ;
    cerr << "\t\t--gss_precision [float]" << endl ;
    cerr << "\t\t\t specify precision in finding optimal value of s using golden section search. default: 1e-5" << endl ;
    cerr << "\t\t--unit_coords" << endl ;
    cerr << "\t\t\t unit for start and stop position in grid and gss search can be defined as chromosome coordinates rather than line in file. default off" << endl ;
    cerr << "\t\t--window [string] [float]" << endl ;
    cerr << "\t\t\t specify size of Markov chain in percent or Morgans. \"p 10\" extends the markov chain 10% of chromosome length on each side of selected site. \"m 0.1\" extends the windows 0.1 Morgan on each side of the selected site." << endl ;
    cerr << "\t\t--traj [int]" << endl ;
    cerr << "\t\t\t change algorithm for generating selection trajectories. 3: 3-point approximation, 4: 4-point approximation, default: forward iteration." << endl ;
    
}

#endif

