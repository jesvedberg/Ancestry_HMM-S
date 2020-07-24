##  A MAKEFILE FOR The Ancestry HMM. TO USE, RUN THE COMMAND "make" VIA COMMAND LINE ##

TCFLAGS = -ltcmalloc 
ARMAFLAGS = -larmadillo

all:
	$(LINK.cc) -std=c++11 -O3 -o ahmm-s ahmms.cpp $(TCFLAGS) $(ARMAFLAGS)


## if you do not have a local install of google perftools, please simply remove the TCFlag link. 
###  $(LINK.cc) -std=c++11 -O3 -o ahmm-s ahmms.cpp $(ARMAFLAGS)

## if you have a local armadillo installation, you will need to provide the directory during compile time and possible also link lblas and lapack 
## our recommendation is to use miniconda3 to do the installation
##	$ conda install -c conda-forge armadillo
## then you will have the appropriate lib and include files in your home directory under subdirectory miniconda3/
## so, replace USERNAME with your unix id on the following line and try this

       #$(LINK.cc) -std=c++11 -O3 -o ahmm-s ahmms.cpp -L /home/USERNAME/miniconda3/lib/ -I /home/USERNAME/miniconda3/include/ $(ARMAFLAGS)
      
## if it builds correctly, you may also need to link the library during runtime
## to do this, add the following line to your ~/.bash_profile or ~/.bashrc
##     export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/USERNAME/miniconda3/lib/

## v0.1 - 2019-07-02 - algorithm working - first alpha
## v0.1.1 - 2019-07-11 - different algorithms for calculating the transtion rates was added
## v0.1.2 - 2019-10-03 - added stochastic trajectory function
## v1.0 - 2020-07-23 - ready to be released
