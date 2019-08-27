
export CXX="clang++"
export CFLAGS=" -I/Users/jlazar/bin/SQuIDS//include -I/usr/local/include -I/usr/local/Cellar/hdf5/1.10.4/lib/../include"
export CXXFLAGS=" -std=c++11"
export LDFLAGS="  -lnuSQuIDS -lpthread -L/Users/jlazar/bin/SQuIDS//lib -lSQuIDS -L/usr/local/lib -lgsl -lgslcblas -lm -L/usr/local/Cellar/hdf5/1.10.4/lib -lhdf5_hl -lhdf5 -L/usr/local/opt/szip/lib -lsz -lz -ldl -lm"

export DYLD_LIBRARY_PATH="/lib:/usr/lib:${HOME}/lib:/usr/local/lib:../lib:/Users/jlazar/bin/SQuIDS//lib:/usr/local/lib:/usr/local/Cellar/hdf5/1.10.4/lib"
