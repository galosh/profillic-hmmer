#!/bin/tcsh

### We begin by assuming that you are presently in the parent dir where you'll like all of the profillic-related subdirs to reside.  eg
# mkdir profillic-dirs
# cd profillic-dirs

# Note that for installing hmmer-3.1, you need a recent version of the autotools (2.63 or above)
autoconf --version
# If it's not high enough, update it.  See the script install-autotools-macosx.sh

## This assumes you already have boost installed.  If not, you'll need to follow the installation instructions at the Boost website.  You'll do something like this:
# Retrieve boost from http://sourceforge.net/projects/boost/files/boost
# Unpack it.  From one dir above it:
# cd boost_1_*
# ./bootstrap.sh 
# sudo bjam install
# cd ..

setenv BOOSTLIB /usr/local/lib
setenv BOOSTINC /usr/local/include
#setenv HMMER2 ~/src/profuse-hmmer

## Download latest seqan (follow instructions at http://trac.seqan.de/wiki/Tutorial/GettingStarted) to get the latest update from the seqan trunk  - here's what I do:
svn co http://svn.seqan.de/seqan/trunk seqan-trunk
## Download the latest hmmer snapshot:
svn checkout https://svn.janelia.org/eddylab/eddys/src/hmmer/trunk hmmer
## Get the rest from git:
git clone git://github.com/galosh/profillic.git 
git clone git://github.com/galosh/prolific.git 
git clone git://github.com/galosh/profuse.git 
git clone git://github.com/pedlefsen/HMMoC-BFloat-Algebra.git 
git clone git://github.com/galosh/profillic-hmmer

cd profillic
ln -s $BOOSTLIB boost-lib
ln -s $BOOSTINC boost-include
ln -s ../prolific
ln -s ../HMMoC-BFloat-Algebra
ln -s ../seqan-trunk
bjam install release
cd ..

cd profuse
ln -s $BOOSTLIB boost-lib
ln -s $BOOSTINC boost-include
ln -s ../prolific
ln -s ../HMMoC-BFloat-Algebra
ln -s ../seqan-trunk
bjam install release
cd ..

cd hmmer
autoconf
./configure
make
cd ..

cd profillic-hmmer
ln -s $BOOSTLIB boost-lib
ln -s $BOOSTINC boost-include
ln -s ../hmmer
ln -s ../seqan-trunk
ln -s ../prolific
ln -s ../HMMoC-BFloat-Algebra
make
cd ..
