BASNAM=`basename $1 .f90`
rm -v $BASNAM.x
gfortran -fbounds-check -Wall $1 -o $BASNAM.x source/libmain.a

