BASNAM=`basename $1 .f90`
rm -v $BASNAM.x
gfortran -fbounds-check -Wall $1 -o $BASNAM.x source/libmain.a
if [ $? -eq 0 ]; then
    echo
    echo "*** Build process successful ***"
else
    echo
    echo "!!! Build process failed !!!"
fi
