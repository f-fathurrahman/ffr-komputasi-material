#!/bin/bash
set -x

BASNAM=`basename $1 .cpp`

PREP1="-DTBTK_MATPLOTLIB_DO_NOT_FORCE_AGG"
PREP2="-DTBTK_USE_OPEN_MP -DTBTK_VERSION_GIT_HASH=\"\""
PREP3="-DTBTK_VERSION_MAJOR=2 -DTBTK_VERSION_MINOR=6 -DTBTK_VERSION_PATCH=4 -DTBTK_WRAP_PRIMITIVE_TYPES=1"
PREPALL="$PREP1 $PREP2 $PREP3"

INC1="-I/usr/include/python3.8"
INC2="-I/home/efefer/mysoftwares/tbtk-2.6.4/include"

LIB1="/home/efefer/mysoftwares/tbtk-2.6.4/lib/TBTK/libTBTK.a"
LIB2="-lpython3.8 -llapack -lblas"

g++ $PREPALL $INC1 $INC2 $1 $LIB1 $LIB2 -o $BASNAM.x

