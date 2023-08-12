#!/bin/bash
set -x

BASNAM=`basename $1 .cpp`

PREP1="-DMyTBTK_MATPLOTLIB_DO_NOT_FORCE_AGG"
PREP2="-DMyTBTK_USE_OPEN_MP -DMyTBTK_VERSION_GIT_HASH=\"\""
PREP3="-DMyTBTK_VERSION_MAJOR=2 -DMyTBTK_VERSION_MINOR=6 -DMyTBTK_VERSION_PATCH=4 -DMyTBTK_WRAP_PRIMITIVE_TYPES=1"
PREPALL="$PREP1 $PREP2 $PREP3"

INC1="-I/usr/include/python3.8"
INC2="-I../include"

#LIB1="/home/efefer/mysoftwares/tbtk-2.6.4/lib/TBTK/libTBTK.a"
LIB1="../MyTBTK_objects/libMyTBTK.a"
LIB2="-lpython3.8 -llapack -lblas"

g++ $PREPALL $INC1 $INC2 $1 $LIB1 $LIB2 -o $BASNAM.x

