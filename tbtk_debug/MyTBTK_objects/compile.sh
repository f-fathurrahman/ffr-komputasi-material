#!/bin/bash
set -x

BASNAM=`basename $1 .cpp`

# Note that macro name also should be changed from TBTK to MyTBTK

PREP1="-DMyTBTK_MATPLOTLIB_DO_NOT_FORCE_AGG "
PREP2="-DMyTBTK_USE_OPEN_MP -DMyTBTK_VERSION_GIT_HASH=\"\""
PREP3="-DMyTBTK_VERSION_MAJOR=2 -DMyTBTK_VERSION_MINOR=6 -DMyTBTK_VERSION_PATCH=4"
PREP4="-DMyTBTK_WRAP_PRIMITIVE_TYPES=1"

PREPALL="$PREP1 $PREP2 $PREP3 $PREP4"

INC1="-I/usr/include/python3.8"
INC2="-I/home/efefer/WORKS/my_github_repos/ffr-komputasi-material/tbtk_debug/include"

g++ -std=c++11 -c $PREPALL $INC1 $INC2 $1 -o $BASNAM.o

