#!/bin/sh

cd $JANITHA_GIT_CODE/mesh-generator/src
make zero
make clean
mv $JANITHA_GIT_CODE/mesh-generator/src/mesh.out $JANITHA_GIT_CODE/bin/
# cd ../../bin
