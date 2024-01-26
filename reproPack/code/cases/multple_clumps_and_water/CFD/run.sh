mkdir 0
cp -f 0_org/* 0/
blockMesh
setFields
decomposePar -force
mpirun -np 12 intFlowIB -parallel > run.log 2>&1
#mpirun -np 12 intFlowTest -parallel > run.log 2>&1
#mpirun -np 12 intFoamIB -parallel > run.log 2>&1
