mkdir 0
cp -f 0_org/* 0/
blockMesh
setFields
decomposePar -force
mpirun -np 24 interFlowIB -parallel > run.log 2>&1
