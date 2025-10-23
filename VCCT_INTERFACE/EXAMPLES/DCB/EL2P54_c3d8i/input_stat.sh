#!/bin/tcsh
rm -f *.txt
rm -f *.npy
rm -f *.npyc
rm -f *.pickle
cp /users0/ndecarva/_VCCT_Interface/CODE_MASTER/*.sh /users0/ndecarva/_VCCT_Interface/VERIFICATION/DCB_Q2D_software_release/
cp /users0/ndecarva/_VCCT_Interface/CODE_MASTER/*.f  /users0/ndecarva/_VCCT_Interface/VERIFICATION/DCB_Q2D_software_release/
cp /users0/ndecarva/_VCCT_Interface/CODE_MASTER/*.py /users0/ndecarva/_VCCT_Interface/VERIFICATION/DCB_Q2D_software_release/
python dir_name.py
chmod 777 *
abaqus j=DCB_CS_el2p54_C3D8I input=DCB_CS_el2p54_C3D8I user=UEL3D_INTERFACE_NVCCT.f cpus=1
