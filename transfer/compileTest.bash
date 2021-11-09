# echo "There is no compilation or test script currently for the general THBI folder."
# echo "Not running THBI compileTest.bash"
echo "Compilation test script for whole THBI folder currently only involves testing HV code execution."
echo "If you are not using H/V, you can skip this step. "
cd ../matlab_to_hv_kernel/example 
matlab -nodisplay -nodesktop -nosplash -batch "run('test_compare_to_MINEOS.m')"
echo "See Figure THBI_ENAM/matlab_to_hv_kernel/example/executionTest.pdf to verify H/V code is working properly. "