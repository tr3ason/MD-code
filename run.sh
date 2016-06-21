cat > orca_job${1} << EOL
! pm3
! engrad
* xyzfile 0 1 geom${1}.xyz
EOL

# RUN ORCA JOB
$HOME/software/orca_3_0_3_linux_x86-64/orca orca_job${1} > /dev/null
