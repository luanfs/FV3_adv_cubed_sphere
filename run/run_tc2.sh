#!/bin/bash


#test case
tc="2"

# 1d advection scheme
hords=(0 8)

# grid types (0-equiedge; 1-equiangular)
gtypes=(2 2 2 2)

# departure point scheme
dps=(1 2 1 2)

# inner adv
iadvs=(1 2 1 2)

# mass fixers
mfs=(0 0 1 1)

# Initial value of N
N=48

#number of grids to be tested (we double the values N for each new grid and divide dt by 2)
Ng=3

cd ..
clear
make


case $tc in
    #------------- advection tests-------------
    1)
    dt=3600
    ;;

    2)
    dt=3600
    ;;

    3)
    dt=1800
    ;;

    4)
    dt=7200
    ;;
esac

# Loop over all schemes
size=${#dps[@]}
h=${#hords[@]}

# Loop on adv scheme
# Perform the loop from 1 to Ng
for ((j=1; j<=$Ng; j++)); do
    for ((k=0; k<=$h-1; k++)); do
        hord=${hords[k]}
        for ((i=0; i<=size-1; i++)); do
            dp=${dps[i]}
            iadv=${iadvs[i]}
            mf=${mfs[i]}
            gtype=${gtypes[i]}

            # Create input.par file
            echo "#Advection test case parameters" > input.par
            echo "#Test case" >> input.par
            echo "$tc" >> input.par
            echo "grid type (0-equiedge; 2-equiangular)" >> input.par
            echo "$gtype" >> input.par
            echo "# N (number of cells)" >> input.par
            echo "$N" >> input.par
            echo "#Time step" >> input.par
            echo "$dt" >> input.par
            echo "#hord scheme" >> input.par
            echo "$hord" >> input.par
            echo "#departure point (1 or 2)" >> input.par
            echo "$dp" >> input.par
            echo "#inner advection (1 or 2)" >> input.par
            echo "$iadv" >> input.par
            echo "#mass fixer (0-none, 1-flux average, 2-div projection)" >> input.par
            echo "$mf" >> input.par
            echo "# Number of plots" >> input.par
            echo "12" >> input.par
            echo "#-----------Description of parameters -------------------------------------" >> input.par
            echo '#Test case' >> input.par
            echo '#  case(2) - gaussian hill at corner with 45 degrees rotated zonal wind' >> input.par
            echo '#  case(3) - two gaussian hills and non divergent deformation wind (nair & lauritzen 2010)' >> input.par
            echo '#  case(4) - two gaussian hills and     divergent deformation wind (nair & lauritzen 2010)' >> input.par

            # Run the executable with input.par
            mv input.par par/input.par

            # Run the code
            ./main
        done
    done
    # Update the values of N and dt for the next iteration
    N=$((N * 2))
    dt=$((dt / 2))
done
