#!/bin/bash

echo "cs.nc"
diff cs.nc org/cs.nc 
echo "lonslats.nc"
diff lonslats.nc org/lonslats.nc 
echo "scrip.nc"
diff jhkim_ne005_unrotated.nc org
