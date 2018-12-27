#!/bin/bash

echo "cs.nc"
diff cs.nc rot/cs.nc 
echo "lonslats.nc"
diff lonslats.nc rot/lonslats.nc 
echo "scrip.nc"
diff jhkim_ne030_unrotated.nc rot/.
