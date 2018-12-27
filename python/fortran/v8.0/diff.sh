#!/bin/bash

echo "diff org_dst/BackgroundMod.F90 dst/."
diff org_dst/BackgroundMod.F90 dst/.
echo "diff org_dst/CostFunMinMod.F90 dst/."
diff org_dst/CostFunMinMod.F90 dst/.
echo "diff org_dst/TestCode.F90 dst/."
diff org_dst/TestCode.F90 dst/.

