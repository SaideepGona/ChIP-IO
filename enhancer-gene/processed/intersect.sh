bedtools multiinter \
-i converted_Astrocyte_EP.txt \
converted_Heart_EP.txt \
converted_Liver_EP.txt \
converted_Skeletal_muscle_EP.txt \
converted_Small_intestine_EP.txt \
-names brain heart liver skel smallint \
-header \
> intersectAHSIL