for w in amber03w  amber03star charmm22st; do 
echo $w
cd $w
python ../test.py
cd ..
done

