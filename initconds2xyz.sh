
au2ang=0.529
#Nwigner=10000 
Nwigner=$1

# create dir to put all the xyz files
xyzdir=xyz
mkdir -p $xyzdir

# read number of atoms
Nat=$(grep Natom initconds | sed 's/Natom//')

# equilibrium.xyz
fname='equilibrium.xyz'
echo $Nat > $xyzdir/$fname
echo " " >> $xyzdir/$fname
head -n $((9+$Nat)) initconds | tail -n $Nat | awk -v SF="$au2ang" '{ print $1, $3*SF, $4*SF, $5*SF }' >> $xyzdir/$fname

# strip away beginning part of initconds
# temporary file ready to loop over indices...
rm tmp
for i in {01..08}
do
	echo " " >> tmp
done
tail -n +$((12+$Nat)) initconds >> tmp

# start loop
for (( i=1; i<=$Nwigner; i++ ))
do
	#echo $i
	fname=$(printf "%06d" $i)_sample.xyz
	echo $Nat > $xyzdir/$fname
	echo " " >> $xyzdir/$fname
	head -n $(( $i*$((10+$Nat)) )) tmp | tail -n $Nat | awk -v SF="$au2ang" '{ print $1, $3*SF, $4*SF, $5*SF }' >> $xyzdir/$fname
done

# cleanup temporary files
rm tmp

