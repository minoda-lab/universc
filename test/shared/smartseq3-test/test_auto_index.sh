R1_file=$1
R2_file=$(echo $R1_file | perl -pne 's/(.*)_R1/$1_R2/' )
I1_file=$(echo $R1_file | perl -pne 's/(.*)_R1/$1_I1/' )
I2_file=$(echo $R1_file | perl -pne 's/(.*)_R1/$1_I2/' )

echo $R1_file
echo $I1_file
echo $R2_file

# copies index 1 to next line (1st to 2nd) and deletes 3rd line
cat $R1_file | sed -E "s/ (.):(.):(.):(.*)\+(.*)$/ \1:\2:\3:\4+\5\n\4/g" | sed "3~5d" > $I1_file
indexlength=$(($(head $I1_file -n 2 | tail -n 1 | wc -c) -1))
echo $indexlength
qualscores=$(seq 1 $indexlength | xargs -I {} printf I)
echo $qualscores
sed -i "4~4s/^.*$/$qualscores/g" $I1_file
# coies index 2 to next line (1st to 2nd) and deletes 3rd line
cat $R1_file | sed -E "s/ (.):(.):(.):(.*)\+(.*)$/ \1:\2:\3:\4+\5\n\5/g" | sed "3~5d" >  $I2_file
index2length=$(($(head $I2_file -n 2 | tail -n 1 | wc -c) -1))
qualscores2=$(seq 1 $index2length | xargs -I {} printf I)
sed -i "4~4s/^.*$/({'$qualscores2'})/g" $I2_file

