#! /bin/sh


# diamond_matrices=("BLOSUM45" "BLOSUM62" "BLOSUM80" "PAM30" "PAM70" "PAM250")
diamond_matrices=("BLOSUM45" "BLOSUM50" "BLOSUM62" "BLOSUM80" "BLOSUM90" "PAM30" "PAM70" "PAM250")

for matrix in "${diamond_matrices[@]}"; do 
	echo 
	command="orthofinder -f fungal_dataset --matrix $matrix"
	echo ">>>>> Run commond $command"
	result=$($command)
	echo
	echo "Command output:"
	echo "$result"
done 
