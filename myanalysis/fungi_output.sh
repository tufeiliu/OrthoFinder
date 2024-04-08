#! /bin/sh

diamond_matrices=("BLOSUM45" "BLOSUM50" "BLOSUM62" "BLOSUM80" "BLOSUM90" "PAM30" "PAM70" "PAM250")
items=("Orthogroups" "Orthologues")

for matrix in "${diamond_matrices[@]}"; do 
    for item in "${items[@]}"; do
        echo 
        command="mkdir -p fungi_output/OrthoFinder/Results_Mar22_$matrix/$item"
        echo ">>>>> Run commond $command"
        result=$($command)
    done
done 


smdir=$(pwd)/$(ls | grep scoring_matrix)
if [ -d "$smdir" ]; then
	for file in "$smdir"/*; do
		if [ -f "$file" ]; then 
            database_filename=$(basename "$file") 
            filename="${database_filename%.*}"

            for item in "${items[@]}"; do
				echo
				custom_command="mkdir -p fungi_output/OrthoFinder/Results_Mar22_$filename/$item"
				echo ">>>>> Run commond $custom_command"
				custom_result=$($custom_command)
            done
		fi
	done
else
	echo "Substitution matrix directory not find!"
fi